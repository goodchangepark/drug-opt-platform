"""Deterministic Stage 4A evidence assembly and medicinal-chemistry strategy ranking.

The engine ranks transformations only.  It never applies a reaction or emits an analog.
"""

from __future__ import annotations

from datetime import datetime, timezone
from statistics import mean

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from sqlalchemy import select

from .activity_models import ActivityMeasurement, ActivityPrediction, AssayDefinition, MatchedMolecularPair
from .admet import ADMETEndpoint, ADMETMeasurement, ADMETPrediction
from .admet_predictor import MODEL_SPECS, comparable_experimental, metabolic_stability_assessment
from .metabolism import MetabolicPredictionRun
from .models import Compound, CompoundVersion, Project


ENGINE_NAME = "Stage 4A deterministic hit-optimization strategy engine"
ENGINE_VERSION = "4A.1.0"

OBJECTIVES = (
    "Improve potency", "Improve metabolic stability", "Improve solubility",
    "Improve permeability", "Reduce CYP inhibition", "Reduce hERG liability",
    "Reduce P-gp inhibition", "Balanced optimization", "Custom",
)

# Numerical weights are deliberately ordinal, not probabilities or endpoint scores.
EVIDENCE_HIERARCHY = (
    {"rank": 1, "type": "Experimental", "weight": 1.00},
    {"rank": 2, "type": "Project-specific validated model/SAR", "weight": 0.85},
    {"rank": 3, "type": "External validated quantitative model", "weight": 0.65},
    {"rank": 4, "type": "External classification model", "weight": 0.45},
    {"rank": 5, "type": "Rule-based hypothesis", "weight": 0.25},
)
EVIDENCE_WEIGHTS = {row["type"]: row["weight"] for row in EVIDENCE_HIERARCHY}
EVIDENCE_WEIGHTS["Calculated"] = 0.70
EVIDENCE_WEIGHTS["Manual override"] = 1.00


def _t(
    identifier, name, reaction, purpose, motif, expected, risk, source, liability_types,
    objectives, version="1.0",
):
    return {
        "id": identifier, "name": name, "reaction_smarts": reaction, "purpose": purpose,
        "applicable_motif": motif, "expected_effect": expected, "possible_risk": risk,
        "source": source, "version": version, "liability_types": liability_types,
        "objectives": objectives,
    }


MMP_SOURCE = "Hussain & Rea, J Chem Inf Model 2010, DOI 10.1021/ci900450m"
BIOISOSTERE_SOURCE = "SwissBioisostere, Nucleic Acids Res 2013, DOI 10.1093/nar/gks1059"
METABOLISM_SOURCE = "Johnson et al., J Med Chem 2020, DOI 10.1021/acs.jmedchem.9b01877"
HERG_SOURCE = "Waring & Johnstone, Bioorg Med Chem Lett 2007, DOI 10.1016/j.bmcl.2006.12.061"

TRANSFORMATION_LIBRARY = (
    _t("MET_F_FLUORINATION", "Fluorination at an oxidizable carbon", "[C;H1,H2,H3:1]>>[C:1]F", "Metabolism", "[C;H1,H2,H3]", "May block or redirect local oxidative metabolism", "May alter potency, lipophilicity, pKa, or create alternate metabolism", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("MET_STERIC_SHIELD", "Steric shielding near a soft spot", "[cH:1]>>[c:1]C", "Metabolism", "[cH]", "May reduce enzyme access to an oxidation-prone site", "Added lipophilicity or steric clash can reduce potency", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("MET_METHYL_REMOVAL", "Remove a labile methyl group", "[N,O,S:1][CH3:2]>>[N,O,S:1]", "Metabolism", "[N,O,S][CH3]", "Removes a potential N/O/S-dealkylation handle", "May alter basicity, solubility, conformation, or target interaction", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("MET_BENZYLIC_CH_REMOVAL", "Remove or contract benzylic CH", "[c:1][CH2:2][*:3]>>[c:1][*:3]", "Metabolism", "[c][CH2][*]", "Eliminates a benzylic C-H oxidation site", "Linker geometry and potency may change substantially", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("MET_HETEROATOM_REPLACEMENT", "Replace oxidizable carbon with nitrogen", "[cH:1]>>[n:1]", "Metabolism", "[cH]", "Can remove an aromatic C-H oxidation site and lower lipophilicity", "Electronic and binding geometry changes may reduce potency", BIOISOSTERE_SOURCE, ["metabolic_stability", "metabolic_soft_spot", "high_lipophilicity"], ["Improve metabolic stability", "Improve solubility"]),
    _t("MET_N_DEALK_BLOCK", "Block N-dealkylation handle", "[N:1][C;H2,H3:2]>>[N:1][C:2](F)", "Metabolism", "[N][C;H2,H3]", "May reduce N-dealkylation at the adjacent carbon", "Basicity and off-target profile can change", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("MET_O_DEALK_BLOCK", "Block O-dealkylation handle", "[O:1][C;H2,H3:2]>>[O:1][C:2](F)", "Metabolism", "[O][C;H2,H3]", "May reduce O-dealkylation at the adjacent carbon", "May redirect metabolism or alter potency", METABOLISM_SOURCE, ["metabolic_stability", "metabolic_soft_spot"], ["Improve metabolic stability"]),
    _t("LIPO_PHENYL_HETEROARYL", "Phenyl to pyridyl/heteroaryl replacement", "[cH:1]>>[n:1]", "Lipophilicity", "c1ccccc1", "Often lowers lipophilicity and may improve solubility", "Regioisomer, pKa, permeability, and potency effects are context-dependent", BIOISOSTERE_SOURCE, ["high_lipophilicity", "low_solubility", "herg", "pgp"], ["Improve solubility", "Reduce hERG liability", "Reduce P-gp inhibition"]),
    _t("LIPO_C_TO_HETERO", "Carbon to heteroatom replacement", "[CH2:1]>>[O:1]", "Lipophilicity", "[CH2]", "Can lower lipophilicity and add a hydrogen-bond acceptor", "Conformation, permeability, and metabolic route may change", BIOISOSTERE_SOURCE, ["high_lipophilicity", "low_solubility", "herg"], ["Improve solubility", "Reduce hERG liability"]),
    _t("LIPO_ALKYL_REDUCTION", "Reduce an alkyl substituent", "[*:1][CH2:2][CH3:3]>>[*:1][CH3:2]", "Lipophilicity", "[*][CH2][CH3]", "May reduce lipophilicity, MW, and nonspecific binding", "Hydrophobic target contacts may be lost", HERG_SOURCE, ["high_lipophilicity", "herg", "pgp"], ["Improve solubility", "Reduce hERG liability", "Reduce P-gp inhibition"]),
    _t("SOL_POLAR_SUBSTITUENT", "Introduce a small polar substituent", "[cH:1]>>[c:1]O", "Solubility", "[cH]", "May improve hydration and lower lipophilicity", "May reduce permeability or introduce conjugation metabolism", BIOISOSTERE_SOURCE, ["low_solubility", "high_lipophilicity"], ["Improve solubility"]),
    _t("SOL_BASIC_CENTER_MOD", "Attenuate a basic amine", "[N;H0;+0:1]>>[N:1]C(=O)", "Solubility/Safety", "[N;H0;+0]", "May reduce basicity-driven hERG risk and change ionization", "Neutralization can reduce solubility or target binding", HERG_SOURCE, ["herg", "high_lipophilicity"], ["Reduce hERG liability", "Balanced optimization"]),
    _t("SOL_AROMATICITY_REDUCTION", "Partially saturate an aromatic ring", "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[CH2:1]1[CH2:2][CH2:3][CH2:4][CH2:5][CH2:6]1", "Solubility", "c1ccccc1", "Increases Fsp3 and may improve solubility", "Large geometry and potency change; only a strategy, not an automatic proposal", BIOISOSTERE_SOURCE, ["low_solubility", "high_aromaticity"], ["Improve solubility", "Balanced optimization"]),
    _t("POT_BIOISOSTERE", "Local bioisostere replacement", "[C:1](=[O:2])[NH:3]>>[S:1](=[O:2])(=O)[NH:3]", "Potency-preserving exploration", "[C](=O)[NH]", "Offers a local replacement with potentially retained interaction geometry", "Bioisosteres are not potency-equivalent without project data", BIOISOSTERE_SOURCE, ["potency", "structural_alert"], ["Improve potency", "Balanced optimization"]),
    _t("POT_RING_BIOISOSTERE", "Ring bioisostere replacement", "c1ccccc1>>n1ccccc1", "Potency-preserving exploration", "c1ccccc1", "Explores ring electronics while retaining approximate topology", "Binding orientation and pKa may change", BIOISOSTERE_SOURCE, ["potency", "high_lipophilicity"], ["Improve potency", "Balanced optimization"]),
    _t("POT_LINKER_REPLACE", "Linker replacement", "[*:1][CH2:2][CH2:3][*:4]>>[*:1][CH2:2][O:3][*:4]", "Potency-preserving exploration", "[*][CH2][CH2][*]", "Can tune geometry, polarity, and flexibility", "May disrupt a protected binding vector", BIOISOSTERE_SOURCE, ["potency", "high_lipophilicity"], ["Improve potency", "Improve solubility", "Balanced optimization"]),
    _t("SAFE_ALERT_REMOVAL", "Remove the matched structural-alert motif", "[*:1][N+:2](=O)[O-:3]>>[*:1][NH2:2]", "Safety", "[N+](=O)[O-]", "Removes a nitro alert motif when that exact alert is present", "Replacement is not guaranteed to retain potency or remove all safety risk", BIOISOSTERE_SOURCE, ["structural_alert", "ames", "dili"], ["Balanced optimization"]),
)


def _confidence_value(value):
    return {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "UNKNOWN": 0, "EXPERIMENTAL": 3}.get(str(value).upper(), 0)


def _activity_summary(db, version_id, assay_id):
    rows = db.scalars(select(ActivityMeasurement).where(
        ActivityMeasurement.version_id == version_id,
        ActivityMeasurement.assay_id == assay_id,
    ).order_by(ActivityMeasurement.created_at)).all()
    if not rows:
        return None
    values = [row.normalized_value_nm for row in rows]
    return {
        "type": "Experimental", "hierarchy": 1, "n": len(values),
        "mean_nm": round(mean(values), 5), "unit": "nM",
        "measurements": [{"id": row.id, "value_nm": row.normalized_value_nm, "source": row.source} for row in rows],
        "confidence": "EXPERIMENTAL",
    }


def _version_project(db, version):
    compound = db.get(Compound, version.compound_row_id)
    return compound, db.get(Project, compound.project_id) if compound else None


def _version_endpoint_evidence(db, version_id, endpoint, endpoint_names):
    """Return the best compatible endpoint evidence for an arbitrary project version."""
    spec = MODEL_SPECS.get(endpoint, {})
    rows = db.scalars(select(ADMETMeasurement).where(
        ADMETMeasurement.version_id == version_id,
    ).order_by(ADMETMeasurement.created_at.desc())).all()
    for measurement in rows:
        normalized, note = comparable_experimental(endpoint, measurement, endpoint_names.get(measurement.endpoint_id, ""))
        if normalized is not None:
            return {
                "value": normalized, "unit": spec.get("unit", measurement.unit),
                "type": "Experimental", "hierarchy": 1, "confidence": "EXPERIMENTAL",
                "measurement_id": measurement.id, "conversion": note,
            }
    predictions = db.scalars(select(ADMETPrediction).where(
        ADMETPrediction.version_id == version_id,
    ).order_by(ADMETPrediction.created_at.desc())).all()
    prediction = next((row for row in predictions if row.model.endpoint_name == endpoint), None)
    if prediction is None:
        return None
    output = prediction.outputs_json or {}
    return {
        "value": prediction.predicted_value, "unit": prediction.unit,
        "classification": output.get("classification"),
        "type": "External classification model" if spec.get("prediction_type") == "binary_classification" else "External validated quantitative model",
        "hierarchy": 4 if spec.get("prediction_type") == "binary_classification" else 3,
        "confidence": prediction.confidence, "applicability_domain": prediction.applicability_domain,
        "prediction_id": prediction.id,
    }


def _mmp_endpoint_effect(endpoint, parent, other):
    if not parent or not other or parent.get("unit") != other.get("unit"):
        return None
    spec = MODEL_SPECS.get(endpoint, {})
    if spec.get("prediction_type") == "binary_classification":
        positive, negative = spec.get("positive_label"), spec.get("negative_label")
        if parent.get("classification") == positive and other.get("classification") == negative:
            direction = "IMPROVED"
        elif parent.get("classification") == negative and other.get("classification") == positive:
            direction = "WORSENED"
        else:
            direction = "UNCHANGED_OR_UNKNOWN"
        return {"endpoint": endpoint, "parent": parent, "other": other, "direction": direction, "delta": None}
    if parent.get("value") is None or other.get("value") is None:
        return None
    delta = float(other["value"]) - float(parent["value"])
    lower_is_better = endpoint.endswith("intrinsic clearance")
    higher_is_better = endpoint in {"Solubility", "Permeability"}
    if not (lower_is_better or higher_is_better):
        direction = "NOT_INTERPRETED"
    else:
        improvement = -delta if lower_is_better else delta
        direction = "IMPROVED" if improvement > 1e-12 else ("WORSENED" if improvement < -1e-12 else "UNCHANGED")
    return {"endpoint": endpoint, "parent": parent, "other": other, "direction": direction, "delta": round(delta, 6)}


def assemble_evidence(db, run):
    version = db.get(CompoundVersion, run.parent_version_id)
    compound, project = _version_project(db, version)
    properties = version.properties_json or {}
    activity = {"selected_assay": None, "experimental": None, "predicted": None, "nearest_compounds": [], "mmp": [], "activity_cliffs": []}
    if run.assay_id:
        assay = db.get(AssayDefinition, run.assay_id)
        activity["selected_assay"] = {
            "id": assay.id, "name": assay.name, "measurement_type": assay.measurement_type,
            "unit": assay.unit, "target": assay.target,
        }
        activity["experimental"] = _activity_summary(db, version.id, assay.id)
        prediction = db.scalar(select(ActivityPrediction).where(
            ActivityPrediction.version_id == version.id,
            ActivityPrediction.assay_id == assay.id,
        ).order_by(ActivityPrediction.created_at.desc()))
        if prediction:
            activity["predicted"] = {
                "type": "Project-specific validated model/SAR", "hierarchy": 2,
                "value_nm": prediction.predicted_value_nm, "unit": "nM",
                "pactivity": prediction.predicted_pactivity, "confidence": prediction.confidence,
                "applicability_domain": prediction.applicability_domain,
                "prediction_type": prediction.prediction_type, "prediction_id": prediction.id,
            }
            activity["nearest_compounds"] = prediction.nearest_neighbors or []
        pairs = db.scalars(select(MatchedMolecularPair).where(
            MatchedMolecularPair.assay_id == assay.id,
            (MatchedMolecularPair.version_a_id == version.id) | (MatchedMolecularPair.version_b_id == version.id),
        ).order_by(MatchedMolecularPair.created_at.desc())).all()
        for pair in pairs:
            other_id = pair.version_b_id if pair.version_a_id == version.id else pair.version_a_id
            other = db.get(CompoundVersion, other_id)
            other_compound = db.get(Compound, other.compound_row_id) if other else None
            row = {
                "pair_id": pair.id, "other_version_id": other_id,
                "other_compound": other_compound.compound_id if other_compound else "Unknown",
                "similarity": pair.similarity, "delta_pactivity": pair.delta_pactivity,
                "parent_is_a": pair.version_a_id == version.id,
                "transformation": pair.transformation_smiles, "is_cliff": pair.is_cliff,
                "type": "Project-specific validated model/SAR", "hierarchy": 2,
                "provenance": pair.provenance_json or {},
            }
            (activity["activity_cliffs"] if pair.is_cliff else activity["mmp"]).append(row)

    endpoints = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project.id))}
    measurements = db.scalars(select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id).order_by(ADMETMeasurement.created_at.desc())).all()
    predictions = db.scalars(select(ADMETPrediction).where(ADMETPrediction.version_id == version.id).order_by(ADMETPrediction.created_at.desc())).all()
    latest = {}
    for prediction in predictions:
        latest.setdefault(prediction.model.endpoint_name, prediction)
    admet = {}
    # Experimental evidence must remain usable even when no model was run or a
    # checkpoint is unavailable.  Scan all registered model endpoint definitions,
    # then add prediction-only endpoints.  comparable_experimental preserves the
    # endpoint/species/unit separation established in Stage 3.
    candidate_endpoints = list(dict.fromkeys([*MODEL_SPECS.keys(), *latest.keys()]))
    for endpoint in candidate_endpoints:
        prediction = latest.get(endpoint)
        spec = MODEL_SPECS.get(endpoint, {})
        compatible = []
        for measurement in measurements:
            normalized, note = comparable_experimental(endpoint, measurement, endpoints.get(measurement.endpoint_id, ""))
            if normalized is not None:
                compatible.append({
                    "measurement_id": measurement.id, "type": "Experimental", "hierarchy": 1,
                    "value": normalized, "unit": spec.get("unit", measurement.unit),
                    "classification": (
                        spec.get("positive_label") if normalized == 1.0 else spec.get("negative_label")
                    ) if spec.get("prediction_type") == "binary_classification" else None,
                    "original_value": measurement.mean_value if measurement.mean_value is not None else measurement.value,
                    "original_unit": measurement.unit, "conversion": note, "confidence": "EXPERIMENTAL",
                    "source": measurement.source, "method": measurement.method,
                    "record_type": "Experimental", "compound_version_id": version.id,
                })
        if not compatible and prediction is None:
            continue
        output = prediction.outputs_json or {} if prediction is not None else {}
        prediction_type = "External classification model" if spec.get("prediction_type") == "binary_classification" else "External validated quantitative model"
        model_provenance = (prediction.model.provenance_json or {}) if prediction is not None else {}
        predicted = ({
            "prediction_id": prediction.id, "type": prediction_type,
            "hierarchy": 4 if prediction_type == "External classification model" else 3,
            "value": prediction.predicted_value, "unit": prediction.unit,
            "classification": output.get("classification"), "probability": output.get("probability"),
            "confidence": prediction.confidence, "applicability_domain": prediction.applicability_domain,
            "assessment": output.get("metabolic_stability_assessment"),
            "model": prediction.model.model_name, "model_version": prediction.model.model_version,
            "record_type": "Predicted", "compound_version_id": version.id,
            "endpoint_definition": output.get("endpoint_definition") or model_provenance.get("endpoint_definition") or spec.get("endpoint_definition"),
            "training_dataset": output.get("training_dataset") or model_provenance.get("training_dataset") or spec.get("training_dataset"),
            "validation": output.get("validation") or model_provenance.get("validation") or spec.get("validation"),
            "license": output.get("license") or model_provenance.get("license") or spec.get("license"),
            "source": output.get("model_source") or model_provenance.get("source") or spec.get("source"),
            "timestamp": prediction.created_at.isoformat(),
        } if prediction is not None else None)
        preferred = compatible[0] if compatible else predicted
        if compatible and endpoint.endswith("intrinsic clearance"):
            preferred = {**compatible[0], "assessment": metabolic_stability_assessment(endpoint, compatible[0]["value"])}
        admet[endpoint] = {"experimental": compatible, "predicted": predicted, "preferred": preferred}
    raw_experimental = [{
        "measurement_id": row.id, "endpoint": endpoints.get(row.endpoint_id, ""),
        "value": row.mean_value if row.mean_value is not None else row.value,
        "unit": row.unit, "species": row.species, "matrix": row.matrix,
        "type": "Experimental", "hierarchy": 1,
    } for row in measurements]

    # Attach endpoint-specific project observations to actual MMP pairs.  This
    # does not infer a causal transformation effect; it records paired evidence
    # so an observed activity-tolerated ADMET direction can outrank generic rules.
    mmp_endpoints = (
        "Solubility", "Permeability", "HLM intrinsic clearance", "RLM intrinsic clearance",
        "hERG liability", "P-gp inhibitor",
    )
    for pair in [*activity["mmp"], *activity["activity_cliffs"]]:
        other_version_id = pair["other_version_id"]
        effects = []
        for endpoint in mmp_endpoints:
            parent_evidence = (admet.get(endpoint) or {}).get("preferred")
            other_evidence = _version_endpoint_evidence(db, other_version_id, endpoint, endpoints)
            effect = _mmp_endpoint_effect(endpoint, parent_evidence, other_evidence)
            if effect:
                effects.append(effect)
        pair["endpoint_effects"] = effects

    metabolism_run = db.scalar(select(MetabolicPredictionRun).where(
        MetabolicPredictionRun.version_id == version.id,
        MetabolicPredictionRun.status == "COMPLETE",
    ).order_by(MetabolicPredictionRun.completed_at.desc()))
    metabolism = {
        "soft_spots": [{
            "id": spot.id, "rank": spot.rank, "atom_index": spot.atom_index,
            "atom_environment": spot.atom_environment, "transformation": spot.transformation,
            "phase": spot.phase, "confidence": spot.confidence,
            "type": "Rule-based hypothesis", "hierarchy": 5,
        } for spot in sorted(metabolism_run.spots, key=lambda row: row.rank)] if metabolism_run else [],
        "metabolite_hypotheses": [{"id": item.id, "rank": item.rank, "transformation": item.transformation, "source_atom": item.source_atom, "confidence": item.confidence} for item in sorted(metabolism_run.metabolites, key=lambda row: row.rank)] if metabolism_run else [],
    }
    return {
        "project": {"id": project.id, "name": project.name},
        "parent": {"compound_id": compound.compound_id, "version_id": version.id, "version_number": version.version_number, "canonical_smiles": version.canonical_smiles},
        "activity": activity,
        "properties": {key: {"value": value, "type": "Calculated", "hierarchy": "Calculated", "engine": "RDKit"} for key, value in properties.items()},
        "structural_alerts": [{"name": row.get("alert_name"), "reason": row.get("reason"), "atom_indices": row.get("matched_atoms", []), "type": "Rule-based hypothesis", "hierarchy": 5} for row in (version.alerts_json or [])],
        "admet": admet, "raw_experimental_admet": raw_experimental,
        "metabolism": metabolism, "evidence_hierarchy": list(EVIDENCE_HIERARCHY),
        "assembled_at": datetime.now(timezone.utc).isoformat(),
    }


def _objective_bonus(objectives, names):
    return 18 if set(objectives).intersection(names) or "Balanced optimization" in objectives else 0


def identify_liabilities(evidence, objectives, constraints, endpoint_weights):
    liabilities = []

    def add(identifier, title, kind, source, severity, rationale, objective_names, endpoint=None, confidence=None, corroboration=None):
        evidence_type = source.get("type", "Rule-based hypothesis")
        weight = EVIDENCE_WEIGHTS.get(evidence_type, 0.25)
        endpoint_weight = float(endpoint_weights.get(endpoint or kind, 1.0))
        low_single = evidence_type == "External classification model" and str(confidence or source.get("confidence")).upper() == "LOW" and not corroboration
        score = (weight * 60 + severity * 20 + _objective_bonus(objectives, objective_names)) * endpoint_weight
        if low_single:
            score = min(score, 42.0)
        liabilities.append({
            "id": identifier, "title": title, "liability_type": kind, "endpoint": endpoint,
            "score": round(score, 2), "severity": severity, "evidence": [source],
            "evidence_type": evidence_type, "confidence": confidence or source.get("confidence", "UNKNOWN"),
            "rationale": rationale,
            "actionability": "SUPPORTING_ONLY" if low_single else "ACTIONABLE",
            "corroboration": corroboration or [],
        })

    activity = evidence["activity"]
    potency_limit = constraints.get("potency_max_nm")
    if potency_limit not in (None, ""):
        source = activity.get("experimental") or activity.get("predicted")
        value = source.get("mean_nm") if source and source.get("type") == "Experimental" else (source or {}).get("value_nm")
        if value is not None and float(value) > float(potency_limit):
            add("LIAB_POTENCY", "Potency constraint not met", "potency", source, min(3, float(value) / float(potency_limit)), f"Preferred activity {value:.4g} nM exceeds the {float(potency_limit):.4g} nM constraint.", ["Improve potency"], endpoint="Activity", confidence=source.get("confidence"))

    props = {key: row.get("value") for key, row in evidence["properties"].items()}
    clogp_limit = float(constraints.get("clogp_max", 4.0))
    if props.get("clogp") is not None and props["clogp"] > clogp_limit:
        add("LIAB_CLOGP", "High cLogP", "high_lipophilicity", {"type": "Calculated", "value": props["clogp"], "unit": "cLogP", "confidence": "HIGH"}, min(3, 1 + (props["clogp"] - clogp_limit) / 2), f"RDKit cLogP {props['clogp']:.2f} exceeds configured maximum {clogp_limit:.2f}.", ["Improve solubility", "Improve metabolic stability", "Reduce hERG liability", "Reduce P-gp inhibition"], endpoint="cLogP", confidence="HIGH")
    if constraints.get("mw_max") not in (None, "") and props.get("molecular_weight") is not None and props["molecular_weight"] > float(constraints["mw_max"]):
        add("LIAB_MW", "Molecular weight constraint exceeded", "high_mw", {"type": "Calculated", "value": props["molecular_weight"], "unit": "Da", "confidence": "HIGH"}, 2, "Calculated molecular weight exceeds the run constraint.", ["Balanced optimization"], endpoint="MW", confidence="HIGH")
    tpsa = props.get("tpsa")
    if tpsa is not None and ((constraints.get("tpsa_min") not in (None, "") and tpsa < float(constraints["tpsa_min"])) or (constraints.get("tpsa_max") not in (None, "") and tpsa > float(constraints["tpsa_max"]))):
        add("LIAB_TPSA", "TPSA outside configured range", "tpsa", {"type": "Calculated", "value": tpsa, "unit": "Å²", "confidence": "HIGH"}, 2, "Calculated TPSA is outside the configured range.", ["Improve permeability", "Improve solubility"], endpoint="TPSA", confidence="HIGH")

    positive_types = {
        "hERG liability": ("herg", "Reduce hERG liability"),
        "P-gp inhibitor": ("pgp", "Reduce P-gp inhibition"),
        "Ames mutagenicity": ("ames", "Balanced optimization"),
        "DILI clinical liability": ("dili", "Balanced optimization"),
    }
    for endpoint, row in evidence["admet"].items():
        preferred = row["preferred"]
        assessment = preferred.get("assessment") or {}
        if endpoint.endswith("intrinsic clearance") and assessment.get("category") == "UNSTABLE":
            add(f"LIAB_{endpoint.split()[0]}_UNSTABLE", f"{endpoint.split()[0]} metabolic instability", "metabolic_stability", preferred, 3, f"{preferred['type']} {endpoint} is classified UNSTABLE using the recorded threshold policy.", ["Improve metabolic stability"], endpoint=endpoint, confidence=preferred.get("confidence"))
        if endpoint == "Solubility" and preferred.get("value") is not None and float(preferred["value"]) < float(constraints.get("logs_min", -4.0)):
            add("LIAB_SOLUBILITY", "Low aqueous solubility", "low_solubility", preferred, 2, "Preferred LogS is below the configured screening threshold; it is not interpreted as pH-specific solubility.", ["Improve solubility"], endpoint=endpoint, confidence=preferred.get("confidence"))
        if endpoint == "Permeability" and preferred.get("value") is not None and float(preferred["value"]) < float(constraints.get("caco2_logpapp_min", -5.5)):
            add("LIAB_PERMEABILITY", "Low Caco-2 permeability", "low_permeability", preferred, 2, "Preferred Caco-2 LogPapp is below the configured screening threshold.", ["Improve permeability"], endpoint=endpoint, confidence=preferred.get("confidence"))
        classification = preferred.get("classification")
        if endpoint in positive_types and classification in {MODEL_SPECS[endpoint].get("positive_label"), "INHIBITOR"}:
            kind, objective = positive_types[endpoint]
            corroboration = []
            if kind in {"herg", "pgp"} and props.get("clogp", 0) > clogp_limit:
                corroboration.append("Calculated high cLogP is directionally consistent but not independent assay confirmation.")
            add(f"LIAB_{kind.upper()}", f"Potential {endpoint} concern", kind, preferred, 2, f"Positive {endpoint} classification. Probability is not a quantitative potency value.", [objective], endpoint=endpoint, confidence=preferred.get("confidence"), corroboration=corroboration)
        if endpoint.startswith("CYP") and endpoint.endswith("inhibitor") and classification == "INHIBITOR":
            add(f"LIAB_{endpoint.replace(' ', '_')}", f"Potential {endpoint} concern", "cyp", preferred, 2, "Positive CYP inhibitor classification; no IC50 is inferred.", ["Reduce CYP inhibition"], endpoint=endpoint, confidence=preferred.get("confidence"))

    for index, alert in enumerate(evidence["structural_alerts"]):
        add(f"LIAB_ALERT_{index}", f"Structural alert: {alert.get('name') or 'matched motif'}", "structural_alert", alert, 2, alert.get("reason") or "RDKit structural-alert rule match.", ["Balanced optimization"], endpoint="Structural alert", confidence="LOW")
    if evidence["metabolism"]["soft_spots"] and any(item["liability_type"] == "metabolic_stability" for item in liabilities):
        spot = evidence["metabolism"]["soft_spots"][0]
        add("LIAB_SOFT_SPOT", f"Metabolic soft spot: {spot['transformation']}", "metabolic_soft_spot", spot, 1.5, "Rule-based atom-level hypothesis supports, but does not prove, the clearance liability.", ["Improve metabolic stability"], endpoint="Metabolic soft spot", confidence=spot["confidence"], corroboration=["Compound-level microsomal instability is also present."])

    liabilities.sort(key=lambda row: (-row["score"], row["id"]))
    for rank, row in enumerate(liabilities, 1):
        row["rank"] = rank
    return liabilities


def _changed_parent_atoms(parent_smiles, other_smiles):
    parent, other = Chem.MolFromSmiles(parent_smiles), Chem.MolFromSmiles(other_smiles)
    if not parent or not other:
        return []
    result = rdFMCS.FindMCS([parent, other], timeout=2, ringMatchesRingOnly=True, completeRingsOnly=True)
    query = Chem.MolFromSmarts(result.smartsString) if result.smartsString else None
    match = parent.GetSubstructMatch(query) if query else ()
    changed = sorted(set(range(parent.GetNumAtoms())) - set(match))
    if not changed and match:
        changed = sorted({atom.GetIdx() for index in match for atom in parent.GetAtomWithIdx(index).GetNeighbors() if atom.GetIdx() not in match})
    return changed


def infer_regions(db, run, evidence):
    version = db.get(CompoundVersion, run.parent_version_id)
    mol = Chem.MolFromSmiles(version.canonical_smiles)
    overrides = run.manual_overrides_json or {}
    protected, modifiable = [], []
    for atoms in overrides.get("protect_atoms", []):
        atom_list = atoms if isinstance(atoms, list) else [atoms]
        protected.append({"id": f"MANUAL_PROTECT_{'_'.join(map(str, atom_list))}", "status": "DO NOT MODIFY", "atom_indices": atom_list, "reason": "Manual override", "risk": "HIGH", "confidence": "HIGH", "source": "Manual override"})
    activity = evidence["activity"]
    if run.assay_id:
        parent_activity = _activity_summary(db, version.id, run.assay_id)
        for cliff in activity["activity_cliffs"]:
            other = db.get(CompoundVersion, cliff["other_version_id"])
            other_activity = _activity_summary(db, other.id, run.assay_id) if other else None
            if not parent_activity or not other_activity or parent_activity["mean_nm"] >= other_activity["mean_nm"]:
                continue
            atoms = _changed_parent_atoms(version.canonical_smiles, other.canonical_smiles)
            if atoms:
                protected.append({"id": f"CLIFF_{cliff['pair_id']}", "status": "HIGH-RISK TO MODIFY", "atom_indices": atoms, "reason": f"Project activity cliff: modifying this region in {cliff['other_compound']} reduced potency by at least {abs(cliff['delta_pactivity']):.2f} pActivity.", "risk": "HIGH", "confidence": "HIGH", "source": "Project experimental SAR/activity cliff"})
    protected_atoms = {atom for row in protected for atom in row["atom_indices"]}
    allow_atoms = {int(atom) for group in overrides.get("allow_atoms", []) for atom in (group if isinstance(group, list) else [group])}
    protected_atoms -= allow_atoms
    protected = [{**row, "atom_indices": [atom for atom in row["atom_indices"] if atom not in allow_atoms]} for row in protected]
    protected = [row for row in protected if row["atom_indices"]]

    for spot in evidence["metabolism"]["soft_spots"][:5]:
        atom = int(spot["atom_index"])
        if atom not in protected_atoms:
            modifiable.append({"id": f"SOFT_SPOT_{spot['id']}", "atom_indices": [atom], "fragment": spot["atom_environment"], "reason": f"Rank {spot['rank']} {spot['transformation']} soft-spot hypothesis", "risk": "MEDIUM", "confidence": spot["confidence"], "source": "Rule-based metabolic soft spot"})
    for index, alert in enumerate(evidence["structural_alerts"]):
        atoms = [int(atom) for atom in alert.get("atom_indices", []) if int(atom) not in protected_atoms]
        if atoms:
            modifiable.append({"id": f"ALERT_{index}", "atom_indices": atoms, "fragment": alert.get("name"), "reason": "Structural-alert motif is a candidate for removal", "risk": "HIGH", "confidence": "LOW", "source": "Rule-based structural alert"})
    clogp = (evidence["properties"].get("clogp") or {}).get("value", 0)
    if clogp and clogp > float((run.constraints_json or {}).get("clogp_max", 4.0)):
        aromatic = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic() and atom.GetSymbol() == "C" and atom.GetIdx() not in protected_atoms][:6]
        if aromatic:
            modifiable.append({"id": "LIPOPHILIC_AROMATIC_REGION", "atom_indices": aromatic, "fragment": "aromatic carbon region", "reason": "Calculated high cLogP; heteroaryl or polarity scan may be considered", "risk": "MEDIUM", "confidence": "MEDIUM", "source": "Calculated property + rule"})
    for pair in activity["mmp"]:
        other = db.get(CompoundVersion, pair["other_version_id"])
        atoms = _changed_parent_atoms(version.canonical_smiles, other.canonical_smiles) if other else []
        atoms = [atom for atom in atoms if atom not in protected_atoms]
        if atoms and abs(pair["delta_pactivity"]) <= 0.3:
            modifiable.append({"id": f"MMP_{pair['pair_id']}", "atom_indices": atoms, "fragment": pair["transformation"], "reason": f"Project MMP {pair['other_compound']} changed this region with <=0.3 pActivity shift", "risk": "LOW", "confidence": "HIGH", "source": "Project experimental MMP"})
    for group in overrides.get("allow_atoms", []):
        atoms = group if isinstance(group, list) else [group]
        modifiable.append({"id": f"MANUAL_ALLOW_{'_'.join(map(str, atoms))}", "atom_indices": atoms, "fragment": "manual selection", "reason": "Manual override allows modification", "risk": "USER_ACCEPTED", "confidence": "HIGH", "source": "Manual override"})
    if not modifiable:
        terminals = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1 and atom.GetDegree() == 1 and atom.GetIdx() not in protected_atoms][:3]
        for atom in terminals:
            modifiable.append({"id": f"EXPOSED_{atom}", "atom_indices": [atom], "fragment": mol.GetAtomWithIdx(atom).GetSymbol(), "reason": "Terminal exposed substituent heuristic; no project SAR support", "risk": "UNKNOWN", "confidence": "LOW", "source": "Rule-based hypothesis"})
    if not protected:
        protected = [{"id": "PROTECTION_UNKNOWN", "status": "UNKNOWN", "atom_indices": [], "reason": "Insufficient project activity-cliff/SAR evidence; no protected pharmacophore is asserted.", "risk": "UNKNOWN", "confidence": "UNKNOWN", "source": "Evidence gap"}]
    return protected, modifiable


def rank_transformations(db, run, evidence, liabilities, protected, modifiable):
    mol = Chem.MolFromSmiles(evidence["parent"]["canonical_smiles"])
    overrides = run.manual_overrides_json or {}
    excluded = set(overrides.get("exclude_transformations", []))
    prioritized = set(overrides.get("prioritize_transformations", []))
    protected_atoms = {atom for row in protected for atom in row.get("atom_indices", [])}
    candidates = []
    for liability in liabilities:
        if liability["actionability"] == "SUPPORTING_ONLY":
            continue
        for rule in TRANSFORMATION_LIBRARY:
            if rule["id"] in excluded or liability["liability_type"] not in rule["liability_types"]:
                continue
            motif = Chem.MolFromSmarts(rule["applicable_motif"])
            matches = list(mol.GetSubstructMatches(motif)) if motif else []
            if not matches:
                continue
            region_atoms = {atom for region in modifiable for atom in region["atom_indices"]}
            match = next((item for item in matches if set(item).intersection(region_atoms)), matches[0])
            overlap = sorted(set(match).intersection(protected_atoms))
            project_mmp = []
            for pair in evidence["activity"]["mmp"]:
                if abs(pair["delta_pactivity"]) > 0.3:
                    continue
                other = db.get(CompoundVersion, pair["other_version_id"])
                changed = _changed_parent_atoms(evidence["parent"]["canonical_smiles"], other.canonical_smiles) if other else []
                if set(changed).intersection(match):
                    project_mmp.append(pair)
            score = liability["score"] * 0.65 + (18 if project_mmp else 0) + (45 if rule["id"] in prioritized else 0)
            score -= 30 if overlap else 0
            confidence = "MEDIUM" if liability["evidence_type"] in {"Experimental", "Project-specific validated model/SAR", "Calculated"} or project_mmp else "LOW"
            candidates.append({
                **{key: rule[key] for key in ("id", "name", "reaction_smarts", "purpose", "applicable_motif", "expected_effect", "possible_risk", "source", "version")},
                "score": round(score, 2), "source_atom_indices": list(match),
                "target_liability": liability["id"], "target_liability_rank": liability["rank"],
                "potency_risk": "HIGH" if overlap else ("LOW" if project_mmp else "MEDIUM"),
                "protected_overlap": overlap,
                "evidence": [f"{liability['evidence_type']}: {liability['title']}"] + ([f"Project MMP tolerated nearby chemistry ({len(project_mmp)} pair(s))"] if project_mmp else []),
                "confidence": confidence, "manual_priority": rule["id"] in prioritized,
                "application_status": "STRATEGY_ONLY — no analog generated",
            })
    for pair in evidence["activity"]["mmp"]:
        if abs(pair["delta_pactivity"]) > 0.3 or "MMP_PROJECT_OBSERVED" in excluded:
            continue
        other = db.get(CompoundVersion, pair["other_version_id"])
        atoms = _changed_parent_atoms(evidence["parent"]["canonical_smiles"], other.canonical_smiles) if other else []
        improved = [row for row in pair.get("endpoint_effects", []) if row["direction"] == "IMPROVED"]
        candidates.append({
            "id": f"MMP_PROJECT_OBSERVED_{pair['pair_id']}", "name": f"Project-observed MMP toward {pair['other_compound']}",
            "reaction_smarts": f"{evidence['parent']['canonical_smiles']}>>{other.canonical_smiles}",
            "purpose": "Potency-preserving project MMP", "applicable_motif": "Project pair",
            "expected_effect": f"Observed project potency shift {pair['delta_pactivity']:+.3f} pActivity" + ("; paired improvement in " + ", ".join(row["endpoint"] for row in improved) if improved else ""),
            "possible_risk": "Project observation may not transfer to other substitutions or endpoints",
            "source": MMP_SOURCE, "version": "project-observed-1", "score": 95.0 + min(15.0, 5.0 * len(improved)),
            "source_atom_indices": atoms, "target_liability": "PROJECT_MMP", "target_liability_rank": 0,
            "potency_risk": "LOW", "protected_overlap": sorted(set(atoms).intersection(protected_atoms)),
            "evidence": ["Project experimental MMP"] + [f"Paired {row['endpoint']} {row['direction']} ({row['parent']['type']} → {row['other']['type']})" for row in improved], "confidence": "HIGH", "manual_priority": False,
            "application_status": "REFERENCE STRATEGY ONLY — existing project compound, no analog generated",
        })
    deduplicated = {}
    for item in candidates:
        key = (item["id"], tuple(item["source_atom_indices"]))
        if key not in deduplicated or item["score"] > deduplicated[key]["score"]:
            deduplicated[key] = item
    ranked = sorted(deduplicated.values(), key=lambda row: (-row["score"], row["id"]))[:15]
    for rank, row in enumerate(ranked, 1):
        row["rank"] = rank
    return ranked


def render_strategy_svg(smiles, protected, modifiable, soft_spots):
    mol = Chem.MolFromSmiles(smiles)
    protected_atoms = {atom for row in protected for atom in row.get("atom_indices", [])}
    modifiable_atoms = {atom for row in modifiable for atom in row.get("atom_indices", [])} - protected_atoms
    soft_atoms = {int(row["atom_index"]) for row in soft_spots} - protected_atoms - modifiable_atoms
    colors = {atom: (0.86, 0.22, 0.22) for atom in protected_atoms}
    colors.update({atom: (0.95, 0.58, 0.12) for atom in modifiable_atoms})
    colors.update({atom: (0.48, 0.25, 0.78) for atom in soft_atoms})
    drawer = rdMolDraw2D.MolDraw2DSVG(650, 430)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=sorted(colors), highlightAtomColors=colors, highlightAtomRadii={atom: 0.42 for atom in colors})
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def analyze_run(db, run):
    evidence = assemble_evidence(db, run)
    liabilities = identify_liabilities(evidence, run.objectives_json or [], run.constraints_json or {}, run.endpoint_weights_json or {})
    protected, modifiable = infer_regions(db, run, evidence)
    transformations = rank_transformations(db, run, evidence, liabilities, protected, modifiable)
    run.evidence_json = evidence
    run.liabilities_json = liabilities
    run.protected_regions_json = protected
    run.modifiable_regions_json = modifiable
    run.transformations_json = transformations
    run.highlighted_svg = render_strategy_svg(evidence["parent"]["canonical_smiles"], protected, modifiable, evidence["metabolism"]["soft_spots"][:3])
    run.engine_name = ENGINE_NAME
    run.engine_version = ENGINE_VERSION
    run.status = "COMPLETE"
    run.message = f"Ranked {len(liabilities)} liabilities and {len(transformations)} transformation strategies; no analog structures generated."
    run.completed_at = datetime.now(timezone.utc)
    return run
