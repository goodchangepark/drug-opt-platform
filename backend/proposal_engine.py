"""Deterministic Stage 4B analog generation, staged rescoring, and ranking.

No language model or random molecular generator is used.  Every structure is a
single or, when enabled, two-step application of an explicit Stage 4A strategy.
"""

from __future__ import annotations

import base64
import math
import pickle
from datetime import datetime, timezone
from functools import lru_cache

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFMCS, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Contrib.SA_Score import sascorer
from sqlalchemy import select

from .activity_models import ActivityMeasurement, AssayDefinition, QSARModel
from .admet_predictor import MODEL_SPECS, model_files_available, predict_endpoint
from .chemistry import ENGINE_VERSION as RDKIT_VERSION, analyze_smiles
from .database import SessionLocal
from .metabolic_soft_spot import ENGINE_VERSION as SOFT_SPOT_VERSION, predict_soft_spots
from .models import Compound, CompoundVersion
from .optimization import OptimizationRun
from .optimization_engine import ENGINE_VERSION as STRATEGY_VERSION, _version_endpoint_evidence
from .proposal import (
    CandidatePredictionSnapshot, CandidateRanking, CandidateRejectionReason,
    CandidateTransformation, OptimizationCandidate, OptimizationProposalRun,
)
from .qsar import (
    DESCRIPTOR_NAMES, FINGERPRINT_CONFIG, applicability, feature_vector,
    fingerprint_and_descriptors, nearest_neighbors, tanimoto_similarity,
    value_from_pactivity,
)


ENGINE_NAME = "Stage 4B deterministic analog proposal engine"
ENGINE_VERSION = "4B.1.0"
RANDOM_SEED = 42  # recorded for reproducibility; generation itself is deterministic.

EXECUTABLE_TRANSFORMATIONS = {
    "MET_F_FLUORINATION": "targeted_atom_addition",
    "MET_STERIC_SHIELD": "targeted_atom_addition",
    "MET_METHYL_REMOVAL": "targeted_atom_deletion",
    "MET_BENZYLIC_CH_REMOVAL": "targeted_linker_contraction",
    "MET_HETEROATOM_REPLACEMENT": "targeted_atom_replacement",
    "MET_N_DEALK_BLOCK": "targeted_alpha_fluorination",
    "MET_O_DEALK_BLOCK": "targeted_alpha_fluorination",
    "LIPO_PHENYL_HETEROARYL": "targeted_atom_replacement",
    "LIPO_C_TO_HETERO": "targeted_atom_replacement",
    "LIPO_ALKYL_REDUCTION": "targeted_atom_deletion",
    "SOL_POLAR_SUBSTITUENT": "targeted_atom_addition",
    "POT_RING_BIOISOSTERE": "targeted_atom_replacement",
    "POT_LINKER_REPLACE": "targeted_atom_replacement",
}
STRATEGY_ONLY_TRANSFORMATIONS = {
    "SOL_BASIC_CENTER_MOD": "No single context-independent acylation/debasicification edit preserves valid valence and target binding.",
    "SOL_AROMATICITY_REDUCTION": "Whole-ring saturation changes geometry/stereochemistry and requires explicit medicinal-chemistry design.",
    "POT_BIOISOSTERE": "Amide-to-sulfonamide atom mapping is not safely context-independent.",
    "SAFE_ALERT_REMOVAL": "Alert remediation must be alert-specific; the generic strategy is not executed automatically.",
}

DEFAULT_WEIGHTS = {
    "Activity": 2.0, "HLM intrinsic clearance": 1.5, "RLM intrinsic clearance": 0.6,
    "Solubility": 1.2, "Permeability": 1.0, "CYP": 0.8, "hERG liability": 1.2,
    "P-gp inhibitor": 0.7, "Structural alerts": 1.0, "Synthetic complexity": 0.7,
    "Similarity": 0.8,
}
OBJECTIVE_WEIGHT_BOOSTS = {
    "Improve potency": {"Activity": 1.8},
    "Improve metabolic stability": {"HLM intrinsic clearance": 1.8, "RLM intrinsic clearance": 1.4},
    "Improve solubility": {"Solubility": 1.8},
    "Improve permeability": {"Permeability": 1.8},
    "Reduce CYP inhibition": {"CYP": 1.8},
    "Reduce hERG liability": {"hERG liability": 1.8},
    "Reduce P-gp inhibition": {"P-gp inhibitor": 1.8},
}
CONFIDENCE_FACTOR = {"EXPERIMENTAL": 1.0, "HIGH": 1.0, "MEDIUM": 0.75, "LOW": 0.5, "UNKNOWN": 0.35, "NOT_AVAILABLE": 0.0}
DOMAIN_FACTOR = {"IN_DOMAIN": 1.0, "IN DOMAIN": 1.0, "BORDERLINE": 0.7, "OUT_OF_DOMAIN": 0.3, "OUT OF DOMAIN": 0.3, "UNKNOWN": 0.5}


def _canonical_product(rw_mol):
    try:
        product = rw_mol.GetMol()
        Chem.SanitizeMol(product)
        if len(Chem.GetMolFrags(product)) != 1:
            return None
        if abs(Chem.GetFormalCharge(product)) > 2:
            return None
        return product
    except Exception:
        return None


def _add_atom(parent, target, atomic_number, bond_type=Chem.BondType.SINGLE):
    rw = Chem.RWMol(parent)
    new_index = rw.AddAtom(Chem.Atom(atomic_number))
    rw.AddBond(int(target), new_index, bond_type)
    return _canonical_product(rw)


def _replace_atom(parent, target, atomic_number, aromatic=False):
    rw = Chem.RWMol(parent)
    atom = rw.GetAtomWithIdx(int(target))
    atom.SetAtomicNum(atomic_number)
    atom.SetFormalCharge(0)
    atom.SetIsAromatic(aromatic)
    atom.SetNoImplicit(False)
    return _canonical_product(rw)


def _remove_atom(parent, target):
    rw = Chem.RWMol(parent)
    rw.RemoveAtom(int(target))
    return _canonical_product(rw)


def _contract_atom(parent, target):
    atom = parent.GetAtomWithIdx(int(target))
    neighbors = [item.GetIdx() for item in atom.GetNeighbors()]
    if len(neighbors) != 2:
        return None
    rw = Chem.RWMol(parent)
    first, second = neighbors
    if rw.GetBondBetweenAtoms(first, second) is None:
        rw.AddBond(first, second, Chem.BondType.SINGLE)
    rw.RemoveAtom(int(target))
    return _canonical_product(rw)


def _candidate_sites(mol, strategy, allowed_atoms):
    motif = Chem.MolFromSmarts(strategy.get("applicable_motif") or "")
    matches = list(mol.GetSubstructMatches(motif)) if motif else []
    preferred = set(strategy.get("source_atom_indices") or [])
    allowed = set(allowed_atoms or [])
    sites = []
    for match in matches:
        if allowed and not set(match).intersection(allowed):
            continue
        if preferred and set(match).intersection(preferred):
            sites.insert(0, match)
        else:
            sites.append(match)
    seen = set()
    return [item for item in sites if not (item in seen or seen.add(item))]


def execute_strategy(smiles, strategy, allowed_atoms=None):
    """Return sanitized, single-edit products with exact source-site provenance."""
    identifier = strategy["id"]
    if identifier.startswith("MMP_PROJECT_OBSERVED_"):
        parts = strategy.get("reaction_smarts", "").split(">>", 1)
        product = Chem.MolFromSmiles(parts[1]) if len(parts) == 2 else None
        return [(product, strategy.get("source_atom_indices") or [])] if product else []
    if identifier not in EXECUTABLE_TRANSFORMATIONS:
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    products = []
    for match in _candidate_sites(mol, strategy, allowed_atoms):
        product, changed = None, []
        if identifier == "MET_F_FLUORINATION":
            target = next((idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetTotalNumHs() > 0), None)
            product, changed = (_add_atom(mol, target, 9), [target]) if target is not None else (None, [])
        elif identifier == "MET_STERIC_SHIELD":
            target = next((idx for idx in match if mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetTotalNumHs() > 0), None)
            product, changed = (_add_atom(mol, target, 6), [target]) if target is not None else (None, [])
        elif identifier == "MET_METHYL_REMOVAL":
            target = next((idx for idx in reversed(match) if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetDegree() == 1), None)
            product, changed = (_remove_atom(mol, target), [target]) if target is not None else (None, [])
        elif identifier == "MET_BENZYLIC_CH_REMOVAL":
            target = match[1] if len(match) >= 3 else None
            product, changed = (_contract_atom(mol, target), [target]) if target is not None else (None, [])
        elif identifier in {"MET_HETEROATOM_REPLACEMENT", "LIPO_PHENYL_HETEROARYL", "POT_RING_BIOISOSTERE"}:
            target = next((idx for idx in match if mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetSymbol() == "C" and mol.GetAtomWithIdx(idx).GetTotalNumHs() > 0), None)
            product, changed = (_replace_atom(mol, target, 7, True), [target]) if target is not None else (None, [])
        elif identifier in {"MET_N_DEALK_BLOCK", "MET_O_DEALK_BLOCK"}:
            target = match[1] if len(match) > 1 and mol.GetAtomWithIdx(match[1]).GetTotalNumHs() > 0 else None
            product, changed = (_add_atom(mol, target, 9), [target]) if target is not None else (None, [])
        elif identifier in {"LIPO_C_TO_HETERO", "POT_LINKER_REPLACE"}:
            search = match if identifier == "LIPO_C_TO_HETERO" else match[1:]
            target = next((idx for idx in search if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetDegree() == 2 and not mol.GetAtomWithIdx(idx).GetIsAromatic()), None)
            product, changed = (_replace_atom(mol, target, 8, False), [target]) if target is not None else (None, [])
        elif identifier == "LIPO_ALKYL_REDUCTION":
            target = match[-1] if mol.GetAtomWithIdx(match[-1]).GetAtomicNum() == 6 and mol.GetAtomWithIdx(match[-1]).GetDegree() == 1 else None
            product, changed = (_remove_atom(mol, target), [target]) if target is not None else (None, [])
        elif identifier == "SOL_POLAR_SUBSTITUENT":
            target = next((idx for idx in match if mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetTotalNumHs() > 0), None)
            product, changed = (_add_atom(mol, target, 8), [target]) if target is not None else (None, [])
        if product is not None:
            products.append((product, changed))
    unique = {}
    for product, changed in products:
        canonical = Chem.MolToSmiles(product, canonical=True, isomericSmiles=False)
        unique.setdefault(canonical, (product, changed))
    return list(unique.values())


def _mcs_mapping(parent, candidate):
    result = rdFMCS.FindMCS(
        [parent, candidate], timeout=3, ringMatchesRingOnly=True,
        completeRingsOnly=True, matchChiralTag=False,
    )
    query = Chem.MolFromSmarts(result.smartsString) if result.smartsString else None
    return (parent.GetSubstructMatch(query), candidate.GetSubstructMatch(query)) if query else ((), ())


def structure_change(parent, candidate):
    parent_match, candidate_match = _mcs_mapping(parent, candidate)
    changed_parent_set = set(range(parent.GetNumAtoms())) - set(parent_match)
    changed_candidate_set = set(range(candidate.GetNumAtoms())) - set(candidate_match)
    parent_to_candidate = dict(zip(parent_match, candidate_match))
    candidate_to_parent = {candidate_index: parent_index for parent_index, candidate_index in parent_to_candidate.items()}
    # MCS can retain every parent atom for an atom addition.  Attribute each
    # added candidate atom to its mapped parent neighbour so protected-region
    # checks see the actual attachment site.
    for candidate_index in list(changed_candidate_set):
        for neighbor in candidate.GetAtomWithIdx(candidate_index).GetNeighbors():
            if neighbor.GetIdx() in candidate_to_parent:
                changed_parent_set.add(candidate_to_parent[neighbor.GetIdx()])
    for parent_index, candidate_index in parent_to_candidate.items():
        parent_atom, candidate_atom = parent.GetAtomWithIdx(parent_index), candidate.GetAtomWithIdx(candidate_index)
        if parent_atom.GetAtomicNum() != candidate_atom.GetAtomicNum() or parent_atom.GetFormalCharge() != candidate_atom.GetFormalCharge() or parent_atom.GetDegree() != candidate_atom.GetDegree():
            changed_parent_set.add(parent_index)
    changed_parent = sorted(changed_parent_set)
    changed_candidate = sorted(changed_candidate_set)
    parent_coverage = len(parent_match) / max(parent.GetNumHeavyAtoms(), 1)
    return changed_parent, changed_candidate, parent_coverage, parent_to_candidate


def stereochemistry_preserved(parent, candidate, mapping):
    Chem.AssignStereochemistry(parent, force=True, cleanIt=True)
    Chem.AssignStereochemistry(candidate, force=True, cleanIt=True)
    # Only encoded stereochemistry can be preserved. Unspecified tetrahedral
    # geometry is not silently promoted to a stereochemical assertion.
    parent_centers = dict(Chem.FindMolChiralCenters(parent, includeUnassigned=False, useLegacyImplementation=False))
    candidate_centers = dict(Chem.FindMolChiralCenters(candidate, includeUnassigned=False, useLegacyImplementation=False))
    for parent_index, label in parent_centers.items():
        candidate_index = mapping.get(parent_index)
        if candidate_index is None or candidate_index not in candidate_centers:
            return False
        candidate_label = candidate_centers[candidate_index]
        parent_neighbors = sorted(atom.GetAtomicNum() for atom in parent.GetAtomWithIdx(parent_index).GetNeighbors())
        candidate_neighbors = sorted(atom.GetAtomicNum() for atom in candidate.GetAtomWithIdx(candidate_index).GetNeighbors())
        # If substituent priorities did not change, an R/S change is a true
        # inversion. If priorities changed, R/S may change without geometry loss.
        if parent_neighbors == candidate_neighbors and candidate_label != label:
            return False
    return True


def _draw_difference(mol, atoms, color):
    drawer = rdMolDraw2D.MolDraw2DSVG(430, 300)
    colors = {idx: color for idx in atoms}
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=atoms, highlightAtomColors=colors, highlightAtomRadii={idx: 0.42 for idx in atoms})
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def chemical_validation(parent_smiles, candidate_smiles, protected_atoms):
    parent, candidate = Chem.MolFromSmiles(parent_smiles), Chem.MolFromSmiles(candidate_smiles)
    if parent is None or candidate is None:
        return {"valid": False, "code": "INVALID_VALENCE_OR_SANITIZATION", "detail": "RDKit sanitization/valence/aromaticity validation failed."}
    if len(Chem.GetMolFrags(candidate)) != 1:
        return {"valid": False, "code": "FRAGMENTED_STRUCTURE", "detail": "Candidate contains more than one disconnected fragment."}
    canonical = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
    parent_canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=False)
    candidate_isomeric = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
    parent_isomeric = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    if canonical == parent_canonical and candidate_isomeric == parent_isomeric:
        return {"valid": False, "code": "PARENT_DUPLICATE", "detail": "Transformation did not change the canonical structure."}
    changed_parent, changed_candidate, coverage, mapping = structure_change(parent, candidate)
    protected_overlap = sorted(set(changed_parent).intersection(protected_atoms))
    if protected_overlap:
        return {"valid": False, "code": "PROTECTED_REGION_MODIFIED", "detail": f"Changed protected parent atoms {protected_overlap}.", "changed_parent_atoms": changed_parent}
    if not stereochemistry_preserved(parent, candidate, mapping):
        return {"valid": False, "code": "STEREOCHEMISTRY_LOSS", "detail": "A mapped parent stereocenter was lost or inverted."}
    if abs(Chem.GetFormalCharge(candidate)) > 2:
        return {"valid": False, "code": "CHARGE_SANITY", "detail": "Absolute formal charge exceeds the configured sanity limit of 2."}
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    similarity = float(DataStructs.TanimotoSimilarity(generator.GetFingerprint(parent), generator.GetFingerprint(candidate)))
    return {
        "valid": True, "canonical_smiles": canonical,
        "isomeric_smiles": Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True),
        "inchikey": Chem.MolToInchiKey(candidate), "similarity": similarity, "mcs_coverage": coverage,
        "changed_parent_atoms": changed_parent, "changed_candidate_atoms": changed_candidate,
        "mapping": mapping,
        "parent_difference_svg": _draw_difference(parent, changed_parent, (0.82, 0.18, 0.18)),
        "candidate_difference_svg": _draw_difference(candidate, changed_candidate, (0.08, 0.55, 0.72)),
    }


def synthetic_feasibility(parent_smiles, candidate_smiles):
    parent, candidate = Chem.MolFromSmiles(parent_smiles), Chem.MolFromSmiles(candidate_smiles)
    score = float(sascorer.calculateScore(candidate))
    parent_complexity, complexity = Descriptors.BertzCT(parent), Descriptors.BertzCT(candidate)
    parent_chiral = len(Chem.FindMolChiralCenters(parent, includeUnassigned=True))
    candidate_chiral = len(Chem.FindMolChiralCenters(candidate, includeUnassigned=True))
    ring_complexity = sum(1 for ring in candidate.GetRingInfo().AtomRings() if len(ring) >= 7)
    if score <= 3.0 and ring_complexity == 0 and candidate_chiral <= parent_chiral:
        classification = "LOW SYNTHETIC COMPLEXITY"
    elif score <= 5.0 and ring_complexity <= 1 and candidate_chiral <= parent_chiral + 1:
        classification = "MODERATE SYNTHETIC COMPLEXITY"
    else:
        classification = "HIGH SYNTHETIC COMPLEXITY"
    return {
        "classification": classification, "sa_score": round(score, 4),
        "method": "Ertl-Schuffenhauer SA score surrogate + deterministic complexity checks",
        "not_synthesis_success_probability": True, "large_ring_count": ring_complexity,
        "stereocenter_change": candidate_chiral - parent_chiral,
        "bertz_complexity": round(float(complexity), 3),
        "delta_bertz_complexity": round(float(complexity - parent_complexity), 3),
    }


def _activity_summary(db, version_id, assay_id):
    rows = db.scalars(select(ActivityMeasurement).where(
        ActivityMeasurement.version_id == version_id, ActivityMeasurement.assay_id == assay_id,
    )).all()
    if not rows:
        return None
    value = sum(row.normalized_value_nm for row in rows) / len(rows)
    return {"record_type": "Experimental", "value_nm": value, "pactivity": -math.log10(value * 1e-9), "unit": "nM", "confidence": "EXPERIMENTAL", "applicability_domain": "IN_DOMAIN", "n": len(rows)}


def predict_candidate_activity(db, project_id, assay_id, smiles, existing_version_id=None):
    if not assay_id:
        return {"status": "MODEL_UNAVAILABLE", "reason": "No assay selected", "record_type": "Unavailable"}
    if existing_version_id:
        experimental = _activity_summary(db, existing_version_id, assay_id)
        if experimental:
            return {"status": "COMPLETE", **experimental, "prediction_preserved": False}
    assay = db.get(AssayDefinition, assay_id)
    target_mol, target_fp, target_desc, _ = fingerprint_and_descriptors(smiles)
    dataset = {"rows": [], "fingerprints": [], "descriptors": []}
    compounds = db.scalars(select(Compound).where(Compound.project_id == project_id)).all()
    for compound in compounds:
        version = next((row for row in compound.versions if row.version_number == compound.current_version), None)
        if not version:
            continue
        summary = _activity_summary(db, version.id, assay_id)
        if not summary:
            continue
        _, fingerprint, descriptors, _ = fingerprint_and_descriptors(version.canonical_smiles)
        dataset["rows"].append({"compound_id": compound.compound_id, "version_id": version.id, "activity_nm": summary["value_nm"], "pactivity": summary["pactivity"]})
        dataset["fingerprints"].append(fingerprint)
        dataset["descriptors"].append([descriptors[name] for name in DESCRIPTOR_NAMES])
    neighbors = nearest_neighbors(target_fp, dataset)
    domain, confidence, max_similarity, outside = applicability(
        neighbors, target_desc,
        {"descriptors": np.array(dataset["descriptors"]) if dataset["descriptors"] else np.empty((0, len(DESCRIPTOR_NAMES)))}
    )
    model_row = db.scalar(select(QSARModel).where(QSARModel.assay_id == assay_id).order_by(QSARModel.created_at.desc()))
    if model_row and len(dataset["rows"]) >= 15:
        fitted = pickle.loads(base64.b64decode(model_row.pickle_data))
        predicted_p = float(np.asarray(fitted["model"].predict(np.vstack([feature_vector(target_fp, target_desc)])))[0])
        prediction_type, uncertainty = f"QSAR {fitted['name']}", None
        model_version = model_row.model_uid
    elif len(dataset["rows"]) >= 5:
        selected = neighbors[:min(5, len(neighbors))]
        weights = np.array([row["similarity"] ** 4 for row in selected])
        predicted_p = float(np.average(np.array([row["pactivity"] for row in selected]), weights=weights))
        prediction_type = "Similarity nearest neighbor"
        uncertainty = float(np.std([row["pactivity"] for row in selected]) / max(len(selected), 1) ** 0.5) if len(selected) > 1 else 0.75
        model_version = f"RDKit-{RDKIT_VERSION}-project-similarity"
    else:
        return {"status": "MODEL_UNAVAILABLE", "reason": "Fewer than five project experimental compounds", "record_type": "Unavailable", "nearest_neighbors": neighbors, "applicability_domain": domain, "confidence": "UNKNOWN"}
    if domain == "OUT OF DOMAIN":
        confidence = "LOW"
    return {
        "status": "COMPLETE", "record_type": "Predicted", "prediction_type": prediction_type,
        "value_nm": value_from_pactivity(predicted_p), "pactivity": predicted_p, "unit": "nM",
        "confidence": confidence, "applicability_domain": domain, "nearest_neighbors": neighbors,
        "uncertainty": uncertainty, "model_version": model_version, "max_similarity": max_similarity,
        "descriptor_outside_training_space": outside,
    }


@lru_cache(maxsize=4096)
def _cached_admet(smiles, endpoint, model_version):
    return predict_endpoint(smiles, endpoint)


def rescore_admet(db, candidate, project_id):
    from .admet import ADMETEndpoint

    results = {}
    endpoint_names = {row.id: row.name for row in db.scalars(select(ADMETEndpoint).where(ADMETEndpoint.project_id == project_id))}
    for endpoint, spec in MODEL_SPECS.items():
        if candidate.existing_version_id:
            existing = _version_endpoint_evidence(db, candidate.existing_version_id, endpoint, endpoint_names)
            if existing and existing.get("type") == "Experimental":
                result = {"status": "COMPLETE", "record_type": "Experimental", "predicted_value": existing["value"], "unit": existing["unit"], "confidence": "EXPERIMENTAL", "applicability_domain": "IN_DOMAIN", **existing}
            else:
                result = None
        else:
            result = None
        if result is None:
            available, reason = model_files_available(endpoint)
            if not available:
                results[endpoint] = {"status": "MODEL_UNAVAILABLE", "reason": reason, "scored": False}
                continue
            try:
                result = dict(_cached_admet(candidate.canonical_smiles, endpoint, spec["model_version"]))
                result["record_type"] = "Predicted"
            except Exception as exc:
                results[endpoint] = {"status": "MODEL_UNAVAILABLE", "reason": f"Inference failed: {exc}", "scored": False}
                continue
        result["scored"] = result.get("status") == "COMPLETE"
        results[endpoint] = result
        if result["scored"]:
            db.add(CandidatePredictionSnapshot(
                candidate_id=candidate.id, stage="Stage 3", endpoint=endpoint,
                record_type=result.get("record_type", "Predicted"),
                value_json={key: value for key, value in result.items() if key not in {"applicability_domain_details"}},
                unit=result.get("unit", spec.get("unit", "")), model_name=spec["display_name"],
                model_version=spec["model_version"], confidence=result.get("confidence", "UNKNOWN"),
                applicability_domain=(result.get("applicability_domain") or {}).get("classification", "UNKNOWN") if isinstance(result.get("applicability_domain"), dict) else result.get("applicability_domain", "UNKNOWN"),
                provenance_json={"endpoint_definition": spec["endpoint_definition"], "dataset": spec["training_dataset"], "validation": spec["validation"], "license": spec["license"], "source": spec["source"], "compound_version_type": "proposal_snapshot"},
            ))
    return results


def _preferred_value(row):
    if not row:
        return None
    value = row.get("predicted_value", row.get("value"))
    return float(value) if value is not None else None


def _domain(result):
    domain = result.get("applicability_domain", "UNKNOWN") if result else "UNKNOWN"
    return domain.get("classification", "UNKNOWN") if isinstance(domain, dict) else domain


def _confidence(result):
    return str((result or {}).get("confidence", "UNKNOWN")).upper()


def _quality(endpoint, result):
    value = _preferred_value(result)
    if value is None:
        return None
    if endpoint == "Solubility":
        return max(0.0, min(1.0, (value + 6.0) / 5.0))
    if endpoint == "Permeability":
        return max(0.0, min(1.0, (value + 7.0) / 3.0))
    if endpoint.endswith("intrinsic clearance"):
        return max(0.0, min(1.0, (3.0 - value) / 3.0))
    if MODEL_SPECS.get(endpoint, {}).get("prediction_type") == "binary_classification":
        return max(0.0, min(1.0, 1.0 - value))
    return None


def _reject(db, candidate, code, detail, stage, evidence_type="Calculated"):
    candidate.status, candidate.rejection_stage = "REJECTED", stage
    db.add(CandidateRejectionReason(candidate=candidate, code=code, detail=detail, stage=stage, hard_constraint=True, evidence_type=evidence_type))


def _cheap_constraints(db, candidate, parent_analysis, constraints, protected_atoms, seen):
    validation = chemical_validation(parent_analysis["identity"]["canonical_smiles"], candidate.canonical_smiles, protected_atoms)
    if not validation["valid"]:
        _reject(db, candidate, validation["code"], validation["detail"], "FILTERING")
        return False
    candidate.canonical_smiles, candidate.isomeric_smiles, candidate.inchikey = validation["canonical_smiles"], validation["isomeric_smiles"], validation["inchikey"]
    candidate.parent_similarity, candidate.mcs_coverage = validation["similarity"], validation["mcs_coverage"]
    candidate.changed_parent_atoms_json, candidate.changed_candidate_atoms_json = validation["changed_parent_atoms"], validation["changed_candidate_atoms"]
    candidate.parent_difference_svg, candidate.candidate_difference_svg = validation["parent_difference_svg"], validation["candidate_difference_svg"]
    if candidate.canonical_smiles in seen:
        retained = seen.get(candidate.canonical_smiles) if isinstance(seen, dict) else None
        if retained is not None:
            for transformation in candidate.transformations:
                db.add(CandidateTransformation(
                    candidate=retained, sequence_number=transformation.sequence_number,
                    transformation_id=transformation.transformation_id, name=transformation.name,
                    reaction_smarts=transformation.reaction_smarts,
                    transformation_version=transformation.transformation_version, source=transformation.source,
                    source_atom_indices_json=transformation.source_atom_indices_json,
                    changed_parent_atoms_json=transformation.changed_parent_atoms_json,
                    execution_status="DUPLICATE_HYPOTHESIS_MERGED",
                    provenance_json={**(transformation.provenance_json or {}), "duplicate_candidate_number": candidate.candidate_number},
                ))
            retained.why_generated = (retained.why_generated + " · Alternative equivalent hypothesis: " + candidate.hypothesis).strip(" ·")
        _reject(db, candidate, "DUPLICATE_ANALOG", "Canonical structure duplicates an earlier generated candidate.", "FILTERING")
        return False
    if isinstance(seen, dict): seen[candidate.canonical_smiles] = candidate
    else: seen.add(candidate.canonical_smiles)
    analysis = analyze_smiles(candidate.isomeric_smiles)
    candidate.stage1_json = analysis
    candidate.structure_svg = analysis["svg"]
    properties, parent_properties = analysis["properties"], parent_analysis["properties"]
    keys = ("molecular_weight", "clogp", "tpsa", "hbd", "hba", "rotatable_bonds", "fraction_csp3", "qed")
    candidate.property_delta_json = {key: round(float(properties[key]) - float(parent_properties[key]), 5) for key in keys if isinstance(properties.get(key), (int, float)) and isinstance(parent_properties.get(key), (int, float))}
    candidate.synthetic_feasibility_json = synthetic_feasibility(parent_analysis["identity"]["canonical_smiles"], candidate.canonical_smiles)
    gates = (
        ("molecular_weight", "mw_max", "EXCESSIVE_MOLECULAR_WEIGHT", "MW"),
        ("clogp", "clogp_max", "CLOGP_CONSTRAINT", "cLogP"),
    )
    for property_name, constraint_name, code, label in gates:
        if constraints.get(constraint_name) is not None and properties[property_name] > float(constraints[constraint_name]):
            _reject(db, candidate, code, f"{label} {properties[property_name]} exceeds hard maximum {constraints[constraint_name]}.", "FILTERING")
            return False
    if constraints.get("tpsa_min") is not None and properties["tpsa"] < float(constraints["tpsa_min"]):
        _reject(db, candidate, "TPSA_CONSTRAINT", f"TPSA {properties['tpsa']} is below {constraints['tpsa_min']} Å².", "FILTERING"); return False
    if constraints.get("tpsa_max") is not None and properties["tpsa"] > float(constraints["tpsa_max"]):
        _reject(db, candidate, "TPSA_CONSTRAINT", f"TPSA {properties['tpsa']} exceeds {constraints['tpsa_max']} Å².", "FILTERING"); return False
    minimum_similarity = float(constraints.get("similarity_min", 0.0))
    if candidate.parent_similarity < minimum_similarity:
        _reject(db, candidate, "LOW_PARENT_SIMILARITY", f"Morgan/Tanimoto {candidate.parent_similarity:.3f} is below {minimum_similarity:.3f}.", "FILTERING"); return False
    parent_alerts = {(row["alert_set"], row["alert_name"]) for row in parent_analysis.get("alerts", [])}
    new_alerts = [row for row in analysis.get("alerts", []) if (row["alert_set"], row["alert_name"]) not in parent_alerts]
    if constraints.get("no_new_structural_alert", True) and new_alerts:
        _reject(db, candidate, "NEW_STRUCTURAL_ALERT", "New alerts: " + ", ".join(row["alert_name"] for row in new_alerts), "FILTERING", "Rule-based hypothesis"); return False
    candidate.status = "FILTERED"
    return True


def _activity_gate(db, candidate, optimization_run, parent_activity, constraints=None):
    candidate.activity_json = predict_candidate_activity(db, candidate.project_id, optimization_run.assay_id, candidate.canonical_smiles, candidate.existing_version_id)
    activity = candidate.activity_json
    if activity.get("status") == "COMPLETE":
        db.add(CandidatePredictionSnapshot(
            candidate_id=candidate.id, stage="Stage 2", endpoint="Activity", record_type=activity["record_type"],
            value_json=activity, unit="nM", model_name=activity.get("prediction_type", "Project experimental activity"),
            model_version=activity.get("model_version", "project-experimental"), confidence=activity.get("confidence", "UNKNOWN"),
            applicability_domain=activity.get("applicability_domain", "UNKNOWN"),
            provenance_json={"source": "Project assay/QSAR/SAR/MMP", "assay_id": optimization_run.assay_id, "experimental_priority": True},
        ))
    constraints = constraints or optimization_run.constraints_json or {}
    fold = constraints.get("do_not_worsen_fold")
    potency_max = constraints.get("potency_max_nm")
    if potency_max and activity.get("status") == "COMPLETE" and activity.get("value_nm") is not None and _domain(activity) != "OUT OF DOMAIN":
        if float(activity["value_nm"]) > float(potency_max):
            _reject(db, candidate, "POTENCY_CONSTRAINT", f"Candidate {activity['value_nm']:.4g} nM exceeds absolute potency constraint {float(potency_max):.4g} nM.", "PREDICTING", activity.get("record_type", "Predicted"))
            return False
    if fold and parent_activity and activity.get("status") == "COMPLETE" and activity.get("value_nm") is not None:
        limit = float(parent_activity["value_nm"]) * float(fold)
        if float(activity["value_nm"]) > limit and _domain(activity) != "OUT OF DOMAIN":
            _reject(db, candidate, "ACTIVITY_WORSENED", f"Candidate {activity['value_nm']:.4g} nM exceeds parent-based {float(fold):.3g}-fold limit {limit:.4g} nM.", "PREDICTING", activity.get("record_type", "Predicted"))
            return False
    return True


def _post_prediction_constraints(db, candidate, optimization_run, constraints=None):
    constraints = constraints or optimization_run.constraints_json or {}
    for endpoint, constraint_name, code in (
        ("Solubility", "logs_min", "SOLUBILITY_CONSTRAINT"),
        ("Permeability", "caco2_logpapp_min", "PERMEABILITY_CONSTRAINT"),
    ):
        threshold, result = constraints.get(constraint_name), (candidate.admet_json or {}).get(endpoint, {})
        value = _preferred_value(result)
        if threshold is not None and value is not None and result.get("status") == "COMPLETE" and _domain(result) != "OUT_OF_DOMAIN" and _confidence(result) != "LOW":
            if value < float(threshold):
                _reject(db, candidate, code, f"{endpoint} {value:.4g} is below configured minimum {float(threshold):.4g}.", "PREDICTING", result.get("record_type", "Predicted"))
                return False
    # Installed hERG predictions are intentionally LOW confidence in Stage 3F.
    # The do-not-increase setting therefore contributes a score/risk penalty but
    # cannot become an automatic hard rejection unless comparable evidence is
    # at least MEDIUM confidence (or experimental).
    if constraints.get("herg_do_not_increase"):
        parent = ((optimization_run.evidence_json or {}).get("admet", {}).get("hERG liability") or {}).get("preferred")
        result = (candidate.admet_json or {}).get("hERG liability", {})
        parent_value, candidate_value = (parent or {}).get("value"), _preferred_value(result)
        if all(value is not None for value in (parent_value, candidate_value)) and _confidence(result) in {"EXPERIMENTAL", "HIGH", "MEDIUM"} and candidate_value > float(parent_value):
            _reject(db, candidate, "HERG_LIABILITY_INCREASE", "Comparable hERG liability increased under sufficient-confidence evidence.", "PREDICTING", result.get("record_type", "Predicted"))
            return False
    return True


def _parent_activity(optimization_run):
    activity = (optimization_run.evidence_json or {}).get("activity", {})
    if activity.get("experimental"):
        return {"record_type": "Experimental", "value_nm": activity["experimental"]["mean_nm"], "confidence": "EXPERIMENTAL"}
    if activity.get("predicted"):
        return {"record_type": "Predicted", "value_nm": activity["predicted"]["value_nm"], "confidence": activity["predicted"]["confidence"]}
    return None


def _soft_spot_change(parent_spots, candidate_result):
    candidate_spots = candidate_result.get("spots", [])
    parent_primary = parent_spots[0] if parent_spots else None
    candidate_primary = candidate_spots[0] if candidate_spots else None
    parent_type = parent_primary.get("transformation") if parent_primary else None
    candidate_top3 = [row["transformation"] for row in candidate_spots[:3]]
    return {
        "parent_primary": parent_primary, "candidate_primary": candidate_primary,
        "parent_primary_absent_from_candidate_top3": bool(parent_type and parent_type not in candidate_top3),
        "new_primary_liability": bool(candidate_primary and candidate_primary.get("transformation") != parent_type),
        "interpretation": "Rule-prior rank comparison; atom indices are molecule-local and are not directly equated.",
    }


def _candidate_objectives(candidate, optimization_run):
    activity = candidate.activity_json or {}
    activity_quality = None
    if activity.get("status") == "COMPLETE" and activity.get("pactivity") is not None:
        activity_quality = max(0.0, min(1.0, (float(activity["pactivity"]) - 4.0) / 6.0))
    admet = candidate.admet_json or {}
    values = {
        "Activity": activity_quality,
        "HLM intrinsic clearance": _quality("HLM intrinsic clearance", admet.get("HLM intrinsic clearance")),
        "RLM intrinsic clearance": _quality("RLM intrinsic clearance", admet.get("RLM intrinsic clearance")),
        "Solubility": _quality("Solubility", admet.get("Solubility")),
        "Permeability": _quality("Permeability", admet.get("Permeability")),
        "hERG liability": _quality("hERG liability", admet.get("hERG liability")),
        "P-gp inhibitor": _quality("P-gp inhibitor", admet.get("P-gp inhibitor")),
        "Structural alerts": 1.0 if not candidate.stage1_json.get("alerts") else max(0.0, 1.0 - 0.2 * len(candidate.stage1_json["alerts"])),
        "Synthetic complexity": max(0.0, min(1.0, (10.0 - candidate.synthetic_feasibility_json.get("sa_score", 10.0)) / 9.0)),
        "Similarity": candidate.parent_similarity,
    }
    cyp = [_quality(endpoint, result) for endpoint, result in admet.items() if endpoint.startswith("CYP") and endpoint.endswith("inhibitor")]
    values["CYP"] = sum(value for value in cyp if value is not None) / len([value for value in cyp if value is not None]) if any(value is not None for value in cyp) else None
    factors = {}
    factors["Activity"] = CONFIDENCE_FACTOR.get(_confidence(activity), 0.35) * DOMAIN_FACTOR.get(_domain(activity), 0.5)
    for endpoint in ("HLM intrinsic clearance", "RLM intrinsic clearance", "Solubility", "Permeability", "hERG liability", "P-gp inhibitor"):
        result = admet.get(endpoint, {})
        factors[endpoint] = CONFIDENCE_FACTOR.get(_confidence(result), 0.35) * DOMAIN_FACTOR.get(_domain(result), 0.5)
    cyp_factors = [CONFIDENCE_FACTOR.get(_confidence(result), 0.35) * DOMAIN_FACTOR.get(_domain(result), 0.5) for endpoint, result in admet.items() if endpoint.startswith("CYP") and endpoint.endswith("inhibitor") and result.get("scored")]
    factors["CYP"] = sum(cyp_factors) / len(cyp_factors) if cyp_factors else 0.0
    factors.update({"Structural alerts": 0.8, "Synthetic complexity": 0.8, "Similarity": 1.0})
    return values, factors


def pareto_fronts(candidates):
    remaining = list(candidates)
    front = 1
    while remaining:
        current = []
        for candidate in remaining:
            vector = candidate.objective_vector_json.get("values", {})
            dominated = False
            for other in remaining:
                if other.id == candidate.id:
                    continue
                other_vector = other.objective_vector_json.get("values", {})
                keys = [key for key, value in vector.items() if value is not None and other_vector.get(key) is not None]
                if keys and all(other_vector[key] >= vector[key] for key in keys) and any(other_vector[key] > vector[key] for key in keys):
                    dominated = True; break
            if not dominated:
                current.append(candidate)
        for candidate in current:
            candidate.pareto_front = front
        remaining = [candidate for candidate in remaining if candidate not in current]
        front += 1


def rank_candidates(db, proposal_run, optimization_run):
    candidates = [row for row in proposal_run.candidates if row.status not in {"REJECTED", "FAILED"} and row.user_decision != "REJECTED"]
    weights = dict(DEFAULT_WEIGHTS)
    for objective in optimization_run.objectives_json or []:
        for endpoint, boost in OBJECTIVE_WEIGHT_BOOSTS.get(objective, {}).items():
            weights[endpoint] = weights.get(endpoint, 1.0) * boost
    weights.update(proposal_run.endpoint_weights_json or optimization_run.endpoint_weights_json or {})
    for candidate in candidates:
        values, factors = _candidate_objectives(candidate, optimization_run)
        usable = {key: value for key, value in values.items() if value is not None and factors.get(key, 0) > 0 and float(weights.get(key, 1.0)) > 0}
        denominator = sum(float(weights.get(key, 1.0)) for key in usable) or 1.0
        objective_term = sum(float(weights.get(key, 1.0)) * value * factors[key] for key, value in usable.items()) / denominator
        structural_penalty = 0.04 * len(candidate.stage1_json.get("alerts", []))
        synthetic_penalty = 0.08 if candidate.synthetic_feasibility_json.get("classification") == "HIGH SYNTHETIC COMPLEXITY" else (0.03 if candidate.synthetic_feasibility_json.get("classification") == "MODERATE SYNTHETIC COMPLEXITY" else 0.0)
        uncertainty_penalty = 0.12 if candidate.activity_json.get("applicability_domain") in {"OUT OF DOMAIN", "OUT_OF_DOMAIN"} else 0.0
        score = 100.0 * max(0.0, objective_term - structural_penalty - synthetic_penalty - uncertainty_penalty)
        candidate.objective_vector_json = {"values": values, "evidence_domain_factors": factors, "weights": weights}
        candidate.ranking_score = round(score, 4)
        candidate.confidence = "HIGH" if min([factor for factor in factors.values() if factor > 0] or [0]) >= 0.75 else ("MEDIUM" if np.mean([factor for factor in factors.values() if factor > 0] or [0]) >= 0.6 else "LOW")
        candidate.applicability_domain = "OUT_OF_DOMAIN" if any(_domain(row) in {"OUT OF DOMAIN", "OUT_OF_DOMAIN"} for row in [candidate.activity_json, *candidate.admet_json.values()] if isinstance(row, dict) and row.get("status") == "COMPLETE") else ("BORDERLINE" if any(_domain(row) == "BORDERLINE" for row in [candidate.activity_json, *candidate.admet_json.values()] if isinstance(row, dict)) else "IN_DOMAIN")
        risks = []
        if candidate.applicability_domain == "OUT_OF_DOMAIN": risks.append("One or more predictions are OUT_OF_DOMAIN")
        if candidate.synthetic_feasibility_json.get("classification") == "HIGH SYNTHETIC COMPLEXITY": risks.append("High synthetic-complexity surrogate")
        for endpoint in ("hERG liability", "Ames mutagenicity", "DILI clinical liability"):
            row = candidate.admet_json.get(endpoint, {})
            if row.get("classification") == MODEL_SPECS.get(endpoint, {}).get("positive_label"):
                risks.append(f"Potential {endpoint} ({row.get('confidence', 'UNKNOWN')} confidence)")
        if candidate.soft_spot_change_json.get("new_primary_liability"):
            risks.append("New primary metabolic soft-spot hypothesis")
        candidate.main_risk = "; ".join(risks) or "No hard constraint violation; prediction uncertainty remains"
    pareto_fronts(candidates)
    ordered = sorted(candidates, key=lambda row: (row.pareto_front or 999, -(row.ranking_score or 0), -row.parent_similarity, row.id))
    selected, used_hypotheses, selected_fps = [], set(), []
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for candidate in ordered:
        transformation_ids = tuple(row.transformation_id for row in candidate.transformations)
        fp = generator.GetFingerprint(Chem.MolFromSmiles(candidate.canonical_smiles))
        max_selected_similarity = max((float(DataStructs.TanimotoSimilarity(fp, other)) for other in selected_fps), default=0.0)
        unique_hypothesis = transformation_ids not in used_hypotheses
        candidate.information_value = "HIGH" if unique_hypothesis and (candidate.pareto_front == 1 or max_selected_similarity < 0.75) else ("MEDIUM" if unique_hypothesis or max_selected_similarity < 0.85 else "LOW")
        if len(selected) < 10 and (unique_hypothesis or max_selected_similarity < 0.82 or len(ordered) <= 10):
            selected.append(candidate); selected_fps.append(fp); used_hypotheses.add(transformation_ids)
    if len(selected) < min(10, len(ordered)):
        for candidate in ordered:
            if candidate not in selected:
                selected.append(candidate)
                if len(selected) == min(10, len(ordered)): break
    for rank, candidate in enumerate(ordered, 1):
        candidate.selected_top10 = candidate in selected or candidate.user_decision == "PROMOTED"
        candidate.status = "TOP_10" if candidate.selected_top10 else "ACCEPTED"
        db.add(CandidateRanking(
            candidate_id=candidate.id, rank=rank, score=candidate.ranking_score or 0,
            pareto_front=candidate.pareto_front or 999, score_formula_version=ENGINE_VERSION,
            score_breakdown_json={
                "formula": "100 × max(0, Σ(weight × normalized objective × evidence-confidence × domain-factor) / Σ(usable weights) − structural risk − synthetic complexity − OOD penalty)",
                "objective_vector": candidate.objective_vector_json,
            },
            diversity_json={"selection": "greedy transformation-hypothesis diversity then Morgan/Tanimoto diversity", "information_value": candidate.information_value},
            selected_top10=candidate.selected_top10,
        ))
    return ordered, selected


def _new_candidate(db, proposal_run, optimization_run, number, product, strategy, changed, sequence=1, user_added=False, existing_version_id=None):
    canonical = Chem.MolToSmiles(product, canonical=True, isomericSmiles=False)
    candidate = OptimizationCandidate(
        proposal_run_id=proposal_run.id, project_id=proposal_run.project_id,
        optimization_run_id=optimization_run.id, parent_version_id=proposal_run.parent_version_id,
        existing_version_id=existing_version_id, candidate_number=number,
        canonical_smiles=canonical, isomeric_smiles=Chem.MolToSmiles(product, canonical=True, isomericSmiles=True),
        inchikey=Chem.MolToInchiKey(product), generation_priority=1 if strategy["id"].startswith("MMP_PROJECT") else (3 if strategy.get("target_liability", "").startswith("LIAB_SOFT") else 4),
        generation_source="Project MMP-derived transformation" if strategy["id"].startswith("MMP_PROJECT") else "Stage 4A ranked curated transformation",
        hypothesis=strategy["name"], why_generated=" · ".join(strategy.get("evidence") or [strategy["purpose"]]),
        expected_benefit=strategy.get("expected_effect", ""), status="GENERATED", user_added=user_added,
    )
    db.add(candidate); db.flush()
    db.add(CandidateTransformation(
        candidate=candidate, sequence_number=sequence, transformation_id=strategy["id"], name=strategy["name"],
        reaction_smarts=strategy.get("reaction_smarts", ""), transformation_version=strategy.get("version", ""),
        source=strategy.get("source", ""), source_atom_indices_json=list(changed), changed_parent_atoms_json=list(changed),
        execution_status="USER_DEFINED" if user_added else "EXECUTED",
        provenance_json={"strategy_rank": strategy.get("rank"), "target_liability": strategy.get("target_liability"), "application_status": strategy.get("application_status"), "engine": ENGINE_NAME, "engine_version": ENGINE_VERSION},
    ))
    return candidate


def _generate_candidates(db, proposal_run, optimization_run):
    parent = db.get(CompoundVersion, proposal_run.parent_version_id)
    parent_mol = Chem.MolFromSmiles(parent.isomeric_smiles)
    protected = {int(atom) for row in (optimization_run.protected_regions_json or []) for atom in row.get("atom_indices", [])}
    allowed = {int(atom) for row in (optimization_run.modifiable_regions_json or []) for atom in row.get("atom_indices", [])}
    max_raw = int((proposal_run.settings_json or {}).get("max_raw_candidates", 120))
    strategies = optimization_run.transformations_json or []
    generated, number = [], 1
    for strategy in strategies:
        if strategy["id"] in STRATEGY_ONLY_TRANSFORMATIONS:
            continue
        # A ranked Stage 4A source site is itself an explicit strategy decision;
        # combine it with separately stored modifiable/soft-spot regions.  Other
        # motif matches remain excluded, so the molecule is never scanned and
        # edited indiscriminately.
        scoped_atoms = allowed.union(int(atom) for atom in strategy.get("source_atom_indices", []))
        for product, changed in execute_strategy(parent.isomeric_smiles, strategy, scoped_atoms):
            if product is None or set(changed).intersection(protected):
                continue
            existing_version = None
            if strategy["id"].startswith("MMP_PROJECT_OBSERVED_"):
                try:
                    pair_id = int(strategy["id"].rsplit("_", 1)[1])
                    pair = next(row for row in (optimization_run.evidence_json or {}).get("activity", {}).get("mmp", []) if row["pair_id"] == pair_id)
                    existing_version = pair["other_version_id"]
                except (ValueError, StopIteration, KeyError):
                    existing_version = None
            generated.append(_new_candidate(db, proposal_run, optimization_run, number, product, strategy, changed, existing_version_id=existing_version)); number += 1
            if len(generated) >= max_raw: break
        if len(generated) >= max_raw: break
    if (proposal_run.settings_json or {}).get("allow_double_transforms", True) and len(generated) >= 2:
        seed_rows = generated[:min(10, len(generated))]
        executable = [row for row in strategies if row["id"] in EXECUTABLE_TRANSFORMATIONS]
        for seed in seed_rows:
            first_id = seed.transformations[0].transformation_id
            first_liability = (seed.transformations[0].provenance_json or {}).get("target_liability")
            for strategy in executable:
                if strategy["id"] == first_id or not first_liability or strategy.get("target_liability") == first_liability:
                    continue
                for product, changed in execute_strategy(seed.isomeric_smiles, strategy, None)[:2]:
                    validation = chemical_validation(parent.isomeric_smiles, Chem.MolToSmiles(product, isomericSmiles=True), protected)
                    if not validation.get("valid") or validation["similarity"] < float((optimization_run.constraints_json or {}).get("similarity_min", 0.0)):
                        continue
                    candidate = _new_candidate(db, proposal_run, optimization_run, number, product, strategy, changed, sequence=2)
                    first = seed.transformations[0]
                    db.add(CandidateTransformation(
                        candidate=candidate, sequence_number=1, transformation_id=first.transformation_id,
                        name=first.name, reaction_smarts=first.reaction_smarts,
                        transformation_version=first.transformation_version, source=first.source,
                        source_atom_indices_json=first.source_atom_indices_json,
                        changed_parent_atoms_json=first.changed_parent_atoms_json,
                        execution_status="EXECUTED", provenance_json=first.provenance_json,
                    ))
                    candidate.hypothesis = first.name + " + " + strategy["name"]
                    candidate.why_generated = "Two independent Stage 4A strategies; retained only after parent-level protection/similarity validation."
                    generated.append(candidate); number += 1
                    if len(generated) >= max_raw: break
                if len(generated) >= max_raw: break
            if len(generated) >= max_raw: break
    return generated


def _set_status(db, proposal_run, status, message):
    proposal_run.status, proposal_run.stage_message = status, message
    db.commit(); db.refresh(proposal_run)


def execute_proposal_run(run_id, session=None):
    """Execute a persisted job; individual candidate failures never abort the run."""
    own_session = session is None
    db = session or SessionLocal()
    try:
        proposal_run = db.get(OptimizationProposalRun, run_id)
        if not proposal_run:
            return
        optimization_run = db.get(OptimizationRun, proposal_run.optimization_run_id)
        parent_version = db.get(CompoundVersion, proposal_run.parent_version_id)
        proposal_run.started_at = datetime.now(timezone.utc)
        _set_status(db, proposal_run, "GENERATING", "Applying executable Stage 4A transformations deterministically.")
        candidates = _generate_candidates(db, proposal_run, optimization_run)
        proposal_run.raw_candidate_count = len(candidates)
        db.commit()
        _set_status(db, proposal_run, "FILTERING", "Running sanitization, structural preservation, Stage 1, similarity, and cheap hard gates.")
        parent_analysis = analyze_smiles(parent_version.isomeric_smiles)
        protected = {int(atom) for row in (optimization_run.protected_regions_json or []) for atom in row.get("atom_indices", [])}
        constraints = {"no_new_structural_alert": True, **(optimization_run.constraints_json or {}), **(proposal_run.hard_constraints_json or {})}
        seen, survivors = {}, []
        for candidate in candidates:
            try:
                if _cheap_constraints(db, candidate, parent_analysis, constraints, protected, seen):
                    survivors.append(candidate)
            except Exception as exc:
                candidate.status = "FAILED"; candidate.main_risk = f"Candidate-local Stage 1 failure: {exc}"
        db.commit()
        _set_status(db, proposal_run, "PREDICTING", "Running project activity, available Stage 3 endpoints, and soft-spot reanalysis for filtered candidates.")
        parent_activity = _parent_activity(optimization_run)
        rescored = []
        parent_spots = (optimization_run.evidence_json or {}).get("metabolism", {}).get("soft_spots", [])
        for candidate in survivors:
            try:
                if not _activity_gate(db, candidate, optimization_run, parent_activity, constraints):
                    continue
                candidate.admet_json = rescore_admet(db, candidate, proposal_run.project_id)
                if not _post_prediction_constraints(db, candidate, optimization_run, constraints):
                    db.commit(); continue
                candidate.soft_spot_json = predict_soft_spots(candidate.canonical_smiles, max_spots=12)
                candidate.soft_spot_change_json = _soft_spot_change(parent_spots, candidate.soft_spot_json)
                candidate.status = "RESCORED"; rescored.append(candidate)
                db.commit()
            except Exception as exc:
                candidate.status = "FAILED"; candidate.main_risk = f"Candidate-local prediction failure: {exc}"; db.commit()
        _set_status(db, proposal_run, "RANKING", "Computing confidence/domain-aware objectives, Pareto fronts, diversity, and information value.")
        ordered, selected = rank_candidates(db, proposal_run, optimization_run)
        proposal_run.accepted_count = len(ordered)
        proposal_run.rejected_count = len([row for row in proposal_run.candidates if row.status in {"REJECTED", "FAILED"}])
        proposal_run.top_count = len(selected)
        proposal_run.model_versions_json = {
            "Stage 1": f"RDKit {RDKIT_VERSION}", "Stage 2": "latest project assay model or similarity",
            "Stage 3": {endpoint: spec["model_version"] for endpoint, spec in MODEL_SPECS.items()},
            "Soft spot": SOFT_SPOT_VERSION,
        }
        proposal_run.summary_json = {
            "engine": ENGINE_NAME, "engine_version": ENGINE_VERSION,
            "executable_transformations": sorted(EXECUTABLE_TRANSFORMATIONS),
            "strategy_only_transformations": STRATEGY_ONLY_TRANSFORMATIONS,
            "generation_policy": "Ranked project MMP → SAR/soft-spot supported → curated transformations; single change first, at most two changes",
            "score_formula": "100 × max(0, weighted confidence/domain-adjusted objective mean − structural − synthetic-complexity − OOD penalties)",
            "pareto_method": "Non-dominated sorting over available normalized objectives; unavailable endpoints omitted",
            "diversity_method": "Greedy unique transformation hypothesis followed by Morgan/Tanimoto diversity",
            "llm_used": False, "pk_run": False,
        }
        proposal_run.status = "COMPLETED" if ordered else "COMPLETED"
        proposal_run.stage_message = f"Generated {len(candidates)} raw candidates; {len(ordered)} accepted, {proposal_run.rejected_count} rejected, {len(selected)} selected."
        if not ordered:
            proposal_run.summary_json["no_valid_analog"] = True
        proposal_run.completed_at = datetime.now(timezone.utc)
        db.commit()
    except Exception as exc:
        db.rollback()
        proposal_run = db.get(OptimizationProposalRun, run_id)
        if proposal_run:
            proposal_run.status = "FAILED"; proposal_run.stage_message = str(exc); proposal_run.completed_at = datetime.now(timezone.utc); db.commit()
    finally:
        if own_session:
            db.close()


def process_user_candidate(db, proposal_run, optimization_run, smiles, reason="User-added analog"):
    product = Chem.MolFromSmiles(smiles)
    if product is None:
        raise ValueError("User analog failed RDKit sanitization")
    strategy = {"id": "USER_DEFINED", "name": "User-added analog", "reaction_smarts": "", "version": ENGINE_VERSION, "source": "User", "purpose": "Manual medicinal chemistry proposal", "expected_effect": reason, "evidence": [reason]}
    number = max([row.candidate_number for row in proposal_run.candidates] or [0]) + 1
    candidate = _new_candidate(db, proposal_run, optimization_run, number, product, strategy, [], user_added=True)
    db.commit()
    parent_version = db.get(CompoundVersion, proposal_run.parent_version_id)
    parent_analysis = analyze_smiles(parent_version.isomeric_smiles)
    protected = {int(atom) for row in (optimization_run.protected_regions_json or []) for atom in row.get("atom_indices", [])}
    constraints = {"no_new_structural_alert": True, **(optimization_run.constraints_json or {}), **(proposal_run.hard_constraints_json or {})}
    seen = {row.canonical_smiles: row for row in proposal_run.candidates if row.id != candidate.id and row.status not in {"REJECTED", "FAILED"}}
    if _cheap_constraints(db, candidate, parent_analysis, constraints, protected, seen) and _activity_gate(db, candidate, optimization_run, _parent_activity(optimization_run), constraints):
        candidate.admet_json = rescore_admet(db, candidate, proposal_run.project_id)
        if _post_prediction_constraints(db, candidate, optimization_run, constraints):
            candidate.soft_spot_json = predict_soft_spots(candidate.canonical_smiles, max_spots=12)
            candidate.soft_spot_change_json = _soft_spot_change((optimization_run.evidence_json or {}).get("metabolism", {}).get("soft_spots", []), candidate.soft_spot_json)
            candidate.status = "RESCORED"
    rank_candidates(db, proposal_run, optimization_run)
    proposal_run.accepted_count = len([row for row in proposal_run.candidates if row.status not in {"REJECTED", "FAILED"}])
    proposal_run.rejected_count = len([row for row in proposal_run.candidates if row.status in {"REJECTED", "FAILED"}])
    proposal_run.top_count = len([row for row in proposal_run.candidates if row.selected_top10])
    db.commit(); db.refresh(candidate)
    return candidate
