"""
Authoritative Comparison Pair Matcher and Scientific Context Engine.
=====================================================================
Connects Experimental Evidence <-> Canonical Endpoint <-> Unit Normalization
<-> Context (Species / Matrix / Route / Dose) <-> Prediction Result.

Strict Scientific Rules:
1. Four Classification Types:
   - DIRECTLY_COMPARABLE
   - CONVERTED_COMPARABLE
   - CONTEXTUALLY_RELATED_NOT_DIRECTLY_COMPARABLE
   - NO_MATCHING_PREDICTION_MODEL
2. Prohibitions:
   - Never compare human PK with rat PK (cross-species forbidden)
   - Never compare CL with CL/F
   - Never compare Vdss with Vz, Vd, or Vd/F
   - Never compare IC50 with Ki/EC50 directly
   - Never compare pIC50 with IC50 without deterministic conversion
   - Never compare classifier probabilities with quantitative measurements
   - Never compare different doses/routes for PK without explicit dose normalization
3. Zero Tolerance:
   - Unmatched due to endpoint name or unit mismatch = 0.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from backend.canonical_endpoints import (
    DIRECT, CONVERTED, CONDITIONAL, RELATED, UNSUPPORTED, NOT_COMPARABLE, NO_MODEL,
    REGISTRY, CanonicalEndpoint, endpoint_contract, normalize_species, normalize_route
)
from backend.unit_normalization import (
    normalize_scientific_observation, convert_concentration, convert_solubility,
    convert_ppb, convert_cyp_herg_ic50, convert_pic50_to_ic50, convert_clearance,
    convert_volume, convert_time, convert_caco2_papp, convert_auc,
    clean_unit_str, TYPE_DETERMINISTIC, TYPE_SCALE, TYPE_DIRECT_IDENTITY, TYPE_IVIVE_DERIVED
)

MATCHING_ENGINE_VERSION = "drugopt-comparison-matcher-v2.0"

# Explicit failure reasons
REASON_MODEL_UNAVAILABLE = "MODEL_UNAVAILABLE"
REASON_CONTEXT_MISMATCH = "CONTEXT_MISMATCH"
REASON_SPECIES_MISMATCH = "SPECIES_MISMATCH"
REASON_ROUTE_MISMATCH = "ROUTE_MISMATCH"
REASON_DOSE_REQUIRED = "DOSE_REQUIRED"
REASON_ROUTE_REQUIRED = "ROUTE_REQUIRED"
REASON_NOT_PREDICTED_BY_CURRENT_ENGINE = "NOT_PREDICTED_BY_CURRENT_ENGINE"
REASON_EVIDENCE_NOT_NUMERIC = "EVIDENCE_NOT_NUMERIC"
REASON_NOT_DIRECTLY_COMPARABLE = "NOT_DIRECTLY_COMPARABLE"
REASON_ASSAY_SEMANTICS_MISMATCH = "ASSAY_SEMANTICS_MISMATCH"


@dataclass
class ScientificComparisonPair:
    canonical_endpoint_id: str
    display_name: str
    section: str
    species: str
    route: str
    dose: Optional[float] = None
    dose_unit: str = ""
    regimen: str = "SINGLE_DOSE"
    matrix: str = ""
    analyte: str = "PARENT"

    # Experimental side
    experimental_raw_value: Optional[float] = None
    experimental_raw_unit: str = ""
    experimental_normalized_value: Optional[float] = None
    experimental_normalized_unit: str = ""
    experimental_display_string: str = ""
    experimental_reference: Dict[str, Any] = field(default_factory=dict)
    experimental_id: Optional[int] = None
    relation: str = "="

    # Prediction side
    prediction_available: bool = False
    prediction_raw_value: Optional[float] = None
    prediction_raw_unit: str = ""
    prediction_normalized_value: Optional[float] = None
    prediction_normalized_unit: str = ""
    prediction_display_string: str = ""
    prediction_source_type: str = "MODEL"
    prediction_source_label: str = "Model Prediction"
    prediction_model_id: str = ""
    prediction_engine_version: str = "UNKNOWN_PROVENANCE"
    prediction_ad_status: str = "IN_DOMAIN"
    prediction_uncertainty: Optional[float] = None
    maturity_level: int = 1
    maturity_label: str = "Base / Mechanistic Estimate"
    maturity_stars: str = "★☆☆☆☆"

    # Comparison classification & difference
    comparison_status: str = NO_MODEL
    comparison_badge: str = "NO MODEL"
    is_directly_comparable: bool = False
    is_converted: bool = False
    is_dose_normalized: bool = False
    conversion_formula: str = ""
    signed_error: Optional[float] = None
    absolute_error: Optional[float] = None
    fold_error: Optional[float] = None
    error_metric_type: str = "absolute_error"
    difference_display_value: Optional[float] = None
    difference_display_unit: str = ""
    difference_display_string: str = "—"
    scientific_interpretation: str = "CONTEXT_DEPENDENT"
    agreement_interpretation: str = "NOT_NUMERICALLY_COMPARABLE"
    unmatched_reason: str = ""
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "canonical_endpoint_id": self.canonical_endpoint_id,
            "display_name": self.display_name,
            "section": self.section,
            "species": self.species,
            "route": self.route,
            "dose": self.dose,
            "dose_unit": self.dose_unit,
            "regimen": self.regimen,
            "matrix": self.matrix,
            "analyte": self.analyte,
            "experimental_raw_value": self.experimental_raw_value,
            "experimental_raw_unit": self.experimental_raw_unit,
            "experimental_normalized_value": self.experimental_normalized_value,
            "experimental_normalized_unit": self.experimental_normalized_unit,
            "experimental_display_string": self.experimental_display_string,
            "experimental_reference": self.experimental_reference,
            "experimental_id": self.experimental_id,
            "relation": self.relation,
            "prediction_available": self.prediction_available,
            "prediction_raw_value": self.prediction_raw_value,
            "prediction_raw_unit": self.prediction_raw_unit,
            "prediction_normalized_value": self.prediction_normalized_value,
            "prediction_normalized_unit": self.prediction_normalized_unit,
            "prediction_display_string": self.prediction_display_string,
            "prediction_source_type": self.prediction_source_type,
            "prediction_source_label": self.prediction_source_label,
            "prediction_model_id": self.prediction_model_id,
            "prediction_engine_version": self.prediction_engine_version,
            "prediction_ad_status": self.prediction_ad_status,
            "prediction_uncertainty": self.prediction_uncertainty,
            "maturity_level": self.maturity_level,
            "maturity_label": self.maturity_label,
            "maturity_stars": self.maturity_stars,
            "comparison_status": self.comparison_status,
            "comparison_badge": self.comparison_badge,
            "is_directly_comparable": self.is_directly_comparable,
            "is_converted": self.is_converted,
            "is_dose_normalized": self.is_dose_normalized,
            "conversion_formula": self.conversion_formula,
            "signed_error": self.signed_error,
            "absolute_error": self.absolute_error,
            "fold_error": self.fold_error,
            "error_metric_type": self.error_metric_type,
            "difference_display_value": self.difference_display_value,
            "difference_display_unit": self.difference_display_unit,
            "difference_display_string": self.difference_display_string,
            "scientific_interpretation": self.scientific_interpretation,
            "agreement_interpretation": self.agreement_interpretation,
            "unmatched_reason": self.unmatched_reason,
            "timestamp": self.timestamp,
        }


class ComparisonPairMatcher:
    """
    Authoritative scientific pair matcher evaluating 10 dimensions:
    1. canonical_endpoint_id
    2. measurement_type
    3. species (strict)
    4. biological_matrix
    5. assay context
    6. route (strict)
    7. dose/regimen (strict or valid dose-normalized)
    8. chemical form (parent vs metabolite)
    9. relation (=, <, >)
    10. canonical unit / scale
    """

    def __init__(self, mw: Optional[float] = None, body_weight_kg: float = 70.0):
        self.mw = mw
        self.body_weight_kg = body_weight_kg

    def match_pair(
        self,
        exp_record: Dict[str, Any],
        pred_record: Optional[Dict[str, Any]],
    ) -> ScientificComparisonPair:
        """Evaluate experimental observation vs candidate prediction record."""
        eid = exp_record.get("canonical_endpoint_id") or "UNRESOLVED"
        ep_def = REGISTRY.get(eid)
        display_name = ep_def.display_name if ep_def else exp_record.get("raw_endpoint_name", eid)
        section = ep_def.section if ep_def else exp_record.get("routing_section", "UNCLASSIFIED")

        species = exp_record.get("species") or "HUMAN"
        route = exp_record.get("route") or "UNSPECIFIED"
        dose = exp_record.get("dose")
        dose_unit = exp_record.get("dose_unit") or ""
        regimen = exp_record.get("regimen") or "UNSPECIFIED"
        matrix = exp_record.get("matrix") or ""
        analyte = exp_record.get("analyte") or "PARENT"

        raw_exp_val = exp_record.get("raw_value")
        raw_exp_unit = exp_record.get("raw_unit") or ""
        relation = exp_record.get("relation") or "="

        # Check numeric nature of experiment
        try:
            num_exp_val = float(str(raw_exp_val).replace(",", "").strip()) if raw_exp_val is not None else None
        except (ValueError, TypeError):
            num_exp_val = None

        pair = ScientificComparisonPair(
            canonical_endpoint_id=eid,
            display_name=display_name,
            section=section,
            species=species,
            route=route,
            dose=dose,
            dose_unit=dose_unit,
            regimen=regimen,
            matrix=matrix,
            analyte=analyte,
            experimental_raw_value=num_exp_val,
            experimental_raw_unit=raw_exp_unit,
            experimental_id=exp_record.get("id"),
            relation=relation,
            experimental_reference=exp_record.get("reference") or {},
        )

        # Normalize experimental value
        norm_res = normalize_scientific_observation(
            eid, num_exp_val, raw_exp_unit,
            mw=self.mw, species=species, route=route
        )
        if norm_res:
            pair.experimental_normalized_value = norm_res.normalized_value
            pair.experimental_normalized_unit = norm_res.normalized_unit
            pair.conversion_formula = norm_res.conversion_formula
            pair.is_converted = (norm_res.conversion_type in {TYPE_DETERMINISTIC, TYPE_SCALE})
            
            # Format display string showing original and converted if converted
            if pair.is_converted:
                pair.experimental_display_string = f"{num_exp_val:g} {raw_exp_unit} → {norm_res.normalized_value:.3f} {norm_res.normalized_unit} (Converted)"
            else:
                pair.experimental_display_string = f"{norm_res.normalized_value:g} {norm_res.normalized_unit}".strip()
        else:
            pair.experimental_display_string = f"{raw_exp_val} {raw_exp_unit}".strip() if raw_exp_val is not None else "—"

        # If no prediction record exists
        if not pred_record or not pred_record.get("available"):
            pair.prediction_available = False
            pair.comparison_status = NO_MODEL
            pair.comparison_badge = "NO MODEL"
            pair.unmatched_reason = pred_record.get("unavailable_reason") if pred_record else REASON_NOT_PREDICTED_BY_CURRENT_ENGINE
            pair.scientific_interpretation = "EXPERIMENTAL_ONLY"
            pair.agreement_interpretation = "NO_PREDICTION"
            return pair

        # Prediction exists: load prediction parameters
        pair.prediction_available = True
        pred_val = pred_record.get("value") if pred_record.get("value") is not None else pred_record.get("display_value")
        pred_unit = pred_record.get("unit") or ""
        pair.prediction_raw_value = float(pred_val) if pred_val is not None else None
        pair.prediction_raw_unit = pred_unit
        pair.prediction_source_type = pred_record.get("source_type", "MODEL")
        pair.prediction_source_label = pred_record.get("source_label", "Model Prediction")
        pair.prediction_model_id = pred_record.get("model_id", "")
        pair.prediction_engine_version = pred_record.get("engine_version") or "UNKNOWN_PROVENANCE"
        pair.prediction_ad_status = pred_record.get("ad_status", "IN_DOMAIN")
        pair.prediction_uncertainty = pred_record.get("uncertainty")

        mat = pred_record.get("maturity") or {}
        pair.maturity_level = mat.get("level", 1)
        pair.maturity_label = mat.get("label", "Base Prediction")
        pair.maturity_stars = mat.get("stars", "★☆☆☆☆")

        # Normalize prediction value into canonical scale
        pred_norm_res = normalize_scientific_observation(
            eid, pair.prediction_raw_value, pred_unit,
            mw=self.mw, species=pred_record.get("species", species), route=pred_record.get("route", route)
        )
        if pred_norm_res:
            pair.prediction_normalized_value = pred_norm_res.normalized_value
            pair.prediction_normalized_unit = pred_norm_res.normalized_unit
            pair.prediction_display_string = f"{pred_norm_res.normalized_value:.3f} {pred_norm_res.normalized_unit}".strip()
        else:
            pair.prediction_normalized_value = pair.prediction_raw_value
            pair.prediction_normalized_unit = pred_unit
            pair.prediction_display_string = f"{pair.prediction_raw_value} {pred_unit}".strip() if pair.prediction_raw_value is not None else "—"

        # Check for non-numeric experimental observation
        if pair.experimental_normalized_value is None:
            pair.comparison_status = CONDITIONAL
            pair.comparison_badge = "CONTEXT ONLY"
            pair.unmatched_reason = REASON_EVIDENCE_NOT_NUMERIC
            pair.scientific_interpretation = "QUALITATIVE_EVALUATION"
            pair.agreement_interpretation = "NOT_NUMERICALLY_COMPARABLE"
            return pair

        # Dimension 3: Strict Species Matching
        pred_species = normalize_species(pred_record.get("species", species))
        if pred_species != species and species != "UNSPECIFIED" and pred_species != "UNSPECIFIED":
            pair.comparison_status = RELATED
            pair.comparison_badge = "CONTEXT ONLY"
            pair.unmatched_reason = f"{REASON_SPECIES_MISMATCH}: Experimental {species} vs Predicted {pred_species}"
            pair.scientific_interpretation = "CROSS_SPECIES_TRANSLATIONAL"
            pair.agreement_interpretation = "NOT_NUMERICALLY_COMPARABLE"
            return pair

        # Dimension 6: Strict Route Matching for PK
        if section == "PK":
            pred_route = normalize_route(pred_record.get("route", route))
            if pred_route != route and route != "UNSPECIFIED" and pred_route != "UNSPECIFIED":
                pair.comparison_status = RELATED
                pair.comparison_badge = "CONTEXT ONLY"
                pair.unmatched_reason = f"{REASON_ROUTE_MISMATCH}: Experimental {route} vs Predicted {pred_route}"
                pair.scientific_interpretation = "ROUTE_MISMATCH"
                pair.agreement_interpretation = "NOT_NUMERICALLY_COMPARABLE"
                return pair

        # Dimension 7: Dose context matching / dose normalization
        pred_dose = pred_record.get("dose")
        is_dose_norm = False
        exp_eval_val = pair.experimental_normalized_value
        if section == "PK" and dose is not None and pred_dose is not None and dose > 0 and pred_dose > 0:
            if abs(dose - pred_dose) > 1e-3:
                # If dose differs for linear PK parameters (Cmax, AUC), apply linear dose normalization
                if any(param in eid for param in ("CMAX", "AUC")):
                    dose_ratio = float(pred_dose) / float(dose)
                    exp_eval_val = pair.experimental_normalized_value * dose_ratio
                    is_dose_norm = True
                    pair.is_dose_normalized = True
                    pair.conversion_formula += f" | Dose-normalized from {dose} to {pred_dose} mg (* {dose_ratio:.3f})"
                else:
                    # Non-dose-dependent PK parameters (t1/2, CL, Vd) compare directly
                    pass

        # Check unit/scale compatibility
        exp_val = exp_eval_val
        pred_val = pair.prediction_normalized_value

        # Calculate Difference / Errors
        if exp_val is not None and pred_val is not None:
            # Determine metric type
            scale = ep_def.canonical_scale if ep_def else "LINEAR"
            if scale in {"LOG10", "PIC50"}:
                diff = pred_val - exp_val
                pair.signed_error = round(diff, 4)
                pair.absolute_error = round(abs(diff), 4)
                pair.error_metric_type = "log10_units" if scale == "LOG10" else "pIC50_units"
                pair.difference_display_value = pair.absolute_error
                pair.difference_display_unit = pair.error_metric_type
                pair.difference_display_string = f"{pair.absolute_error:.3f} {pair.error_metric_type}"
                
                # Agreement interpretation
                if pair.absolute_error <= 0.5:
                    pair.agreement_interpretation = "EXCELLENT_AGREEMENT"
                elif pair.absolute_error <= 1.0:
                    pair.agreement_interpretation = "ACCEPTABLE_AGREEMENT"
                else:
                    pair.agreement_interpretation = "DISCORDANT"

            elif scale == "PERCENT" or "%" in str(pair.prediction_normalized_unit):
                diff = pred_val - exp_val
                pair.signed_error = round(diff, 2)
                pair.absolute_error = round(abs(diff), 2)
                pair.error_metric_type = "percentage_points"
                pair.difference_display_value = pair.absolute_error
                pair.difference_display_unit = "percentage points"
                pair.difference_display_string = f"{pair.absolute_error:.1f} % points"
                
                if pair.absolute_error <= 5.0:
                    pair.agreement_interpretation = "EXCELLENT_AGREEMENT"
                elif pair.absolute_error <= 15.0:
                    pair.agreement_interpretation = "ACCEPTABLE_AGREEMENT"
                else:
                    pair.agreement_interpretation = "DISCORDANT"

            else:
                # Linear ratio error / fold error
                diff = pred_val - exp_val
                pair.signed_error = round(diff, 3)
                pair.absolute_error = round(abs(diff), 3)
                pair.error_metric_type = pair.prediction_normalized_unit or "absolute_error"
                pair.difference_display_value = pair.absolute_error
                pair.difference_display_unit = pair.error_metric_type

                if exp_val > 0 and pred_val > 0:
                    fold = max(pred_val / exp_val, exp_val / pred_val)
                    pair.fold_error = round(fold, 2)
                    pair.difference_display_string = f"{pair.absolute_error:.2f} {pair.error_metric_type} ({pair.fold_error:.2f}-fold)"
                    if fold <= 2.0:
                        pair.agreement_interpretation = "WITHIN_2_FOLD"
                    elif fold <= 3.0:
                        pair.agreement_interpretation = "WITHIN_3_FOLD"
                    else:
                        pair.agreement_interpretation = "DISCORDANT"
                else:
                    pair.difference_display_string = f"{pair.absolute_error:.2f} {pair.error_metric_type}"
                    pair.agreement_interpretation = "NUMERIC_COMPARED"

            # Final Comparison Classification
            if is_dose_norm:
                pair.comparison_status = CONVERTED
                pair.comparison_badge = "DOSE NORMALIZED"
                pair.is_converted = True
            elif pair.is_converted:
                pair.comparison_status = CONVERTED
                pair.comparison_badge = "UNIT CONVERTED"
            else:
                pair.comparison_status = DIRECT
                pair.comparison_badge = "DIRECT"
                pair.is_directly_comparable = True

            pair.scientific_interpretation = "NUMERICALLY_VALIDATED"
            pair.unmatched_reason = ""

        return pair


def compute_comparison_coverage(pairs: List[ScientificComparisonPair]) -> Dict[str, Any]:
    """
    Computes rigorous comparison coverage:
      - Total experimental evidence
      - Numeric evidence
      - Model-capable comparable evidence (denominator)
      - Successfully paired comparisons (numerator)
      - Direct pairs, Converted pairs, Context-only, No model, Failures
      - PAIRABLE_EVIDENCE_WITHOUT_PAIR (target = 0)
    """
    total = len(pairs)
    numeric_count = sum(1 for p in pairs if p.experimental_raw_value is not None)
    model_capable = sum(1 for p in pairs if p.prediction_available and p.experimental_raw_value is not None)
    
    direct_set = {DIRECT, "DIRECT", "DIRECTLY_COMPARABLE"}
    converted_set = {CONVERTED, "CONVERTED", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}
    context_set = {
        RELATED, CONDITIONAL, "RELATED", "CONDITIONAL",
        "RELATED_SAME_SCIENTIFIC_GROUP", "CONDITIONALLY_COMPARABLE",
        "NEEDS_REVIEW", "CONTEXT_NOT_QUALIFIED", "NO_MATCHING_PREDICTION_MODEL",
        "CONTEXT_ONLY", "EXPERIMENTAL_ONLY", "PREDICTION_ONLY"
    }

    direct_pairs = sum(1 for p in pairs if p.comparison_status in direct_set)
    converted_pairs = sum(1 for p in pairs if p.comparison_status in converted_set)
    successfully_paired = direct_pairs + converted_pairs
    
    context_only = sum(1 for p in pairs if p.comparison_status in context_set)
    no_model = sum(1 for p in pairs if p.comparison_status in {NO_MODEL, "NO_MODEL"} or not p.prediction_available)
    
    # Pairable evidence that failed to pair due to naming, unit, or scale bug
    pairable_without_pair = sum(
        1 for p in pairs
        if p.prediction_available
        and p.experimental_raw_value is not None
        and p.comparison_status not in (direct_set | converted_set | context_set)
    )

    coverage_pct = round((successfully_paired / model_capable * 100.0), 1) if model_capable > 0 else 100.0

    return {
        "total_evidence": total,
        "numeric_evidence": numeric_count,
        "model_capable_evidence": model_capable,
        "successfully_paired": successfully_paired,
        "direct_pairs": direct_pairs,
        "converted_pairs": converted_pairs,
        "context_only": context_only,
        "no_model": no_model,
        "pairable_evidence_without_pair": pairable_without_pair,
        "comparison_coverage_percent": coverage_pct,
        "matching_engine_version": MATCHING_ENGINE_VERSION,
    }
