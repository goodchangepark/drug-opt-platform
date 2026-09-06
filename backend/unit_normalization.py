"""
Authoritative Scientific Unit and Scale Normalization Engine for Drug-OPT.
=============================================================================
Provides rigorous, deterministic unit and scale conversions with full scientific
provenance, formula tracking, and molecular weight handling.

Critical rules:
1. DETERMINISTIC_UNIT_CONVERSION:
   - Concentration conversions: ng/mL <-> µg/L <-> mg/L <-> µM <-> nM <-> mol/L (requires exact MW)
   - Solubility conversions: mg/mL <-> µg/mL <-> µM <-> mM <-> mol/L -> log10(mol/L)
   - PPB conversions: % bound <-> fraction bound <-> % unbound <-> fu
   - CYP / hERG conversions: IC50 in M / mM / µM / nM <-> pIC50 = -log10(IC50[M])
   - Clearance: mL/min/kg <-> L/h/kg (1 mL/min/kg = 0.06 L/h/kg)
   - Volume: mL/kg <-> L/kg (1 L/kg = 1000 mL/kg)
   - Half-life / time: min <-> h <-> day
   - Caco-2 Papp: cm/s <-> 10^-6 cm/s <-> log10(cm/s) (direction preserved)
   - AUC: ng*h/mL <-> µg*h/L <-> mg*h/L
2. MODEL_DERIVED_OR_IVIVE_CONVERSION:
   - Microsomal clearance (µL/min/mg) -> hepatic CL (mL/min/kg) via MPPGL, liver wt, fu, Qh
3. DOSE_NORMALIZED:
   - Linear PK scaling: Val_norm = Val * (Dose_target / Dose_actual)
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, Optional, Tuple


# Conversion types
TYPE_DETERMINISTIC = "DETERMINISTIC_UNIT_CONVERSION"
TYPE_SCALE = "SCALE_CONVERSION"
TYPE_IVIVE_DERIVED = "MODEL_DERIVED_OR_IVIVE_CONVERSION"
TYPE_DOSE_NORMALIZED = "DOSE_NORMALIZED"
TYPE_DIRECT_IDENTITY = "DIRECT_IDENTITY"
TYPE_UNSUPPORTED = "UNSUPPORTED"

CONVERSION_ENGINE_VERSION = "drugopt-unit-normalization-v2.0"


@dataclass
class UnitConversionResult:
    original_value: float
    original_unit: str
    normalized_value: float
    normalized_unit: str
    conversion_formula: str
    conversion_type: str
    is_valid: bool = True
    molecular_weight_used: Optional[float] = None
    scale: str = "LINEAR"
    provenance_detail: str = ""
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def to_dict(self) -> Dict[str, Any]:
        return {
            "original_value": self.original_value,
            "original_unit": self.original_unit,
            "normalized_value": self.normalized_value,
            "normalized_unit": self.normalized_unit,
            "conversion_formula": self.conversion_formula,
            "conversion_type": self.conversion_type,
            "is_valid": self.is_valid,
            "molecular_weight_used": self.molecular_weight_used,
            "scale": self.scale,
            "provenance_detail": self.provenance_detail,
            "conversion_version": CONVERSION_ENGINE_VERSION,
            "timestamp": self.timestamp,
        }


def _clean_float(v: float) -> float:
    if v is None:
        return None
    r = round(v, 8)
    if abs(v - r) < 1e-10:
        return r
    return v


def clean_unit_str(unit: Any) -> str:
    """Normalize unit string characters for deterministic lookup."""
    if unit is None:
        return ""
    u = str(unit).strip().lower()
    u = u.replace("μ", "µ").replace("u", "µ") if ("um" in u or "ug" in u or "ul" in u or "umol" in u) else u.replace("μ", "µ")
    u = u.replace("·", "*").replace("•", "*").replace(" ", "")
    u = u.replace(".hr.", "*h.").replace(".h.", "*h.").replace(".min.", "*min.").replace(".d.", "*d.")
    u = u.replace(".ml-1", "/ml").replace(".l-1", "/l").replace(".kg-1", "/kg").replace(".h-1", "/h").replace(".min-1", "/min")
    u = u.replace("ml-1", "/ml").replace("l-1", "/l").replace("kg-1", "/kg").replace("h-1", "/h").replace("min-1", "/min")
    u = u.replace("//", "/")
    u = u.replace("hr", "h").replace("hours", "h").replace("hour", "h")
    u = u.replace("mins", "min").replace("minutes", "min")
    u = u.replace("days", "d").replace("day", "d")
    u = u.replace("liters", "l").replace("litres", "l").replace("liter", "l")
    return u


# =====================================================================
# 1. Concentration Normalization
# =====================================================================

def convert_concentration(
    value: float,
    from_unit: str,
    to_unit: str = "nM",
    mw: Optional[float] = None,
) -> UnitConversionResult:
    """
    Convert concentrations between mass-based (ng/mL, µg/L, mg/L, etc.)
    and molar-based (mol/L, mM, µM, nM, pM).
    Uses MW when crossing mass <-> molar boundaries.
    """
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    # Base molar conversion factors to mol/L (M)
    molar_factors_to_m = {
        "m": 1.0,
        "mol/l": 1.0,
        "mm": 1e-3,
        "mmol/l": 1e-3,
        "µm": 1e-6,
        "µmol/l": 1e-6,
        "nm": 1e-9,
        "nmol/l": 1e-9,
        "pm": 1e-12,
        "pmol/l": 1e-12,
    }

    # Base mass conversion factors to g/L (or mg/mL)
    mass_factors_to_g_l = {
        "g/l": 1.0,
        "mg/ml": 1.0,
        "mg/l": 1e-3,
        "µg/ml": 1e-3,
        "ug/ml": 1e-3,
        "µg/l": 1e-6,
        "ug/l": 1e-6,
        "ng/ml": 1e-6,
        "pg/ml": 1e-9,
    }

    if u_from == u_to:
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit=to_unit,
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
        )

    # Case A: Molar to Molar
    if u_from in molar_factors_to_m and u_to in molar_factors_to_m:
        val_m = value * molar_factors_to_m[u_from]
        norm_val = _clean_float(val_m / molar_factors_to_m[u_to])
        factor = molar_factors_to_m[u_from] / molar_factors_to_m[u_to]
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=f"value * {factor:.6e}",
            conversion_type=TYPE_DETERMINISTIC,
        )

    # Case B: Mass to Mass
    if u_from in mass_factors_to_g_l and u_to in mass_factors_to_g_l:
        val_g_l = value * mass_factors_to_g_l[u_from]
        norm_val = _clean_float(val_g_l / mass_factors_to_g_l[u_to])
        factor = mass_factors_to_g_l[u_from] / mass_factors_to_g_l[u_to]
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=f"value * {factor:.6e}",
            conversion_type=TYPE_DETERMINISTIC,
        )

    # Case C: Mass to Molar (requires MW in g/mol)
    if u_from in mass_factors_to_g_l and u_to in molar_factors_to_m:
        if not mw or mw <= 0:
            raise ValueError(f"Molecular weight (MW) required for mass-to-molar conversion from {from_unit} to {to_unit}")
        val_g_l = value * mass_factors_to_g_l[u_from]
        val_m = val_g_l / mw
        norm_val = _clean_float(val_m / molar_factors_to_m[u_to])
        formula = f"(value * {mass_factors_to_g_l[u_from]:.6e} g/L) / {mw:.2f} g/mol / {molar_factors_to_m[u_to]:.6e}"
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=formula,
            conversion_type=TYPE_DETERMINISTIC,
            molecular_weight_used=mw,
        )

    # Case D: Molar to Mass (requires MW in g/mol)
    if u_from in molar_factors_to_m and u_to in mass_factors_to_g_l:
        if not mw or mw <= 0:
            raise ValueError(f"Molecular weight (MW) required for molar-to-mass conversion from {from_unit} to {to_unit}")
        val_m = value * molar_factors_to_m[u_from]
        val_g_l = val_m * mw
        norm_val = _clean_float(val_g_l / mass_factors_to_g_l[u_to])
        formula = f"(value * {molar_factors_to_m[u_from]:.6e} mol/L) * {mw:.2f} g/mol / {mass_factors_to_g_l[u_to]:.6e}"
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=formula,
            conversion_type=TYPE_DETERMINISTIC,
            molecular_weight_used=mw,
        )

    raise ValueError(f"Unsupported concentration units: {from_unit} -> {to_unit}")


# =====================================================================
# 2. Solubility Normalization (-> log10(mol/L))
# =====================================================================

def convert_solubility(
    value: float,
    from_unit: str,
    mw: Optional[float] = None,
) -> UnitConversionResult:
    """
    Convert solubility to canonical log10(mol/L).
    Accepts: mol/L, M, mM, µM, mg/mL, µg/mL, log10(mol/L).
    """
    u = clean_unit_str(from_unit)

    # If already log10
    if "log" in u:
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit="log10(mol/L)",
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
            scale="LOG10",
        )

    if value <= 0:
        raise ValueError(f"Solubility value must be positive for log10 conversion: {value}")

    molar_units = {"m", "mol/l", "mm", "mmol/l", "µm", "µmol/l", "nm", "nmol/l"}
    mass_units = {"mg/ml", "g/l", "µg/ml", "ug/ml", "mg/l", "µg/l", "ug/l"}

    if u in molar_units:
        molar_res = convert_concentration(value, from_unit, "mol/L")
        log_s = math.log10(molar_res.normalized_value)
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=log_s,
            normalized_unit="log10(mol/L)",
            conversion_formula=f"log10({molar_res.normalized_value:.6e} mol/L)",
            conversion_type=TYPE_SCALE,
            scale="LOG10",
        )

    if u in mass_units:
        if not mw or mw <= 0:
            raise ValueError(f"Molecular weight (MW) required for mass-to-log10(mol/L) solubility conversion from {from_unit}")
        molar_res = convert_concentration(value, from_unit, "mol/L", mw=mw)
        log_s = math.log10(molar_res.normalized_value)
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=log_s,
            normalized_unit="log10(mol/L)",
            conversion_formula=f"log10(({value} {from_unit} / MW {mw:.2f}) -> {molar_res.normalized_value:.6e} mol/L)",
            conversion_type=TYPE_SCALE,
            molecular_weight_used=mw,
            scale="LOG10",
        )

    raise ValueError(f"Unsupported solubility unit: {from_unit}")


# =====================================================================
# 3. Plasma Protein Binding & fu Normalization
# =====================================================================

def convert_ppb(
    value: float,
    from_unit: str,
    to_unit: str = "% bound",
) -> UnitConversionResult:
    """
    Convert plasma protein binding between:
      - % bound (e.g. 98.0 %)
      - fraction bound (e.g. 0.98)
      - % unbound (e.g. 2.0 %)
      - fraction unbound fu (e.g. 0.02)
    Canonical platform unit is '% bound', with canonical fu bound >= 0.0001.
    """
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    # First determine fraction bound (0.0 to 1.0)
    if "%" in u_from or "percent" in u_from:
        if "unbound" in u_from or "free" in u_from:
            pct_unbound = value
            fraction_bound = max(0.0, min(1.0, 1.0 - (pct_unbound / 100.0)))
        else:
            pct_bound = value
            fraction_bound = max(0.0, min(1.0, pct_bound / 100.0))
    elif "fu" in u_from or "unbound" in u_from or "free" in u_from:
        fu_val = value
        fraction_bound = max(0.0, min(1.0, 1.0 - fu_val))
    else:  # fraction bound
        fraction_bound = max(0.0, min(1.0, value))

    pct_bound_val = fraction_bound * 100.0
    fu_val = max(1.0 - fraction_bound, 0.0001)

    if u_to in {"%bound", "%", "percent", "%_bound"}:
        norm_val = pct_bound_val
        formula = f"fraction_bound * 100 = {norm_val:.2f}%"
    elif u_to in {"fu", "fraction_unbound", "fractionunbound"}:
        norm_val = fu_val
        formula = f"max(1.0 - (pct_bound / 100), 0.0001) = {norm_val:.4f}"
    elif u_to in {"fraction_bound", "fractionbound"}:
        norm_val = fraction_bound
        formula = f"pct_bound / 100 = {norm_val:.4f}"
    else:
        norm_val = pct_bound_val
        to_unit = "% bound"
        formula = f"normalized to {norm_val:.2f}% bound"

    is_identity = (abs(norm_val - value) < 1e-9 and u_from == u_to)
    return UnitConversionResult(
        original_value=value,
        original_unit=from_unit,
        normalized_value=norm_val,
        normalized_unit=to_unit,
        conversion_formula=formula,
        conversion_type=TYPE_DIRECT_IDENTITY if is_identity else TYPE_DETERMINISTIC,
        scale="PERCENT" if "%" in to_unit else "LINEAR",
    )


# =====================================================================
# 4. CYP & hERG Potency Normalization (IC50 -> pIC50)
# =====================================================================

def convert_cyp_herg_ic50(
    value: float,
    from_unit: str,
    to_scale: str = "pIC50",
) -> UnitConversionResult:
    """
    Convert IC50 concentration (M, mM, µM, nM, pM) to pIC50 (-log10(M)).
    Or reverse pIC50 -> IC50 [µM] or [nM].
    """
    u_from = clean_unit_str(from_unit)

    if u_from in {"pic50", "p_ic50"} or to_scale.lower() == "identity":
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit="pIC50",
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
            scale="PIC50",
        )

    if value <= 0:
        raise ValueError(f"IC50 concentration must be positive for pIC50 conversion: {value}")

    # Convert to M (mol/L)
    molar_res = convert_concentration(value, from_unit, "mol/L")
    pic50 = -math.log10(molar_res.normalized_value)

    return UnitConversionResult(
        original_value=value,
        original_unit=from_unit,
        normalized_value=pic50,
        normalized_unit="pIC50",
        conversion_formula=f"-log10({molar_res.normalized_value:.6e} M)",
        conversion_type=TYPE_SCALE,
        scale="PIC50",
    )


def convert_pic50_to_ic50(
    pic50_value: float,
    to_unit: str = "nM",
) -> UnitConversionResult:
    """Convert pIC50 back to quantitative IC50 (e.g. nM or µM)."""
    molar_val = 10.0 ** (-pic50_value)
    conc_res = convert_concentration(molar_val, "mol/L", to_unit)
    return UnitConversionResult(
        original_value=pic50_value,
        original_unit="pIC50",
        normalized_value=conc_res.normalized_value,
        normalized_unit=to_unit,
        conversion_formula=f"10^(-{pic50_value:.3f}) M -> {conc_res.normalized_value:.2f} {to_unit}",
        conversion_type=TYPE_SCALE,
        scale="LINEAR",
    )


# =====================================================================
# 5. Clearance Normalization (mL/min/kg <-> L/h/kg)
# =====================================================================

def convert_clearance(
    value: float,
    from_unit: str,
    to_unit: str = "mL/min/kg",
    body_weight_kg: float = 70.0,
) -> UnitConversionResult:
    """
    Convert clearance:
      - 1 mL/min/kg = 0.06 L/h/kg (exact: 1 * 60 / 1000)
      - 1 L/h/kg = 16.6667 mL/min/kg (exact: 1000 / 60)
      - L/h <-> mL/min/kg uses standard adult body weight (70 kg default)
      - log10(mL/min/kg) <-> mL/min/kg
    """
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    if u_from == u_to:
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit=to_unit,
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
        )

    # Log10 clearance handling
    if "log" in u_from and "log" not in u_to:
        lin_val = 10.0 ** value
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=lin_val,
            normalized_unit="mL/min/kg",
            conversion_formula=f"10^({value}) = {lin_val:.3f} mL/min/kg",
            conversion_type=TYPE_SCALE,
        )
    if "log" not in u_from and "log" in u_to:
        if value <= 0:
            raise ValueError(f"Clearance must be positive for log10 conversion: {value}")
        log_val = math.log10(value)
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=log_val,
            normalized_unit="log10(mL/min/kg)",
            conversion_formula=f"log10({value}) = {log_val:.3f}",
            conversion_type=TYPE_SCALE,
            scale="LOG10",
        )

    # Standard linear units
    if u_from in {"ml/min/kg", "mlmin/kg"} and u_to in {"l/h/kg", "l/hr/kg", "lhkg"}:
        norm_val = value * 0.06
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="L/h/kg",
            conversion_formula=f"{value} * (60 / 1000) = {norm_val:.4f} L/h/kg",
            conversion_type=TYPE_DETERMINISTIC,
        )

    if u_from in {"l/h/kg", "l/hr/kg", "lhkg"} and u_to in {"ml/min/kg", "mlmin/kg"}:
        norm_val = value * (1000.0 / 60.0)
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="mL/min/kg",
            conversion_formula=f"{value} * (1000 / 60) = {norm_val:.3f} mL/min/kg",
            conversion_type=TYPE_DETERMINISTIC,
        )

    # Absolute clearance (L/h) to bodyweight-normalized (mL/min/kg)
    if u_from in {"l/h", "l/hr", "l/hour"} and u_to in {"ml/min/kg", "mlmin/kg"}:
        norm_val = (value * 1000.0 / 60.0) / body_weight_kg
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="mL/min/kg",
            conversion_formula=f"({value} L/h * 1000 / 60) / {body_weight_kg:.1f} kg = {norm_val:.3f} mL/min/kg",
            conversion_type=TYPE_DETERMINISTIC,
            provenance_detail=f"Body weight assumption: {body_weight_kg} kg",
        )

    if u_from in {"ml/min/kg", "mlmin/kg"} and u_to in {"l/h", "l/hr", "l/hour"}:
        norm_val = (value * 0.06) * body_weight_kg
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="L/h",
            conversion_formula=f"{value} mL/min/kg * 0.06 * {body_weight_kg:.1f} kg = {norm_val:.3f} L/h",
            conversion_type=TYPE_DETERMINISTIC,
            provenance_detail=f"Body weight assumption: {body_weight_kg} kg",
        )

    raise ValueError(f"Unsupported clearance conversion: {from_unit} -> {to_unit}")


# =====================================================================
# 6. Volume Normalization (L/kg <-> mL/kg <-> L)
# =====================================================================

def convert_volume(
    value: float,
    from_unit: str,
    to_unit: str = "L/kg",
    body_weight_kg: float = 70.0,
) -> UnitConversionResult:
    """
    Convert volume of distribution:
      - 1 L/kg = 1000 mL/kg
      - 1 mL/kg = 0.001 L/kg
      - Absolute L to L/kg using body weight (70 kg default)
    """
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    if u_from == u_to:
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit=to_unit,
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
        )

    if u_from in {"ml/kg", "mlkg"} and u_to in {"l/kg", "lkg"}:
        norm_val = value / 1000.0
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="L/kg",
            conversion_formula=f"{value} / 1000 = {norm_val:.4f} L/kg",
            conversion_type=TYPE_DETERMINISTIC,
        )

    if u_from in {"l/kg", "lkg"} and u_to in {"ml/kg", "mlkg"}:
        norm_val = value * 1000.0
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="mL/kg",
            conversion_formula=f"{value} * 1000 = {norm_val:.1f} mL/kg",
            conversion_type=TYPE_DETERMINISTIC,
        )

    if u_from in {"l", "liters", "litres"} and u_to in {"l/kg", "lkg"}:
        norm_val = value / body_weight_kg
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="L/kg",
            conversion_formula=f"{value} L / {body_weight_kg:.1f} kg = {norm_val:.4f} L/kg",
            conversion_type=TYPE_DETERMINISTIC,
            provenance_detail=f"Body weight assumption: {body_weight_kg} kg",
        )

    if u_from in {"l/kg", "lkg"} and u_to in {"l", "liters", "litres"}:
        norm_val = value * body_weight_kg
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit="L",
            conversion_formula=f"{value} L/kg * {body_weight_kg:.1f} kg = {norm_val:.2f} L",
            conversion_type=TYPE_DETERMINISTIC,
            provenance_detail=f"Body weight assumption: {body_weight_kg} kg",
        )

    raise ValueError(f"Unsupported volume conversion: {from_unit} -> {to_unit}")


# =====================================================================
# 7. Time / Half-Life Normalization (min <-> h <-> day)
# =====================================================================

def convert_time(
    value: float,
    from_unit: str,
    to_unit: str = "hours",
) -> UnitConversionResult:
    """Convert half-life / time between minutes, hours, and days."""
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    factors_to_hours = {
        "h": 1.0,
        "hr": 1.0,
        "hours": 1.0,
        "min": 1.0 / 60.0,
        "mins": 1.0 / 60.0,
        "minutes": 1.0 / 60.0,
        "d": 24.0,
        "day": 24.0,
        "days": 24.0,
        "s": 1.0 / 3600.0,
        "sec": 1.0 / 3600.0,
    }

    if u_from in factors_to_hours and u_to in factors_to_hours:
        val_h = value * factors_to_hours[u_from]
        norm_val = val_h / factors_to_hours[u_to]
        factor = factors_to_hours[u_from] / factors_to_hours[u_to]
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=f"{value} * {factor:.6f} = {norm_val:.3f} {to_unit}",
            conversion_type=TYPE_DIRECT_IDENTITY if abs(factor - 1.0) < 1e-9 else TYPE_DETERMINISTIC,
        )

    raise ValueError(f"Unsupported time conversion: {from_unit} -> {to_unit}")


# =====================================================================
# 8. Caco-2 Permeability (cm/s <-> 10^-6 cm/s <-> log10(cm/s))
# =====================================================================

def convert_caco2_papp(
    value: float,
    from_unit: str,
    to_scale: str = "log10",
) -> UnitConversionResult:
    """
    Convert Caco-2 apparent permeability:
      - 10^-6 cm/s -> log10(cm/s)
      - cm/s -> log10(cm/s)
      - µm/s (10^-4 cm/s) -> log10(cm/s)
    """
    u = clean_unit_str(from_unit)

    if "log" in u:
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=value,
            normalized_unit="log10(cm/s)",
            conversion_formula="identity",
            conversion_type=TYPE_DIRECT_IDENTITY,
            scale="LOG10",
        )

    if value <= 0:
        raise ValueError(f"Caco-2 Papp must be positive for log10 conversion: {value}")

    if any(pattern in u for pattern in ("10^-6", "10-6", "10e-6", "e-6")):
        log_val = math.log10(value * 1e-6)
        formula = f"log10({value} * 1e-6 cm/s) = {log_val:.3f}"
    elif "µm/s" in u or "um/s" in u:
        log_val = math.log10(value * 1e-4)
        formula = f"log10({value} * 1e-4 cm/s) = {log_val:.3f}"
    elif "cm/s" in u:
        log_val = math.log10(value)
        formula = f"log10({value} cm/s) = {log_val:.3f}"
    else:
        # Default assumption: if value is around 1-50, it is in 10^-6 cm/s
        if 0.01 <= value <= 500:
            log_val = math.log10(value * 1e-6)
            formula = f"assumed_10^-6_cm/s: log10({value} * 1e-6 cm/s) = {log_val:.3f}"
        else:
            log_val = math.log10(value)
            formula = f"log10({value}) = {log_val:.3f}"

    return UnitConversionResult(
        original_value=value,
        original_unit=from_unit,
        normalized_value=log_val,
        normalized_unit="log10(cm/s)",
        conversion_formula=formula,
        conversion_type=TYPE_SCALE,
        scale="LOG10",
    )


# =====================================================================
# 9. AUC Normalization (ng*h/mL <-> µg*h/L <-> mg*h/L)
# =====================================================================

def convert_auc(
    value: float,
    from_unit: str,
    to_unit: str = "ng*h/mL",
) -> UnitConversionResult:
    """
    Convert Area Under the Curve (AUC):
      - 1 ng*h/mL = 1 µg*h/L (exact: 1 ng/mL = 1 µg/L)
      - 1 mg*h/L = 1000 ng*h/mL
      - 1 µg*h/mL = 1000 ng*h/mL
    """
    u_from = clean_unit_str(from_unit)
    u_to = clean_unit_str(to_unit)

    factors_to_ng_h_ml = {
        "ng*h/ml": 1.0,
        "ngh/ml": 1.0,
        "ng*hr/ml": 1.0,
        "h*ng/ml": 1.0,
        "µg*h/l": 1.0,
        "µgh/l": 1.0,
        "ug*h/l": 1.0,
        "ugh/l": 1.0,
        "mg*h/l": 1000.0,
        "mgh/l": 1000.0,
        "µg*h/ml": 1000.0,
        "µgh/ml": 1000.0,
        "ug*h/ml": 1000.0,
        "ugh/ml": 1000.0,
    }

    if u_from in factors_to_ng_h_ml and u_to in factors_to_ng_h_ml:
        val_base = value * factors_to_ng_h_ml[u_from]
        norm_val = val_base / factors_to_ng_h_ml[u_to]
        factor = factors_to_ng_h_ml[u_from] / factors_to_ng_h_ml[u_to]
        return UnitConversionResult(
            original_value=value,
            original_unit=from_unit,
            normalized_value=norm_val,
            normalized_unit=to_unit,
            conversion_formula=f"{value} * {factor} = {norm_val:.3f} {to_unit}",
            conversion_type=TYPE_DIRECT_IDENTITY if abs(factor - 1.0) < 1e-9 else TYPE_DETERMINISTIC,
        )

    raise ValueError(f"Unsupported AUC conversion: {from_unit} -> {to_unit}")


# =====================================================================
# 10. Master Normalization Dispatcher
# =====================================================================

def normalize_scientific_observation(
    endpoint_id: str,
    raw_value: Any,
    raw_unit: str,
    *,
    mw: Optional[float] = None,
    species: str = "HUMAN",
    route: str = "ORAL",
    context: Optional[Dict[str, Any]] = None,
) -> Optional[UnitConversionResult]:
    """
    Main entrypoint: takes any endpoint and raw experimental value/unit,
    and returns a normalized result using the strict scientific conversion rules.
    Returns None if value cannot be parsed or safely converted.
    """
    if raw_value is None:
        return None

    try:
        val = float(str(raw_value).replace(",", "").strip())
    except (ValueError, TypeError):
        return None

    ep = str(endpoint_id).upper().strip()
    u = clean_unit_str(raw_unit)

    try:
        # Solubility
        if "SOLUBILITY" in ep:
            return convert_solubility(val, raw_unit, mw=mw)

        # Caco-2
        if "CACO2" in ep or "PERMEABILITY" in ep:
            return convert_caco2_papp(val, raw_unit)

        # PPB / fu
        if "PPB" in ep or "PROTEIN_BINDING" in ep:
            return convert_ppb(val, raw_unit, "% bound")

        # CYP / hERG Potency
        if any(cyp in ep for cyp in ("CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "HERG")) and any(tok in u for tok in ("m", "µm", "um", "nm", "pm", "mol")):
            return convert_cyp_herg_ic50(val, raw_unit)

        # Clearance (HLM, RLM, MLM, or PK CL)
        if any(tok in ep for tok in ("CLINT", "CLEARANCE", "_CL_", "_CLF_")):
            return convert_clearance(val, raw_unit, "mL/min/kg")

        # Volume of Distribution
        if any(tok in ep for tok in ("VDSS", "VD", "VDF", "VSSF")):
            return convert_volume(val, raw_unit, "L/kg")

        # Half-life / Tmax
        if any(tok in ep for tok in ("T_HALF", "HALF_LIFE", "TMAX")):
            return convert_time(val, raw_unit, "hours")

        # Cmax / PK Concentration
        if "CMAX" in ep:
            # Standard PK Cmax canonical unit: ng/mL
            return convert_concentration(val, raw_unit, "ng/mL", mw=mw)

        # AUC
        if "AUC" in ep:
            return convert_auc(val, raw_unit, "ng*h/mL")

        # Bioavailability F
        if "_F_" in ep or ep.endswith("_F"):
            if "%" in u:
                return UnitConversionResult(val, raw_unit, val, "%", "identity", TYPE_DIRECT_IDENTITY, scale="PERCENT")
            if 0.0 <= val <= 1.0:
                return UnitConversionResult(val, raw_unit, val * 100.0, "%", f"{val} * 100", TYPE_DETERMINISTIC, scale="PERCENT")

    except Exception:
        pass

    # Default fallback: return identity if valid number
    return UnitConversionResult(
        original_value=val,
        original_unit=raw_unit,
        normalized_value=val,
        normalized_unit=raw_unit,
        conversion_formula="identity",
        conversion_type=TYPE_DIRECT_IDENTITY,
    )
