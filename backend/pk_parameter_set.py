"""
PK Critical Parameter Set & 1-Compartment Disposition Engine.
============================================================

Scientific Architecture:
1. Canonical fu / PPB:
   - fu = max(0.0001, min(1.0, (100.0 - PPB) / 100.0))
   - Strictly bounded at fu >= 0.0001 to prevent unphysical unbound clearance.
2. Well-Stirred Hepatic IVIVE:
   - Standard physiological parameters:
     * Qh = 20.7 mL/min/kg (1449 mL/min for 70 kg human)
     * Liver weight = 25.5 g liver / kg body weight
     * Microsomal protein (MPPGL) = 45.0 mg MSP / g liver
     * Scaling factor = 25.5 * 45.0 = 1147.5 mg MSP / kg body weight
     * Blood-to-plasma ratio Rb = 1.0 (default)
   - CLh,blood = (Qh * fu * CLint,scaled) / (Qh + fu * CLint,scaled / Rb)
   - CLh,plasma = CLh,blood * Rb
   - Extraction ratio Eh = CLh,blood / Qh
   - Hepatic bioavailability Fh = 1.0 - Eh
3. Route-Aware Absorption:
   - Caco-2 Papp (log10 cm/s) mapped to intestinal fraction absorbed (Fa)
   - Gut availability (Fg)
   - Total oral bioavailability F = Fa * Fg * Fh
   - Absorption rate constant ka (1/h)
4. 1-Compartment Disposition & Simulation:
   - IV bolus: C(t) = (Dose / Vd) * exp(-ke * t)
   - Oral first-order: C(t) = (F * Dose * ka) / (Vd * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
   - Non-compartmental parameters: C0 / Cmax, Tmax, AUC0-inf, t1/2, CL (L/h/70kg), Vss (L/70kg)
5. Analytical & Monte Carlo Uncertainty Propagation:
   - 90% confidence intervals (5th - 95th percentiles) across disposition parameters.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import numpy as np

# Standard Human Physiological Constants (Davies & Morris 1993, Ring et al. 2011, Barter et al. 2007)
HUMAN_QH_ML_MIN_KG = 20.7142857
HUMAN_LIVER_WT_G_KG = 25.5
HUMAN_MPPGL_MG_G = 45.0
HUMAN_MICROSOMAL_SCALING_FACTOR = HUMAN_LIVER_WT_G_KG * HUMAN_MPPGL_MG_G  # 1147.5 mg MSP / kg
HUMAN_DEFAULT_BW_KG = 70.0
HUMAN_DEFAULT_RB = 1.0

# Physical Lower Bound for fu (fraction unbound)
MIN_FU_BOUND = 0.0001
MAX_FU_BOUND = 1.0000


@dataclass
class PKParameterSet:
    """Canonical representation of in vivo pharmacokinetic disposition parameters."""
    compound_name: str
    smiles: str
    species: str
    route: str
    dose_mg: float
    body_weight_kg: float
    
    # Primary Disposition Parameters
    cl_plasma_ml_min_kg: float
    cl_plasma_l_hr_70kg: float
    vdss_l_kg: float
    vdss_l_70kg: float
    half_life_hr: float
    ke_hr_inv: float
    
    # Hepatic & Absorption Components
    cl_int_hlm_ml_min_kg: float
    cl_int_scaled_ml_min_kg: float
    fu_plasma: float
    ppb_percent: float
    extraction_ratio_eh: float
    fh_hepatic_availability: float
    fa_fraction_absorbed: float
    fg_gut_availability: float
    f_oral_est: float
    ka_hr_inv: float
    
    # Simulated Non-Compartmental Parameters
    cmax_ng_ml: float
    tmax_hr: float
    auc_inf_ng_hr_ml: float
    c0_ng_ml: Optional[float] = None
    
    # Status & Metadata
    status: str = "COMPLETE"
    confidence: str = "MEDIUM"
    assumptions: List[str] = field(default_factory=list)
    confidence_interval_90: Dict[str, Dict[str, float]] = field(default_factory=dict)
    provenance: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def canonical_fu_from_ppb(ppb_percent: Optional[float]) -> float:
    """
    Derive canonical fraction unbound (fu) from plasma protein binding rate (% bound).
    Strictly applies physical lower bound 0.0001 (0.01% free) to prevent unbounded clearance.
    """
    if ppb_percent is None or math.isnan(ppb_percent):
        return 0.10  # Standard fallback: 10% unbound
    fu = (100.0 - float(ppb_percent)) / 100.0
    return max(MIN_FU_BOUND, min(MAX_FU_BOUND, fu))


def calculate_well_stirred_clearance(
    cl_int_input: float,
    fu: float,
    rb: float = HUMAN_DEFAULT_RB,
    is_raw_microsomal: bool = False,
    qh_ml_min_kg: float = HUMAN_QH_ML_MIN_KG,
    liver_wt_g_kg: float = HUMAN_LIVER_WT_G_KG,
    mppgl_mg_g: float = HUMAN_MPPGL_MG_G,
) -> Dict[str, float]:
    """
    Well-stirred liver model clearance IVIVE.
    CLh = (Qh * fu * CLint,scaled) / (Qh + fu * CLint,scaled / Rb)
    """
    if is_raw_microsomal:
        # Input in uL/min/mg MSP -> scale by (liver_wt * mppgl) / 1000 = 1.1475 mL/min/kg
        scaling_factor = (liver_wt_g_kg * mppgl_mg_g) / 1000.0
        cl_int_scaled = cl_int_input * scaling_factor
    else:
        # Already scaled to mL/min/kg
        cl_int_scaled = cl_int_input

    cl_int_scaled = max(0.0, cl_int_scaled)
    fu = max(MIN_FU_BOUND, min(MAX_FU_BOUND, fu))
    rb = max(0.2, min(5.0, rb))

    # Unbound intrinsic clearance in blood
    cl_u_b = (fu * cl_int_scaled) / rb
    
    # Well-stirred clearance in blood
    cl_h_blood = (qh_ml_min_kg * cl_u_b) / (qh_ml_min_kg + cl_u_b)
    
    # Plasma clearance = Blood clearance * Rb
    cl_h_plasma = cl_h_blood * rb
    
    # Extraction ratio
    eh = cl_h_blood / qh_ml_min_kg
    eh = max(0.0, min(0.9999, eh))
    
    # Hepatic bioavailability
    fh = 1.0 - eh

    return {
        "cl_int_scaled_ml_min_kg": round(cl_int_scaled, 3),
        "cl_h_blood_ml_min_kg": round(cl_h_blood, 3),
        "cl_h_plasma_ml_min_kg": round(cl_h_plasma, 3),
        "cl_plasma_l_hr_70kg": round(cl_h_plasma * 60.0 * 70.0 / 1000.0, 3),
        "extraction_ratio_eh": round(eh, 4),
        "hepatic_availability_fh": round(fh, 4),
    }


def estimate_oral_absorption_components(
    caco2_log10_cm_s: Optional[float] = None,
    solubility_logs: Optional[float] = None,
    dose_mg: float = 100.0,
) -> Tuple[float, float, float]:
    """
    Estimate Fa (fraction absorbed), Fg (gut availability), and ka (absorption rate constant 1/h).
    """
    # 1. Permeability component Fa from Caco-2
    if caco2_log10_cm_s is not None:
        papp = 10.0 ** caco2_log10_cm_s
        if caco2_log10_cm_s >= -5.0:  # High permeability (> 10 * 10^-6 cm/s)
            fa = 0.95
            ka = 1.20
        elif caco2_log10_cm_s >= -5.5:  # Moderate-high (3.16 to 10 * 10^-6 cm/s)
            fa = 0.80
            ka = 0.85
        elif caco2_log10_cm_s >= -6.0:  # Moderate (1 to 3.16 * 10^-6 cm/s)
            fa = 0.60
            ka = 0.55
        elif caco2_log10_cm_s >= -6.5:  # Moderate-low
            fa = 0.40
            ka = 0.35
        else:  # Low permeability (< 0.3 * 10^-6 cm/s)
            fa = 0.20
            ka = 0.25
    else:
        fa = 0.80
        ka = 0.80

    # 2. Solubility modulation (dissolution rate / Maximum Absorbable Dose)
    if solubility_logs is not None:
        # logS in log10(mol/L)
        if solubility_logs < -4.5:  # Poorly soluble (< 30 uM)
            fa = min(fa, 0.50)
            ka = min(ka, 0.40)
        elif solubility_logs < -3.5:  # Moderate solubility
            fa = min(fa, 0.75)
            ka = min(ka, 0.70)

    # 3. Gut availability Fg (first-pass gut metabolism, e.g. CYP3A4 in enterocytes)
    # Default high unless extensive gut metabolism is known
    fg = 0.90

    return round(fa, 3), round(fg, 3), round(ka, 3)


def simulate_one_compartment_disposition(
    dose_mg: float,
    route: str,
    cl_plasma_ml_min_kg: float,
    vdss_l_kg: float,
    f_oral: float = 1.0,
    ka_hr_inv: float = 1.0,
    body_weight_kg: float = HUMAN_DEFAULT_BW_KG,
    time_points_hr: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """
    Simulate 1-compartment IV bolus or Oral first-order disposition profile.
    Units:
      Dose: mg
      CL: mL/min/kg -> L/h total
      Vd: L/kg -> L total
      Concentration: ng/mL (ug/L)
      AUC: ng*h/mL
      t1/2: h
    """
    route = route.upper()
    total_vd_l = vdss_l_kg * body_weight_kg
    total_cl_l_hr = (cl_plasma_ml_min_kg * 60.0 / 1000.0) * body_weight_kg
    
    if total_vd_l <= 0.0 or total_cl_l_hr <= 0.0:
        raise ValueError("Vd and CL must be positive")

    ke = total_cl_l_hr / total_vd_l
    half_life = math.log(2.0) / ke

    if time_points_hr is None:
        # Generate default time course across 5 half-lives
        max_t = min(72.0, max(12.0, 5.0 * half_life))
        time_points_hr = list(np.linspace(0.0, max_t, 100))

    concentrations = []
    dose_ug = dose_mg * 1000.0

    if route == "IV":
        c0 = dose_ug / total_vd_l
        cmax = c0
        tmax = 0.0
        auc_inf = dose_ug / total_cl_l_hr

        for t in time_points_hr:
            c = c0 * math.exp(-ke * t)
            concentrations.append(round(c, 2))
            
        sim_result = {
            "route": "IV",
            "c0_ng_ml": round(c0, 2),
            "cmax_ng_ml": round(cmax, 2),
            "tmax_hr": 0.0,
            "auc_inf_ng_hr_ml": round(auc_inf, 2),
            "half_life_hr": round(half_life, 2),
            "ke_hr_inv": round(ke, 4),
            "time_course": [{"time_hr": round(t, 2), "conc_ng_ml": c} for t, c in zip(time_points_hr, concentrations)],
        }
    else:  # Oral
        c0 = 0.0
        f = max(0.01, min(1.0, f_oral))
        ka = max(0.05, ka_hr_inv)
        auc_inf = (f * dose_ug) / total_cl_l_hr

        if abs(ka - ke) < 1e-5:
            tmax = 1.0 / ke
            cmax = (f * dose_ug / total_vd_l) * (ke * tmax) * math.exp(-ke * tmax)
        else:
            tmax = math.log(ka / ke) / (ka - ke)
            cmax = (f * dose_ug * ka) / (total_vd_l * (ka - ke)) * (math.exp(-ke * tmax) - math.exp(-ka * tmax))

        cmax = max(0.0, cmax)
        tmax = max(0.0, tmax)

        for t in time_points_hr:
            if abs(ka - ke) < 1e-5:
                c = (f * dose_ug / total_vd_l) * (ke * t) * math.exp(-ke * t)
            else:
                c = (f * dose_ug * ka) / (total_vd_l * (ka - ke)) * (math.exp(-ke * t) - math.exp(-ka * t))
            concentrations.append(round(max(0.0, c), 2))

        sim_result = {
            "route": "ORAL",
            "c0_ng_ml": 0.0,
            "cmax_ng_ml": round(cmax, 2),
            "tmax_hr": round(tmax, 2),
            "auc_inf_ng_hr_ml": round(auc_inf, 2),
            "half_life_hr": round(half_life, 2),
            "ke_hr_inv": round(ke, 4),
            "time_course": [{"time_hr": round(t, 2), "conc_ng_ml": c} for t, c in zip(time_points_hr, concentrations)],
        }

    return sim_result


def propagate_pk_uncertainty(
    cl_mean: float,
    cl_rel_std: float = 0.30,
    vd_mean: float = 1.0,
    vd_rel_std: float = 0.35,
    f_mean: float = 0.70,
    f_rel_std: float = 0.25,
    ka_mean: float = 1.0,
    dose_mg: float = 100.0,
    route: str = "ORAL",
    n_draws: int = 1000,
    seed: int = 42,
) -> Dict[str, Dict[str, float]]:
    """
    Monte Carlo uncertainty propagation sampling from parameter distributions.
    Emits 5th, 50th, and 95th percentiles (90% confidence interval).
    """
    rng = np.random.default_rng(seed)
    
    # Log-normal distribution for clearance & volume
    cl_draws = rng.lognormal(mean=math.log(max(0.1, cl_mean)), sigma=cl_rel_std, size=n_draws)
    vd_draws = rng.lognormal(mean=math.log(max(0.05, vd_mean)), sigma=vd_rel_std, size=n_draws)
    
    # Truncated normal / beta for oral bioavailability F
    f_draws = np.clip(rng.normal(loc=f_mean, scale=f_mean * f_rel_std, size=n_draws), 0.05, 0.99)
    
    ke_draws = (cl_draws * 60.0 / 1000.0) / vd_draws
    thalf_draws = np.log(2.0) / ke_draws

    dose_ug = dose_mg * 1000.0
    bw = HUMAN_DEFAULT_BW_KG

    if route.upper() == "IV":
        auc_draws = dose_ug / (cl_draws * bw * 60.0 / 1000.0)
        cmax_draws = dose_ug / (vd_draws * bw)
    else:
        auc_draws = (f_draws * dose_ug) / (cl_draws * bw * 60.0 / 1000.0)
        # Approximate Cmax draws
        cmax_draws = []
        for f_i, vd_i, ke_i in zip(f_draws, vd_draws, ke_draws):
            ka = ka_mean
            if abs(ka - ke_i) < 1e-5:
                tm = 1.0 / ke_i
                cm = (f_i * dose_ug / (vd_i * bw)) * (ke_i * tm) * math.exp(-ke_i * tm)
            else:
                tm = math.log(ka / ke_i) / (ka - ke_i)
                cm = (f_i * dose_ug * ka) / (vd_i * bw * (ka - ke_i)) * (math.exp(-ke_i * tm) - math.exp(-ka * tm))
            cmax_draws.append(max(0.0, cm))
        cmax_draws = np.array(cmax_draws)

    def p_stats(arr: np.ndarray) -> Dict[str, float]:
        return {
            "p05": round(float(np.percentile(arr, 5)), 2),
            "p50": round(float(np.percentile(arr, 50)), 2),
            "p95": round(float(np.percentile(arr, 95)), 2),
        }

    return {
        "cl_plasma_ml_min_kg": p_stats(cl_draws),
        "vdss_l_kg": p_stats(vd_draws),
        "half_life_hr": p_stats(thalf_draws),
        "auc_inf_ng_hr_ml": p_stats(auc_draws),
        "cmax_ng_ml": p_stats(cmax_draws),
    }


def build_pk_parameter_set(
    smiles: str,
    compound_name: str = "Candidate Compound",
    dose_mg: float = 100.0,
    route: str = "ORAL",
    body_weight_kg: float = HUMAN_DEFAULT_BW_KG,
    cl_int_hlm_ml_min_kg: Optional[float] = None,
    ppb_percent: Optional[float] = None,
    vdss_l_kg: Optional[float] = None,
    caco2_log10_cm_s: Optional[float] = None,
    solubility_logs: Optional[float] = None,
) -> PKParameterSet:
    """
    Build a complete PKParameterSet from available experimental or predicted parameters.
    """
    # 1. Canonical fu from PPB
    if ppb_percent is None:
        ppb_percent = 90.0  # default assumption
        ppb_assumed = True
    else:
        ppb_assumed = False
    fu = canonical_fu_from_ppb(ppb_percent)

    # 2. HLM intrinsic clearance
    if cl_int_hlm_ml_min_kg is None:
        cl_int_hlm = 15.0  # default moderate clearance (mL/min/kg)
        clint_assumed = True
    else:
        cl_int_hlm = float(cl_int_hlm_ml_min_kg)
        clint_assumed = False

    ivive_res = calculate_well_stirred_clearance(cl_int_hlm, fu, rb=HUMAN_DEFAULT_RB)

    # 3. Volume of Distribution (Vdss)
    if vdss_l_kg is None:
        # Predict Vdss via Lombardo ionization model from smiles
        from backend.ionization import analyze_ionization, IonizationClass
        ion = analyze_ionization(smiles)
        ion_class = ion.get("ionization_class", IonizationClass.NEUTRAL)
        clogp = ion.get("clogp", 1.5)
        logd74 = ion.get("physiological_state_7_4", {}).get("estimated_logd74", clogp)
        eff_lipo = max(-1.5, min(4.5, float(logd74 if logd74 is not None else clogp)))

        if ion_class == IonizationClass.ACID:
            vd_est = 0.08 + 0.15 * fu + 0.05 * fu * (10.0 ** (0.2 * eff_lipo))
            vd_est = max(0.05, min(1.5, vd_est))
        elif ion_class == IonizationClass.BASE:
            vd_est = 0.6 + 0.4 * fu + 0.30 * fu * (10.0 ** (0.35 * eff_lipo))
            vd_est = max(0.2, min(30.0, vd_est))
        elif ion_class in (IonizationClass.ZWITTERION_POSSIBLE, IonizationClass.AMPHOLYTE):
            vd_est = 0.3 + 0.2 * fu + 0.10 * fu * (10.0 ** (0.2 * eff_lipo))
            vd_est = max(0.1, min(5.0, vd_est))
        else:
            vd_est = 0.6 + 0.4 * fu + 0.15 * fu * (10.0 ** (0.3 * eff_lipo))
            vd_est = max(0.1, min(20.0, vd_est))
        vdss_val = round(vd_est, 3)
        vdss_assumed = True
    else:
        vdss_val = float(vdss_l_kg)
        vdss_assumed = False

    # 4. Oral Absorption
    fa, fg, ka = estimate_oral_absorption_components(caco2_log10_cm_s, solubility_logs, dose_mg)
    f_oral = round(fa * fg * ivive_res["hepatic_availability_fh"], 3)
    if route.upper() == "IV":
        f_disp = 1.0
    else:
        f_disp = f_oral

    # 5. 1-Compartment Simulation
    cl_plasma = ivive_res["cl_h_plasma_ml_min_kg"]
    sim = simulate_one_compartment_disposition(
        dose_mg=dose_mg,
        route=route,
        cl_plasma_ml_min_kg=cl_plasma,
        vdss_l_kg=vdss_val,
        f_oral=f_disp,
        ka_hr_inv=ka,
        body_weight_kg=body_weight_kg,
    )

    # 6. Uncertainty Propagation
    ci_90 = propagate_pk_uncertainty(
        cl_mean=cl_plasma,
        vd_mean=vdss_val,
        f_mean=f_disp,
        ka_mean=ka,
        dose_mg=dose_mg,
        route=route,
    )

    assumptions = [
        "Well-stirred hepatic venous equilibrium disposition (Qh = 20.7 mL/min/kg, 1147.5 mg MSP/kg)",
        "Linear one-compartment pharmacokinetic elimination and distribution",
        f"Blood-to-plasma concentration ratio assumed Rb = {HUMAN_DEFAULT_RB}",
    ]
    if ppb_assumed:
        assumptions.append("Plasma protein binding fallback applied (PPB = 90.0%, fu = 0.10)")
    if clint_assumed:
        assumptions.append("HLM intrinsic clearance fallback applied (CLint = 15.0 mL/min/kg)")
    if vdss_assumed:
        assumptions.append(f"Steady-state volume Vdss estimated from Lombardo physicochemical model ({vdss_val} L/kg)")

    return PKParameterSet(
        compound_name=compound_name,
        smiles=smiles,
        species="Human",
        route=route.upper(),
        dose_mg=dose_mg,
        body_weight_kg=body_weight_kg,
        cl_plasma_ml_min_kg=cl_plasma,
        cl_plasma_l_hr_70kg=ivive_res["cl_plasma_l_hr_70kg"],
        vdss_l_kg=vdss_val,
        vdss_l_70kg=round(vdss_val * body_weight_kg, 2),
        half_life_hr=sim["half_life_hr"],
        ke_hr_inv=sim["ke_hr_inv"],
        cl_int_hlm_ml_min_kg=cl_int_hlm,
        cl_int_scaled_ml_min_kg=ivive_res["cl_int_scaled_ml_min_kg"],
        fu_plasma=fu,
        ppb_percent=ppb_percent,
        extraction_ratio_eh=ivive_res["extraction_ratio_eh"],
        fh_hepatic_availability=ivive_res["hepatic_availability_fh"],
        fa_fraction_absorbed=fa,
        fg_gut_availability=fg,
        f_oral_est=f_oral,
        ka_hr_inv=ka,
        cmax_ng_ml=sim["cmax_ng_ml"],
        tmax_hr=sim["tmax_hr"],
        auc_inf_ng_hr_ml=sim["auc_inf_ng_hr_ml"],
        c0_ng_ml=sim.get("c0_ng_ml"),
        status="COMPLETE",
        confidence="HIGH" if not (ppb_assumed or clint_assumed or vdss_assumed) else "MEDIUM",
        assumptions=assumptions,
        confidence_interval_90=ci_90,
        provenance={
            "engine": "Drug-OPT PK Critical Parameter Engine v1.0",
            "ivive_method": "WELL_STIRRED_HEPATIC",
            "disposition_model": "ONE_COMPARTMENT_FIRST_ORDER",
            "qh_ml_min_kg": HUMAN_QH_ML_MIN_KG,
            "liver_wt_g_kg": HUMAN_LIVER_WT_G_KG,
            "mppgl_mg_g": HUMAN_MPPGL_MG_G,
        }
    )
