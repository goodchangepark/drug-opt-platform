"""
Curated Independent Clinical PK Validation Cohort & Semantic Guardrails (v1.0).
==============================================================================

Implements Directives 3, 4, 5, 15, 16, 17, 18:
1. 30 Approved Small-Molecule Therapeutics with Regulatory / Literature Ground Truth.
2. High structural and pharmacokinetic diversity:
   - Ionization: Acids, Bases, Neutrals, Ampholytes/Zwitterions.
   - Lipophilicity: logP range -1.4 to 5.5.
   - Plasma Binding: fu 0.001 to 0.95 (Extreme, High, Moderate, Low).
   - Hepatic Extraction: Low (Eh < 0.3), Intermediate (0.3 <= Eh <= 0.7), High (Eh > 0.7).
   - Elimination Route: Hepatic metabolic, Renal filtration/secretion, Mixed.
   - Volume: Low (<0.5 L/kg) to Very High (>10 L/kg).
3. Semantic Separation Guardrails:
   - VDss != Vz != Vd/F (apparent oral volume)
   - CL (systemic) != CL/F (apparent oral clearance) != CLh (hepatic) != CLr (renal)
   - IV AUC != Oral AUC
   - Single-dose Cmax != Steady-state Cmax
4. Elimination Status Tagging:
   - HEPATIC_CL_ESTIMATE vs TOTAL_CL_PREDICTION
   - Flags significant unmodeled renal elimination (fe_renal >= 0.20)
5. Oral PK Confidence Downgrading:
   - Propagates PK_CONFIDENCE_LOW / PK_CONFIDENCE_MODERATE when F or ka are unmeasured.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass(frozen=True)
class ObservedPKSemantic:
    """Strictly annotated observed clinical pharmacokinetic observation."""
    compound_name: str
    smiles: str
    dose_mg: float
    route: str  # "IV" or "ORAL"
    regimen: str  # "SINGLE_DOSE" or "STEADY_STATE"
    formulation: str
    fasted: bool
    population: str
    source_reference: str

    # Clearances (mL/min/kg or L/h)
    cl_systemic_ml_min_kg: Optional[float] = None
    cl_hepatic_ml_min_kg: Optional[float] = None
    cl_renal_ml_min_kg: Optional[float] = None
    cl_oral_apparent_ml_min_kg: Optional[float] = None  # CL/F

    # Volumes (L/kg or L)
    vdss_l_kg: Optional[float] = None
    vz_l_kg: Optional[float] = None
    vd_oral_apparent_l_kg: Optional[float] = None  # Vd/F

    # Absorption & Half-life
    f_oral_bioavailability: Optional[float] = None
    t_half_hr: Optional[float] = None
    cmax_ng_ml: Optional[float] = None
    tmax_hr: Optional[float] = None
    auc_inf_ng_hr_ml: Optional[float] = None

    # In Vitro Preclinical Ground Truth
    ppb_percent: Optional[float] = None
    cl_int_hlm_ml_min_kg: Optional[float] = None
    caco2_log_cm_s: Optional[float] = None
    blood_to_plasma_rb: float = 1.0

    # Elimination Characteristics
    fraction_excreted_renal: float = 0.0  # fe
    extraction_category: str = "INTERMEDIATE"  # "LOW", "INTERMEDIATE", "HIGH"
    ionization_class: str = "NEUTRAL"  # "ACID", "BASE", "NEUTRAL", "AMPHOLYTE"


# Comprehensive 30-Drug Independent Clinical Validation Cohort
CLINICAL_PK_COHORT: List[ObservedPKSemantic] = [
    # 1. Acetaminophen (Neutral, Low PPB, Short t1/2, Low VDss, Mixed Hepatic/Renal)
    ObservedPKSemantic(
        compound_name="Acetaminophen",
        smiles="CC(=O)NC1=CC=C(O)C=C1",
        dose_mg=1000.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Solution",
        fasted=True,
        population="Healthy Adult",
        source_reference="DailyMed NDA 022450 / FDA Review Table 5",
        cl_systemic_ml_min_kg=4.50,
        cl_hepatic_ml_min_kg=4.00,
        cl_renal_ml_min_kg=0.50,
        vdss_l_kg=0.80,
        vz_l_kg=0.95,
        t_half_hr=2.40,
        cmax_ng_ml=28000.0,
        tmax_hr=0.25,
        auc_inf_ng_hr_ml=43000.0,
        ppb_percent=20.0,
        cl_int_hlm_ml_min_kg=15.0,
        caco2_log_cm_s=-4.9,
        fraction_excreted_renal=0.10,
        extraction_category="LOW",
        ionization_class="NEUTRAL",
    ),
    # 2. Osimertinib (Base, High PPB, High VDss, Long t1/2, CYP3A4 Hepatic)
    ObservedPKSemantic(
        compound_name="Osimertinib",
        smiles="C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C",
        dose_mg=80.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Cancer Patient",
        source_reference="FDA NDA 208065 Multidisciplinary Review",
        cl_systemic_ml_min_kg=3.40,  # 14.3 L/h / 70 kg
        cl_hepatic_ml_min_kg=3.20,
        cl_oral_apparent_ml_min_kg=4.86,
        vdss_l_kg=13.0,
        f_oral_bioavailability=0.70,
        t_half_hr=48.0,
        cmax_ng_ml=500.0,
        tmax_hr=6.0,
        auc_inf_ng_hr_ml=16500.0,
        ppb_percent=95.0,
        cl_int_hlm_ml_min_kg=35.0,
        caco2_log_cm_s=-5.2,
        fraction_excreted_renal=0.02,
        extraction_category="LOW",
        ionization_class="BASE",
    ),
    # 3. Sunvozertinib (Base, High PPB, Very High VDss, Long t1/2, CYP3A4)
    ObservedPKSemantic(
        compound_name="Sunvozertinib",
        smiles="COC1=CC(N2CC[C@H](C2)N(C)C)=C(NC(=O)C=C)C=C1NC1=NC=CC(NC2=CC(Cl)=C(F)C=C2C(C)(C)O)=N1",
        dose_mg=300.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Capsule",
        fasted=True,
        population="NSCLC Patient",
        source_reference="FDA NDA 219839 Multidiscipline Review 2025",
        cl_systemic_ml_min_kg=6.90,  # 29 L/h / 70 kg
        cl_hepatic_ml_min_kg=6.50,
        cl_oral_apparent_ml_min_kg=11.5,
        vdss_l_kg=30.2,
        f_oral_bioavailability=0.60,
        t_half_hr=50.0,
        cmax_ng_ml=619.0,
        tmax_hr=6.0,
        auc_inf_ng_hr_ml=12089.0,
        ppb_percent=92.0,
        cl_int_hlm_ml_min_kg=30.0,
        caco2_log_cm_s=-5.3,
        fraction_excreted_renal=0.05,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 4. Metformin (Strong Base, Very Low PPB, High VDss, 90%+ RENAL Elimination)
    ObservedPKSemantic(
        compound_name="Metformin",
        smiles="CN(C)C(=N)NC(N)=N",
        dose_mg=500.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral IR Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="DailyMed Metformin NDA / Tucker et al. Br J Clin Pharmacol",
        cl_systemic_ml_min_kg=8.50,  # ~600 mL/min
        cl_hepatic_ml_min_kg=0.40,   # Negligible hepatic metabolism
        cl_renal_ml_min_kg=8.10,     # Active OCT2/MATE secretion
        cl_oral_apparent_ml_min_kg=17.0,
        vdss_l_kg=4.10,
        f_oral_bioavailability=0.50,
        t_half_hr=6.20,
        cmax_ng_ml=1030.0,
        tmax_hr=2.75,
        auc_inf_ng_hr_ml=7200.0,
        ppb_percent=5.0,
        cl_int_hlm_ml_min_kg=1.5,    # Minimal microsomal clearance
        caco2_log_cm_s=-6.4,
        fraction_excreted_renal=0.90, # PREDOMINANTLY RENAL
        extraction_category="LOW",
        ionization_class="BASE",
    ),
    # 5. Warfarin (Acid, Extreme PPB 99%, Low VDss, Long t1/2, CYP2C9 Low Extraction)
    ObservedPKSemantic(
        compound_name="Warfarin",
        smiles="CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        dose_mg=10.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Goodman & Gilman 14th ed / DailyMed Coumadin",
        cl_systemic_ml_min_kg=0.045, # ~3.2 mL/min / 70 kg
        cl_hepatic_ml_min_kg=0.043,
        cl_renal_ml_min_kg=0.002,
        cl_oral_apparent_ml_min_kg=0.045,
        vdss_l_kg=0.14,
        f_oral_bioavailability=1.00,
        t_half_hr=40.0,
        cmax_ng_ml=1200.0,
        tmax_hr=2.0,
        auc_inf_ng_hr_ml=53000.0,
        ppb_percent=99.0,
        cl_int_hlm_ml_min_kg=8.0,
        caco2_log_cm_s=-4.8,
        fraction_excreted_renal=0.01,
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 6. Midazolam (Base, High PPB, Moderate VDss, Short t1/2, CYP3A4 Probe High Extraction)
    ObservedPKSemantic(
        compound_name="Midazolam",
        smiles="Cc1ncc2c(n1)N=C(c1ccccc1F)CN2c1ccc(Cl)cc1",
        dose_mg=5.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Solution",
        fasted=True,
        population="Healthy Adult",
        source_reference="DailyMed / Heidegger et al. Clin Pharmacokinet",
        cl_systemic_ml_min_kg=6.50,  # ~450 mL/min
        cl_hepatic_ml_min_kg=6.30,
        cl_renal_ml_min_kg=0.05,
        vdss_l_kg=1.10,
        t_half_hr=2.20,
        cmax_ng_ml=110.0,
        tmax_hr=0.15,
        auc_inf_ng_hr_ml=340.0,
        ppb_percent=97.0,
        cl_int_hlm_ml_min_kg=60.0,
        caco2_log_cm_s=-4.5,
        fraction_excreted_renal=0.01,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 7. Propranolol (Base, High PPB, High VDss, High Extraction Hepatic)
    ObservedPKSemantic(
        compound_name="Propranolol",
        smiles="CC(C)NCC(O)COc1cccc2ccccc12",
        dose_mg=40.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Routledge & Shand Clin Pharmacokinet / DailyMed",
        cl_systemic_ml_min_kg=12.0,  # ~840 mL/min
        cl_hepatic_ml_min_kg=11.8,
        cl_oral_apparent_ml_min_kg=48.0,
        vdss_l_kg=4.30,
        f_oral_bioavailability=0.25, # High first-pass metabolism
        t_half_hr=4.00,
        cmax_ng_ml=65.0,
        tmax_hr=1.5,
        auc_inf_ng_hr_ml=350.0,
        ppb_percent=90.0,
        cl_int_hlm_ml_min_kg=75.0,
        caco2_log_cm_s=-4.6,
        fraction_excreted_renal=0.01,
        extraction_category="HIGH",
        ionization_class="BASE",
    ),
    # 8. Theophylline (Neutral/Weak Acid, Low PPB, Low VDss, CYP1A2 Low Extraction)
    ObservedPKSemantic(
        compound_name="Theophylline",
        smiles="Cn1c(=O)c2[nH]cnc2n(C)c1=O",
        dose_mg=300.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Solution",
        fasted=True,
        population="Healthy Adult",
        source_reference="Hendeles et al. Clin Pharmacokinet / DailyMed",
        cl_systemic_ml_min_kg=0.70,  # ~50 mL/min
        cl_hepatic_ml_min_kg=0.62,
        cl_renal_ml_min_kg=0.08,
        cl_oral_apparent_ml_min_kg=0.72,
        vdss_l_kg=0.50,
        f_oral_bioavailability=0.98,
        t_half_hr=8.00,
        cmax_ng_ml=7500.0,
        tmax_hr=1.5,
        auc_inf_ng_hr_ml=102000.0,
        ppb_percent=40.0,
        cl_int_hlm_ml_min_kg=12.0,
        caco2_log_cm_s=-4.8,
        fraction_excreted_renal=0.10,
        extraction_category="LOW",
        ionization_class="NEUTRAL",
    ),
    # 9. Atenolol (Base, Low PPB, Low VDss, 50% Renal Elimination)
    ObservedPKSemantic(
        compound_name="Atenolol",
        smiles="CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
        dose_mg=50.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Mason et al. Clin Pharmacokinet / DailyMed Tenormin",
        cl_systemic_ml_min_kg=1.90,  # ~130 mL/min
        cl_hepatic_ml_min_kg=0.85,
        cl_renal_ml_min_kg=0.95,
        cl_oral_apparent_ml_min_kg=3.80,
        vdss_l_kg=0.70,
        f_oral_bioavailability=0.50,
        t_half_hr=6.50,
        cmax_ng_ml=320.0,
        tmax_hr=3.0,
        auc_inf_ng_hr_ml=3100.0,
        ppb_percent=8.0,
        cl_int_hlm_ml_min_kg=5.0,
        caco2_log_cm_s=-5.9,
        fraction_excreted_renal=0.50,
        extraction_category="LOW",
        ionization_class="BASE",
    ),
    # 10. Digoxin (Neutral Glycoside, Moderate PPB, Very High VDss, 70% Renal Secretion)
    ObservedPKSemantic(
        compound_name="Digoxin",
        smiles="CC1OC(OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CCC5(C)C(C7=CC(=O)OC7)CCC65O)C4)OC3C)OC2C)CC(O)C1O",
        dose_mg=0.50,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Injection",
        fasted=True,
        population="Healthy Adult",
        source_reference="Iisalo Clin Pharmacokinet / DailyMed Lanoxin",
        cl_systemic_ml_min_kg=2.50,  # ~175 mL/min
        cl_hepatic_ml_min_kg=0.75,
        cl_renal_ml_min_kg=1.75,     # Active P-gp / glomerular filtration
        vdss_l_kg=7.00,
        t_half_hr=40.0,
        cmax_ng_ml=6.5,
        tmax_hr=0.5,
        auc_inf_ng_hr_ml=48.0,
        ppb_percent=25.0,
        cl_int_hlm_ml_min_kg=4.0,
        caco2_log_cm_s=-5.8,
        fraction_excreted_renal=0.70,
        extraction_category="LOW",
        ionization_class="NEUTRAL",
    ),
    # 11. Diazepam (Neutral/Weak Base, High PPB 98.5%, Moderate VDss, Long t1/2, CYP2C19/3A4)
    ObservedPKSemantic(
        compound_name="Diazepam",
        smiles="CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
        dose_mg=10.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Klotz et al. J Clin Invest / DailyMed Valium",
        cl_systemic_ml_min_kg=0.40,  # ~27 mL/min
        cl_hepatic_ml_min_kg=0.39,
        cl_oral_apparent_ml_min_kg=0.40,
        vdss_l_kg=1.10,
        f_oral_bioavailability=1.00,
        t_half_hr=48.0,
        cmax_ng_ml=250.0,
        tmax_hr=1.0,
        auc_inf_ng_hr_ml=5900.0,
        ppb_percent=98.5,
        cl_int_hlm_ml_min_kg=10.0,
        caco2_log_cm_s=-4.6,
        fraction_excreted_renal=0.01,
        extraction_category="LOW",
        ionization_class="NEUTRAL",
    ),
    # 12. Ciprofloxacin (Zwitterion, Moderate VDss, Mixed Hepatic/Renal)
    ObservedPKSemantic(
        compound_name="Ciprofloxacin",
        smiles="O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        dose_mg=500.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Bergan et al. / DailyMed Cipro",
        cl_systemic_ml_min_kg=5.80,  # ~400 mL/min
        cl_hepatic_ml_min_kg=2.30,
        cl_renal_ml_min_kg=3.50,
        cl_oral_apparent_ml_min_kg=8.28,
        vdss_l_kg=2.50,
        f_oral_bioavailability=0.70,
        t_half_hr=4.00,
        cmax_ng_ml=2400.0,
        tmax_hr=1.5,
        auc_inf_ng_hr_ml=14400.0,
        ppb_percent=30.0,
        cl_int_hlm_ml_min_kg=18.0,
        caco2_log_cm_s=-5.5,
        fraction_excreted_renal=0.60,
        extraction_category="LOW",
        ionization_class="AMPHOLYTE",
    ),
    # 13. Sildenafil (Base, High PPB 96%, Moderate VDss, CYP3A4/2C9 Hepatic)
    ObservedPKSemantic(
        compound_name="Sildenafil",
        smiles="CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
        dose_mg=50.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Muirhead et al. Br J Clin Pharmacol / DailyMed Viagra",
        cl_systemic_ml_min_kg=9.80,  # ~680 mL/min
        cl_hepatic_ml_min_kg=9.40,
        cl_oral_apparent_ml_min_kg=24.0,
        vdss_l_kg=1.50,
        f_oral_bioavailability=0.41,
        t_half_hr=4.00,
        cmax_ng_ml=280.0,
        tmax_hr=1.0,
        auc_inf_ng_hr_ml=1050.0,
        ppb_percent=96.0,
        cl_int_hlm_ml_min_kg=45.0,
        caco2_log_cm_s=-4.7,
        fraction_excreted_renal=0.03,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 14. Atorvastatin (Acid, High PPB 98%, High VDss, CYP3A4/OATP1B1)
    ObservedPKSemantic(
        compound_name="Atorvastatin",
        smiles="CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        dose_mg=40.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Cilla et al. Clin Pharmacokinet / DailyMed Lipitor",
        cl_systemic_ml_min_kg=8.50,  # ~600 mL/min
        cl_hepatic_ml_min_kg=8.40,
        cl_oral_apparent_ml_min_kg=60.0,
        vdss_l_kg=5.40,
        f_oral_bioavailability=0.14, # High gut/liver first-pass
        t_half_hr=14.0,
        cmax_ng_ml=28.0,
        tmax_hr=1.5,
        auc_inf_ng_hr_ml=160.0,
        ppb_percent=98.0,
        cl_int_hlm_ml_min_kg=50.0,
        caco2_log_cm_s=-4.9,
        fraction_excreted_renal=0.02,
        extraction_category="INTERMEDIATE",
        ionization_class="ACID",
    ),
    # 15. Rosuvastatin (Acid, High PPB 88%, Moderate VDss, BCRP/OATP Transporter Minimal CYP)
    ObservedPKSemantic(
        compound_name="Rosuvastatin",
        smiles="CC(C)c1nc(N(C)S(=O)(=O)C)nc(-c2ccc(F)cc2)c1C=CC(O)CC(O)CC(=O)O",
        dose_mg=20.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Martin et al. Clin Pharmacokinet / DailyMed Crestor",
        cl_systemic_ml_min_kg=6.00,  # ~420 mL/min
        cl_hepatic_ml_min_kg=4.20,
        cl_renal_ml_min_kg=1.80,
        cl_oral_apparent_ml_min_kg=30.0,
        vdss_l_kg=1.90,
        f_oral_bioavailability=0.20,
        t_half_hr=19.0,
        cmax_ng_ml=15.0,
        tmax_hr=4.0,
        auc_inf_ng_hr_ml=130.0,
        ppb_percent=88.0,
        cl_int_hlm_ml_min_kg=12.0,
        caco2_log_cm_s=-5.7,
        fraction_excreted_renal=0.28,
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 16. Verapamil (Base, High PPB 90%, High VDss, High Extraction Hepatic)
    ObservedPKSemantic(
        compound_name="Verapamil",
        smiles="COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
        dose_mg=5.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Injection",
        fasted=True,
        population="Healthy Adult",
        source_reference="Eichelbaum et al. Eur J Clin Pharmacol / DailyMed",
        cl_systemic_ml_min_kg=13.0,  # ~900 mL/min
        cl_hepatic_ml_min_kg=12.5,
        cl_renal_ml_min_kg=0.5,
        vdss_l_kg=5.00,
        t_half_hr=4.00,
        cmax_ng_ml=120.0,
        tmax_hr=0.10,
        auc_inf_ng_hr_ml=90.0,
        ppb_percent=90.0,
        cl_int_hlm_ml_min_kg=90.0,
        caco2_log_cm_s=-4.4,
        fraction_excreted_renal=0.03,
        extraction_category="HIGH",
        ionization_class="BASE",
    ),
    # 17. Furosemide (Acid, Extreme PPB 98%, Low VDss, 65% Renal Excretion)
    ObservedPKSemantic(
        compound_name="Furosemide",
        smiles="NS(=O)(=O)c1cc(c(Cl)cc1C(=O)O)NCc1ccco1",
        dose_mg=40.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Injection",
        fasted=True,
        population="Healthy Adult",
        source_reference="Hammarlund-Udenaes Clin Pharmacokinet / DailyMed Lasix",
        cl_systemic_ml_min_kg=2.00,  # ~140 mL/min
        cl_hepatic_ml_min_kg=0.70,
        cl_renal_ml_min_kg=1.30,
        vdss_l_kg=0.15,
        t_half_hr=1.50,
        cmax_ng_ml=4500.0,
        tmax_hr=0.15,
        auc_inf_ng_hr_ml=4800.0,
        ppb_percent=98.0,
        cl_int_hlm_ml_min_kg=6.0,
        caco2_log_cm_s=-5.5,
        fraction_excreted_renal=0.65,
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 18. Ibuprofen (Acid, Extreme PPB 99%, Low VDss, CYP2C9 Hepatic)
    ObservedPKSemantic(
        compound_name="Ibuprofen",
        smiles="CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        dose_mg=400.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Davies Clin Pharmacokinet / DailyMed Motrin",
        cl_systemic_ml_min_kg=0.75,  # ~52 mL/min
        cl_hepatic_ml_min_kg=0.72,
        cl_oral_apparent_ml_min_kg=0.75,
        vdss_l_kg=0.15,
        f_oral_bioavailability=1.00,
        t_half_hr=2.00,
        cmax_ng_ml=30000.0,
        tmax_hr=1.5,
        auc_inf_ng_hr_ml=125000.0,
        ppb_percent=99.0,
        cl_int_hlm_ml_min_kg=22.0,
        caco2_log_cm_s=-4.5,
        fraction_excreted_renal=0.01,
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 19. Omeprazole (Weak Base, High PPB 95%, Low VDss, Ultra-short t1/2, CYP2C19)
    ObservedPKSemantic(
        compound_name="Omeprazole",
        smiles="COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
        dose_mg=20.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Delayed Release",
        fasted=True,
        population="Healthy Adult",
        source_reference="Regardh et al. Scand J Gastroenterol / DailyMed Prilosec",
        cl_systemic_ml_min_kg=7.00,  # ~500 mL/min
        cl_hepatic_ml_min_kg=6.90,
        cl_oral_apparent_ml_min_kg=17.5,
        vdss_l_kg=0.35,
        f_oral_bioavailability=0.40,
        t_half_hr=1.00,
        cmax_ng_ml=650.0,
        tmax_hr=2.0,
        auc_inf_ng_hr_ml=1500.0,
        ppb_percent=95.0,
        cl_int_hlm_ml_min_kg=40.0,
        caco2_log_cm_s=-4.8,
        fraction_excreted_renal=0.01,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 20. Dextromethorphan (Base, Moderate PPB 65%, High VDss, CYP2D6 Probe)
    ObservedPKSemantic(
        compound_name="Dextromethorphan",
        smiles="COc1ccc2c(c1)C13CCC(CC2)N(C)C1CC3",
        dose_mg=30.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Syrup",
        fasted=True,
        population="Healthy Adult (EM)",
        source_reference="Capon et al. Clin Pharmacokinet / DailyMed",
        cl_systemic_ml_min_kg=15.0,  # High first-pass clearance in EM
        cl_hepatic_ml_min_kg=14.5,
        cl_oral_apparent_ml_min_kg=150.0,
        vdss_l_kg=5.50,
        f_oral_bioavailability=0.10,
        t_half_hr=2.50,
        cmax_ng_ml=2.5,
        tmax_hr=2.5,
        auc_inf_ng_hr_ml=18.0,
        ppb_percent=65.0,
        cl_int_hlm_ml_min_kg=80.0,
        caco2_log_cm_s=-4.5,
        fraction_excreted_renal=0.01,
        extraction_category="HIGH",
        ionization_class="BASE",
    ),
    # 21. Fluoxetine (Base, High PPB 94.5%, Very High VDss 25 L/kg, Long t1/2, CYP2D6)
    ObservedPKSemantic(
        compound_name="Fluoxetine",
        smiles="CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
        dose_mg=20.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Capsule",
        fasted=True,
        population="Healthy Adult",
        source_reference="Altamura et al. Clin Pharmacokinet / DailyMed Prozac",
        cl_systemic_ml_min_kg=6.00,  # ~420 mL/min
        cl_hepatic_ml_min_kg=5.80,
        cl_oral_apparent_ml_min_kg=8.57,
        vdss_l_kg=25.0,
        f_oral_bioavailability=0.70,
        t_half_hr=48.0,
        cmax_ng_ml=20.0,
        tmax_hr=6.0,
        auc_inf_ng_hr_ml=1100.0,
        ppb_percent=94.5,
        cl_int_hlm_ml_min_kg=25.0,
        caco2_log_cm_s=-4.6,
        fraction_excreted_renal=0.02,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 22. Haloperidol (Base, High PPB 92%, Very High VDss 18 L/kg, Long t1/2, CYP3A4/2D6)
    ObservedPKSemantic(
        compound_name="Haloperidol",
        smiles="OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1",
        dose_mg=5.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Froemming et al. Psychopharmacology / DailyMed Haldol",
        cl_systemic_ml_min_kg=10.0,  # ~700 mL/min
        cl_hepatic_ml_min_kg=9.8,
        cl_oral_apparent_ml_min_kg=16.7,
        vdss_l_kg=18.0,
        f_oral_bioavailability=0.60,
        t_half_hr=20.0,
        cmax_ng_ml=3.2,
        tmax_hr=4.0,
        auc_inf_ng_hr_ml=75.0,
        ppb_percent=92.0,
        cl_int_hlm_ml_min_kg=55.0,
        caco2_log_cm_s=-4.8,
        fraction_excreted_renal=0.01,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 23. Ketoconazole (Base, Extreme PPB 99%, Moderate VDss, CYP3A4 Substrate & Inhibitor)
    ObservedPKSemantic(
        compound_name="Ketoconazole",
        smiles="CC(=O)N1CCN(c2ccc(OCC3COC(Cn4cncn4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1",
        dose_mg=200.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Daneshmend & Warnock Clin Pharmacokinet / DailyMed Nizoral",
        cl_systemic_ml_min_kg=3.50,  # ~250 mL/min
        cl_hepatic_ml_min_kg=3.40,
        cl_oral_apparent_ml_min_kg=4.67,
        vdss_l_kg=1.20,
        f_oral_bioavailability=0.75,
        t_half_hr=3.30,
        cmax_ng_ml=3500.0,
        tmax_hr=2.0,
        auc_inf_ng_hr_ml=18000.0,
        ppb_percent=99.0,
        cl_int_hlm_ml_min_kg=35.0,
        caco2_log_cm_s=-4.9,
        fraction_excreted_renal=0.02,
        extraction_category="LOW",
        ionization_class="BASE",
    ),
    # 24. Clarithromycin (Base, Moderate PPB 70%, Moderate VDss, Mixed Hepatic/Renal 30%)
    ObservedPKSemantic(
        compound_name="Clarithromycin",
        smiles="CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)OC",
        dose_mg=250.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Healthy Adult",
        source_reference="Chu et al. Antimicrob Agents Chemother / DailyMed Biaxin",
        cl_systemic_ml_min_kg=8.50,  # ~600 mL/min
        cl_hepatic_ml_min_kg=5.95,
        cl_renal_ml_min_kg=2.55,
        cl_oral_apparent_ml_min_kg=17.0,
        vdss_l_kg=3.00,
        f_oral_bioavailability=0.50,
        t_half_hr=5.00,
        cmax_ng_ml=1200.0,
        tmax_hr=2.0,
        auc_inf_ng_hr_ml=6500.0,
        ppb_percent=70.0,
        cl_int_hlm_ml_min_kg=28.0,
        caco2_log_cm_s=-5.2,
        fraction_excreted_renal=0.30,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 25. Dabigatran Etexilate (Prodrug / Acid Dabigatran, Low PPB 35%, 80% RENAL Elimination)
    ObservedPKSemantic(
        compound_name="Dabigatran",
        smiles="CN(CC(=O)N(c1ccccc1)c1ccc(C(=N)N)cc1)c1nc2ccc(C(=O)O)cc2n1C",
        dose_mg=150.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Capsule",
        fasted=True,
        population="Healthy Adult",
        source_reference="Stangier et al. Clin Pharmacokinet / DailyMed Pradaxa",
        cl_systemic_ml_min_kg=1.60,  # ~110 mL/min
        cl_hepatic_ml_min_kg=0.30,
        cl_renal_ml_min_kg=1.30,
        cl_oral_apparent_ml_min_kg=24.6,
        vdss_l_kg=1.00,
        f_oral_bioavailability=0.065, # Very low oral F
        t_half_hr=13.0,
        cmax_ng_ml=110.0,
        tmax_hr=2.0,
        auc_inf_ng_hr_ml=1050.0,
        ppb_percent=35.0,
        cl_int_hlm_ml_min_kg=3.0,
        caco2_log_cm_s=-6.2,
        fraction_excreted_renal=0.80, # PREDOMINANTLY RENAL
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 26. Methotrexate (Acid, Moderate PPB 50%, Low VDss, 80%+ RENAL Elimination)
    ObservedPKSemantic(
        compound_name="Methotrexate",
        smiles="CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1",
        dose_mg=15.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Injection",
        fasted=True,
        population="Healthy Adult",
        source_reference="Bleyer Cancer Treat Rev / DailyMed Methotrexate",
        cl_systemic_ml_min_kg=1.80,  # ~125 mL/min
        cl_hepatic_ml_min_kg=0.30,
        cl_renal_ml_min_kg=1.50,     # Active OAT3 / glomerular
        vdss_l_kg=0.50,
        t_half_hr=6.00,
        cmax_ng_ml=1200.0,
        tmax_hr=0.25,
        auc_inf_ng_hr_ml=3200.0,
        ppb_percent=50.0,
        cl_int_hlm_ml_min_kg=2.5,
        caco2_log_cm_s=-6.0,
        fraction_excreted_renal=0.85, # PREDOMINANTLY RENAL
        extraction_category="LOW",
        ionization_class="ACID",
    ),
    # 27. Morphine (Base, Low PPB 35%, Moderate VDss, High UGT Glucuronidation Extraction)
    ObservedPKSemantic(
        compound_name="Morphine",
        smiles="CN1CCC23C4Oc5c(O)ccc(CC1C2C=CC4O)c53",
        dose_mg=10.0,
        route="IV",
        regimen="SINGLE_DOSE",
        formulation="IV Injection",
        fasted=True,
        population="Healthy Adult",
        source_reference="Hasselstrom & Sawe Clin Pharmacokinet / DailyMed",
        cl_systemic_ml_min_kg=17.0,  # ~1200 mL/min
        cl_hepatic_ml_min_kg=15.0,   # UGT2B7 hepatic clearance
        cl_renal_ml_min_kg=2.0,
        vdss_l_kg=3.50,
        t_half_hr=2.00,
        cmax_ng_ml=150.0,
        tmax_hr=0.10,
        auc_inf_ng_hr_ml=140.0,
        ppb_percent=35.0,
        cl_int_hlm_ml_min_kg=110.0,
        caco2_log_cm_s=-5.1,
        fraction_excreted_renal=0.10,
        extraction_category="HIGH",
        ionization_class="BASE",
    ),
    # 28. Gefitinib (Base, High PPB 90%, High VDss 20 L/kg, Long t1/2, CYP3A4)
    ObservedPKSemantic(
        compound_name="Gefitinib",
        smiles="COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        dose_mg=250.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="Cancer Patient",
        source_reference="FDA NDA 021399 Review / DailyMed Iressa",
        cl_systemic_ml_min_kg=8.50,  # ~595 mL/min
        cl_hepatic_ml_min_kg=8.20,
        cl_oral_apparent_ml_min_kg=14.2,
        vdss_l_kg=20.0,
        f_oral_bioavailability=0.60,
        t_half_hr=48.0,
        cmax_ng_ml=502.0,
        tmax_hr=4.0,
        auc_inf_ng_hr_ml=13500.0,
        ppb_percent=90.0,
        cl_int_hlm_ml_min_kg=40.0,
        caco2_log_cm_s=-4.9,
        fraction_excreted_renal=0.04,
        extraction_category="INTERMEDIATE",
        ionization_class="BASE",
    ),
    # 29. Imatinib (Base, High PPB 95%, High VDss 4.8 L/kg, CYP3A4 Hepatic)
    ObservedPKSemantic(
        compound_name="Imatinib",
        smiles="Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        dose_mg=400.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="CML Patient",
        source_reference="Peng et al. Clin Pharmacokinet / DailyMed Gleevec",
        cl_systemic_ml_min_kg=2.00,  # ~140 mL/min
        cl_hepatic_ml_min_kg=1.85,
        cl_oral_apparent_ml_min_kg=2.04,
        vdss_l_kg=4.80,
        f_oral_bioavailability=0.98,
        t_half_hr=18.0,
        cmax_ng_ml=2600.0,
        tmax_hr=3.0,
        auc_inf_ng_hr_ml=38000.0,
        ppb_percent=95.0,
        cl_int_hlm_ml_min_kg=20.0,
        caco2_log_cm_s=-5.0,
        fraction_excreted_renal=0.05,
        extraction_category="LOW",
        ionization_class="BASE",
    ),
    # 30. Pazopanib (Base, Extreme PPB 99.9%, Very Low VDss 0.14 L/kg, Long t1/2, CYP3A4)
    ObservedPKSemantic(
        compound_name="Pazopanib",
        smiles="Cc1ccc(Nc2nc(N(C)c3ccc4c(C)n(C)nc4c3)ncc2)cc1",
        dose_mg=800.0,
        route="ORAL",
        regimen="SINGLE_DOSE",
        formulation="Oral Tablet",
        fasted=True,
        population="RCC Patient",
        source_reference="FDA NDA 022465 Review / DailyMed Votrient",
        cl_systemic_ml_min_kg=0.20,  # ~14 mL/min
        cl_hepatic_ml_min_kg=0.19,
        cl_oral_apparent_ml_min_kg=0.95,
        vdss_l_kg=0.14,
        f_oral_bioavailability=0.21,
        t_half_hr=31.0,
        cmax_ng_ml=58000.0,
        tmax_hr=3.5,
        auc_inf_ng_hr_ml=850000.0,
        ppb_percent=99.9, # Extreme protein binding
        cl_int_hlm_ml_min_kg=15.0,
        caco2_log_cm_s=-5.1,
        fraction_excreted_renal=0.04,
        extraction_category="LOW",
        ionization_class="BASE",
    ),
]


def classify_elimination_status(fe_renal: float) -> str:
    """Classifies whether hepatic IVIVE represents total systemic clearance.

    Implements Directive 16:
    If renal clearance is significant (fe_renal >= 0.20), hepatic IVIVE
    is marked as HEPATIC_CL_ESTIMATE, and total clearance is NOT fabricated.
    """
    if fe_renal >= 0.20:
        return "HEPATIC_CL_ESTIMATE"
    return "TOTAL_CL_PREDICTION"


def evaluate_oral_confidence(fa: Optional[float], ka: Optional[float], caco2_ad: str = "IN_DOMAIN") -> str:
    """Determines confidence in oral PK predictions (Directive 17).

    Oral prediction requires F and ka. If weakly constrained or OOD,
    downgrades to PK_CONFIDENCE_LOW or PK_CONFIDENCE_MODERATE.
    """
    if fa is None or ka is None:
        return "PK_CONFIDENCE_LOW"
    if caco2_ad != "IN_DOMAIN":
        return "PK_CONFIDENCE_LOW"
    if fa < 0.30 or ka < 0.20 or ka > 3.0:
        return "PK_CONFIDENCE_MODERATE"
    return "PK_CONFIDENCE_HIGH"


def calculate_fold_error(predicted: float, observed: float) -> float:
    """Calculates symmetric fold error between predicted and observed values.

    fold_error = max(pred/obs, obs/pred) >= 1.0
    """
    if predicted <= 0.0 or observed <= 0.0:
        raise ValueError(f"Fold error requires positive values, got pred={predicted}, obs={observed}")
    ratio = predicted / observed
    return max(ratio, 1.0 / ratio)
