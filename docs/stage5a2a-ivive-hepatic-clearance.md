# Stage 5A-2A — IVIVE Hepatic Clearance Foundation

## Scope boundary

This stage estimates hepatic clearance only. It does not predict renal or other
non-hepatic clearance, total systemic clearance, Vd, absorption, absolute oral
bioavailability, Cmax/AUC, concentration-time profiles, PBPK, efficacy, or PK/PD.

## Data architecture and isolation

- `ivive_input_sets`: immutable, project- and CompoundVersion-scoped IVIVE inputs.
- `physiological_parameter_sets`: species/parameter/value/unit/reference/version.
- `physiological_parameter_overrides`: immutable project-scoped, study-specific overrides.
- `pk_ivive_method_registry`: PK/IVIVE equations, assumptions, version, and references;
  separate from the ADMET AI model registry.
- `ivive_runs`: immutable input/physiology snapshot, equations, parameter-set version,
  outputs, assumptions, warnings, confidence ceiling, timestamp, and SHA-256 input hash.

All input discovery filters by the exact `CompoundVersion`, project, and species.
Existing ADMET measurements and predictions remain unchanged.

## Input selection policy

1. Experimental hepatocyte Clint
2. Experimental microsomal Clint
3. Project-calibrated validated quantitative Clint
4. External quantitative prediction
5. No calculation

Experimental values precede predictions. Qualitative or classification-only results
are excluded and are never converted into a quantitative intrinsic clearance.

Plasma binding uses experimental PPB/fu,p before predicted PPB. A PPB measurement is
normalized only when its unit explicitly denotes percent/fraction bound or unbound.

## Unit and scaling engine

Canonical whole-body clearance unit: `mL/min/kg`.

Supported raw intrinsic-clearance units:

- microsomes: `µL/min/mg protein`, `mL/min/mg`
- hepatocytes: `µL/min/10^6 cells`, `mL/min/10^6 cells`

Supported scaled units:

- `mL/min/kg`
- `L/h/kg`
- `log10(mL/min/kg)` for a declared pre-scaled Clint

Microsomal scaling:

`Clint (mL/min/kg) = raw Clint (µL/min/mg) × MPPGL (mg/g) × liver weight (g/kg) / 1000`

Hepatocyte scaling:

`Clint (mL/min/kg) = raw Clint (µL/min/10^6 cells) × hepatocellularity (10^6 cells/g) × liver weight (g/kg) / 1000`

`RAW_MICROSOMAL`, `RAW_HEPATOCYTE`, and `PRESCALED_CLINT` are hard types. A
type/unit mismatch fails. A pre-scaled OpenADMET output is linearized as
`10^(log10 Clint)` and has `scaling_count = 0`; MPPGL and liver weight are not applied.

## Species physiology v1.0

The coherent five-species defaults reproduce Table 2 of Ring et al. (2011). Hepatic
blood flow was converted from whole-animal flow using the table's standard body
weight; relative liver weight was converted from percent body weight to g/kg.

| Species | Qh (mL/min/kg) | Liver (g/kg) | MPPGL (mg/g) | Hepatocellularity (10^6/g) | Default BW (kg) |
|---|---:|---:|---:|---:|---:|
| Mouse | 120 | 54.9 | 47 | 128 | 0.02 |
| Rat | 67.6 | 36.6 | 47 | 128 | 0.25 |
| Dog | 30.9 | 32.9 | 58 | 187.5 | 10 |
| Monkey | 43.6 | 24.8 | 32 | 99 | 5 |
| Human | 20.714 | 25.7 | 32 | 99 | 70 |

Parameter-set version: `PHRMA-CPCDC-2011-v1.0`.

References retained on every row:

- Ring BJ et al. J Pharm Sci. 2011;100:4090–4110. DOI `10.1002/jps.22552`.
- Davies B, Morris T. Pharm Res. 1993;10:1093–1095. DOI `10.1023/A:1018943613122`.
- Brown RP et al. Toxicol Ind Health. 1997;13:407–484. DOI `10.1177/074823379701300401`.
- Barter ZE et al. Curr Drug Metab. 2007;8:33–45. DOI `10.2174/138920007779315053`.

Overrides are project-scoped immutable records. The canonical normalized value and
the raw value/unit, source, notes, confidence, timestamp, and default parameter-set
version are all retained. UI labels defaults and overrides as `DEFAULT PHYSIOLOGY`
and `USER OVERRIDE`.

## Blood and plasma binding

`fu,p` is the unbound fraction in plasma. For an experimental blood/plasma ratio:

`fu,b = fu,p / (B/P)`

The result must be in `(0, 1]`. If B/P is missing, the compound-specific value is not
invented. The current policy permits the explicit plasma-based approximation
`fu,b ≈ fu,p`, adds a warning, labels the binding basis `PLASMA_APPROXIMATION`, and
caps result confidence at MEDIUM. If PPB/fu,p is missing, calculation is unavailable.

## Well-stirred hepatic model

`CLh = Qh × fu,b × Clint / (Qh + fu,b × Clint)`

`Eh = CLh / Qh`

`Fh = 1 − Eh`

Extraction classification is Low for `Eh < 0.3`, Intermediate for `0.3 ≤ Eh ≤ 0.7`,
and High for `Eh > 0.7`. Threshold provenance is retained from the public clearance
scaling literature (e.g. Lombardo et al., Pharm Res. 2014;31:3499–3511,
PMCID `PMC4181675`). `Fh` is always displayed as **Predicted Hepatic Availability
(Fh)** and is not absolute oral bioavailability F.

The method registry is extensible by method key. Only `WELL_STIRRED` is active in
this stage; parallel-tube and dispersion methods have not been implemented.

## Experimental PK comparison

The latest valid, species-matched IV NCA `CL` for the exact CompoundVersion is shown
as **Observed IV systemic CL** beside **Predicted hepatic CL**. The difference, fold
error reference, and `CLh / observed CL` ratio are descriptive only. The ratio is
called **Estimated hepatic contribution**, and values over 100% generate a warning.
No ADMET endpoint feedback or model weighting is updated from this comparison.

## Confidence and missing-data safety

The run confidence is the minimum confidence across selected Clint, PPB, B/P when
present, and required physiology. It can never exceed its weakest required input.
A missing-B/P approximation caps confidence at MEDIUM. An unavailable calculation
has `NOT_AVAILABLE` confidence and still stores an auditable run snapshot.

The output always states:

- `Hepatic Clearance Estimate`
- `Non-hepatic Clearance: Not modeled`
- `Predicted Total Clearance: null / not generated`

## Public and synthetic numerical validation

Machine-readable cases are in `validation/stage5a2a_ivive_validation.json`.

- Li et al. (2025 public HepaRG IVIVE): diclofenac hepatocyte Clint 261
  µL/min/10^6 cells scales to about 664 mL/min/kg with the consensus human scalar.
  With fu,p 0.005 and no B/P it gives 2.87 mL/min/kg, close to the published 2.91.
  The small difference is expected because the paper includes `fu,inc = 0.963`, which
  this stage does not silently infer.
- Obach (1999 public 29-drug HLM table): prescaled Clint examples 2.1 mL/min/kg
  (diphenhydramine) and 189 mL/min/kg (diclofenac) exercise low/high unbound
  extraction boundaries without re-scaling.
- Exact synthetic cases verify low and high extraction, unit conversions, MPPGL and
  hepatocellularity scaling, log-linearization, and double-scaling prevention.

The diclofenac literature calculation-reproduction fold difference is about 1.01.
This is not an observed systemic-CL validation metric. For project data, the engine
computes fold error and can aggregate average absolute fold error, percent within
2-fold, and percent within 3-fold, while preserving the limitation that observed IV
CL is total systemic clearance and predicted CLh is hepatic only.

Published large-dataset context further documents systematic IVIVE underprediction:
Bowman et al. (2020, PMCID `PMC7325626`) reported 42.2% of human microsome and 30.7%
of human hepatocyte predictions within twofold in the Wood et al. datasets. These
literature results are context, not claimed as validation performance of this local
implementation.

