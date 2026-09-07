"""Create the authoritative endpoint completion board without changing models."""
from __future__ import annotations
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend.prediction_engine_registry import CURRENT_ENGINE_ID  # noqa: E402

PROD = CURRENT_ENGINE_ID


def row(endpoint, semantics, model, candidate, dev, locked, independent, metric, current, target, ad, independence, status, blocker, action, experiment=None):
    return {"endpoint": endpoint, "canonical_semantics": semantics, "production_model": model,
            "candidate_model": candidate, "development_N": dev, "locked_validation_N": locked,
            "independent_validation_N": independent, "metric": metric, "current_value": current,
            "target_value": target, "AD_status": ad, "source_independence": independence,
            "status": status, "blocker": blocker, "last_experiment": experiment or blocker, "next_action": action}


def main():
    rows = [
        row("HUMAN_FU_PLASMA", "human plasma fraction unbound", "existing PPB route", "Extra Trees + RDKit descriptors", 232, 0, 0, "AAFE", 1.317, "AAFE <=2.0; within2 >=70%", "CV only", "DrugBank/reference only", "EXTERNAL_ACCESS_BLOCKED", "v6.2 PPB qualification", "Acquire independent human plasma fu cohort (>=50)."),
        row("HUMAN_TOTAL_IV_CL", "human IV total systemic clearance", "none", "Random Forest + RDKit descriptors", 748, 7, 0, "AAFE", 1.700, "AAFE <=2.5; within2 >=50%", "same-source only", "PKSmart only", "EXTERNAL_ACCESS_BLOCKED", "v6.2 total-IV-CL qualification", "Acquire independent human IV total-CL cohort (>=30)."),
        row("HUMAN_VDSS", "human IV true Vss/Vdss", "mechanistic consensus", "none", 50, 0, 0, "AAFE", None, "AAFE <=2.5; within2 >=50%", "strict N=50", "DrugBank/reference only", "EXTERNAL_ACCESS_BLOCKED", "v6.2 strict VDss audit", "Acquire independent record-level true Vss data (>=100)."),
        row("PKA_ACID", "experimental acidic macro/site pKa", "ionization_smarts_rules_v1", "none", 0, 0, 0, "MAE pKa", None, "MAE <=0.7", "not established", "no qualified experimental set", "EXTERNAL_ACCESS_BLOCKED", "v6.4 pKa closure", "Acquire licensed experimental acid pKa dataset."),
        row("PKA_BASE", "experimental basic macro/site pKa", "ionization_smarts_rules_v1", "none", 0, 0, 0, "MAE pKa", None, "MAE <=0.7", "not established", "no qualified experimental set", "EXTERNAL_ACCESS_BLOCKED", "v6.4 pKa closure", "Acquire licensed experimental base pKa dataset."),
        row("LOGD_7_4", "experimental logD at explicit pH 7.4", "henderson_hasselbalch_logd_v1", "none", 0, 0, 0, "MAE logD", None, "MAE <=0.7", "not established", "license unresolved", "EXTERNAL_ACCESS_BLOCKED", "v6.4 logD closure", "Resolve license and obtain pH-specific dataset."),
        row("HUMAN_HEPATOCYTE_CLINT", "human hepatocyte intrinsic clearance", "none", "none", 0, 0, 0, "AAFE", None, "AAFE <=2.5", "none", "no qualified checkpoint", "EXTERNAL_ACCESS_BLOCKED", "v6.2 evidence closure", "Acquire context-rich human hepatocyte dataset/checkpoint."),
        row("HUMAN_RENAL_CL", "human renal clearance", "routing/readiness only", "none", 0, 0, 0, "AAFE", None, "AAFE <=2.5", "none", "insufficient quantitative evidence", "EXTERNAL_ACCESS_BLOCKED", "v5.3 renal readiness", "Acquire CLr/fe/urinary recovery with IV context."),
        row("HUMAN_HEPATIC_CL", "human hepatic clearance", "HLM IVIVE", "none", 30, 0, 0, "AAFE", 3.295, "AAFE <=2.5", "known limitations", "assisted/limited", "EXTERNAL_ACCESS_BLOCKED", "v5.6 prediction-only chain; no independent compatible cohort", "Acquire independent human hepatic/nonrenal clearance with prediction-only fu/Clint context."),
        row("QUANT_CYP2C19", "quantitative human CYP2C19 inhibition/substrate endpoint", "none", "none", 0, 0, 0, "endpoint-specific", None, "independent quantitative validation", "none", "no qualified model", "EXTERNAL_ACCESS_BLOCKED", "v5.8 quantitative transporter/CYP audit", "Acquire exact assay-compatible quantitative data."),
        row("QUANT_PGP", "quantitative human P-gp endpoint", "none", "none", 0, 0, 0, "endpoint-specific", None, "independent quantitative validation", "none", "no qualified model", "EXTERNAL_ACCESS_BLOCKED", "v5.8 quantitative transporter/CYP audit", "Acquire dose-response quantitative P-gp data."),
        row("QUANT_BCRP", "quantitative human BCRP endpoint", "none", "none", 0, 0, 0, "endpoint-specific", None, "independent quantitative validation", "none", "no qualified model", "EXTERNAL_ACCESS_BLOCKED", "v5.8 quantitative transporter/CYP audit", "Acquire dose-response quantitative BCRP data."),
        row("HUMAN_F", "context-aware absolute bioavailability", "mechanistic context route", "none", 0, 0, 0, "absolute F error", None, "MAE <=0.15", "none", "heterogeneous clinical context", "EXTERNAL_ACCESS_BLOCKED", "v5.4 evidence closure", "Acquire paired IV/PO studies with dose/formulation."),
        row("KA", "context-aware absorption rate constant", "none", "none", 0, 0, 0, "endpoint-specific", None, "independent continuous validation", "none", "no homogeneous dataset", "EXTERNAL_ACCESS_BLOCKED", "v5.2 PK context resolver", "Acquire formulation/time-course absorption data."),
        row("IV_HALF_LIFE", "human IV half-life from predicted total CL + VDss", "assisted PK equations", "none", 0, 7, 0, "AAFE", None, "AAFE <=2.5; within2 >=50%", "depends on CL/VDss", "upstream independent gaps", "EXTERNAL_ACCESS_BLOCKED", "v5.6 IV half-life research; CL/VDss external qualification required", "Acquire independent human IV CL + true VDss + half-life records."),
        row("ORAL_AUC", "human oral AUC with explicit F/route/context", "frozen oral PK route", "none", 0, 0, 0, "AAFE", 3.282, "AAFE <=2.5; within2 >=50%", "upstream dependent", "assisted baseline", "EXTERNAL_ACCESS_BLOCKED", "v5.3 PK benchmark; independent CL/F/ka/Vd evidence absent", "Acquire paired IV/PO clinical PK with explicit F, ka/context, CL and Vd semantics."),
        row("ORAL_CMAX", "human oral Cmax with explicit F/ka/route/context", "frozen oral PK route", "none", 0, 0, 0, "AAFE", 2.540, "AAFE <=2.5; within2 >=50%", "upstream dependent", "assisted baseline", "EXTERNAL_ACCESS_BLOCKED", "v5.3 PK benchmark; independent F/ka/context evidence absent", "Acquire paired oral PK time-course records with formulation, F and ka-compatible context."),
        row("SOLUBILITY", "qualified aqueous solubility", "production v3.3.2 route", "none", 0, 0, 0, "production validation", "stable", "no regression", "established", "qualified production lineage", "DONE", "v3.3.2 production", "Regression anchor; preserve."),
        row("CACO2_PAPP_AB", "Caco-2 A to B Papp", "production v3.3.2 route", "none", 0, 0, 0, "production validation", "stable", "no regression", "established", "qualified production lineage", "DONE", "v3.3.2 production", "Regression anchor; preserve."),
        row("HLM_CLINT", "human liver microsomal intrinsic clearance", "production v3.3.2 route", "none", 179, 0, 0, "production validation", "stable", "no regression", "established", "qualified production lineage", "DONE", "v3.3.2 production", "Keep semantics distinct from hepatic/total CL."),
        row("CYP3A4", "human CYP3A4 assay-defined endpoint", "production v3.3.2 route", "none", 0, 0, 0, "production validation", "stable", "no regression", "established", "qualified production lineage", "DONE", "v3.3.2 production", "Regression anchor; preserve."),
        row("HERG", "human hERG assay-defined liability", "production v3.3.2 route", "none", 0, 0, 0, "production validation", "stable", "no regression", "established", "qualified production lineage", "DONE", "v3.3.2 production", "Regression anchor; preserve."),
    ]
    artifact = {"artifact": "ENDPOINT_COMPLETION_BOARD_V70", "created_at": str(date.today()), "production_engine_version": PROD, "reference_library_n": 1000, "rows": rows, "status_policy": ["NOT_STARTED", "DATA_BUILDING", "MODEL_DEVELOPMENT", "VALIDATION", "IMPROVING", "PRODUCTION_REVIEW", "DONE", "EXTERNAL_ACCESS_BLOCKED"], "note": "EXTERNAL_ACCESS_BLOCKED is used only where completion requires licensed/credentialed data unavailable in this environment; it is not a claim that the endpoint is scientifically complete."}
    (ROOT / "validation/endpoint_completion_board.json").write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"production_engine_version": PROD, "endpoints": len(rows), "done": sum(x["status"] == "DONE" for x in rows), "external_access_blocked": sum(x["status"] == "EXTERNAL_ACCESS_BLOCKED" for x in rows)}, indent=2))


if __name__ == "__main__": main()
