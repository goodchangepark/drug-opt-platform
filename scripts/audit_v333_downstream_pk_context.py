"""Persist honest downstream PK readiness states for v3.3.3."""
import json
from datetime import datetime, timezone
from pathlib import Path

rows = []
for endpoint in ("IV_HALF_LIFE", "ORAL_F", "ORAL_KA", "ORAL_AUC", "ORAL_CMAX", "ORAL_TMAX"):
    rows.append({
        "endpoint": endpoint,
        "status": "CONTEXT_REQUIRED" if endpoint != "IV_HALF_LIFE" else "MODEL_UNAVAILABLE",
        "prediction_mode": "FULL_PREDICTION",
        "experimental_validation_N": 0,
        "reason": "Requires species/context-qualified upstream CL/Vd/F/ka and frozen OOF values; no generic structure-only value is emitted.",
    })
out = {
    "artifact": "V333_DOWNSTREAM_PK_CONTEXT_AUDIT",
    "generated_at": datetime.now(timezone.utc).isoformat(),
    "engine_id": "drugopt-prediction-engine-v3@3.3.3",
    "rows": rows,
    "safety": "No experimental upstream target values were propagated into predictions.",
}
Path("validation/v333_downstream_pk_context_audit.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
print("wrote validation/v333_downstream_pk_context_audit.json")
