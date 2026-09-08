"""Keep the v3.3.3 campaign restartable after release validation."""
import json
from datetime import datetime, timezone
from pathlib import Path

p = Path("validation/v333_completion_campaign.json")
s = json.loads(p.read_text())
s["status"] = "IN_PROGRESS"
s["phase"] = "NUMERICAL_OPTIMIZATION"
s["current_task"] = "ppb_strategy_ladder_v80"
s["exact_next_action"] = "Evaluate train-fold similarity-weighted kNN log-fu and descriptor-space local regression on the frozen V72 benchmark; inverse-transform before scoring."
s["endpoint_campaign"]["dynamic_queue"] = ["PPB_KNN_LOGFU_V80", "PPB_DESCRIPTOR_LOCAL_V80", "HEPATIC_CL_DEPENDENCY_REEVALUATION", "VDSS_HYBRID_REVIEW", "DOWNSTREAM_PK_CONTEXT_AUDIT"]
s["endpoint_campaign"]["next_strategy"] = "Leakage-safe similarity-weighted kNN/local regression candidates (V80); do not rerun rejected V74-V76 families."
s["updated_at"] = datetime.now(timezone.utc).isoformat()
p.write_text(json.dumps(s, indent=2) + "\n")
print("queued PPB V80 continuation")
