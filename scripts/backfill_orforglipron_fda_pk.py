#!/usr/bin/env python3
"""Accept the bounded set of source-qualified Orforglipron FDA PK facts.

This does not promote scraped candidate rows wholesale.  Each row below has
an explicit endpoint and study context, and is idempotently linked to the
current Orforglipron compound version.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import select

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.database import SessionLocal
from backend.models import Compound, ExternalExperimentalEvidence

FDA_LABEL = "https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/220934s003lbl.pdf"
FDA_REVIEW = "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2026/220934Orig1s000MultidisciplineR.pdf"
ROWS = (
    {
        "endpoint": "HUMAN_PK_F_ORAL", "raw_endpoint": "absolute oral bioavailability",
        "value": "77", "unit": "%", "route": "ORAL", "dose": 0.8,
        "regimen": "SINGLE_DOSE", "source_id": "FDA-FOUNDAYO-LABEL-12.3-F-0.8MG",
        "source_url": FDA_LABEL,
        "reference": "FOUNDAYO US Prescribing Information (2026), Clinical Pharmacology 12.3: geometric mean absolute bioavailability 77% after 0.8 mg.",
        "conditions": "Human parent-drug absolute oral bioavailability; 0.8 mg oral dose; product tablet; fasted/fed status not stated.",
    },
    {
        "endpoint": "HUMAN_PK_VD_IV", "raw_endpoint": "steady-state volume of distribution",
        "value": "285", "unit": "L", "route": "IV", "dose": None,
        "regimen": "UNSPECIFIED", "source_id": "FDA-FOUNDAYO-LABEL-12.3-VDSS",
        "source_url": FDA_LABEL,
        "reference": "FOUNDAYO US Prescribing Information (2026), Clinical Pharmacology 12.3: mean steady-state volume of distribution approximately 285 L following intravenous dosing in healthy subjects.",
        "conditions": "Human parent-drug Vdss following intravenous dosing in healthy subjects; dose and regimen not stated.",
    },
    {
        "endpoint": "HUMAN_PK_CL_UNSPECIFIED", "raw_endpoint": "systemic clearance",
        "value": "7.15", "unit": "L/h", "route": "UNSPECIFIED", "dose": None,
        "regimen": "UNSPECIFIED", "source_id": "FDA-FOUNDAYO-LABEL-12.3-CL",
        "source_url": FDA_LABEL,
        "reference": "FOUNDAYO US Prescribing Information (2026), Clinical Pharmacology 12.3: mean systemic clearance 7.15 L/hour.",
        "conditions": "Human parent-drug systemic clearance. The source does not establish an IV or oral CL/F measurement context; route intentionally retained as unspecified.",
    },
    {
        "endpoint": "HUMAN_PK_CMAX_UNSPECIFIED", "raw_endpoint": "Cmax",
        "value": "149", "unit": "ng/mL", "route": "UNSPECIFIED", "dose": 36.0,
        "regimen": "MULTIPLE_DOSE", "source_id": "FDA-NDA220934-36MG-CMAXSS",
        "source_url": FDA_REVIEW,
        "reference": "FDA NDA 220934 multidisciplinary review: 149 ng/mL Cmax at 36 mg human exposure; source does not supply sufficient formulation/route details for numeric pairing.",
        "conditions": "Human parent-drug Cmax context: 36 mg, multiple-dose exposure. Formulation and route not asserted beyond the source text; display-only, not numeric-pairable with an unrelated scenario.",
    },
)


def main() -> int:
    with SessionLocal() as db:
        compound = db.scalar(select(Compound).where(Compound.name == "Orforglipron"))
        if not compound:
            raise RuntimeError("Orforglipron compound not found")
        version = next((v for v in compound.versions if v.version_number == compound.current_version), None)
        if not version:
            raise RuntimeError("Orforglipron current version not found")
        accepted = []
        for source in ROWS:
            key = hashlib.sha256((str(version.id) + "|" + source["source_id"]).encode()).hexdigest()
            row = db.scalar(select(ExternalExperimentalEvidence).where(ExternalExperimentalEvidence.provenance_key == key))
            context = {
                "conditions": source["conditions"], "species": "HUMAN", "route": source["route"],
                "dose": source["dose"], "dose_unit": "mg" if source["dose"] is not None else "",
                "regimen": source["regimen"], "analyte": "PARENT", "matrix": "PLASMA",
                "formulation": "UNSPECIFIED", "pk_context_version": "drugopt-pk-context-v5.4-curated",
                "context_completeness": "PARTIAL" if source["route"] == "UNSPECIFIED" else "QUALIFIED",
            }
            values = dict(compound_version_id=version.id, provenance_key=key, cas_number=compound.cas_number or "",
                raw_endpoint_name=source["raw_endpoint"], raw_value=source["value"], raw_relation="=", raw_unit=source["unit"],
                assay_type="CLINICAL_PK", assay_conditions_json=context, species="HUMAN", source_database="FDA / Regulatory",
                source_record_id=source["source_id"], source_document_id="NDA220934", reference_text=source["reference"],
                source_url=source["source_url"], identity_match_status="EXACT_STRUCTURE_MATCH", endpoint_match_status="CANONICAL_MATCH",
                mapping_status="CANONICAL_PK_CURATED", evidence_origin="EXPERIMENTAL_EXTERNAL", canonical_endpoint_id=source["endpoint"],
                normalized_value=source["value"], normalized_unit=source["unit"], normalization_rule="FDA source unit retained", normalization_version="drugopt-unit-normalization-v5.4",
                comparability_status="DIRECTLY_COMPARABLE" if source["route"] != "UNSPECIFIED" else "CONTEXT_REQUIRED",
                source_quality_class="A", duplicate_status="DISTINCT_MEASUREMENT", evidence_state="EXTERNAL_IMPORTED",
                accepted_at=datetime.now(timezone.utc), qualification_status="QUALIFIED", routing_section="PK", routing_reason="FDA primary clinical PK source",
                lifecycle_status="ACTIVE")
            if row:
                for key_name, value in values.items():
                    setattr(row, key_name, value)
            else:
                row = ExternalExperimentalEvidence(**values)
                db.add(row)
            accepted.append({"endpoint": source["endpoint"], "value": source["value"], "unit": source["unit"], "source": source["source_id"]})
        db.commit()
    artifact = {"artifact": "orforglipron_fda_pk_backfill_v1", "created_at": datetime.now(timezone.utc).isoformat(),
                "compound": "Orforglipron", "alias": "LY3502970", "accepted_observations": accepted,
                "root_cause": "Qualified FDA PK evidence existed as auto-qualified candidates and was omitted from the accepted-evidence PK display stream."}
    Path("validation/orforglipron_pk_backfill_v1.json").write_text(json.dumps(artifact, indent=2) + "\n")
    print(json.dumps({"accepted": len(accepted)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
