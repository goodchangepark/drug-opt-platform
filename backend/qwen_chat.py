"""Qwen3.5 9B Local Section & Comparison Chat Service for Drug-OPT."""

from __future__ import annotations
import json
import logging
import urllib.request
import urllib.error
from typing import Any, List, Dict, Optional
from sqlalchemy.orm import Session
from sqlalchemy import select

from backend.models import (
    Compound, CompoundVersion, Project, ExternalExperimentalEvidence
)
from backend.activity_models import ActivityMeasurement, ActivityPrediction
from backend.admet import ADMETMeasurement, ADMETPrediction
from backend.pk import PKStudy, PKNCAResult
from backend.ivive import IVIVERun
from backend.metabolism import MetabolicSoftSpot

logger = logging.getLogger("drugopt.qwen_chat")

OLLAMA_URL = "http://localhost:11434"
QWEN_MODEL = "qwen3.5:9b"

SYSTEM_PROMPT = """당신은 Drug-OPT에서 저분자 신약개발과 medicinal chemistry를 지원하는 전문 AI assistant이다.

기본 규칙:
1. 사용자가 다른 언어를 명시하지 않는 한 모든 답변은 한국어로 작성한다.
2. 제공된 Drug-OPT의 구조, 실험값, prediction 결과를 최우선 근거로 사용한다.
3. 데이터에 없는 숫자, 실험 결과, 구조적 특징을 만들어내지 않는다.
4. 정보가 없으면 "현재 데이터에서는 확인할 수 없습니다."라고 명확히 표현한다.
5. 답변은 핵심 위주로 간결하게 작성하며 최대 3가지 주요 제안으로 제한한다.

Medicinal chemistry 규칙:
6. CYP inhibition, metabolic stability, permeability, BBB/CNS penetration,
   PPB, P-gp, hERG, potency를 서로 구분하고 trade-off를 함께 고려한다.
7. CYP inhibition과 metabolic stability를 동일한 현상으로 취급하지 않는다.
8. metabolic soft-spot blocking이 CYP inhibition 감소를 보장한다고 표현하지 않는다.
9. 구조 또는 SMILES가 제공되지 않은 경우 C-2, para-position 등 특정 위치나
   특정 substituent가 실제 존재한다고 가정하거나 만들어내지 않는다.
10. 구조 정보가 부족하면 "해당 기능기가 존재한다면"과 같은 조건부 표현을 사용한다.
11. 실제 구조가 제공되면 구조에 존재하는 functional group을 중심으로 제안한다.
12. SMILES만으로 특정 원자 위치를 확신하기 어려운 경우 위치를 임의로 단정하지 않는다.
13. P-gp, hERG, BBB 효과는 예측값이나 구조적 근거가 없으면 단정하지 않는다.
14. CNS 최적화에서는 cLogP/logD, TPSA, HBD, pKa, P-gp 및 unbound fraction의
    균형을 고려한다.
15. fluorination, heteroaryl replacement, pKa 조절, bioisostere replacement,
    soft-spot blocking 등의 전략은 적용 가능한 구조적 근거가 있을 때만 제안한다.
16. 근거가 제한적이면 사실처럼 단정하지 말고 가능성 또는 가설이라고 표현한다.

비교 질문에서는 각 물질의 장점, liability 및 중요한 차이를 명확하게 설명한다.
답변은 원칙적으로 3~6개의 짧은 문장 또는 3개의 간단한 항목 이내로 작성한다.
"""

def _call_qwen(prompt: str, timeout: int = 75) -> str:
    """Send prompt to local Ollama Qwen3.5 9B using native chat API."""
    payload = {
        "model": QWEN_MODEL,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ],
        "stream": False,
        "think": False,
        "keep_alive": "30m",
        "options": {
            "num_ctx": 4096,
            "num_predict": 300,
            "temperature": 0.2,
        },
    }

    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        f"{OLLAMA_URL}/api/chat",
        data=data,
        headers={"Content-Type": "application/json"},
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = json.loads(resp.read().decode("utf-8"))
            message = body.get("message") or {}
            answer = message.get("content", "").strip()
            return answer or "현재 데이터에서는 확인할 수 없습니다."
    except urllib.error.URLError as exc:
        logger.error(f"Ollama connection error: {exc}")
        return f"로컬 Qwen3.5 9B 서비스 연결에 실패했습니다: {exc}"
    except Exception as exc:
        logger.error(f"Qwen generation error: {exc}")
        return f"답변 생성 중 오류가 발생했습니다: {exc}"

def build_compound_section_context(
    db: Session,
    compound: Compound,
    version: Optional[CompoundVersion],
    section: str,
    workspace_data: Optional[Dict[str, Any]] = None,
) -> str:
    """Build compact, structured text context for a specific compound and section."""
    lines = [
        f"Compound: {compound.name} (CAS: {compound.cas_number or 'N/A'})",
        f"Section: {section.upper()}",
    ]
    if not version:
        lines.append("Status: Draft compound (No calculated version yet)")
        return "\n".join(lines)

    lines.append(f"Structure SMILES: {version.canonical_smiles or version.original_smiles}")
    props = version.properties_json or {}
    calc = version.calculation_json or {}

    sec = section.lower()

    if sec in ("overview", "properties"):
        lines.append("Physicochemical Properties (RDKit Calculated):")
        if props.get("molecular_weight"):
            lines.append(f"- Molecular Weight (MW): {props['molecular_weight']:.1f} g/mol")
        if props.get("clogp") is not None:
            lines.append(f"- cLogP: {props['clogp']:.2f}")
        if props.get("tpsa") is not None:
            lines.append(f"- TPSA: {props['tpsa']:.1f} Å²")
        if props.get("hbd") is not None:
            lines.append(f"- HBD: {props['hbd']}, HBA: {props.get('hba', 'N/A')}")
        if props.get("rotatable_bonds") is not None:
            lines.append(f"- Rotatable Bonds: {props['rotatable_bonds']}")
        if props.get("fraction_csp3") is not None:
            lines.append(f"- Fraction Csp3 (Fsp3): {props['fraction_csp3']:.2f}")
        if props.get("qed") is not None:
            lines.append(f"- QED Drug-Likeness: {props['qed']:.2f}")
        rules = calc.get("rules", {})
        if rules:
            rule_summary = ", ".join([f"{k}: {v.get('result')}" for k, v in rules.items()])
            lines.append(f"- Drug-likeness Rules: {rule_summary}")
        alerts = version.alerts_json or []
        if alerts:
            alert_names = ", ".join([a.get("alert_name", "") for a in alerts])
            lines.append(f"- Structural Alerts: {alert_names}")
        else:
            lines.append("- Structural Alerts: None detected")

    elif sec == "activity":
        lines.append("Activity Measurements & Predictions:")
        measurements = list(db.scalars(
            select(ActivityMeasurement).where(ActivityMeasurement.version_id == version.id)
        ))
        if measurements:
            for m in measurements:
                lines.append(f"- [Exp] Assay: {m.assay_name} | Type: {m.measurement_type} | Value: {m.value} {m.unit or ''} | Target: {m.target or 'N/A'}")
        else:
            lines.append("- [Exp] No internal activity measurements recorded.")

        predictions = list(db.scalars(
            select(ActivityPrediction).where(ActivityPrediction.version_id == version.id)
        ))
        if predictions:
            for p in predictions:
                lines.append(f"- [Pred] Model: {p.model_name} | Predicted Value: {p.predicted_value:.2f} {p.unit or ''} (Confidence: {p.confidence or 'N/A'})")

        # Qualified external activity evidence
        ext_evidence = list(db.scalars(
            select(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id == version.id
            )
        ))
        act_ev = [row for row in ext_evidence if "IC50" in (row.canonical_endpoint_id or row.raw_endpoint_name or "") or "EC50" in (row.canonical_endpoint_id or row.raw_endpoint_name or "") or "Ki" in (row.canonical_endpoint_id or row.raw_endpoint_name or "")]
        if act_ev:
            lines.append("Qualified External Activity Evidence:")
            for e in act_ev[:8]:
                lines.append(f"- [External Exp] Endpoint: {e.canonical_endpoint_id or e.raw_endpoint_name} | Value: {e.normalized_value or e.raw_value} {e.normalized_unit or e.raw_unit or ''} | Source: {e.source_database}")

    elif sec == "admet":
        lines.append("ADMET Endpoint Profile (Experimental & Predictions):")
        measurements = list(db.scalars(
            select(ADMETMeasurement).where(ADMETMeasurement.version_id == version.id)
        ))
        if measurements:
            for m in measurements:
                lines.append(f"- [Exp] Endpoint: {m.endpoint} | Value: {m.value} {m.unit or ''}")
        predictions = list(db.scalars(
            select(ADMETPrediction).where(ADMETPrediction.version_id == version.id)
        ))
        if predictions:
            for p in predictions:
                lines.append(f"- [Pred] Endpoint: {p.endpoint} | Predicted: {p.predicted_value} {p.unit or ''} | Assessment: {p.assessment or 'N/A'}")

        # Qualified external ADMET evidence
        ext_evidence = list(db.scalars(
            select(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id == version.id
            )
        ))
        admet_ev = [row for row in ext_evidence if any(k in (row.canonical_endpoint_id or row.raw_endpoint_name or "").upper() for k in ("CACO", "SOLUBILITY", "PPB", "HERG", "AMES", "DILI", "BBB", "LOGD"))]
        if admet_ev:
            lines.append("Qualified External ADMET Evidence:")
            for e in admet_ev[:8]:
                lines.append(f"- [External Exp] Endpoint: {e.canonical_endpoint_id or e.raw_endpoint_name} | Value: {e.normalized_value or e.raw_value} {e.normalized_unit or e.raw_unit or ''} | Source: {e.source_database}")

    elif sec == "metabolism":
        lines.append("Metabolism & Clearance Profile:")
        meta_preds = list(db.scalars(
            select(ADMETPrediction).where(
                ADMETPrediction.version_id == version.id,
                ADMETPrediction.endpoint.in_(["HLM", "RLM", "MLM", "CYP3A4 Inh", "CYP2D6 Inh", "CYP2C9 Inh", "P-gp Inh"])
            )
        ))
        if meta_preds:
            for p in meta_preds:
                lines.append(f"- [Pred] {p.endpoint}: {p.predicted_value} {p.unit or ''} ({p.assessment or 'N/A'})")

        soft_spots = list(db.scalars(
            select(MetabolicSoftSpot).where(MetabolicSoftSpot.version_id == version.id)
        ))
        if soft_spots:
            spots_txt = ", ".join([f"Atom {s.atom_index} ({s.cyp_isoform or 'CYP'}: {s.vulnerability_score})" for s in soft_spots[:5]])
            lines.append(f"- Soft spots: {spots_txt}")

        ext_evidence = list(db.scalars(
            select(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id == version.id
            )
        ))
        met_ev = [row for row in ext_evidence if any(k in (row.canonical_endpoint_id or row.raw_endpoint_name or "").upper() for k in ("CYP", "HLM", "RLM", "MLM", "P-GP", "BCRP", "METABOL"))]
        if met_ev:
            lines.append("Qualified Metabolism Evidence:")
            for e in met_ev[:8]:
                lines.append(f"- [External Exp] Endpoint: {e.canonical_endpoint_id or e.raw_endpoint_name} | Value: {e.normalized_value or e.raw_value} {e.normalized_unit or e.raw_unit or ''} | Source: {e.source_database}")

    elif sec == "pk":
        lines.append("Pharmacokinetics & Translational Profile:")
        studies = list(db.scalars(
            select(PKStudy).where(PKStudy.version_id == version.id)
        ))
        if studies:
            for st in studies:
                nca = db.scalars(select(PKNCAResult).where(PKNCAResult.pk_study_id == st.id)).first()
                if nca:
                    lines.append(f"- [Study: {st.species} {st.route} {st.dose or ''} {st.dose_unit or ''}] CL: {nca.cl or nca.cl_obs} {nca.cl_unit or 'L/h'}, Vd: {nca.vz or nca.vz_obs} {nca.vz_unit or 'L'}, t1/2: {nca.terminal_half_life or 'N/A'} h, AUC: {nca.auclast or 'N/A'}, Cmax: {nca.cmax or 'N/A'}")
        else:
            lines.append("- No preclinical in vivo study rows.")

        ivive_runs = list(db.scalars(
            select(IVIVERun).where(IVIVERun.version_id == version.id)
        ))
        if ivive_runs:
            for iv in ivive_runs:
                calc_cl = (iv.outputs_json or {}).get("calculated_clearance") or iv.outputs_json or {}
                lines.append(f"- [IVIVE {iv.species or 'Human'}] Hepatic CL: {calc_cl.get('cl_blood', 'N/A')} mL/min/kg | Method: {iv.parameter_set_version or 'IVIVE'}")

        ext_evidence = list(db.scalars(
            select(ExternalExperimentalEvidence).where(
                ExternalExperimentalEvidence.compound_version_id == version.id
            )
        ))
        pk_ev = [row for row in ext_evidence if any(k in (row.canonical_endpoint_id or row.raw_endpoint_name or "").upper() for k in ("CLEARANCE", "HALF-LIFE", "HALF_LIFE", "AUC", "CMAX", "VSS", "VD", "F_", "BIOAVAILABILITY")) or "PK" in (row.raw_endpoint_name or "").upper()]
        if pk_ev:
            lines.append("Qualified Clinical / External PK Evidence:")
            for e in pk_ev[:8]:
                lines.append(f"- [External PK] {e.canonical_endpoint_id or e.raw_endpoint_name}: {e.normalized_value or e.raw_value} {e.normalized_unit or e.raw_unit or ''} (Context: {e.species or 'Clinical'})")

    return "\n".join(lines)


def build_comparison_context(
    db: Session,
    project: Project,
    compounds: List[Compound],
    comparison_data: Optional[Dict[str, Any]] = None,
) -> str:
    """Build structured comparison context for selected compounds."""
    lines = [
        f"Project: {project.name} (Target: {project.target or 'N/A'})",
        f"Comparing {len(compounds)} Compounds: {', '.join([c.name for c in compounds])}",
        "Structured Comparison Matrix:",
    ]

    for c in compounds:
        lines.append(f"\n--- Compound: {c.name} ---")
        v = next((v for v in c.versions if v.version_number == c.current_version), c.versions[-1] if c.versions else None)
        if not v:
            lines.append("Status: Draft compound")
            continue
        props = v.properties_json or {}
        lines.append(f"MW: {props.get('molecular_weight', 'N/A')} | cLogP: {props.get('clogp', 'N/A')} | TPSA: {props.get('tpsa', 'N/A')} | QED: {props.get('qed', 'N/A')}")

        # ADMET key endpoints
        preds = list(db.scalars(select(ADMETPrediction).where(ADMETPrediction.version_id == v.id)))
        if preds:
            pred_map = {p.endpoint: f"{p.predicted_value} {p.unit or ''}" for p in preds}
            lines.append(f"ADMET Predictions: {', '.join([f'{k}={val}' for k, val in list(pred_map.items())[:6]])}")

        # PK summary
        studies = list(db.scalars(select(PKStudy).where(PKStudy.version_id == v.id)))
        if studies:
            pk_notes = []
            for st in studies:
                nca = db.scalars(select(PKNCAResult).where(PKNCAResult.pk_study_id == st.id)).first()
                if nca and nca.cl:
                    pk_notes.append(f"{st.species} CL={nca.cl} {nca.cl_unit or ''}, t1/2={nca.terminal_half_life or 'N/A'}h")
            if pk_notes:
                lines.append(f"PK Studies: {'; '.join(pk_notes)}")

    return "\n".join(lines)


def answer_section_question(
    db: Session,
    compound_id: int,
    section: str,
    question: str,
    workspace_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Process a question regarding a single compound in a specific section."""
    compound = db.scalar(select(Compound).where(Compound.id == compound_id))
    if not compound:
        return {"error": "Compound not found", "answer": "해당 물질을 찾을 수 없습니다."}

    version = next((v for v in compound.versions if v.version_number == compound.current_version), compound.versions[-1] if compound.versions else None)
    context_text = build_compound_section_context(db, compound, version, section, workspace_data)

    user_prompt = f"""[Current Structured Context]
{context_text}

Question: {question.strip()}"""
    answer = _call_qwen(user_prompt)
    return {
        "model": QWEN_MODEL,
        "compound_name": compound.name,
        "section": section,
        "answer": answer,
    }


def answer_comparison_question(
    db: Session,
    project_id: int,
    compound_ids: List[int],
    question: str,
    comparison_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Process a comparison question across selected compounds."""
    project = db.scalar(select(Project).where(Project.id == project_id))
    if not project:
        return {"error": "Project not found", "answer": "해당 프로젝트를 찾을 수 없습니다."}

    compounds = list(db.scalars(select(Compound).where(Compound.id.in_(compound_ids), Compound.project_id == project_id)))
    if not compounds:
        return {"error": "No compounds selected", "answer": "비교할 물질이 선택되지 않았습니다."}

    context_text = build_comparison_context(db, project, compounds, comparison_data)

    user_prompt = f"""[Selected Compounds Comparison Context]
{context_text}

Question: {question.strip()}"""
    answer = _call_qwen(user_prompt)
    return {
        "model": QWEN_MODEL,
        "project_name": project.name,
        "compounds": [c.name for c in compounds],
        "answer": answer,
    }
