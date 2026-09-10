"""Fail-closed model artifact authority for Stable Core current predictions.

The engine routing registry describes *which* endpoint route is approved.  This
module additionally binds a selectable current value to files that actually
exist in this checkout.  The stored snapshot hash is a deterministic SHA256 of
the complete registered file bundle; a changed or missing component therefore
invalidates current selection instead of silently inheriting endpoint metadata.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from .prediction_engine_registry import (
    ROUTE_MODEL_UNAVAILABLE,
    get_current_production_routing,
)

ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ModelArtifactRegistration:
    endpoint_id: str
    model_id: str
    model_version: str
    model_route: str
    artifact_paths: tuple[str, ...]
    validation_artifact: str

    def resolved_paths(self) -> tuple[Path, ...]:
        return tuple((ROOT / path).resolve() for path in self.artifact_paths)


# Only routes with an inspectable installed implementation are admitted.  A
# route omitted here remains visible in registry metadata but cannot publish a
# CurrentPredictionSnapshot until its exact executable bundle is registered.
_ARTIFACT_PATHS: dict[str, tuple[str, ...]] = {
    "SOLUBILITY_GENERIC": ("models/admetica/solubility/model_v2_1.pt",),
    "CACO2_PAPP_AB": ("models/admetica/caco2/model_v2_1.pt",),
    "HUMAN_PPB": ("models/admetica/ppbr/model_v2_1.pt",),
    "HLM_CLINT": ("models/openadmet/microsomal_clearance/model.pth",),
    "RLM_CLINT": ("models/openadmet/microsomal_clearance/model.pth",),
    "MLM_CLINT": ("models/openadmet/microsomal_clearance/model.pth",),
    "MW": ("backend/chemistry.py",),
    "CLOGP": ("backend/chemistry.py",),
    "TPSA": ("backend/chemistry.py",),
    "HBD": ("backend/chemistry.py",),
    "HBA": ("backend/chemistry.py",),
    "ROTB": ("backend/chemistry.py",),
    "FSP3": ("backend/chemistry.py",),
    "QED": ("backend/chemistry.py",),
    "FORMAL_CHARGE": ("backend/chemistry.py",),
    "HEAVY_ATOM_COUNT": ("backend/chemistry.py",),
    "METABOLIC_SOFT_SPOTS": ("models/sygma/phase1.txt", "models/sygma/phase2.txt"),
    "METABOLITE_HYPOTHESES": ("models/sygma/phase1.txt", "models/sygma/phase2.txt"),
    "PKA": ("backend/ionization.py",),
    "LOGD_7_4": ("backend/ionization.py",),
    "CYP1A2_INHIBITOR_CLASS": ("models/admetica/cyp/cyp1a2-inhibitor/model_v2_1.pt",),
    "CYP2C9_INHIBITOR_CLASS": ("models/admetica/cyp/cyp2c9-inhibitor/model_v2_1.pt",),
    "CYP2C19_INHIBITOR_CLASS": ("models/admetica/cyp/cyp2c19-inhibitor/model_v2_1.pt",),
    "CYP2D6_INHIBITOR_CLASS": ("models/admetica/cyp/cyp2d6-inhibitor/model_v2_1.pt",),
    "CYP3A4_INHIBITOR_CLASS": ("models/admetica/cyp/cyp3a4-inhibitor/model_v2_1.pt",),
    "CYP2C9_SUBSTRATE": ("models/admetica/cyp/cyp2c9-substrate/model_v2_1.pt",),
    "CYP2D6_SUBSTRATE": ("models/admetica/cyp/cyp2d6-substrate/model_v2_1.pt",),
    "CYP3A4_SUBSTRATE": ("models/admetica/cyp/cyp3a4-substrate/model_v2_1.pt",),
    "PGP_INHIBITION": ("models/admetica/transporter/pgp-inhibitor/model_v2_1.pt",),
    "HERG_CLASS": ("models/admetica/safety/herg/model_v2_1.pt",),
    "AMES_MUTAGENICITY": tuple(f"models/admet_ai/classification/model_{index}.pt" for index in range(5)),
    "DILI_LIABILITY": tuple(f"models/admet_ai/classification/model_{index}.pt" for index in range(5)),
}

VALIDATION_ARTIFACT = "validation/model_maturity_ui_audit.json"


def _route_for(endpoint_id: str) -> dict | None:
    endpoint = str(endpoint_id or "").strip().upper()
    return next(
        (row for row in get_current_production_routing() if row.get("endpoint_id") == endpoint),
        None,
    )


def model_artifact_registration(endpoint_id: str) -> ModelArtifactRegistration | None:
    route = _route_for(endpoint_id)
    endpoint = str(endpoint_id or "").strip().upper()
    paths = _ARTIFACT_PATHS.get(endpoint)
    if not route or route.get("route") == ROUTE_MODEL_UNAVAILABLE or not paths:
        return None
    return ModelArtifactRegistration(
        endpoint_id=endpoint,
        model_id=str(route.get("model_or_ensemble") or "").strip(),
        model_version=str(route.get("model_version_hash") or "").strip(),
        model_route=str(route.get("route") or "").strip(),
        artifact_paths=paths,
        validation_artifact=VALIDATION_ARTIFACT,
    )


def artifact_bundle_sha256(registration: ModelArtifactRegistration) -> str | None:
    digest = hashlib.sha256()
    for relative, path in zip(registration.artifact_paths, registration.resolved_paths()):
        try:
            path.relative_to(ROOT)
        except ValueError:
            return None
        if not path.is_file():
            return None
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def validation_artifact_exists(registration: ModelArtifactRegistration) -> bool:
    path = (ROOT / registration.validation_artifact).resolve()
    try:
        path.relative_to(ROOT)
    except ValueError:
        return False
    return path.is_file()


def validation_registration_error(registration: ModelArtifactRegistration) -> str | None:
    """Require an endpoint validation record for the exact routed model class."""
    path = (ROOT / registration.validation_artifact).resolve()
    if not path.is_file():
        return "VALIDATION_ARTIFACT_MISSING"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "VALIDATION_ARTIFACT_INVALID"
    record = next(
        (row for row in payload.get("rows", []) if row.get("endpoint") == registration.endpoint_id),
        None,
    )
    if record is None:
        return "VALIDATION_RECORD_MISSING"
    if str(record.get("current_model") or "") != registration.model_route:
        return "VALIDATION_MODEL_ROUTE_MISMATCH"
    if str(record.get("registry_status") or "").upper() in {"", "MODEL_UNAVAILABLE", "UNAVAILABLE"}:
        return "VALIDATION_STATUS_UNQUALIFIED"
    return None


def resolve_model_artifact(
    endpoint_id: str,
    model_id: str,
    model_version: str,
    artifact_hash: str,
) -> tuple[ModelArtifactRegistration | None, str | None]:
    registration = model_artifact_registration(endpoint_id)
    if registration is None:
        return None, "MODEL_NOT_REGISTERED"
    if model_id != registration.model_id:
        return None, "MODEL_ID_REGISTRY_MISMATCH"
    if model_version != registration.model_version:
        return None, "MODEL_VERSION_REGISTRY_MISMATCH"
    validation_error = validation_registration_error(registration)
    if validation_error:
        return None, validation_error
    actual_hash = artifact_bundle_sha256(registration)
    if actual_hash is None:
        return None, "MODEL_ARTIFACT_MISSING"
    if artifact_hash != actual_hash:
        return None, "MODEL_ARTIFACT_HASH_MISMATCH"
    return registration, None
