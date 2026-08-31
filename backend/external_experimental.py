"""CAS-gated public experimental-evidence lookup.

Only a CAS is sent to public services.  Returned public structures are
standardized and compared locally; no private SMILES are transmitted.
"""
from __future__ import annotations

import json
import re
import time
from datetime import datetime, timezone
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from rdkit import Chem

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
_CACHE: dict[str, tuple[float, dict]] = {}
_TTL_SECONDS = 900


def cas_status(value: str | None) -> str:
    value = (value or "").strip()
    if not value:
        return "EMPTY"
    if not CAS_RE.fullmatch(value):
        return "INVALID_FORMAT"
    digits = value.replace("-", "")
    check = sum(int(digit) * (index + 1) for index, digit in enumerate(reversed(digits[:-1]))) % 10
    return "VALID" if check == int(digits[-1]) else "INVALID_CHECKSUM"


def valid_cas(value: str | None) -> bool:
    return cas_status(value) == "VALID"


def _get_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "Drug-OPT/1.0 public-evidence lookup (CAS-only)"})
    try:
        with urlopen(request, timeout=12) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        if exc.code == 404:
            return {"_not_found": True}
        if exc.code == 429:
            return {"_error": "SOURCE_RATE_LIMITED"}
        return {"_error": f"SOURCE_HTTP_{exc.code}"}
    except (URLError, TimeoutError, ValueError):
        return {"_error": "SOURCE_TEMPORARILY_UNAVAILABLE"}


def _identity_status(public_smiles: str, local_inchikey: str) -> dict:
    molecule = Chem.MolFromSmiles(public_smiles or "")
    if not molecule:
        return {"status": "STRUCTURE_MISMATCH", "public_inchikey": ""}
    public_inchikey = Chem.MolToInchiKey(molecule)
    if public_inchikey == local_inchikey:
        return {"status": "EXACT_STRUCTURE_MATCH", "public_inchikey": public_inchikey}
    if public_inchikey[:14] == (local_inchikey or "")[:14]:
        return {"status": "CONNECTIVITY_MATCH_STEREO_UNCERTAIN", "public_inchikey": public_inchikey}
    return {"status": "STRUCTURE_MISMATCH", "public_inchikey": public_inchikey}


def lookup(cas_number: str, local_inchikey: str) -> dict:
    """Resolve CAS through official PubChem PUG REST and return preview-only evidence."""
    cached = _CACHE.get(cas_number)
    if cached and time.monotonic() - cached[0] < _TTL_SECONDS:
        identity = dict(cached[1])
    else:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + quote(cas_number) + "/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
        body = _get_json(url)
        if body.get("_not_found"):
            return {"status": "CAS_NOT_FOUND", "cas_number": cas_number, "records": []}
        if body.get("_error"):
            return {"status": body["_error"], "cas_number": cas_number, "records": []}
        row = (body.get("PropertyTable", {}).get("Properties") or [{}])[0]
        identity = {"cid": row.get("CID"), "canonical_smiles": row.get("ConnectivitySMILES") or row.get("CanonicalSMILES", ""),
                    "isomeric_smiles": row.get("SMILES") or row.get("IsomericSMILES", ""), "inchikey": row.get("InChIKey", "")}
        _CACHE[cas_number] = (time.monotonic(), dict(identity))
    match = _identity_status(identity.get("isomeric_smiles") or identity.get("canonical_smiles", ""), local_inchikey)
    identity.update(match)
    # PubChem calculated properties are deliberately excluded. Activity adapters
    # will add only source-recorded values with a reference in future releases.
    return {"status": "RESULTS_AVAILABLE" if match["status"] == "EXACT_STRUCTURE_MATCH" else match["status"],
            "cas_number": cas_number, "identity": identity, "records": [],
            "source_notice": "PubChem computed properties are not experimental evidence and are not shown as importable values.",
            "retrieved_at": datetime.now(timezone.utc).isoformat()}
