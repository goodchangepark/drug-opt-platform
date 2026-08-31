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
    records = []
    if match["status"] == "EXACT_STRUCTURE_MATCH" and identity.get("cid"):
        records.extend(_pubchem_experimental_annotations(identity["cid"], cas_number))
        records.extend(_chembl_activities(identity.get("public_inchikey", ""), cas_number))
    return {"status": "RESULTS_AVAILABLE" if match["status"] == "EXACT_STRUCTURE_MATCH" else match["status"],
            "cas_number": cas_number, "identity": identity, "records": records,
            "source_notice": "PubChem computed properties are excluded. Activity values retain source assay context and are not silently mapped to project assays.",
            "retrieved_at": datetime.now(timezone.utc).isoformat()}


def _walk(value):
    if isinstance(value, dict):
        yield value
        for item in value.values():
            yield from _walk(item)
    elif isinstance(value, list):
        for item in value:
            yield from _walk(item)


def _pubchem_experimental_annotations(cid: int, cas_number: str) -> list[dict]:
    body = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Experimental%20Properties")
    if body.get("_error") or body.get("_not_found"):
        return []
    records = []
    for section in _walk(body):
        info = section.get("Information") if isinstance(section, dict) else None
        heading = section.get("TOCHeading", "") if isinstance(section, dict) else ""
        if not isinstance(info, list):
            continue
        for item in info:
            value = item.get("Value", {}).get("StringWithMarkup", [])
            raw_value = "; ".join(x.get("String", "") for x in value if x.get("String"))
            refs = item.get("Reference", []) or []
            reference = "; ".join((r.get("SourceName") or r.get("URL") or "") for r in refs if (r.get("SourceName") or r.get("URL")))
            if not raw_value:
                continue
            records.append({"source": "PubChem", "source_record_id": f"CID:{cid}:{item.get('Name','Experimental property')}",
                            "endpoint": item.get("Name") or heading or "Experimental property", "value": raw_value,
                            "unit": "", "relation": "=", "conditions": item.get("Description", ""),
                            "reference": reference or "REFERENCE_UNRESOLVED", "source_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                            "reference_status": "REFERENCE_RESOLVED" if reference else "REFERENCE_UNRESOLVED",
                            "identity_match_status": "EXACT_STRUCTURE_MATCH", "endpoint_match_status": "ASSAY_CONTEXT_REQUIRED",
                            "import_eligible": False, "evidence_origin": "EXPERIMENTAL_EXTERNAL"})
    return records


def _chembl_activities(inchikey: str, cas_number: str) -> list[dict]:
    if not inchikey:
        return []
    molecules = _get_json("https://www.ebi.ac.uk/chembl/api/data/molecule.json?format=json&limit=5&molecule_structures__standard_inchi_key=" + quote(inchikey))
    molecule_rows = molecules.get("molecules", []) if isinstance(molecules, dict) else []
    if not molecule_rows:
        return []
    chembl_id = molecule_rows[0].get("molecule_chembl_id")
    if not chembl_id:
        return []
    data = _get_json("https://www.ebi.ac.uk/chembl/api/data/activity.json?format=json&limit=100&molecule_chembl_id=" + quote(chembl_id))
    records = []
    for row in data.get("activities", []) if isinstance(data, dict) else []:
        value, unit, kind = row.get("standard_value"), row.get("standard_units"), row.get("standard_type")
        if value in (None, "") or kind not in {"IC50", "EC50", "Ki", "Kd"}:
            continue
        activity_id = str(row.get("activity_id", ""))
        assay_id, document_id = str(row.get("assay_chembl_id", "")), str(row.get("document_chembl_id", ""))
        reference = "ChEMBL " + activity_id + (" · " + document_id if document_id else "")
        records.append({"source": "ChEMBL", "source_record_id": activity_id, "endpoint": kind, "value": str(value), "unit": unit or "",
                        "relation": row.get("standard_relation") or "=", "target": row.get("target_chembl_id", ""),
                        "assay_id": assay_id, "document_id": document_id, "conditions": row.get("assay_description") or "",
                        "species": row.get("target_organism") or "", "reference": reference if activity_id else "REFERENCE_UNRESOLVED",
                        "source_url": "https://www.ebi.ac.uk/chembl/explore/activities/" + activity_id,
                        "reference_status": "REFERENCE_RESOLVED" if activity_id else "REFERENCE_UNRESOLVED",
                        "identity_match_status": "EXACT_STRUCTURE_MATCH", "endpoint_match_status": "ASSAY_CONTEXT_REQUIRED",
                        "import_eligible": False, "evidence_origin": "EXPERIMENTAL_EXTERNAL"})
    return records
