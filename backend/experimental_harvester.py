"""Explicit, public-identifier-only experimental evidence harvesting.

Adapters return preview records only.  They never receive a local SMILES and
never promote literature text or predicted public values to experimental data.
"""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from typing import Any, Protocol
from urllib.parse import quote

from .external_experimental import _get_json, _identity_status

QUALITY_A = "A"
QUALITY_B = "B"
QUALITY_C = "C"
QUALITY_D = "D"


@dataclass
class PublicIdentity:
    cas: str = ""
    name: str = ""
    synonyms: list[str] = field(default_factory=list)
    pubchem_cid: str = ""
    inchikey: str = ""
    chembl_id: str = ""
    dtxsid: str = ""
    doi: str = ""
    pmid: str = ""
    identity_status: str = "UNRESOLVED"

    def public_query_terms(self) -> list[str]:
        return [x for x in (self.cas, self.name, self.pubchem_cid, self.chembl_id, self.dtxsid, *self.synonyms) if x]

    def to_dict(self) -> dict:
        return vars(self) | {"synonyms": self.synonyms}


class EvidenceSource(Protocol):
    name: str
    def status(self) -> str: ...
    def harvest(self, identity: PublicIdentity) -> list[dict]: ...


def _record(source: str, source_record_id: str, endpoint: str, value: Any = "", unit: str = "", **extra) -> dict:
    result = {"source": source, "source_record_id": str(source_record_id), "endpoint": endpoint,
              "value": str(value), "unit": unit, "relation": "=", "evidence_origin": "EXPERIMENTAL_EXTERNAL",
              "source_quality_class": QUALITY_B, "reference_status": "REFERENCE_UNRESOLVED",
              "identity_match_status": "EXACT_STRUCTURE_MATCH", "endpoint_match_status": "ASSAY_CONTEXT_REQUIRED",
              "import_eligible": False}
    result.update(extra)
    return result


def evidence_fingerprint(identity: PublicIdentity, row: dict) -> str:
    """Cross-source provenance fingerprint: preserves, rather than collapses, observations."""
    payload = "|".join(str(row.get(key, "")).strip().lower() for key in (
        "endpoint", "measurement_type", "target", "value", "relation", "unit", "species", "assay_id", "doi", "pmid",
    ))
    return hashlib.sha256(f"{identity.inchikey.lower()}|{payload}".encode()).hexdigest()


def deduplicate(identity: PublicIdentity, records: list[dict]) -> list[dict]:
    seen: dict[str, dict] = {}
    result = []
    for row in records:
        key = evidence_fingerprint(identity, row)
        if key in seen:
            row["duplicate_status"] = "SAME_MEASUREMENT"
            row["duplicate_of"] = seen[key].get("source_record_id", "")
        else:
            row["duplicate_status"] = "DISTINCT_MEASUREMENT"
            seen[key] = row
        row["provenance_fingerprint"] = key
        result.append(row)
    return result


def resolve_public_identity(*, cas: str = "", name: str = "", pubchem_cid: str = "", chembl_id: str = "", dtxsid: str = "", local_inchikey: str = "") -> PublicIdentity:
    """Resolve only identifiers explicitly supplied by the user or local public metadata."""
    identity = PublicIdentity(cas=cas.strip(), name=name.strip(), pubchem_cid=str(pubchem_cid or "").strip(), chembl_id=chembl_id.strip(), dtxsid=dtxsid.strip())
    query = identity.pubchem_cid or identity.cas or identity.name
    if not query:
        return identity
    route = f"cid/{quote(query)}" if identity.pubchem_cid else f"name/{quote(query)}"
    body = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{route}/property/CanonicalSMILES,IsomericSMILES,InChIKey,Title/JSON")
    row = (body.get("PropertyTable", {}).get("Properties") or [{}])[0]
    identity.pubchem_cid = str(row.get("CID") or identity.pubchem_cid)
    identity.inchikey = row.get("InChIKey", "")
    if row.get("Title") and not identity.name: identity.name = row["Title"]
    public_smiles = row.get("SMILES") or row.get("IsomericSMILES") or row.get("ConnectivitySMILES") or row.get("CanonicalSMILES", "")
    if local_inchikey and public_smiles:
        identity.identity_status = _identity_status(public_smiles, local_inchikey)["status"]
    elif identity.inchikey:
        identity.identity_status = "PUBLIC_IDENTIFIER_RESOLVED"
    return identity


class PubChemPUGViewAdapter:
    name = "PubChem PUG View"
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        if not identity.pubchem_cid: return []
        body = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{quote(identity.pubchem_cid)}/JSON?heading=Experimental%20Properties")
        rows = []
        def walk(v):
            if isinstance(v, dict):
                yield v
                for child in v.values(): yield from walk(child)
            elif isinstance(v, list):
                for child in v: yield from walk(child)
        for section in walk(body):
            for info in section.get("Information", []) if isinstance(section, dict) else []:
                strings = info.get("Value", {}).get("StringWithMarkup", [])
                value = "; ".join(x.get("String", "") for x in strings if isinstance(x, dict) and x.get("String"))
                refs = info.get("Reference", []) or []
                reference = "; ".join((x if isinstance(x, str) else x.get("SourceName") or x.get("URL") or "") for x in refs if isinstance(x, (str, dict)))
                if value:
                    rows.append(_record(self.name, f"CID:{identity.pubchem_cid}:{info.get('Name','property')}", info.get("Name") or section.get("TOCHeading", "Experimental property"), value, reference=reference or "REFERENCE_UNRESOLVED", reference_status="REFERENCE_RESOLVED" if reference else "REFERENCE_UNRESOLVED", source_url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{identity.pubchem_cid}"))
        return rows


class PubChemBioAssayAdapter:
    name = "PubChem BioAssay"
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        if not identity.pubchem_cid: return []
        body = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{quote(identity.pubchem_cid)}/assaysummary/JSON")
        rows = []
        for item in body.get("Table", {}).get("Row", []) or body.get("assays", []) or []:
            aid = item.get("AID") or item.get("aid") or ""
            outcome = item.get("Activity Outcome") or item.get("outcome") or ""
            rows.append(_record(self.name, f"AID:{aid}:CID:{identity.pubchem_cid}", item.get("Assay Name") or item.get("name") or "BioAssay", item.get("Value") or item.get("value") or outcome, item.get("Unit") or item.get("unit") or "", assay_id=f"AID:{aid}", target=item.get("Target") or item.get("target") or "", activity_outcome=outcome, conditions=item.get("Description") or "", source_url=f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}", reference=f"PubChem AID {aid}", reference_status="REFERENCE_RESOLVED_SOURCE_RECORD"))
        return rows


class ChEMBLAdapter:
    name = "ChEMBL"
    ALLOWED_ASSAY_TYPES = {"B", "F", "A", "T", "P"}
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        if not identity.inchikey and not identity.chembl_id: return []
        chembl_id = identity.chembl_id
        if not chembl_id:
            molecules = _get_json("https://www.ebi.ac.uk/chembl/api/data/molecule.json?format=json&limit=1&molecule_structures__standard_inchi_key=" + quote(identity.inchikey))
            chembl_id = ((molecules.get("molecules") or [{}])[0]).get("molecule_chembl_id", "")
        if not chembl_id: return []
        identity.chembl_id = chembl_id
        data = _get_json("https://www.ebi.ac.uk/chembl/api/data/activity.json?format=json&limit=100&molecule_chembl_id=" + quote(chembl_id))
        rows = []
        for item in data.get("activities", []) or []:
            assay_type = item.get("assay_type", "")
            value, kind = item.get("standard_value"), item.get("standard_type")
            if assay_type not in self.ALLOWED_ASSAY_TYPES or value in (None, "") or not kind: continue
            activity_id = item.get("activity_id", "")
            doc = item.get("document_chembl_id", "")
            rows.append(_record(self.name, activity_id, kind, value, item.get("standard_units") or "", relation=item.get("standard_relation") or "=", measurement_type=kind, assay_type=assay_type, assay_id=item.get("assay_chembl_id", ""), document_id=doc, target=item.get("target_chembl_id", ""), species=item.get("target_organism", ""), conditions=item.get("assay_description", ""), reference=f"ChEMBL {activity_id}" + (f" · {doc}" if doc else ""), reference_status="REFERENCE_RESOLVED_SOURCE_RECORD", source_url=f"https://www.ebi.ac.uk/chembl/explore/activities/{activity_id}"))
        return rows


class CompToxAdapter:
    name = "CompTox"
    def status(self): return "CONFIGURED" if os.getenv("COMPTOX_API_KEY") else "NOT_CONFIGURED"
    def harvest(self, identity):
        if self.status() != "CONFIGURED": return []
        # API-key implementation point: never calls CompTox without the user-provided key.
        return []


class BindingDBAdapter:
    name = "BindingDB"
    def status(self): return "ADAPTER_READY"
    def harvest(self, identity):
        if not identity.inchikey: return []
        # BindingDB's public endpoint availability varies; fail closed without a response.
        data = _get_json("https://bindingdb.org/axis2/services/BDBService/getLigandsByInchiKey?inchiKey=" + quote(identity.inchikey))
        rows = []
        for item in data.get("records", []) if isinstance(data, dict) else []:
            kind = item.get("measurement_type", "")
            if kind not in {"IC50", "Ki", "Kd"}: continue
            rows.append(_record(self.name, item.get("id", ""), kind, item.get("value", ""), item.get("unit", ""), target=item.get("target", ""), uniprot=item.get("uniprot", ""), reference=item.get("reference", "BindingDB"), reference_status="REFERENCE_RESOLVED_SOURCE_RECORD"))
        return rows


class EuropePMCAdapter:
    name = "Europe PMC"
    ENDPOINTS = ("solubility", "Caco-2", "plasma protein binding", "microsomal stability", "CYP3A4", "hERG", "pKa", "logD")
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        term = identity.cas or identity.name
        if not term: return []
        query = f'("{term}") AND ({" OR ".join(self.ENDPOINTS)})'
        data = _get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=25&query=" + quote(query))
        rows = []
        for item in data.get("resultList", {}).get("result", []) or []:
            pmid = str(item.get("pmid") or "")
            doi = item.get("doi") or ""
            rows.append(_record(self.name, pmid or item.get("id", ""), "Literature candidate", "", "", pmid=pmid, pmcid=item.get("pmcid", ""), doi=doi, publication_title=item.get("title", ""), journal=item.get("journalTitle", ""), publication_year=item.get("pubYear", ""), abstract=item.get("abstractText", ""), oa_fulltext=bool(item.get("isOpenAccess") == "Y"), reference=(f"PMID: {pmid}" if pmid else f"DOI: {doi}"), reference_status="REFERENCE_RESOLVED_PMID" if pmid else ("REFERENCE_RESOLVED_DOI" if doi else "REFERENCE_UNRESOLVED"), source_quality_class=QUALITY_C, evidence_status="LITERATURE_CANDIDATE"))
        return rows


class PKDBAdapter:
    name = "PK-DB"
    def status(self): return "ADAPTER_READY"
    def harvest(self, identity): return []  # optional public PK adapter; never mixes with in-vitro clearance


def configured_adapters() -> list[EvidenceSource]:
    return [PubChemPUGViewAdapter(), PubChemBioAssayAdapter(), ChEMBLAdapter(), CompToxAdapter(), BindingDBAdapter(), EuropePMCAdapter(), PKDBAdapter()]


def harvest_public_evidence(identity: PublicIdentity, enabled_sources: set[str] | None = None) -> dict:
    enabled = enabled_sources or {adapter.name for adapter in configured_adapters()}
    statuses, rows = {}, []
    for adapter in configured_adapters():
        statuses[adapter.name] = adapter.status()
        if adapter.name in enabled and adapter.status() in {"CONFIGURED", "ADAPTER_READY"}:
            rows.extend(adapter.harvest(identity))
    return {"identity": identity.to_dict(), "sources": statuses, "records": deduplicate(identity, rows)}
