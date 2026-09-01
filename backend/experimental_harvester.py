"""Explicit, public-identifier-only experimental evidence harvesting.

Adapters return preview records only.  They never receive a local SMILES and
never promote literature text or predicted public values to experimental data.
"""
from __future__ import annotations

import hashlib
import io
import os
import re
import subprocess
import zipfile
from urllib.request import Request, urlopen
from dataclasses import dataclass, field
from typing import Any, Protocol
from urllib.parse import quote

from .external_experimental import _get_json, _identity_status
from .experimental_endpoint_aliases import VERSION as ENDPOINT_ALIAS_VERSION, classify_experimental_endpoint
from .qualification_contract import QUALIFICATION_VERSION as QUALIFICATION_CONTRACT_VERSION, aggregate_qualification

QUALITY_A = "A"
QUALITY_B = "B"
QUALITY_C = "C"
QUALITY_D = "D"
HARVESTER_SEARCH_VERSION = "experimental-harvester-v3.1"
DOCUMENT_PARSER_VERSION = "regulatory-supplement-parser-v1"
# v4 names every qualification stage explicitly.  The raw evidence remains
# compatible with older persisted runs and is requalified on read/search.
QUALIFICATION_VERSION = "drugopt-experimental-qualification-v4"


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
    unii: str = ""
    approval: dict = field(default_factory=dict)
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
    classification = classify_experimental_endpoint(label=endpoint, assay_type=result.get("assay_type", ""), description=result.get("conditions", ""), species=result.get("species", ""), cell_line=result.get("cell_line", ""), unit=unit)
    result.update({"evidence_category": result.get("evidence_category", classification["category"]), "canonical_endpoint_candidate": classification["endpoint"],
                   "context_qualified": result.get("context_qualified", classification["qualified"]), "qualification_reason": classification["reason"],
                   "endpoint_alias_version": ENDPOINT_ALIAS_VERSION})
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
    # PubChem is a public identity registry.  Pull aliases after CID
    # resolution so a development code/brand entered by the user expands to
    # the same verified public identity rather than becoming a special case.
    if identity.pubchem_cid:
        synonym_body = _get_json(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{quote(identity.pubchem_cid)}/synonyms/JSON")
        aliases = ((synonym_body.get("InformationList", {}).get("Information") or [{}])[0].get("Synonym") or [])
        identity.synonyms = list(dict.fromkeys(str(x).strip() for x in aliases if str(x).strip()))[:80]
        for alias in identity.synonyms:
            if not identity.chembl_id and re.fullmatch(r"CHEMBL\d+", alias, re.I): identity.chembl_id = alias.upper()
            if not identity.dtxsid and re.fullmatch(r"DTXSID\d+", alias, re.I): identity.dtxsid = alias.upper()
            if not identity.unii and (re.fullmatch(r"[A-Z0-9]{10}", alias, re.I) or alias.upper().startswith("UNII-")):
                identity.unii = alias.removeprefix("UNII-")
    return identity


def _get_document_bytes(url: str, max_bytes: int = 30_000_000) -> tuple[bytes, str]:
    """Fetch a bounded public document in memory only (never persisted)."""
    try:
        req = Request(url, headers={"User-Agent": "Drug-OPT/1.0 public-evidence harvester"})
        with urlopen(req, timeout=8) as response:
            body = response.read(max_bytes + 1)
            content_type = response.headers.get("Content-Type", "")
        return (b"", content_type) if len(body) > max_bytes else (body, content_type)
    except Exception:
        return b"", ""


def _get_document_text(url: str, max_bytes: int = 30_000_000) -> str:
    """Bounded public document extraction; raw files are never persisted."""
    body, content_type = _get_document_bytes(url, max_bytes)
    if not body:
        return ""
    if "pdf" in content_type.lower() or url.lower().endswith(".pdf"):
        try:
            proc = subprocess.run(["pdftotext", "-layout", "-", "-"], input=body, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=12, check=False)
            return proc.stdout.decode("utf-8", errors="replace") if proc.returncode == 0 else ""
        except (OSError, subprocess.TimeoutExpired):
            return ""
    return body.decode("utf-8", errors="replace")


def _supplement_text(url: str) -> tuple[str, str, str]:
    """Return text, file type, and deterministic parse state for public supplements."""
    body, content_type = _get_document_bytes(url, max_bytes=12_000_000)
    suffix = url.rsplit("?", 1)[0].rsplit(".", 1)[-1].lower() if "." in url.rsplit("?", 1)[0] else ""
    file_type = "PDF" if suffix == "pdf" or "pdf" in content_type.lower() else ("CSV" if suffix == "csv" else ("XLSX" if suffix == "xlsx" else ("DOCX" if suffix == "docx" else "OTHER")))
    if not body:
        return "", file_type, "SUPPLEMENT_DOWNLOAD_FAILED"
    if file_type == "PDF":
        try:
            proc = subprocess.run(["pdftotext", "-layout", "-", "-"], input=body, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=12, check=False)
            text = proc.stdout.decode("utf-8", errors="replace") if proc.returncode == 0 else ""
            return text, file_type, "SUPPLEMENT_PARSED" if text.strip() else "SUPPLEMENT_TEXT_EXTRACTION_FAILED"
        except (OSError, subprocess.TimeoutExpired):
            return "", file_type, "SUPPLEMENT_TEXT_EXTRACTION_FAILED"
    if file_type == "CSV":
        return body.decode("utf-8", errors="replace"), file_type, "SUPPLEMENT_PARSED"
    if file_type in {"XLSX", "DOCX"}:
        try:
            with zipfile.ZipFile(io.BytesIO(body)) as archive:
                xml = "\n".join(archive.read(name).decode("utf-8", errors="replace") for name in archive.namelist() if name.endswith(".xml"))
            text = re.sub(r"<[^>]+>", " ", xml)
            return text, file_type, "SUPPLEMENT_PARSED" if text.strip() else "SUPPLEMENT_TEXT_EXTRACTION_FAILED"
        except (zipfile.BadZipFile, KeyError):
            return "", file_type, "SUPPLEMENT_UNSUPPORTED_FORMAT"
    return "", file_type, "SUPPLEMENT_UNSUPPORTED_FORMAT"


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
            label = item.get("Assay Name") or item.get("name") or "BioAssay"
            rows.append(_record(self.name, f"AID:{aid}:CID:{identity.pubchem_cid}", label, item.get("Value") or item.get("value") or outcome, item.get("Unit") or item.get("unit") or "", assay_id=f"AID:{aid}", target=item.get("Target") or item.get("target") or "", activity_outcome=outcome, conditions=item.get("Description") or "", source_url=f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}", reference=f"PubChem AID {aid}", reference_status="REFERENCE_RESOLVED_SOURCE_RECORD"))
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
    ENDPOINTS = ("solubility", "Caco-2", "plasma protein binding", "microsomal stability", "intrinsic clearance", "metabolism", "CYP3A4", "hERG", "pKa", "logD", "pharmacokinetics", "Cmax", "AUC", "half-life", "bioavailability")
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        codes = [x for x in identity.synonyms if re.fullmatch(r"[A-Za-z]{2,}\d+[A-Za-z0-9-]*", x or "") and not x.upper().startswith(("CHEMBL", "DTXSID", "UNII"))]
        aliases = list(dict.fromkeys([identity.name, *(codes[:1])]))
        if not aliases: return []
        found = []
        for term in aliases:
            query = f'"{term}" AND ({" OR ".join(self.ENDPOINTS)})'
            found.extend((_get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=100&query=" + quote(query)).get("resultList", {}).get("result") or []))
        seen_publication, results = set(), []
        for item in found:
            key = item.get("pmid") or item.get("doi") or item.get("id")
            if key and key not in seen_publication:
                seen_publication.add(key); results.append(item)
        rows = []
        supplement_lookups = 0
        for item in results[:100]:
            pmid = str(item.get("pmid") or "")
            doi = item.get("doi") or ""
            pmcid = item.get("pmcid", "")
            rows.append(_record(self.name, pmid or item.get("id", ""), "Literature candidate", "", "", pmid=pmid, pmcid=pmcid, doi=doi, publication_title=item.get("title", ""), journal=item.get("journalTitle", ""), publication_year=item.get("pubYear", ""), abstract=item.get("abstractText", ""), oa_fulltext=bool(item.get("isOpenAccess") == "Y"), reference=(f"PMID: {pmid}" if pmid else f"DOI: {doi}"), reference_status="REFERENCE_RESOLVED_PMID" if pmid else ("REFERENCE_RESOLVED_DOI" if doi else "REFERENCE_UNRESOLVED"), source_quality_class=QUALITY_C, record_status="LITERATURE_CANDIDATE"))
            title_text = re.sub(r"<[^>]+>", "", str(item.get("title", ""))).lower()
            is_identity_paper = any(alias and alias.lower() in (title_text + " " + str(item.get("abstractText", "")).lower()) for alias in [identity.name, *codes])
            is_primary_identity_paper = bool(identity.name and title_text.startswith(identity.name.lower()))
            # PMCID is the authoritative availability signal for the PMC
            # full-text endpoint; Europe PMC's search-level OA flag can be
            # absent for author manuscripts that nevertheless expose JATS
            # supplementary links.
            if pmcid and supplement_lookups < 3 and is_primary_identity_paper:
                supplement_lookups += 1
                # Europe PMC's fullTextXML includes JATS supplementary-media
                # links.  It is more reliable than an HTML ``?report=xml``
                # rendition and gives an auditable parent publication.
                xml = _get_document_text(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{quote(pmcid)}/fullTextXML", max_bytes=5_000_000)
                hrefs = sorted(set(re.findall(r"xlink:href=['\"]([^'\"]+\.(?:pdf|csv|xlsx|docx))", xml, re.I)))[:12]
                for href in hrefs:
                    # PMC mirrors occasionally omit author-hosted files.  We
                    # still record the discovered source and its exact parse
                    # state; inaccessible material is never presented as data.
                    url = href if href.startswith("http") else f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/bin/{href.split('/')[-1]}"
                    source_id = f"{pmcid}:SUP:{len(rows)}"
                    base = dict(evidence_category="LITERATURE", context_qualified=False, parent_pmcid=pmcid, parent_doi=doi,
                                parent_pmid=pmid, publication_title=item.get("title", ""), supplement_filename=href.split("/")[-1],
                                supplement_url=url, reference=f"PMCID: {pmcid}" + (f" · DOI: {doi}" if doi else ""),
                                reference_status="REFERENCE_RESOLVED_PMID" if pmid else "REFERENCE_RESOLVED_SOURCE_RECORD",
                                source_quality_class="A2", source_url=url)
                    rows.append(_record(self.name, source_id, "Supplementary file", "", "", record_status="SUPPLEMENTARY_DISCOVERED", **base))
                    text, file_type, parse_state = _supplement_text(url)
                    rows.append(_record(self.name, source_id + ":PARSE", "Supplementary file", "", "", record_status=parse_state, supplement_file_type=file_type, **base))
                    if parse_state != "SUPPLEMENT_PARSED":
                        continue
                    for line_no, line in enumerate(text.splitlines(), 1):
                        context = line.strip()
                        if not context or not RegulatoryAdapter._terms.search(context) or not RegulatoryAdapter._value.search(context):
                            continue
                        endpoint, value, unit = RegulatoryAdapter.endpoint_value_from_context(context)
                        if not endpoint or not value or not unit:
                            continue
                        category = RegulatoryAdapter.category_for_endpoint(endpoint)
                        rows.append(_record(self.name, source_id + f":L{line_no}", endpoint, value, unit,
                                            evidence_category=category, record_status="LITERATURE_NUMERIC_CANDIDATE",
                                            page_or_line=f"line {line_no}", raw_context=context, conditions=context,
                                            supplement_file_type=file_type, **base))
        return rows


class PKDBAdapter:
    name = "PK-DB"
    def status(self): return "CONFIGURED"
    def harvest(self, identity):
        term = identity.name or identity.cas
        if not term: return []
        data = _get_json("https://pk-db.com/api/v1/studies/?substance=" + quote(term.lower()))
        studies = data.get("data", {}).get("data", []) if isinstance(data, dict) else []
        rows = []
        for study in studies[:25]:
            ref = study.get("reference") or {}
            reference = (f"PMID: {ref.get('pmid')}" if ref.get("pmid") else (f"DOI: {ref.get('doi')}" if ref.get("doi") else f"PK-DB {study.get('sid','')}"))
            rows.append(_record(self.name, str(study.get("sid") or study.get("pk") or "study"), "PK study available", "", "", evidence_category="PK", context_qualified=False, record_status="PK_STUDY_CANDIDATE", study_id=study.get("sid", ""), study=study.get("name", ""), population="Context in PK-DB study", reference=reference, reference_status="REFERENCE_RESOLVED_PMID" if ref.get("pmid") else ("REFERENCE_RESOLVED_DOI" if ref.get("doi") else "REFERENCE_RESOLVED_SOURCE_RECORD"), source_url="https://pk-db.com/api/v1/studies/" + quote(str(study.get("sid", ""))) + "/", conditions="PK-DB public study; inspect dose, route, population and output context before import"))
        return rows


class RegulatoryAdapter:
    """Official Drugs@FDA dossier discovery plus bounded label/review extraction."""
    name = "FDA / Regulatory"
    _terms = re.compile(r"\b(Cmax|Tmax|AUC|half[- ]life|clearance|CL/F|volume of distribution|Vd/F|bioavailability|protein binding|Caco-?2|permeability|solubility|hepatocyte|microsomal|CYP\d[A-Z0-9]+|P-?gp|BCRP|metabolite|feces|urine|excret)\b", re.I)
    _value = re.compile(r"(?:<|>|≤|≥)?\s*\d+(?:[.,]\d+)?\s*(?:ng/mL|µg/mL|mg/L|µM|nM|%|h(?:ours?)?|min(?:utes)?|mL/min(?:/kg)?|µL/min(?:/mg|/10\s*\^?6\s*cells)?|L(?:/kg)?|×\s*10\s*[−-]?\s*\d+\s*cm/s)", re.I)
    def status(self): return "CONFIGURED"

    @classmethod
    def endpoint_value_from_context(cls, context: str) -> tuple[str, str, str]:
        """Choose the measurement label nearest the value, not the first word.

        FDA table rows regularly include a distant Cmax reference beside a
        plasma-binding percentage.  Selecting the first term caused a PPB
        value to masquerade as Cmax.  This deterministic rule is generic and
        retains the full source context for review.
        """
        numeric = cls._value.search(context)
        if not numeric:
            return "", "", ""
        raw_numeric = numeric.group(0).strip()
        numeric_match = re.match(r"([<>≤≥]?\s*\d+(?:[.,]\d+)?)\s*(.*)", raw_numeric)
        value = numeric_match.group(1).replace(" ", "") if numeric_match else raw_numeric
        unit = numeric_match.group(2).strip() if numeric_match else ""
        terms = list(cls._terms.finditer(context))
        # A percentage in a plasma-protein context is PPB even when a nearby
        # explanatory Cmax mentions a different value/units.
        if "%" in unit and re.search(r"plasma protein|protein binding|fraction unbound", context, re.I):
            return "protein binding", value, unit
        endpoint = min(terms, key=lambda m: abs(m.start() - numeric.start())).group(1) if terms else ""
        return endpoint, value, unit

    @staticmethod
    def category_for_endpoint(endpoint: str) -> str:
        return "PK" if endpoint.lower() in {"cmax", "tmax", "auc", "half-life", "clearance", "cl/f", "volume of distribution", "vd/f", "bioavailability"} else ("METABOLISM" if re.search(r"cyp|hepatocyte|microsomal|metabolite|feces|urine|excret", endpoint, re.I) else "ADMET")
    def harvest(self, identity):
        short_aliases = [x for x in identity.synonyms if x and len(x) <= 16 and re.fullmatch(r"[\w\-]+", x) and not re.fullmatch(r"\d{2,}-\d{2}-\d", x)]
        aliases = list(dict.fromkeys([identity.name, identity.unii, *short_aliases]))[:12]
        applications = []
        for alias in aliases:
            if not alias: continue
            data = _get_json("https://api.fda.gov/drug/drugsfda.json?limit=10&search=products.brand_name:%22" + quote(alias) + "%22")
            applications.extend(data.get("results", []) or [])
            data = _get_json("https://api.fda.gov/drug/drugsfda.json?limit=10&search=products.active_ingredients.name:%22" + quote(alias) + "%22")
            applications.extend(data.get("results", []) or [])
            if applications: break
        unique_apps = {str(app.get("application_number")): app for app in applications if app.get("application_number")}
        rows=[]
        for app_number, app in unique_apps.items():
            identity.approval["FDA"] = {"status": "APPROVED", "application_number": app_number}
            for submission in app.get("submissions", []) or []:
                for doc in submission.get("application_docs", []) or []:
                    url, doc_type = doc.get("url", ""), doc.get("type", "Other")
                    if url.startswith("http:"): url = "https:" + url[5:]
                    document_class = "FDA_LABEL" if doc_type.lower() == "label" else ("FDA_MULTIDISCIPLINARY_REVIEW" if doc_type.lower() == "review" else "FDA_OTHER")
                    doc_id = str(doc.get("id") or hashlib.sha256(url.encode()).hexdigest()[:12])
                    rows.append(_record(self.name, f"{app_number}:DOC:{doc_id}", "Regulatory document", "", "", evidence_category="REGULATORY", record_status="DOCUMENT_DISCOVERED", context_qualified=False, application_id=app_number, document_type=document_class, document_date=doc.get("date", ""), reference=f"Drugs@FDA {app_number} · {doc_type}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_quality_class="A1", source_url=url))
                    text = _get_document_text(url)
                    parse_state = "DOCUMENT_PARSED" if text else "DOCUMENT_DOWNLOAD_OR_TEXT_EXTRACTION_FAILED"
                    # A document with text but no qualified measurement terms
                    # is still transparently accounted for.
                    extracted = self._extract_document(app_number, document_class, doc_id, url, text)
                    if text and not extracted:
                        parse_state = "DOCUMENT_NO_RELEVANT_SECTION"
                    rows.append(_record(self.name, f"{app_number}:DOC:{doc_id}:PARSE", "Regulatory document", "", "", evidence_category="REGULATORY", record_status=parse_state, context_qualified=False, application_id=app_number, document_type=document_class, document_date=doc.get("date", ""), reference=f"Drugs@FDA {app_number} · {doc_type}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_quality_class="A1", source_url=url))
                    rows.extend(extracted)
                # The Drugs@FDA review TOC deterministically exposes public
                # review classes for an approved application.
                if submission.get("submission_type") == "ORIG" and submission.get("submission_status") == "AP":
                    app_digits = re.sub(r"\D", "", app_number)
                    year = str(submission.get("submission_status_date", ""))[:4]
                    if len(app_digits) == 6 and year:
                        base = f"https://www.accessdata.fda.gov/drugsatfda_docs/nda/{year}/{app_digits}Orig{submission.get('submission_number','1')}s000"
                        for suffix, kind in (("MultidisciplineR.pdf", "FDA_MULTIDISCIPLINARY_REVIEW"),):
                            url = base + suffix
                            text = _get_document_text(url)
                            if text:
                                rows.append(_record(self.name, f"{app_number}:{suffix}", "Regulatory review", "", "", evidence_category="REGULATORY", record_status="DOCUMENT_PARSED", context_qualified=False, application_id=app_number, document_type=kind, reference=f"Drugs@FDA {app_number} · {kind}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_quality_class="A1", source_url=url))
                                rows.extend(self._extract_document(app_number, kind, suffix, url, text))
        return rows

    def _extract_document(self, app: str, doc_type: str, doc_id: str, url: str, text: str) -> list[dict]:
        rows=[]
        if not text: return rows
        lines = text.splitlines()
        for line_no, line in enumerate(lines, 1):
            # Tables often wrap a row label and its value over adjacent lines.
            context = " ".join(lines[line_no - 1:line_no + 2]).strip()
            if not self._terms.search(context) or not self._value.search(context): continue
            endpoint, value, unit = self.endpoint_value_from_context(context)
            if not endpoint or not value or not unit:
                continue
            category = self.category_for_endpoint(endpoint)
            species = "Human" if re.search(r"\bhuman|patients?|healthy subjects?\b", context, re.I) else ""
            rows.append(_record(self.name, f"{app}:{doc_id}:L{line_no}", endpoint, value, unit, evidence_category=category, context_qualified=False, record_status="REGULATORY_CANDIDATE", application_id=app, document_type=doc_type, page_or_line=f"line {line_no}", raw_context=context, species=species, reference=f"Drugs@FDA {app} · {doc_type}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_quality_class="A1", source_url=url, conditions=context))
            if len(rows) >= 400: break
        return rows


class NMPAAdapter:
    """Official NMPA English-site approval discovery (no third-party index)."""
    name = "NMPA / Regulatory"

    def status(self): return "CONFIGURED"

    def harvest(self, identity):
        aliases = list(dict.fromkeys([identity.name, *identity.synonyms]))
        rows = []
        seen = set()
        for alias in aliases[:12]:
            if not alias or alias in seen:
                continue
            seen.add(alias)
            data = _get_json("https://english.nmpa.gov.cn/dataservice/api/search?index=2%40NMPA&keywords=" + quote(alias) + "&limit=20&page=1")
            for item in data.get("data", []) or []:
                title = re.sub(r"<[^>]+>", "", str(item.get("title", "")))
                abstract = re.sub(r"<[^>]+>", "", str(item.get("abstractdesc", "")))
                text = f"{title} {abstract}".lower()
                # Do not infer approval from a mere search hit.  The official
                # article itself must explicitly state the marketing decision.
                if "approved" not in text or "marketing" not in text:
                    continue
                url = item.get("pubUrl", "")
                key = str(item.get("id") or url)
                if not key or key in seen:
                    continue
                seen.add(key)
                identity.approval["NMPA"] = {"status": "APPROVED", "approval_date": str(item.get("pubTime", ""))[:10], "source_id": key}
                rows.append(_record(self.name, f"NMPA:{key}", "NMPA approval", "", "", evidence_category="REGULATORY",
                                    record_status="NMPA_APPROVAL_CONFIRMED_DOCUMENT_NOT_PUBLICLY_ACCESSIBLE", context_qualified=False,
                                    approval_date=str(item.get("pubTime", ""))[:10], approval_title=title, source_quality_class="A1",
                                    reference=f"NMPA official notice {key}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_url=url,
                                    conditions=abstract))
            if rows:
                break
        return rows


def configured_adapters() -> list[EvidenceSource]:
    return [PubChemPUGViewAdapter(), PubChemBioAssayAdapter(), ChEMBLAdapter(), CompToxAdapter(), BindingDBAdapter(), EuropePMCAdapter(), PKDBAdapter(), RegulatoryAdapter(), NMPAAdapter()]


def harvest_public_evidence(identity: PublicIdentity, enabled_sources: set[str] | None = None) -> dict:
    enabled = enabled_sources or {adapter.name for adapter in configured_adapters()}
    statuses, rows = {}, []
    for adapter in configured_adapters():
        statuses[adapter.name] = adapter.status()
        if adapter.name in enabled and adapter.status() in {"CONFIGURED", "ADAPTER_READY"}:
            try:
                rows.extend(adapter.harvest(identity))
            except Exception as exc:
                # A source outage or parser defect must not discard successful
                # evidence from other adapters or be reported as zero data.
                statuses[adapter.name] = "ERROR"
                statuses[f"{adapter.name} error"] = f"{exc.__class__.__name__}: {exc}"
    unique = deduplicate(identity, rows)
    reference_qualified = [r for r in unique if str(r.get("reference_status", "")).startswith("REFERENCE_RESOLVED")]
    comparable = [r for r in unique if r.get("comparability_status") in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}]
    related = [r for r in unique if r.get("comparability_status") == "RELATED_NOT_SAME_ENDPOINT"]
    literature = [r for r in unique if r.get("record_status") == "LITERATURE_CANDIDATE" or r.get("source") == "Europe PMC"]
    categories = {category: sum(1 for row in unique if row.get("evidence_category") == category) for category in ("ACTIVITY", "ADMET", "METABOLISM", "PK", "TOXICITY", "REGULATORY", "LITERATURE")}
    documents = {kind: sum(1 for row in unique if row.get("record_status") == kind) for kind in (
        "DOCUMENT_DISCOVERED", "DOCUMENT_PARSED", "DOCUMENT_NO_RELEVANT_SECTION", "DOCUMENT_DOWNLOAD_OR_TEXT_EXTRACTION_FAILED",
        "SUPPLEMENTARY_DISCOVERED", "SUPPLEMENT_PARSED", "SUPPLEMENT_DOWNLOAD_FAILED", "SUPPLEMENT_TEXT_EXTRACTION_FAILED",
        "SUPPLEMENT_UNSUPPORTED_FORMAT", "LITERATURE_CANDIDATE", "LITERATURE_NUMERIC_CANDIDATE",
    )}
    source_counts = {key: {"found": 0 if value in {"NOT_CONFIGURED", "ADAPTER_READY"} else sum(1 for r in rows if r.get("source") == key), "unique": 0, "numeric": 0, "identity_qualified": 0, "reference_qualified": 0, "endpoint_qualified": 0, "context_qualified": 0, "prediction_pairable": 0, "direct": 0, "conditional": 0, "related": 0, "ready_to_import": 0, "adaptation_eligible": 0} for key, value in statuses.items()}
    qualification = aggregate_qualification(unique, raw_source_counts={key: value["found"] for key, value in source_counts.items()})
    for key, value in qualification["sources"].items():
        source_counts[key] = value
    return {"identity": identity.to_dict(), "sources": statuses, "records": unique,
            "summary": {"sources_searched": sum(1 for v in statuses.values() if v not in {"NOT_CONFIGURED", "ADAPTER_READY"}),
                        "harvester_search_version": HARVESTER_SEARCH_VERSION, "document_parser_version": DOCUMENT_PARSER_VERSION,
                        "qualification_version": QUALIFICATION_CONTRACT_VERSION, "raw_records": len(rows), "raw_found": len(rows), "unique_records": len(unique),
                        "reference_qualified": len(reference_qualified), "directly_comparable": len(comparable),
                        "related_evidence": len(related), "literature_candidates": len(literature),
                        "duplicates_removed": len(rows) - len(unique), "imported": 0, "categories": categories, "documents": documents},
            "source_counts": source_counts}
