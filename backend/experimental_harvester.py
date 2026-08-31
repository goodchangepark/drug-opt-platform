"""Explicit, public-identifier-only experimental evidence harvesting.

Adapters return preview records only.  They never receive a local SMILES and
never promote literature text or predicted public values to experimental data.
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess
from urllib.request import Request, urlopen
from dataclasses import dataclass, field
from typing import Any, Protocol
from urllib.parse import quote

from .external_experimental import _get_json, _identity_status
from .experimental_endpoint_aliases import VERSION as ENDPOINT_ALIAS_VERSION, classify_experimental_endpoint

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


def _get_document_text(url: str, max_bytes: int = 30_000_000) -> str:
    """Bounded public document extraction; raw files are never persisted."""
    try:
        req = Request(url, headers={"User-Agent": "Drug-OPT/1.0 public-evidence harvester"})
        with urlopen(req, timeout=8) as response:
            body = response.read(max_bytes + 1)
            content_type = response.headers.get("Content-Type", "")
        if len(body) > max_bytes: return ""
        if "pdf" in content_type.lower() or url.lower().endswith(".pdf"):
            proc = subprocess.run(["pdftotext", "-layout", "-", "-"], input=body, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=12, check=False)
            return proc.stdout.decode("utf-8", errors="replace") if proc.returncode == 0 else ""
        return body.decode("utf-8", errors="replace")
    except Exception:
        return ""


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
            if pmcid and item.get("isOpenAccess") == "Y" and ("selective egfr inhibitor" in str(item.get("title", "")).lower() or supplement_lookups < 1):
                supplement_lookups += 1
                xml = _get_document_text(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{quote(pmcid)}/?report=xml", max_bytes=5_000_000)
                for href in sorted(set(re.findall(r"(?:href|xlink:href)=['\"]([^'\"]+(?:supp|table|data)[^'\"]*)", xml, re.I)))[:12]:
                    url = href if href.startswith("http") else f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/bin/{href.split('/')[-1]}"
                    rows.append(_record(self.name, f"{pmcid}:SUP:{len(rows)}", "Supplementary file", "", "", evidence_category="LITERATURE", context_qualified=False, record_status="SUPPLEMENTARY_DISCOVERED", parent_pmcid=pmcid, parent_doi=doi, supplement_url=url, reference=f"PMCID: {pmcid}" + (f" · DOI: {doi}" if doi else ""), reference_status="REFERENCE_RESOLVED_PMID" if pmid else "REFERENCE_RESOLVED_SOURCE_RECORD", source_quality_class="A2", source_url=url))
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
                    rows.extend(self._extract_document(app_number, document_class, doc_id, url, text))
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
            endpoint = self._terms.search(context).group(1)
            value = self._value.search(context).group(0)
            category = "PK" if endpoint.lower() in {"cmax", "tmax", "auc", "half-life", "clearance", "cl/f", "volume of distribution", "vd/f", "bioavailability"} else ("METABOLISM" if re.search(r"cyp|hepatocyte|microsomal|metabolite|feces|urine|excret", endpoint, re.I) else "ADMET")
            rows.append(_record(self.name, f"{app}:{doc_id}:L{line_no}", endpoint, value, "", evidence_category=category, context_qualified=False, record_status="REGULATORY_CANDIDATE", application_id=app, document_type=doc_type, page_or_line=f"line {line_no}", raw_context=context, reference=f"Drugs@FDA {app} · {doc_type}", reference_status="REFERENCE_RESOLVED_REGULATORY", source_quality_class="A1", source_url=url, conditions="Review raw regulatory context before normalization"))
            if len(rows) >= 400: break
        return rows


def configured_adapters() -> list[EvidenceSource]:
    return [PubChemPUGViewAdapter(), PubChemBioAssayAdapter(), ChEMBLAdapter(), CompToxAdapter(), BindingDBAdapter(), EuropePMCAdapter(), PKDBAdapter(), RegulatoryAdapter()]


def harvest_public_evidence(identity: PublicIdentity, enabled_sources: set[str] | None = None) -> dict:
    enabled = enabled_sources or {adapter.name for adapter in configured_adapters()}
    statuses, rows = {}, []
    for adapter in configured_adapters():
        statuses[adapter.name] = adapter.status()
        if adapter.name in enabled and adapter.status() in {"CONFIGURED", "ADAPTER_READY"}:
            rows.extend(adapter.harvest(identity))
    unique = deduplicate(identity, rows)
    reference_qualified = [r for r in unique if str(r.get("reference_status", "")).startswith("REFERENCE_RESOLVED")]
    comparable = [r for r in unique if r.get("comparability_status") in {"DIRECTLY_COMPARABLE", "COMPARABLE_AFTER_DETERMINISTIC_CONVERSION"}]
    related = [r for r in unique if r.get("comparability_status") == "RELATED_NOT_SAME_ENDPOINT"]
    literature = [r for r in unique if r.get("record_status") == "LITERATURE_CANDIDATE" or r.get("source") == "Europe PMC"]
    categories = {category: sum(1 for row in unique if row.get("evidence_category") == category) for category in ("ACTIVITY", "ADMET", "METABOLISM", "PK", "TOXICITY", "REGULATORY", "LITERATURE")}
    documents = {kind: sum(1 for row in unique if row.get("record_status") == kind) for kind in ("DOCUMENT_DISCOVERED", "DOCUMENT_PARSED", "SUPPLEMENTARY_DISCOVERED", "LITERATURE_CANDIDATE")}
    source_counts = {key: {"found": 0 if value in {"NOT_CONFIGURED", "ADAPTER_READY"} else sum(1 for r in rows if r.get("source") == key), "qualified": sum(1 for r in unique if r.get("source") == key and r.get("context_qualified") and str(r.get("reference_status", "")).startswith("REFERENCE_RESOLVED"))} for key, value in statuses.items()}
    return {"identity": identity.to_dict(), "sources": statuses, "records": unique,
            "summary": {"sources_searched": sum(1 for v in statuses.values() if v not in {"NOT_CONFIGURED", "ADAPTER_READY"}),
                        "raw_records": len(rows), "unique_records": len(unique),
                        "reference_qualified": len(reference_qualified), "directly_comparable": len(comparable),
                        "related_evidence": len(related), "literature_candidates": len(literature),
                        "duplicates_removed": len(rows) - len(unique), "imported": 0, "categories": categories, "documents": documents},
            "source_counts": source_counts}
