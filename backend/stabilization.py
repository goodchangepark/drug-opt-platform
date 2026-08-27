"""Conservative classification helpers for stabilization project cleanup."""

from __future__ import annotations

import re


CONFIRMED_PATTERNS = (
    re.compile(r"^__STABILIZATION_E2E_TEMP__$", re.I),
    re.compile(r"^E2E(?:\s|$)", re.I),
    re.compile(r"^Stage\s*[345](?:[A-Z0-9-]*)?\s+(?:Test|Browser|Validation|Acceptance|Missing|Route|Freeze|F\s+Precedence)", re.I),
    re.compile(r"^(?:TEST|Test|test)$"),
    re.compile(r"^(?:TEMP|Temporary)(?:\s|$)", re.I),
    re.compile(r"^(?:Browser|Synthetic)\s+(?:Test|E2E|Validation)", re.I),
)


def classify_project(project: dict) -> tuple[str, str]:
    """Classify only strong development markers; uncertainty always survives."""
    name = str(project.get("project_name") or project.get("name") or "").strip()
    target = str(project.get("target") or "").strip()
    description = str(project.get("description") or "").strip()
    compound_count = int(project.get("compound_count") or 0)
    for pattern in CONFIRMED_PATTERNS:
        if pattern.search(name):
            return "CONFIRMED_TEST", f"Project name matches anchored development marker: {pattern.pattern}"
    combined = f"{target} {description}".lower()
    if ("temporary" in combined or "browser acceptance" in combined or "automated e2e" in combined) and re.search(r"\b(test|e2e|validation|acceptance)\b", combined):
        return "CONFIRMED_TEST", "Target/description explicitly identifies a temporary automated validation workflow"
    if name and target and (description or compound_count > 0):
        return "KEEP", "No development marker; project has substantive research metadata or compounds"
    return "AMBIGUOUS", "Insufficient positive evidence to identify either development or genuine research use"

