#!/usr/bin/env python3
"""Browser acceptance against the running production Drug-OPT application."""

from __future__ import annotations

import json
import os
from pathlib import Path

from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException, StaleElementReferenceException
from selenium.webdriver.support.ui import WebDriverWait

from scripts.prediction_first_browser_e2e import DriverService, api, click, driver, open_project_compound


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "mobocertinib_real_predict_browser_e2e"
OUT.mkdir(parents=True, exist_ok=True)


def overview_row(browser, endpoint: str):
    return browser.find_element(
        By.CSS_SELECTOR,
        f"#overview-developability-results .developability-summary-row[data-endpoint='{endpoint}']",
    )


def value(browser, endpoint: str) -> str:
    for _ in range(10):
        try:
            return overview_row(browser, endpoint).find_element(By.TAG_NAME, "strong").text.strip()
        except (NoSuchElementException, StaleElementReferenceException):
            continue
    raise AssertionError(f"stale value row: {endpoint}")


def status(browser, endpoint: str) -> str:
    for _ in range(10):
        try:
            return overview_row(browser, endpoint).find_element(By.CLASS_NAME, "developability-status").text.strip()
        except (NoSuchElementException, StaleElementReferenceException):
            continue
    raise AssertionError(f"stale status row: {endpoint}")


def wait_for_row(wait, endpoint: str):
    return wait.until(lambda d: d.find_elements(
        By.CSS_SELECTOR,
        f"#overview-developability-results .developability-summary-row[data-endpoint='{endpoint}']",
    ))


def wait_value(wait, endpoint: str, expected: str) -> None:
    def check(browser):
        try:
            return value(browser, endpoint) == expected
        except (NoSuchElementException, StaleElementReferenceException, AssertionError):
            return False
    assert wait.until(check), f"unexpected UI value: {endpoint}"


def wait_status(wait, endpoint: str, expected: str) -> None:
    def check(browser):
        try:
            return status(browser, endpoint) == expected
        except (NoSuchElementException, StaleElementReferenceException, AssertionError):
            return False
    assert wait.until(check), f"unexpected UI status: {endpoint}"


def main() -> None:
    projects = api("/api/projects")
    project = next(row for row in projects if row["name"] == "EGFR")
    detail = api(f"/api/projects/{project['id']}")
    mobo = next(row for row in detail["compounds"] if row["name"] == "Mobocertinib")
    other = next(row for row in detail["compounds"] if row["row_id"] != mobo["row_id"])
    version_id = mobo["version"]["id"]
    profile = api(f"/api/compound-versions/{version_id}/developability-profile")
    profile_rows = {row["query_endpoint"]: row for row in profile["availability_catalog"]}
    expected = {
        "SOLUBILITY_GENERIC": "-5.6371 log10(mol/L)",
        "CACO2_PAPP_AB": "-4.9117 log10(cm/s)",
        "HUMAN_PPB": "92.7 % bound",
        "HLM_CLINT": "1.5253 log10(mL/min/kg)",
        "CYP3A4_INHIBITION": "7.1835 pIC50",
        "HERG_LIABILITY": "7.4388 pIC50",
    }
    result = {
        "contract": "MobocertinibRealPredictProductionBrowserE2E/1",
        "base_url": os.environ.get("DRUGOPT_E2E_BASE_URL", "http://127.0.0.1:8765"),
        "project_id": project["id"],
        "compound_row_id": mobo["row_id"],
        "compound_version_id": version_id,
        "production_data_modified_by_explicit_predict": True,
    }
    with DriverService():
        browser = driver(1440, 900)
        wait = WebDriverWait(browser, 180)
        try:
            browser.get(result["base_url"])
            wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(browser, wait, "EGFR", "Mobocertinib")
            click(browser, "▶ PREDICT")
            wait.until(lambda d: "Current Prediction: ✓ Persisted" in d.page_source)
            wait.until(lambda d: d.find_elements(By.ID, "overview-developability-results"))
            for endpoint, expected_value in expected.items():
                wait_for_row(wait, endpoint)
                wait_value(wait, endpoint, expected_value)
                wait_status(wait, endpoint, "PREDICTED")
            for endpoint in ("HIA", "AMES_MUTAGENICITY", "HUMAN_PK_AUC_ORAL"):
                wait_for_row(wait, endpoint)
            wait_status(wait, "HIA", "MODEL_NOT_REGISTERED")
            wait_status(wait, "AMES_MUTAGENICITY", "MODEL_UNAVAILABLE")
            wait_status(wait, "HUMAN_PK_AUC_ORAL", "CONTEXT_REQUIRED")
            browser.save_screenshot(str(OUT / "production_overview.png"))
            result["overview_predict"] = "PASS"

            browser.refresh()
            wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(browser, wait, "EGFR", "Mobocertinib")
            wait.until(lambda d: d.find_elements(By.ID, "overview-developability-results"))
            for endpoint, expected_value in expected.items():
                wait_for_row(wait, endpoint)
                wait_value(wait, endpoint, expected_value)
            result["reload_persistence"] = "PASS"

            # The production UI Compare flow: return to the project list,
            # select Mobocertinib plus one real EGFR compound, then open Compare.
            click(browser, "Back to Compounds")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, ".compound-row"))
            rows = browser.find_elements(By.CSS_SELECTOR, ".compound-row")
            mobo_row = next(row for row in rows if row.get_attribute("data-compound-id") == mobo["compound_id"])
            other_row = next(row for row in rows if row.get_attribute("data-compound-id") != mobo["compound_id"])
            for row in (mobo_row, other_row):
                checkbox = row.find_element(By.CSS_SELECTOR, "input.compound-select")
                browser.execute_script("arguments[0].click()", checkbox)
            click(browser, "Compare Selected")
            wait.until(lambda d: "Selected Compound Comparison" in d.page_source)
            assert "Mobocertinib" in browser.page_source
            result["compare_ui"] = "PASS"
            browser.save_screenshot(str(OUT / "production_compare.png"))
        finally:
            browser.quit()

    after = api(f"/api/compound-versions/{version_id}/developability-profile")
    after_rows = {row["query_endpoint"]: row for row in after["availability_catalog"]}
    comparison = api(f"/api/projects/{project['id']}/compare?ids={mobo['row_id']},{other['row_id']}")
    compare_row = next(row for row in comparison["compounds"] if row["row_id"] == mobo["row_id"])
    parity = {}
    for label, endpoint in {"Solubility": "SOLUBILITY_GENERIC", "Caco-2": "CACO2_PAPP_AB", "PPB": "HUMAN_PPB", "HLM": "HLM_CLINT"}.items():
        parity[label] = {
            "snapshot_id_equal": compare_row["prediction_snapshot_ids"][label] == after_rows[endpoint]["prediction"]["snapshot_id"],
            "value_equal": compare_row[label] == after_rows[endpoint]["prediction"]["value"],
            "snapshot_id": after_rows[endpoint]["prediction"]["snapshot_id"],
            "value": after_rows[endpoint]["prediction"]["value"],
        }
        assert parity[label]["snapshot_id_equal"] and parity[label]["value_equal"], label
    result["compare_snapshot_value_parity"] = parity
    result["status"] = "PASS"
    (OUT / "result.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
