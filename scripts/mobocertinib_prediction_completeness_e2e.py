#!/usr/bin/env python3
"""Desktop/mobile Mobocertinib Predict All acceptance on an isolated DB."""

from __future__ import annotations

import json
import os
from pathlib import Path

from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

from scripts.prediction_first_browser_e2e import (
    DriverService,
    api,
    click,
    click_tab,
    driver,
    endpoint_cell,
    open_project_compound,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation" / "e2e_mobocertinib_prediction_completeness"
OUT.mkdir(parents=True, exist_ok=True)


def assert_value(browser, endpoint: str) -> None:
    assert endpoint_cell(browser, endpoint, ".prediction-value") not in {"", "—"}, endpoint
    assert endpoint_cell(browser, endpoint, ".developability-status") in {"PREDICTED", "EXPERIMENTAL_AVAILABLE"}, endpoint


def overview_row(browser, endpoint: str):
    section = browser.find_element(By.ID, "overview-developability-results")
    return section.find_element(By.CSS_SELECTOR, f".developability-summary-row[data-endpoint='{endpoint}']")


def assert_overview_order(browser) -> None:
    properties = browser.find_element(By.CSS_SELECTOR, ".overview-properties")
    results = browser.find_element(By.ID, "overview-developability-results")
    follows = browser.execute_script(
        "return Boolean(arguments[0].compareDocumentPosition(arguments[1]) & Node.DOCUMENT_POSITION_FOLLOWING)",
        properties,
        results,
    )
    assert follows, "calculated compound Properties must precede core prediction results"


def main() -> None:
    projects = api("/api/projects")
    project = next(row for row in projects if row["name"] == "EGFR")
    detail = api(f"/api/projects/{project['id']}")
    matches = [row for row in detail["compounds"] if row["name"] == "Mobocertinib"]
    assert len(matches) == 1, "must use the one existing EGFR Mobocertinib"
    mobocertinib = matches[0]
    version_id = mobocertinib["version"]["id"]
    before = api(f"/api/compound-versions/{version_id}/developability-profile")
    before_rows = {row["query_endpoint"]: row for row in before["availability_catalog"]}
    assert before_rows["SOLUBILITY_GENERIC"]["prediction"] is None
    assert before_rows["HIA"]["status"] == "MODEL_NOT_REGISTERED"

    result = {
        "contract": "MobocertinibPredictionCompletenessBrowserE2E/1",
        "project_id": project["id"],
        "compound_row_id": mobocertinib["row_id"],
        "compound_version_id": version_id,
        "production_data_modified": False,
    }
    with DriverService():
        desktop = driver(1440, 900)
        wait = WebDriverWait(desktop, 300)
        try:
            desktop.get(os.environ["DRUGOPT_E2E_BASE_URL"])
            wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(desktop, wait, "EGFR", "Mobocertinib")
            click(desktop, "▶ PREDICT")
            wait.until(lambda d: "Predicted:" in d.page_source and "Failed: 0" in d.page_source)

            wait.until(lambda d: d.find_elements(By.ID, "overview-developability-results"))
            assert_overview_order(desktop)
            expected_overview = {
                "SOLUBILITY_GENERIC": "-5.6371 log10(mol/L)",
                "CACO2_PAPP_AB": "-4.9117 log10(cm/s)",
                "HUMAN_PPB": "92.7 % bound",
                "HLM_CLINT": "1.5253 log10(mL/min/kg)",
                "CYP3A4_INHIBITION": "7.1835 pIC50",
                "HERG_LIABILITY": "7.4388 pIC50",
            }
            for endpoint, expected in expected_overview.items():
                row = overview_row(desktop, endpoint)
                assert row.find_element(By.TAG_NAME, "strong").text.strip() == expected, endpoint
                assert row.find_element(By.CLASS_NAME, "developability-status").text.strip() == "PREDICTED", endpoint
            assert overview_row(desktop, "HIA").find_element(By.CLASS_NAME, "developability-status").text.strip() == "MODEL_NOT_REGISTERED"
            assert overview_row(desktop, "AMES_MUTAGENICITY").find_element(By.CLASS_NAME, "developability-status").text.strip() == "MODEL_UNAVAILABLE"
            assert overview_row(desktop, "HUMAN_PK_AUC_ORAL").find_element(By.CLASS_NAME, "developability-status").text.strip() == "CONTEXT_REQUIRED"
            desktop.execute_script("arguments[0].scrollIntoView({block:'start'})", desktop.find_element(By.CSS_SELECTOR, ".overview-properties"))
            desktop.save_screenshot(str(OUT / "desktop_1440x900_mobocertinib_overview.png"))
            result["overview_properties_then_core_results"] = "PASS"

            click_tab(desktop, "ADMET")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']"))
            for endpoint in ("SOLUBILITY_GENERIC", "CACO2_PAPP_AB", "HUMAN_PPB", "HERG_LIABILITY", "HERG_CLASS", "DILI_LIABILITY"):
                assert_value(desktop, endpoint)
            assert endpoint_cell(desktop, "HIA", ".developability-status") == "MODEL_NOT_REGISTERED"
            assert endpoint_cell(desktop, "BBB_PENETRATION", ".developability-status") == "MODEL_NOT_REGISTERED"
            assert endpoint_cell(desktop, "AMES_MUTAGENICITY", ".developability-status") == "MODEL_UNAVAILABLE"
            assert endpoint_cell(desktop, "PAMPA_PERMEABILITY", ".developability-status") == "MODEL_UNAVAILABLE"

            click_tab(desktop, "METABOLISM")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='HLM_CLINT']"))
            for endpoint in ("HLM_CLINT", "RLM_CLINT", "MLM_CLINT", "CYP1A2_INHIBITION", "CYP2C9_INHIBITION", "CYP2D6_INHIBITION", "CYP3A4_INHIBITION"):
                assert_value(desktop, endpoint)
            assert endpoint_cell(desktop, "METABOLIC_SOFT_SPOTS", ".developability-status") == "MECHANISTIC_ONLY"
            assert "Predicted Metabolites" in desktop.page_source

            click_tab(desktop, "PK")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='HUMAN_PK_AUC_ORAL']"))
            assert endpoint_cell(desktop, "HUMAN_PK_AUC_ORAL", ".developability-status") == "CONTEXT_REQUIRED"
            assert endpoint_cell(desktop, "HUMAN_PK_CL_IV", ".developability-status") == "CURRENT_DATA_CEILING"
            desktop.save_screenshot(str(OUT / "desktop_1440x900_mobocertinib_pk.png"))

            desktop.refresh()
            wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(desktop, wait, "EGFR", "Mobocertinib")
            click_tab(desktop, "ADMET")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='SOLUBILITY_GENERIC']"))
            assert_value(desktop, "SOLUBILITY_GENERIC")
            result["desktop_predict_reload"] = "PASS"
        finally:
            desktop.quit()

        mobile = driver(390, 844)
        mobile_wait = WebDriverWait(mobile, 180)
        try:
            mobile.get(os.environ["DRUGOPT_E2E_BASE_URL"])
            mobile_wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(mobile, mobile_wait, "EGFR", "Mobocertinib")
            mobile_wait.until(lambda d: d.find_elements(By.ID, "overview-developability-results"))
            assert_overview_order(mobile)
            assert overview_row(mobile, "SOLUBILITY_GENERIC").is_displayed()
            assert overview_row(mobile, "HLM_CLINT").is_displayed()
            overview_overflow = mobile.execute_script("return document.documentElement.scrollWidth-document.documentElement.clientWidth")
            assert overview_overflow <= 1, f"mobile Overview horizontal overflow: {overview_overflow}px"
            mobile.execute_script("arguments[0].scrollIntoView({block:'start'})", mobile.find_element(By.CSS_SELECTOR, ".overview-properties"))
            mobile.save_screenshot(str(OUT / "mobile_390x844_mobocertinib_overview.png"))
            result["mobile_overview_properties_then_core_results"] = "PASS"
            click_tab(mobile, "ADMET")
            mobile_wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']"))
            assert_value(mobile, "CACO2_PAPP_AB")
            row = mobile.find_element(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']")
            for label in ("Endpoint", "Prediction", "Experimental", "Status"):
                assert row.find_element(By.CSS_SELECTOR, f"td[data-label='{label}']").is_displayed()
            overflow = mobile.execute_script("return document.documentElement.scrollWidth-document.documentElement.clientWidth")
            assert overflow <= 1, f"mobile horizontal overflow: {overflow}px"
            mobile.save_screenshot(str(OUT / "mobile_390x844_mobocertinib_admet.png"))
            result["mobile_prediction_first_no_overflow"] = "PASS"
        finally:
            mobile.quit()

    after = api(f"/api/compound-versions/{version_id}/developability-profile")
    after_rows = {row["query_endpoint"]: row for row in after["availability_catalog"]}
    assert after["availability_summary"]["AVAILABLE_CURRENT"] == 32
    other_id = next(row["row_id"] for row in detail["compounds"] if row["row_id"] != mobocertinib["row_id"] and row.get("version"))
    compared = api(f"/api/projects/{project['id']}/compare?ids={mobocertinib['row_id']},{other_id}")
    mobo_compare = next(row for row in compared["compounds"] if row["row_id"] == mobocertinib["row_id"])
    checks = {
        "Solubility": "SOLUBILITY_GENERIC", "Caco-2": "CACO2_PAPP_AB", "PPB": "HUMAN_PPB",
        "HLM": "HLM_CLINT", "CYP3A4 pIC50": "CYP3A4_INHIBITION", "hERG pIC50": "HERG_LIABILITY",
    }
    for label, endpoint in checks.items():
        assert mobo_compare["prediction_snapshot_ids"][label] == after_rows[endpoint]["prediction"]["snapshot_id"]
        assert mobo_compare[label] == after_rows[endpoint]["prediction"]["value"]
    result.update({
        "available_current": 32,
        "compare_same_snapshot_and_value": "PASS",
        "status": "PASS",
    })
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
