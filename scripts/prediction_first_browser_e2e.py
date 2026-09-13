#!/usr/bin/env python3
"""Desktop/mobile acceptance for the prediction-first developability UX.

Run only through ``run_isolated_e2e.py``.  This scenario intentionally invokes
Predict and creates test evidence, but every mutation is confined to the
private E2E database copy.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
import urllib.request
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.support.ui import WebDriverWait


ROOT = Path(__file__).resolve().parents[1]
BASE_URL = os.environ["DRUGOPT_E2E_BASE_URL"].rstrip("/")
CHROMEDRIVER = "/snap/bin/chromium.chromedriver"
OUT = ROOT / "validation" / "e2e_prediction_first"
OUT.mkdir(parents=True, exist_ok=True)


def api(path: str, method: str = "GET", payload=None):
    data = None if payload is None else json.dumps(payload).encode()
    request = urllib.request.Request(
        BASE_URL + path, data=data, method=method,
        headers={"Content-Type": "application/json", "X-DrugOPT-Caller": "E2E"},
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        return json.loads(response.read())


class DriverService:
    def __enter__(self):
        self.process = subprocess.Popen(
            [CHROMEDRIVER, "--port=9517", "--allowed-ips=127.0.0.1"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        for _ in range(40):
            try:
                urllib.request.urlopen("http://127.0.0.1:9517/status", timeout=1)
                return self
            except Exception:
                time.sleep(.25)
        raise RuntimeError("ChromeDriver did not start")

    def __exit__(self, *_):
        self.process.terminate()
        try:
            self.process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.process.kill()


def driver(width: int, height: int):
    options = Options()
    for argument in ("--headless", "--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"):
        options.add_argument(argument)
    options.add_argument(f"--window-size={width},{height}")
    return webdriver.Remote("http://127.0.0.1:9517", options=options)


def click(driver, text: str):
    matches = driver.find_elements(By.XPATH, f"//button[normalize-space()='{text}']")
    visible = next((row for row in matches if row.is_displayed() and row.is_enabled()), None)
    if not visible:
        raise AssertionError(f"visible button {text!r} not found")
    driver.execute_script("arguments[0].click()", visible)


def click_tab(driver, label: str):
    for button in driver.find_elements(By.CSS_SELECTOR, ".detail-tabs button"):
        label_nodes = button.find_elements(By.TAG_NAME, "span")
        if label_nodes and label_nodes[0].text.strip() == label:
            driver.execute_script("arguments[0].click()", button)
            return
    raise AssertionError(f"tab {label!r} not found")


def open_project_compound(driver, wait, project_name: str, compound_name: str):
    if not any(row.is_displayed() for row in driver.find_elements(By.XPATH, "//button[normalize-space()='Projects']")):
        menu = next((row for row in driver.find_elements(By.CSS_SELECTOR, ".menu-toggle") if row.is_displayed()), None)
        if menu:
            driver.execute_script("arguments[0].click()", menu)
    click(driver, "Projects")
    wait.until(lambda d: project_name in d.page_source)
    click(driver, project_name)
    wait.until(lambda d: compound_name in d.page_source)
    click(driver, compound_name)
    wait.until(lambda d: "Developability Summary" in d.page_source)


def endpoint_cell(driver, endpoint: str, selector: str):
    for _ in range(20):
        try:
            row = driver.find_element(By.CSS_SELECTOR, f"tr[data-endpoint='{endpoint}']")
            return row.find_element(By.CSS_SELECTOR, selector).text.strip()
        except StaleElementReferenceException:
            time.sleep(.1)
    raise AssertionError(f"endpoint row {endpoint} remained stale")


def main():
    # Add the disposable compound to an existing visible project in the
    # private database. New projects are intentionally forced synthetic in
    # E2E mode and therefore hidden from the ordinary portfolio UI.
    project = next(row for row in api("/api/projects") if row["id"] == 324)
    compound = api(f"/api/projects/{project['id']}/compounds", "POST", {
        "compound_id": "NEW-LIKE-001", "name": "New-like Discovery Compound",
        "smiles": "CCOc1ccc(C(=O)N2CCN(C)CC2)cc1", "calculate": False,
    })
    version_id = compound["version"]["id"]
    initial = api(f"/api/compound-versions/{version_id}/developability-profile")
    initial_rows = {row["query_endpoint"]: row for row in initial["availability_catalog"]}
    assert initial_rows["CACO2_PAPP_AB"]["status"] == "ON_DEMAND"
    assert initial_rows["PAMPA_PERMEABILITY"]["status"] == "MODEL_UNAVAILABLE"
    assert initial_rows["HUMAN_PK_AUC_ORAL"]["status"] == "CONTEXT_REQUIRED"

    results = {"contract": "PredictionFirstBrowserE2E/v1", "new_compound": compound["row_id"]}
    with DriverService():
        desktop = driver(1440, 900)
        wait = WebDriverWait(desktop, 300)
        try:
            desktop.get(BASE_URL)
            wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(desktop, wait, project["name"], compound["name"])
            for label in ("PAMPA permeability", "Plasma Stability (PS)", "hERG quantitative inhibition"):
                assert label in desktop.page_source
            click(desktop, "▶ PREDICT")
            wait.until(lambda d: "Predicted:" in d.page_source and "Failed:" in d.page_source)
            click_tab(desktop, "ADMET")
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']"))
            assert endpoint_cell(desktop, "CACO2_PAPP_AB", ".prediction-value") != "—"
            assert endpoint_cell(desktop, "HUMAN_PPB", ".prediction-value") != "—"
            assert endpoint_cell(desktop, "HERG_CLASS", ".prediction-value") != "—"
            assert endpoint_cell(desktop, "HERG_LIABILITY", ".developability-status") == "CURRENT_DATA_CEILING"
            assert endpoint_cell(desktop, "PAMPA_PERMEABILITY", ".developability-status") == "MODEL_UNAVAILABLE"
            assert endpoint_cell(desktop, "PLASMA_STABILITY", ".developability-status") == "MODEL_UNAVAILABLE"

            api(f"/api/projects/{project['id']}/compounds/{compound['row_id']}/experimental", "POST", {
                "canonical_endpoint_id": "RLM_CLINT", "raw_endpoint": "RLM intrinsic clearance",
                "raw_value": "1.70", "raw_unit": "log10(mL/min/kg)", "species": "Rat",
                "matrix": "microsomes", "study_id": "PRED-FIRST-E2E", "measurement_type": "intrinsic clearance",
            })
            desktop.refresh()
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']"))
            assert endpoint_cell(desktop, "CACO2_PAPP_AB", ".prediction-value") != "—"
            click_tab(desktop, "METABOLISM")
            wait.until(lambda d: "Metabolic Stability (MS)" in d.page_source)
            for label in ("HLM (microsomal stability)", "RLM (microsomal stability)", "MLM (microsomal stability)", "CYP Inhibition — Quantitative", "Metabolic Soft Spots", "Predicted Metabolites"):
                assert label in desktop.page_source
            wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='RLM_CLINT']"))
            rlm_rows = desktop.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='RLM_CLINT']")
            assert any(row.find_element(By.CSS_SELECTOR, ".prediction-value").text.strip() != "—" and row.find_element(By.CSS_SELECTOR, ".experimental-value").text.strip() != "—" for row in rlm_rows)

            click_tab(desktop, "PK")
            wait.until(lambda d: "Contextual PK" in d.page_source)
            assert "CONTEXT_REQUIRED" in desktop.page_source
            contextual = desktop.find_element(By.XPATH, "//h3[normalize-space()='Contextual PK']")
            desktop.execute_script("arguments[0].scrollIntoView({block:'start'})", contextual)
            desktop.save_screenshot(str(OUT / "desktop_1440x900_pk.png"))

            profile = api(f"/api/compound-versions/{version_id}/developability-profile")
            caco = next(row for row in profile["groups"]["absorption"] if row["query_endpoint"] == "CACO2_PAPP_AB")
            scientific = api(f"/api/compound-versions/{version_id}/scientific-tabs/admet")
            scientific_caco = next(row for row in scientific["rows"] if row["canonical_endpoint"] == "CACO2_PAPP_AB")
            assert caco["prediction"]["snapshot_id"] == scientific_caco["prediction"]["snapshot_id"]

            open_project_compound(desktop, wait, "GLP-1 (small molecule)", "Orforglipron")
            click_tab(desktop, "PK")
            wait.until(lambda d: "Oral bioavailability (F)" in d.page_source)
            assert endpoint_cell(desktop, "HUMAN_PK_F_ORAL", ".experimental-value") != "—"
            results["orforglipron_experimental_visible"] = True
        finally:
            desktop.quit()

        mobile = driver(390, 844)
        mobile_wait = WebDriverWait(mobile, 120)
        try:
            mobile.get(BASE_URL)
            mobile_wait.until(lambda d: d.find_elements(By.CLASS_NAME, "shell"))
            open_project_compound(mobile, mobile_wait, project["name"], compound["name"])
            click_tab(mobile, "ADMET")
            mobile_wait.until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']"))
            overflow = mobile.execute_script("return document.documentElement.scrollWidth-document.documentElement.clientWidth")
            assert overflow <= 1, f"mobile horizontal overflow: {overflow}px"
            row = mobile.find_element(By.CSS_SELECTOR, "tr[data-endpoint='CACO2_PAPP_AB']")
            for label in ("Endpoint", "Prediction", "Experimental", "Status"):
                assert row.find_element(By.CSS_SELECTOR, f"td[data-label='{label}']").is_displayed()
            mobile.execute_script("arguments[0].scrollIntoView({block:'center'})", row)
            mobile.save_screenshot(str(OUT / "mobile_390x844_admet.png"))
            results["mobile_no_horizontal_overflow"] = True
        finally:
            mobile.quit()

    p300 = api("/api/projects/300?page=1&page_size=1")
    reference_version = p300["compounds"][0]["version"]["id"]
    assert api(f"/api/compound-versions/{reference_version}/developability-profile")["compound_version_id"] == reference_version
    results.update({
        "predict_reload_persistence": True,
        "same_snapshot_api_and_ui": True,
        "project_300_reference_profile": True,
        "status": "PASS",
    })
    (OUT / "result.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
