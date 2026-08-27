#!/usr/bin/env python3
"""Focused Chromium E2E Acceptance Test for Stage 5B-2 Extravascular PK Simulation Engine."""

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

BASE_URL = os.environ.get("STAGE5B2_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
PROJECT_NAME = f"Stage 5B-2 Browser Acceptance {RUN_ID}"


def click_button(driver, text, last=False, timeout=60):
    def match(current):
        rows = []
        for item in current.find_elements(By.TAG_NAME, "button"):
            try:
                if item.is_displayed() and item.text.strip() == text:
                    rows.append(item)
            except StaleElementReferenceException:
                return None
        return (rows[-1] if last else rows[0]) if rows else None

    button = WebDriverWait(driver, timeout).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    driver.execute_script("arguments[0].click();", button)
    return button


def browser_api(driver, method, path, payload=None, timeout=180):
    driver.set_script_timeout(timeout)
    result = driver.execute_async_script(
        """
const done=arguments[arguments.length-1], method=arguments[0], path=arguments[1], payload=arguments[2];
fetch('/api'+path,{method,headers:{'Content-Type':'application/json'},body:payload===null?undefined:JSON.stringify(payload)})
 .then(async response=>{const text=await response.text();if(!response.ok)throw new Error(response.status+' '+text);return text?JSON.parse(text):null})
 .then(data=>done({ok:true,data})).catch(error=>done({ok:false,error:String(error)}));
""",
        method, path, payload,
    )
    if not result["ok"]:
        raise RuntimeError(result["error"])
    return result["data"]


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for arg in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1900,1800"):
        options.add_argument(arg)

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    project_id = None
    checks = []

    try:
        driver.get(BASE_URL + "/?stage5b2-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        # 1. Create temporary project & compound
        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "Extravascular PK Simulation acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after run.",
        })
        project_id = project["id"]

        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "PROP-5B2-E2E", "name": "Propranolol Stage5B2", "smiles": "CC(C)NCC(O)COc1cccc2ccccc12",
            "notes": "Stage 5B-2 browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        # 2. Seed IV study + observations + NCA
        iv_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 2200.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 0.5, "concentration_raw": 950.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 480.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 2.0, "concentration_raw": 160.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 30.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # 3. Seed PO study + observations + NCA
        po_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat PO Gavage PK", "species": "Rat", "route": "PO", "dose": 20.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{po_study['id']}/observations", [
            {"time_raw": 0.25, "concentration_raw": 120.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 0.5, "concentration_raw": 280.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 450.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 2.0, "concentration_raw": 320.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 110.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 6.0, "concentration_raw": 25.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{po_study['id']}/run-nca", {"selection_mode": "AUTO"})
        checks.append({"name": "IV and PO PK NCA Seeding", "status": "PASS"})

        # 4. Assemble PK foundation parameter sets
        browser_api(driver, "POST", f"/compound-versions/{version_id}/assemble-pk", {"species": "Rat", "route": "IV"})
        browser_api(driver, "POST", f"/compound-versions/{version_id}/assemble-pk", {"species": "Rat", "route": "PO"})
        checks.append({"name": "PK Parameter Set Assembly", "status": "PASS"})

        # 5. Navigate via UI: Home -> Project -> Compound -> PK tab
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 90).until(lambda d: PROJECT_NAME in d.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)

        WebDriverWait(driver, 45).until(lambda d: "PROP-5B2-E2E" in d.page_source)
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'PROP-5B2-E2E')]]")
        driver.execute_script("arguments[0].click();", row.find_element(By.XPATH, ".//button[normalize-space()='Open']"))

        WebDriverWait(driver, 45).until(lambda d: "PK" in d.page_source)
        click_button(driver, "PK")
        time.sleep(1.0)
        checks.append({"name": "Compound Detail & PK Tab Navigation", "status": "PASS"})

        # 6. Verify PK Simulation Section & Route Switcher
        WebDriverWait(driver, 45).until(lambda d: "PK SIMULATION — Extravascular" in d.page_source)
        checks.append({"name": "Extravascular PK Simulation DOM Render", "status": "PASS"})

        # Verify Parameter Provenance
        WebDriverWait(driver, 45).until(lambda d: "Parameter Provenance & Evidence Hierarchy" in d.page_source)
        checks.append({"name": "PO Route Parameter Hierarchy Verification", "status": "PASS"})

        # 7. Run Single Dose PO Simulation
        click_button(driver, "RUN PK SIMULATION")
        time.sleep(2.0)

        WebDriverWait(driver, 45).until(lambda d: "CALCULATED PK SIMULATION: PO" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Cmax" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Tmax" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "AUCinf (Analytical)" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Observed vs Simulated Residual Analysis" in d.page_source)
        checks.append({"name": "PO Single Dose Simulation & Experimental Overlay", "status": "PASS"})

        # 8. Test Route Switching to IV, SC, IP
        click_button(driver, "SC (Subcutaneous)")
        time.sleep(1.0)
        WebDriverWait(driver, 45).until(lambda d: "Parameter Provenance & Evidence Hierarchy (SC" in d.page_source)

        click_button(driver, "IP (Intraperitoneal)")
        time.sleep(1.0)
        WebDriverWait(driver, 45).until(lambda d: "Parameter Provenance & Evidence Hierarchy (IP" in d.page_source)

        click_button(driver, "IV (Intravenous)")
        time.sleep(1.0)
        WebDriverWait(driver, 45).until(lambda d: "Parameter Provenance & Evidence Hierarchy (IV" in d.page_source)
        checks.append({"name": "Cross-Route Switcher (PO, SC, IP, IV)", "status": "PASS"})

        # Switch back to PO and test Repeated Dosing
        click_button(driver, "PO (Oral)")
        time.sleep(0.5)

        freq_select = driver.find_element(By.XPATH, "//label[text()='Dosing Frequency']/following-sibling::select")
        for opt in freq_select.find_elements(By.TAG_NAME, "option"):
            if opt.text == "Repeated Dosing":
                opt.click()
                break
        time.sleep(0.5)

        click_button(driver, "RUN PK SIMULATION")
        time.sleep(2.0)

        WebDriverWait(driver, 45).until(lambda d: "Accumulation Ratio (R_acc)" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "LINEAR PK ASSUMPTION" in d.page_source)
        checks.append({"name": "Repeated Dosing Superposition & Accumulation Render", "status": "PASS"})

        # 9. Verify 0 Uncaught JS errors in browser console
        logs = driver.get_log("browser")
        js_errors = [l for l in logs if l.get("level") == "SEVERE" and any(err in l.get("message", "") for err in ("Uncaught", "TypeError", "ReferenceError", "SyntaxError", "RangeError"))]
        if js_errors:
            checks.append({"name": "Zero Uncaught JS Errors Check", "status": f"FAIL: {js_errors}"})
        else:
            checks.append({"name": "Zero Uncaught JS Errors Check", "status": "PASS"})

        # Save Artifacts
        val_dir = Path("/home/xavier/chem/drug-opt-platform/validation")
        val_dir.mkdir(exist_ok=True)
        screenshot_path = val_dir / "stage5b2_browser_e2e.png"
        png_bytes = driver.get_screenshot_as_png()
        with open(str(screenshot_path), "wb") as f:
            f.write(png_bytes)
        print(f"Screenshot written directly: {screenshot_path}, size={len(png_bytes)}")

        results_path = val_dir / "stage5b2_browser_e2e_results.json"
        summary = {
            "stage": "5B-2",
            "run_id": RUN_ID,
            "base_url": BASE_URL,
            "checks": checks,
            "screenshot": str(screenshot_path.resolve()),
            "status": "PASS" if all(c["status"] == "PASS" for c in checks) else "FAIL"
        }
        with open(str(results_path), "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Results JSON saved: {results_path}")
        print(json.dumps(summary, indent=2))

    finally:
        if project_id:
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": PROJECT_NAME})
            except Exception as e:
                print(f"Cleanup error: {e}")
        driver.quit()


if __name__ == "__main__":
    main()
