#!/usr/bin/env python3
"""Focused Chromium E2E Acceptance Test for Stage 5B-3 PK Validation, Cross-Species Scaling & Translational Foundation."""

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

BASE_URL = os.environ.get("STAGE5B3_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
PROJECT_NAME = f"Stage 5B-3 Browser Acceptance {RUN_ID}"


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
        driver.get(BASE_URL + "/?stage5b3-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        # 1. Create temporary project & compound
        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "Translational PK Acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after run.",
        })
        project_id = project["id"]

        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "TRANS-5B3-001", "name": "Translational Lead 5B3", "smiles": "CC(C)NCC(O)COc1ccccc1",
            "notes": "Stage 5B-3 browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        # 2. Seed Multi-Species PK Data: Mouse, Rat, Dog, Monkey
        # Mouse IV
        m_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Mouse IV Bolus PK", "species": "Mouse", "route": "IV", "dose": 10.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 3500.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 0.5, "concentration_raw": 800.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 2.0, "concentration_raw": 80.0, "subject_group_id": "Mouse_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Rat IV
        r_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 2000.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 450.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 40.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Dog IV
        d_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Dog IV Bolus PK", "species": "Dog", "route": "IV", "dose": 2.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/observations", [
            {"time_raw": 0.25, "concentration_raw": 1200.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 2.0, "concentration_raw": 500.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 8.0, "concentration_raw": 80.0, "subject_group_id": "Dog_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Monkey IV
        k_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Monkey IV Bolus PK", "species": "Monkey", "route": "IV", "dose": 3.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/observations", [
            {"time_raw": 0.167, "concentration_raw": 1800.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 1.5, "concentration_raw": 550.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 6.0, "concentration_raw": 60.0, "subject_group_id": "Monkey_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Rat PO
        r_po_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat PO Gavage PK", "species": "Rat", "route": "PO", "dose": 20.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{r_po_study['id']}/observations", [
            {"time_raw": 0.5, "concentration_raw": 300.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 480.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 120.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{r_po_study['id']}/run-nca", {"selection_mode": "AUTO"})

        checks.append({"name": "Multi-Species PK Seeding (4 Species)", "status": "PASS"})

        # 3. Assemble PK Foundation and trigger Translational Profile API
        browser_api(driver, "GET", f"/compound-versions/{version_id}/translational-pk")
        checks.append({"name": "Translational PK Profile API", "status": "PASS"})

        # 4. Navigate via UI: Home -> Project -> Compound -> PK tab
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 90).until(lambda d: PROJECT_NAME in d.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)

        WebDriverWait(driver, 45).until(lambda d: "TRANS-5B3-001" in d.page_source)
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'TRANS-5B3-001')]]")
        driver.execute_script("arguments[0].click();", row.find_element(By.XPATH, ".//button[normalize-space()='Open']"))

        WebDriverWait(driver, 45).until(lambda d: "PK" in d.page_source)
        click_button(driver, "PK")
        time.sleep(1.5)
        checks.append({"name": "PK Tab Navigation", "status": "PASS"})

        # 5. Verify Section 5: Translational PK & Allometry DOM Elements
        WebDriverWait(driver, 45).until(lambda d: "TRANSLATIONAL PK & CROSS-SPECIES ALLOMETRIC SCALING" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Cross-Species In Vivo PK Observation Matrix" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Clearance (CL) Allometry (Animal IV Data)" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Volume of Distribution (Vss) Allometry" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Human PK Translational Comparison & Readiness Scorecard" in d.page_source)
        checks.append({"name": "Translational PK Section DOM Render", "status": "PASS"})

        # 6. Verify Allometry Exponent & Human Extrapolation Render
        WebDriverWait(driver, 45).until(lambda d: "Exponent b (CL)" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Human Extrapolated CL" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Leave-One-Species-Out (LOSO)" in d.page_source)
        checks.append({"name": "Allometry Parameters & LOSO Validation Render", "status": "PASS"})

        # 7. Verify Section 6: PK Validation DOM Elements
        WebDriverWait(driver, 45).until(lambda d: "PK VALIDATION & PREDICTION ERROR METRICS" in d.page_source)
        WebDriverWait(driver, 45).until(lambda d: "Predicted vs Observed" in d.page_source)
        checks.append({"name": "PK Validation Section DOM Render", "status": "PASS"})

        # 8. Verify Zero Uncaught JS Errors in Browser Console
        logs = driver.get_log("browser")
        js_errors = [l for l in logs if l.get("level") == "SEVERE" and any(err in l.get("message", "") for err in ("Uncaught", "TypeError", "ReferenceError", "SyntaxError", "RangeError"))]
        if js_errors:
            checks.append({"name": "Zero Uncaught JS Errors Check", "status": f"FAIL: {js_errors}"})
        else:
            checks.append({"name": "Zero Uncaught JS Errors Check", "status": "PASS"})

        # Save Artifacts directly
        val_dir = Path("/home/xavier/chem/drug-opt-platform/validation")
        val_dir.mkdir(exist_ok=True, parents=True)
        screenshot_path = val_dir / "stage5b3_browser_e2e.png"
        artifact_path = Path("/home/xavier/.gemini/antigravity-cli/brain/8f3aa156-1b0b-44de-b92b-e15a94d0e222/stage5b3_browser_e2e.png")
        artifact_path.parent.mkdir(exist_ok=True, parents=True)

        png_bytes = driver.get_screenshot_as_png()
        with open(str(screenshot_path), "wb") as f:
            f.write(png_bytes)
        with open(str(artifact_path), "wb") as f:
            f.write(png_bytes)
        print(f"Screenshot written: {screenshot_path} and {artifact_path}, size={len(png_bytes)}")

        results_path = val_dir / "stage5b3_browser_e2e_results.json"
        summary = {
            "stage": "5B-3",
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
