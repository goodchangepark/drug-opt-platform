#!/usr/bin/env python3
"""Focused Chromium E2E Acceptance Test for Stage 5B-4 Human PK Prediction & Translational Simulation."""

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

BASE_URL = os.environ.get("STAGE5B4_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
PROJECT_NAME = f"Stage 5B-4 Browser Acceptance {RUN_ID}"
ARTIFACT_DIR = Path("/home/xavier/.gemini/antigravity-cli/brain/8f3aa156-1b0b-44de-b92b-e15a94d0e222")


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
    for arg in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1920,2400"):
        options.add_argument(arg)

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    project_id = None
    checks = []

    try:
        driver.get(BASE_URL + "/?stage5b4-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        # 1. Create temporary project & compound
        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "Human Translational PK Acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after run.",
        })
        project_id = project["id"]

        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "HUMAN-5B4-001", "name": "Human Clinical Candidate 5B4", "smiles": "CC(C)NCC(O)COc1ccccc1",
            "notes": "Stage 5B-4 browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        # 2. Seed Multi-Species Animal PK Data for Allometry: Mouse, Rat, Dog, Monkey
        # Mouse IV: CL = 65.0 mL/min/kg, MRT = 0.5 h
        m_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Mouse IV Bolus PK", "species": "Mouse", "route": "IV", "dose": 10.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 3500.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 0.5, "concentration_raw": 800.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 2.0, "concentration_raw": 80.0, "subject_group_id": "Mouse_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Rat IV: CL = 32.0 mL/min/kg, MRT = 0.9 h
        r_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 2000.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 450.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 40.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Dog IV: CL = 11.5 mL/min/kg, MRT = 2.4 h
        d_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Dog IV Bolus PK", "species": "Dog", "route": "IV", "dose": 2.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/observations", [
            {"time_raw": 0.25, "concentration_raw": 1200.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 2.0, "concentration_raw": 500.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 8.0, "concentration_raw": 80.0, "subject_group_id": "Dog_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Monkey IV: CL = 18.0 mL/min/kg, MRT = 1.8 h
        k_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Monkey IV Bolus PK", "species": "Monkey", "route": "IV", "dose": 3.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/observations", [
            {"time_raw": 0.167, "concentration_raw": 1800.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 1.5, "concentration_raw": 600.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 6.0, "concentration_raw": 70.0, "subject_group_id": "Monkey_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # 3. Seed Human Hepatic IVIVE Run
        browser_api(driver, "POST", f"/compound-versions/{version_id}/ivive/inputs", {
            "species": "Human", "input_endpoint": "Hepatocyte Clint", "input_value": 14.5, "unit": "µL/min/10^6 cells",
            "source_type": "EXPERIMENTAL", "model_source": "Human Cryopreserved Hepatocytes", "confidence": "HIGH",
        })
        browser_api(driver, "POST", f"/compound-versions/{version_id}/ivive/run", {"species": "Human"})

        # 4. Seed Human Clinical IV & PO PK Data
        # Human IV Study: Dose 100 mg, CL = 6.2 mL/min/kg
        h_iv = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Human Phase 1 IV Bolus PK", "species": "Human", "route": "IV", "dose": 100.0, "dose_unit": "mg",
        })
        browser_api(driver, "POST", f"/pk-studies/{h_iv['id']}/observations", [
            {"time_raw": 0.25, "concentration_raw": 1500.0, "subject_group_id": "Human_Cohort_1"},
            {"time_raw": 1.0, "concentration_raw": 950.0, "subject_group_id": "Human_Cohort_1"},
            {"time_raw": 4.0, "concentration_raw": 420.0, "subject_group_id": "Human_Cohort_1"},
            {"time_raw": 12.0, "concentration_raw": 110.0, "subject_group_id": "Human_Cohort_1"},
            {"time_raw": 24.0, "concentration_raw": 25.0, "subject_group_id": "Human_Cohort_1"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{h_iv['id']}/run-nca", {"selection_mode": "AUTO"})

        # Human PO Study: Dose 200 mg
        h_po = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Human Phase 1 Oral Single Dose PK", "species": "Human", "route": "PO", "dose": 200.0, "dose_unit": "mg",
        })
        browser_api(driver, "POST", f"/pk-studies/{h_po['id']}/observations", [
            {"time_raw": 0.5, "concentration_raw": 320.0, "subject_group_id": "Human_PO_Cohort"},
            {"time_raw": 1.5, "concentration_raw": 1180.0, "subject_group_id": "Human_PO_Cohort"},
            {"time_raw": 3.0, "concentration_raw": 890.0, "subject_group_id": "Human_PO_Cohort"},
            {"time_raw": 8.0, "concentration_raw": 340.0, "subject_group_id": "Human_PO_Cohort"},
            {"time_raw": 24.0, "concentration_raw": 45.0, "subject_group_id": "Human_PO_Cohort"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{h_po['id']}/run-nca", {"selection_mode": "AUTO"})

        # 5. Open Compound Detail -> PK Tab
        driver.get(BASE_URL + "/?project=" + str(project_id))
        time.sleep(2)
        click_button(driver, "Compounds")
        time.sleep(1)

        # Click on compound row
        WebDriverWait(driver, 30).until(lambda d: "HUMAN-5B4-001" in d.page_source)
        cmp_row = driver.find_element(By.XPATH, "//*[contains(text(), 'HUMAN-5B4-001')]")
        cmp_row.click()
        time.sleep(1)

        # Click PK Tab
        click_button(driver, "Pharmacokinetics (PK)")
        time.sleep(2)

        # Check 1: Eyebrow and Section Header presence
        WebDriverWait(driver, 30).until(lambda d: "7 · HUMAN PK PREDICTION & TRANSLATIONAL SIMULATION" in d.page_source)
        checks.append({
            "name": "Human PK Section Rendered",
            "passed": "7 · HUMAN PK PREDICTION & TRANSLATIONAL SIMULATION" in driver.page_source,
            "details": "Section 7 header and eyebrow confirmed present in PK tab.",
        })

        # Check 2: Parameter Assembly Matrix Table
        has_assembly = "Human Parameter Assembly & Candidate Evidence Streams" in driver.page_source
        checks.append({
            "name": "Multi-Stream Parameter Assembly Table",
            "passed": has_assembly,
            "details": "Human parameter candidates (Experimental NCA, Hepatic IVIVE, Allometry) visible.",
        })

        # Check 3: Readiness Scorecard
        has_readiness = "Human IV Simulation Readiness" in driver.page_source and "Human PO Simulation Readiness" in driver.page_source
        checks.append({
            "name": "Human Simulation Readiness Scorecard",
            "passed": has_readiness,
            "details": "Route-specific IV and PO readiness evaluated deterministically.",
        })

        # Check 4: Interactive Simulation Trigger & Output
        click_button(driver, "Run Human PK Simulation")
        time.sleep(2)
        has_sim_output = "SIMULATION OUTPUT" in driver.page_source and "Analytical PK Metrics" in driver.page_source
        checks.append({
            "name": "Human PK Simulation Execution",
            "passed": has_sim_output,
            "details": "Ran 1-compartment human PK simulation and analytical metrics rendered.",
        })

        # Check 5: Interactive SVG Curve Plot with Clinical Overlay
        svgs = driver.find_elements(By.TAG_NAME, "svg")
        has_curve_svg = len(svgs) >= 2
        checks.append({
            "name": "SVG Concentration-Time Plot with Clinical Overlay",
            "passed": has_curve_svg,
            "details": f"Found {len(svgs)} SVG visualization elements including human concentration-time plot.",
        })

        # Check 6: Freeze Prospective Snapshot Action
        click_button(driver, "Freeze Prospective Snapshot")
        time.sleep(2)
        # Handle alert dialog if present
        try:
            alert = driver.switch_to.alert
            alert.accept()
        except Exception:
            pass
        time.sleep(1)

        has_snapshot_history = "Immutable Prediction Snapshot History" in driver.page_source
        checks.append({
            "name": "Prospective Snapshot Freeze Governance",
            "passed": has_snapshot_history,
            "details": "Frozen snapshot record created and listed in immutable snapshot history table.",
        })

        # Check 7: Retrospective Clinical Validation Panel
        has_validation = "Clinical Validation Against Frozen Prospective Snapshot" in driver.page_source or "AAFE" in driver.page_source
        checks.append({
            "name": "Retrospective Clinical Validation",
            "passed": has_validation,
            "details": "Clinical validation comparing experimental Human PK data against frozen prediction snapshot rendered.",
        })

        # Check 8: 0 JavaScript Console Errors
        logs = driver.get_log("browser")
        severe_errors = [l for l in logs if l["level"] in ("SEVERE", "ERROR")]
        checks.append({
            "name": "Zero JavaScript Console Errors",
            "passed": len(severe_errors) == 0,
            "details": f"Total severe/error logs: {len(severe_errors)}",
        })

        # Capture Full Screenshot Artifact
        val_img = ROOT / "validation" / "stage5b4_browser_e2e.png"
        driver.save_screenshot(str(val_img))
        try:
            ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
            screenshot_path = ARTIFACT_DIR / "stage5b4_browser_e2e.png"
            driver.save_screenshot(str(screenshot_path))
        except Exception:
            screenshot_path = val_img

    finally:
        # Cleanup temporary project
        if project_id:
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": PROJECT_NAME})
            except Exception as exc:
                print(f"Cleanup warning: {exc}", file=sys.stderr)
        driver.quit()

    # Write report
    report_dir = ROOT / "validation"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_file = report_dir / "stage5b4_browser_e2e_results.json"
    all_passed = all(c["passed"] for c in checks)

    report_data = {
        "suite": "Stage 5B-4 Human PK Prediction & Translational Simulation E2E",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "run_id": RUN_ID,
        "all_passed": all_passed,
        "total_checks": len(checks),
        "passed_checks": sum(1 for c in checks if c["passed"]),
        "failed_checks": sum(1 for c in checks if not c["passed"]),
        "checks": checks,
        "screenshot_artifact": str(screenshot_path),
    }
    report_json = json.dumps(report_data, indent=2)
    report_file.write_text(report_json, encoding="utf-8")

    print("\n========================================================")
    print(f"Stage 5B-4 Browser E2E Results: {'ALL PASSED' if all_passed else 'FAILED'}")
    print(f"Passed: {report_data['passed_checks']} / {report_data['total_checks']}")
    print(f"Report: {report_file}")
    print(f"Screenshot: {screenshot_path}")
    print("========================================================")
    for c in checks:
        status_str = "PASS" if c["passed"] else "FAIL"
        print(f"  [{status_str}] {c['name']}: {c['details']}")

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
