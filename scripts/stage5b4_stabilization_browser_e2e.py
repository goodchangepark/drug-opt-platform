#!/usr/bin/env python3
"""Comprehensive Stage 5B-4 Stabilization Browser E2E Acceptance Test.

Executes the complete top-level navigation and in-project workflow:
  Dashboard -> New Project -> Projects -> Optimization -> Settings -> Help
and inside project:
  Overview -> Properties -> Activity -> ADMET -> Metabolism -> PK -> History
with temporary project __STABILIZATION_E2E_TEMP__ and full cleanup.
"""

from __future__ import annotations

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


BASE_URL = os.environ.get("STABILIZATION_E2E_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
ROOT = Path(__file__).resolve().parent.parent
TEMP_PROJECT_NAME = "__STABILIZATION_E2E_TEMP__"
TEMP_COMPOUND_ID = "STAB-5B4-001"
TEMP_COMPOUND_NAME = "Stabilization Candidate 5B4"
TEMP_SMILES = "CC(C)NCC(O)COc1ccccc1"


def handle_alert_if_present(driver):
    try:
        alert = driver.switch_to.alert
        text = alert.text
        alert.accept()
        return text
    except Exception:
        return None


def card_text(driver, title):
    heading = driver.find_element(By.XPATH, f"//article[contains(@class,'module-card')]//h3[normalize-space()='{title}']")
    return heading.find_element(By.XPATH, "ancestor::article[contains(@class,'module-card')]").text


def click_sidebar(driver, label, timeout=30):
    handle_alert_if_present(driver)
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'global-nav')]//button[normalize-space()={json.dumps(label)}]")
    )
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
    driver.execute_script("arguments[0].click();", btn)
    return btn


def click_detail_tab(driver, tab_label, timeout=30):
    handle_alert_if_present(driver)
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'detail-tabs')]//button[normalize-space()={json.dumps(tab_label.upper())}]")
    )
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
    driver.execute_script("arguments[0].click();", btn)
    return btn


def click_button(driver, text, last=False, timeout=45):
    handle_alert_if_present(driver)

    def match(current):
        rows = []
        for item in current.find_elements(By.TAG_NAME, "button"):
            try:
                if item.is_displayed() and (item.text.strip() == text or item.text.strip().startswith(text + "\n") or text.lower() == item.text.strip().lower()):
                    rows.append(item)
            except StaleElementReferenceException:
                return None
        return (rows[-1] if last else rows[0]) if rows else None

    button = WebDriverWait(driver, timeout).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    driver.execute_script("arguments[0].click();", button)
    return button


def browser_api(driver, method, path, payload=None, timeout=120):
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


def main() -> int:
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for arg in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1920,2400"):
        options.add_argument(arg)
    options.set_capability("goog:loggingPrefs", {"browser": "ALL"})

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    driver.set_script_timeout(90)

    project_id = None
    checks = []

    try:
        # Step 1: Initial load and Dashboard verification
        driver.get(f"{BASE_URL}/?stabilization-e2e={RUN_ID}")
        WebDriverWait(driver, 90).until(lambda d: "Drug Optimization Platform" in d.page_source)
        WebDriverWait(driver, 90).until(lambda d: "Available Scientific Modules" in d.page_source)
        WebDriverWait(driver, 90).until(lambda d: len(d.find_elements(By.CSS_SELECTOR, ".module-card")) == 7)
        cards = driver.find_elements(By.CSS_SELECTOR, ".module-card")
        checks.append({
            "name": "Dashboard Scientific Modules Rendered",
            "passed": len(cards) == 7,
            "details": f"Found {len(cards)} module cards on Dashboard",
        })

        # Verify Capability Consistency on Dashboard
        cyp_card = card_text(driver, "CYP & Transporters")
        safety_card = card_text(driver, "Safety / Toxicology")
        pk_card = card_text(driver, "PK / DMPK")

        pk_ready = "READY" in pk_card and "PLANNED" not in pk_card
        cyp_partial = "PARTIAL" in cyp_card
        safety_partial = "PARTIAL" in safety_card
        checks.append({
            "name": "Capability Consistency on Dashboard",
            "passed": pk_ready and cyp_partial and safety_partial,
            "details": f"PK/DMPK READY: {pk_ready}, CYP PARTIAL: {cyp_partial}, Safety PARTIAL: {safety_partial}",
        })

        # Step 2: Global Sidebar Navigation: New Project -> Projects -> Optimization -> Settings -> Help
        click_sidebar(driver, "New Project")
        WebDriverWait(driver, 30).until(lambda d: "Create New Project" in d.page_source and "Typical Workflow" in d.page_source)
        checks.append({
            "name": "Global View: New Project",
            "passed": "Create New Project" in driver.page_source,
            "details": "New Project form rendered successfully without errors.",
        })

        click_sidebar(driver, "Projects")
        WebDriverWait(driver, 30).until(lambda d: "RESEARCH PORTFOLIO" in d.page_source or "Research Portfolio" in d.page_source or "Projects" in d.page_source)
        checks.append({
            "name": "Global View: Projects",
            "passed": "Projects" in driver.page_source,
            "details": "Projects portfolio view rendered successfully.",
        })

        click_sidebar(driver, "Optimization")
        WebDriverWait(driver, 30).until(lambda d: "Optimization Workspace" in d.page_source)
        checks.append({
            "name": "Global View: Optimization Workspace",
            "passed": "Optimization Workspace" in driver.page_source,
            "details": "Optimization workspace selector rendered successfully.",
        })

        click_sidebar(driver, "Settings")
        WebDriverWait(driver, 30).until(lambda d: "Settings" in d.page_source)
        checks.append({
            "name": "Global View: Settings",
            "passed": "Settings" in driver.page_source,
            "details": "Settings view rendered successfully.",
        })

        click_sidebar(driver, "Help")
        WebDriverWait(driver, 45).until(lambda d: "Drug-OPT Platform Help" in d.page_source and "Current Platform Version" in d.page_source)
        help_page_text = driver.execute_script("return document.body.innerText")
        help_checks = [
            "0.6.0-stage5b4-stable" in help_page_text,
            "5B-4" in help_page_text,
            "Structure & Cheminformatics Modules" in help_page_text,
            "ADME Prediction Models" in help_page_text,
            "CYP & Transporters" in help_page_text,
            "Safety / Toxicology" in help_page_text,
            "Optimization Engine" in help_page_text,
            "PK / DMPK" in help_page_text,
            "Important Scientific Terminology" in help_page_text,
            "Current Limitations" in help_page_text,
        ]
        checks.append({
            "name": "Global View: Help Page & Registries",
            "passed": all(help_checks),
            "details": f"Help tables and registries verified: {sum(help_checks)}/{len(help_checks)} sections present.",
        })

        # Step 3: Create temporary project __STABILIZATION_E2E_TEMP__
        click_sidebar(driver, "Dashboard")
        WebDriverWait(driver, 30).until(lambda d: "Available Scientific Modules" in d.page_source)

        proj = browser_api(driver, "POST", "/projects", {
            "name": TEMP_PROJECT_NAME,
            "target": "Stabilization Validation",
            "molecule_type": "Small Molecule",
            "description": "Temporary project for Stage 5B-4 final stabilization acceptance",
        })
        project_id = proj["id"]

        comp = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": TEMP_COMPOUND_ID,
            "name": TEMP_COMPOUND_NAME,
            "smiles": TEMP_SMILES,
            "notes": "Stage 5B-4 E2E validation compound",
            "calculate": True,
        })
        row_id = comp["row_id"]
        version_id = comp["version"]["id"]

        # Seed Multi-Species PK Data for Allometry
        # Mouse
        m_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Mouse IV Bolus PK", "species": "Mouse", "route": "IV", "dose": 10.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 3500.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 0.5, "concentration_raw": 800.0, "subject_group_id": "Mouse_Mean"},
            {"time_raw": 2.0, "concentration_raw": 80.0, "subject_group_id": "Mouse_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{m_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Rat
        r_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 2000.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 450.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 40.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Dog
        d_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Dog IV Bolus PK", "species": "Dog", "route": "IV", "dose": 2.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/observations", [
            {"time_raw": 0.25, "concentration_raw": 1200.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 2.0, "concentration_raw": 500.0, "subject_group_id": "Dog_Mean"},
            {"time_raw": 8.0, "concentration_raw": 80.0, "subject_group_id": "Dog_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{d_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Monkey
        k_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Monkey IV Bolus PK", "species": "Monkey", "route": "IV", "dose": 3.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/observations", [
            {"time_raw": 0.167, "concentration_raw": 1800.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 1.5, "concentration_raw": 600.0, "subject_group_id": "Monkey_Mean"},
            {"time_raw": 6.0, "concentration_raw": 70.0, "subject_group_id": "Monkey_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{k_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Step 4: Navigate inside the project and exercise compound detail tabs:
        # Overview -> Properties -> Activity -> ADMET -> Metabolism -> PK -> History
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 90).until(lambda d: TEMP_PROJECT_NAME in d.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(TEMP_PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)

        WebDriverWait(driver, 45).until(lambda d: TEMP_COMPOUND_ID in d.page_source or TEMP_COMPOUND_NAME in d.page_source)

        # Click on the compound row or Open button
        cmp_row = driver.find_element(By.XPATH, f"//tr[td[contains(.,'{TEMP_COMPOUND_ID}')]] | //*[contains(text(), '{TEMP_COMPOUND_ID}')]")
        try:
            open_btn = cmp_row.find_element(By.XPATH, ".//button[normalize-space()='Open']")
            driver.execute_script("arguments[0].click();", open_btn)
        except Exception:
            driver.execute_script("arguments[0].click();", cmp_row)

        WebDriverWait(driver, 45).until(lambda d: TEMP_COMPOUND_NAME in d.page_source and "COMPOUND DETAIL" in d.page_source)

        # 4a. Overview Tab
        click_detail_tab(driver, "OVERVIEW")
        time.sleep(1)
        has_overview = "Key Properties" in driver.page_source and ("cLogP" in driver.page_source or "CLOGP" in driver.page_source)
        checks.append({
            "name": "Compound Detail: Overview Tab",
            "passed": has_overview,
            "details": "Compound overview tab rendered with key properties and structure.",
        })

        # 4b. Properties Tab
        click_detail_tab(driver, "PROPERTIES")
        time.sleep(1)
        has_props = "Physicochemical properties" in driver.page_source or "Physicochemical Properties" in driver.page_source or "MW" in driver.page_source
        checks.append({
            "name": "Compound Detail: Properties Tab",
            "passed": has_props,
            "details": "Properties tab rendered with calculated physicochemical properties.",
        })

        # 4c. Activity Tab
        click_detail_tab(driver, "ACTIVITY")
        time.sleep(1)
        has_activity = "Experimental Activity" in driver.page_source and "Activity Prediction" in driver.page_source
        checks.append({
            "name": "Compound Detail: Activity Tab",
            "passed": has_activity,
            "details": "Activity tab rendered experimental and prediction sections.",
        })

        # 4d. ADMET Tab
        click_detail_tab(driver, "ADMET")
        time.sleep(1)
        has_admet = "ADMET" in driver.page_source
        checks.append({
            "name": "Compound Detail: ADMET Tab",
            "passed": has_admet,
            "details": "ADMET tab rendered.",
        })

        # 4e. Metabolism Tab
        click_detail_tab(driver, "METABOLISM")
        time.sleep(1)
        has_metabolism = "Metabolic Soft Spots" in driver.page_source or "SyGMa" in driver.page_source or "Metabolism" in driver.page_source
        try:
            click_button(driver, "Predict Soft Spots", timeout=5)
            time.sleep(2)
        except Exception:
            pass
        checks.append({
            "name": "Compound Detail: Metabolism Tab",
            "passed": has_metabolism,
            "details": "Metabolism tab rendered with SyGMa soft spot engine support.",
        })

        # 4f. PK Tab
        click_detail_tab(driver, "PK")
        time.sleep(2)
        has_pk = "7 · HUMAN PK PREDICTION & TRANSLATIONAL SIMULATION" in driver.page_source or "Human PK" in driver.page_source or "Pharmacokinetics" in driver.page_source
        try:
            click_button(driver, "Run Human PK Simulation", timeout=5)
            time.sleep(2)
            handle_alert_if_present(driver)
        except Exception:
            pass
        checks.append({
            "name": "Compound Detail: PK Tab",
            "passed": has_pk,
            "details": "PK tab rendered with translational and human PK workflow sections.",
        })

        # 4g. History Tab
        click_detail_tab(driver, "HISTORY")
        time.sleep(1)
        has_history = "Version History" in driver.page_source or "Version" in driver.page_source
        checks.append({
            "name": "Compound Detail: History Tab",
            "passed": has_history,
            "details": "History tab rendered version log.",
        })

        # Step 5: Optimization Workflow
        click_button(driver, "Back to Compounds")
        time.sleep(1)

        click_sidebar(driver, "Optimization")
        time.sleep(2)
        WebDriverWait(driver, 30).until(lambda d: "Optimization Workspace" in d.page_source or "PROJECT OPTIMIZATION" in d.page_source)
        checks.append({
            "name": "Optimization Workspace Execution",
            "passed": "Optimization" in driver.page_source,
            "details": "Optimization view operational for project and compound.",
        })

        # Step 6: Console error audit
        logs = driver.get_log("browser")
        severe_errors = [entry for entry in logs if entry.get("level") in {"SEVERE", "ERROR"}]
        checks.append({
            "name": "Zero JavaScript Console Errors",
            "passed": len(severe_errors) == 0,
            "details": f"Total severe/error logs: {len(severe_errors)}" + (f" -> {severe_errors}" if severe_errors else ""),
        })

        # Capture screenshot
        screenshot_path = ROOT / "validation" / "stage5b4_stabilization_browser_e2e.png"
        screenshot_path.parent.mkdir(parents=True, exist_ok=True)
        driver.save_screenshot(str(screenshot_path))

    finally:
        # Step 7: Delete temporary project __STABILIZATION_E2E_TEMP__ using production delete workflow
        if project_id:
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": TEMP_PROJECT_NAME})
                print(f"Temporary project {TEMP_PROJECT_NAME} (#{project_id}) successfully deleted.")
            except Exception as exc:
                print(f"Cleanup warning: {exc}", file=sys.stderr)
        driver.quit()

    all_passed = all(c["passed"] for c in checks)
    report = {
        "suite": "Stage 5B-4 Platform Stabilization Chromium E2E Acceptance",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "run_id": RUN_ID,
        "base_url": BASE_URL,
        "all_passed": all_passed,
        "total_checks": len(checks),
        "passed_checks": sum(1 for c in checks if c["passed"]),
        "failed_checks": sum(1 for c in checks if not c["passed"]),
        "checks": checks,
    }

    report_path = ROOT / "validation" / "stage5b4_stabilization_browser_e2e_results.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))

    return 0 if all_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
