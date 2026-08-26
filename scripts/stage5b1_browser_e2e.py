#!/usr/bin/env python3
"""Focused Chromium E2E Acceptance Test for Stage 5B-1 IV PK Simulation Engine."""

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

BASE_URL = os.environ.get("STAGE5B1_BASE_URL", "http://127.0.0.1:8766")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
PROJECT_NAME = f"Stage 5B-1 Browser Acceptance {RUN_ID}"
ROOT = Path(__file__).parents[1]


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
    result = {"stage": "5B-1", "run_id": RUN_ID, "base_url": BASE_URL, "checks": checks}

    try:
        driver.get(BASE_URL + "/?stage5b1-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        health = browser_api(driver, "GET", "/health")
        if health.get("step") != "5B-1":
            raise AssertionError(f"Unexpected health payload: {health}")
        checks.append({"name": "Health Step 5B-1", "status": "PASS"})

        # 1. Create temporary project & compound
        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "IV Simulation acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after run.",
        })
        project_id = project["id"]

        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "SIM-E2E-001", "name": "SIM-E2E-001", "smiles": "CCOc1ccccc1",
            "notes": "Stage 5B-1 browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        # 2. Seed IV study + observations + NCA
        iv_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 100.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 50.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 12.5, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/run-nca", {"selection_mode": "AUTO"})
        checks.append({"name": "IV NCA Seed", "status": "PASS"})

        # 3. Assemble PK foundation parameter set
        browser_api(driver, "POST", f"/compound-versions/{version_id}/assemble-pk", {"species": "Rat", "route": "IV"})
        checks.append({"name": "PK Parameter Set Assembly", "status": "PASS"})

        # 4. Navigate via UI: Home -> Project -> Compound -> PK tab
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 90).until(lambda d: PROJECT_NAME in d.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)

        WebDriverWait(driver, 45).until(lambda d: "SIM-E2E-001" in d.page_source)
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'SIM-E2E-001')]]")
        driver.execute_script("arguments[0].click();", row.find_element(By.XPATH, ".//button[normalize-space()='Open']"))

        pk_tab = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, "//nav[contains(@class,'detail-tabs')]//button[normalize-space()='PK']")
        )
        driver.execute_script("arguments[0].click();", pk_tab)

        def check_pk_sim(d):
            text = d.execute_script("return document.body ? document.body.innerText : '';")
            return "PK SIMULATION — IV Concentration-Time Engine" in text

        WebDriverWait(driver, 60).until(check_pk_sim)
        checks.append({"name": "PK SIMULATION Section UI Render", "status": "PASS"})

        # 5. Click RUN SIMULATION
        click_button(driver, "RUN SIMULATION")
        time.sleep(1)

        def check_results(d):
            text = d.execute_script("return document.body ? document.body.innerText : '';")
            return "CALCULATED PK SIMULATION" in text and "AUC Numerical Cross-Check" in text

        WebDriverWait(driver, 60).until(check_results)
        checks.append({"name": "IV Bolus Simulation Execution & Results", "status": "PASS"})

        # 6. Capture screenshot
        out_dir = ROOT / "validation"
        out_dir.mkdir(parents=True, exist_ok=True)
        img_path = out_dir / "stage5b1_browser_e2e.png"
        driver.save_screenshot(str(img_path))
        result["screenshot"] = str(img_path)

        # 7. Clean up temporary project via API inside browser
        if project_id:
            browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": PROJECT_NAME})
            result["temporary_project_deleted"] = True
            checks.append({"name": "Cleanup Temp Project", "status": "PASS"})

        result["status"] = "PASS"
    except Exception as exc:
        result["status"] = "FAIL"
        result["error"] = str(exc)
        out_dir = ROOT / "validation"
        out_dir.mkdir(parents=True, exist_ok=True)
        driver.save_screenshot(str(out_dir / "stage5b1_browser_e2e_failure.png"))
    finally:
        driver.quit()
        json_path = ROOT / "validation" / "stage5b1_browser_e2e_results.json"
        json_path.write_text(json.dumps(result, indent=2))
        print(json.dumps(result, indent=2))
        if result["status"] != "PASS":
            sys.exit(1)


if __name__ == "__main__":
    main()
