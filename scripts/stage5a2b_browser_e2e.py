#!/usr/bin/env python3
"""Focused Chromium E2E for Stage 5A-2B PK Parameter Foundation & Route Assembly."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

BASE_URL = os.environ.get("STAGE5A2B_BASE_URL", "http://127.0.0.1:8766")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
PROJECT_NAME = f"Stage 5A-2B Browser Acceptance {RUN_ID}"
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
    result = {"stage": "5A-2B", "run_id": RUN_ID, "base_url": BASE_URL, "checks": checks}

    try:
        driver.get(BASE_URL + "/?stage5a2b-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        health = browser_api(driver, "GET", "/health")
        if health.get("step") != "5A-2B":
            raise AssertionError(f"Unexpected health payload: {health}")

        # 1. Create temporary project & compound
        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "PK Foundation acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after the run.",
        })
        project_id = project["id"]

        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "PKF-E2E-001", "name": "PKF-E2E-001", "smiles": "CCOc1ccccc1",
            "notes": "Stage 5A-2B browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        # 2. Seed IV and PO studies
        iv_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus PK", "species": "Rat", "route": "IV", "dose": 5.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 100.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 1.0, "concentration_raw": 50.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 4.0, "concentration_raw": 12.5, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{iv_study['id']}/run-nca", {"selection_mode": "AUTO"})

        po_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat PO PK", "species": "Rat", "route": "PO", "dose": 10.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{po_study['id']}/observations", [
            {"time_raw": 0.5, "concentration_raw": 30.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 2.0, "concentration_raw": 60.0, "subject_group_id": "Rat_Mean"},
            {"time_raw": 8.0, "concentration_raw": 10.0, "subject_group_id": "Rat_Mean"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{po_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # 3. Seed ADMET / IVIVE measurements (Microsomal Clint, PPB, Caco-2)
        browser_api(driver, "POST", f"/projects/{project_id}/admet/measurements", {
            "version_id": version_id, "endpoint": "RLM intrinsic clearance", "species": "Rat",
            "matrix": "Rat liver microsomes", "value": 15, "unit": "µL/min/mg protein",
            "method": "Microsomal depletion", "source": "E2E fixture",
        })
        browser_api(driver, "POST", f"/projects/{project_id}/admet/measurements", {
            "version_id": version_id, "endpoint": "Plasma Protein Binding (PPB)", "species": "Rat",
            "matrix": "Plasma", "value": 85, "unit": "% bound",
            "method": "Equilibrium dialysis", "source": "E2E fixture",
        })
        browser_api(driver, "POST", f"/projects/{project_id}/admet/measurements", {
            "version_id": version_id, "endpoint": "Caco-2 Permeability", "species": "Human",
            "matrix": "Caco-2 cell monolayer", "value": 8.5, "unit": "10^-6 cm/s",
            "method": "A-to-B permeability", "source": "E2E fixture",
        })

        # Run IVIVE
        browser_api(driver, "POST", f"/compound-versions/{version_id}/ivive/run", {
            "species": "Rat", "method_key": "WELL_STIRRED",
        })
        checks.append("Seeded IV/PO PK studies, ADMET measurements, and IVIVE run.")

        # 4. Fetch PK Foundation via API & verify
        foundation = browser_api(driver, "GET", f"/compound-versions/{version_id}/pk-foundation?species=Rat")
        if not foundation.get("route_parameter_sets"):
            raise AssertionError(f"Invalid foundation response: {foundation}")
        checks.append("Verified PK Foundation API endpoint returns route-aware parameter sets.")

        # 5. Open Compound Detail > PK tab via UI navigation
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 90).until(lambda d: PROJECT_NAME in d.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)

        WebDriverWait(driver, 45).until(lambda d: "PKF-E2E-001" in d.page_source)
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'PKF-E2E-001')]]")
        driver.execute_script("arguments[0].click();", row.find_element(By.XPATH, ".//button[normalize-space()='Open']"))

        pk_tab = WebDriverWait(driver, 45).until(
            lambda d: d.find_element(By.XPATH, "//nav[contains(@class,'detail-tabs')]//button[normalize-space()='PK']")
        )
        driver.execute_script("arguments[0].click();", pk_tab)

        def check_pk_foundation(d):
            text = d.execute_script("return document.body ? document.body.innerText : '';")
            return "PK Parameter Foundation" in text or "5 · PK PARAMETER FOUNDATION" in text

        WebDriverWait(driver, 60).until(check_pk_foundation)
        checks.append("Compound Detail > PK rendered PK Parameter Foundation & Route Assembly section.")

        out_path = ROOT / "validation" / "stage5a2b_browser_e2e.png"
        out_path.parent.mkdir(exist_ok=True)
        driver.save_screenshot(str(out_path))
        result["screenshot"] = str(out_path)

        if project_id:
            browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": PROJECT_NAME})
            result["temporary_project_deleted"] = True
            checks.append("Cleaned up temporary E2E project.")

        result["status"] = "PASS"
        result["timestamp"] = datetime.now(timezone.utc).isoformat()
        res_file = ROOT / "validation" / "stage5a2b_browser_e2e_results.json"
        res_file.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result, indent=2))

    except Exception as exc:
        result["status"] = "FAIL"
        result["error"] = str(exc)
        try:
            result["body_text"] = driver.execute_script("return document.body ? document.body.innerText : '';")
            result["browser_logs"] = driver.get_log("browser")
        except Exception:
            pass
        fail_path = ROOT / "validation" / "stage5a2b_browser_e2e_failure.png"
        fail_path.parent.mkdir(exist_ok=True)
        driver.save_screenshot(str(fail_path))
        result["screenshot"] = str(fail_path)
        print(json.dumps(result, indent=2))
        raise
    finally:
        driver.quit()


if __name__ == "__main__":
    main()
