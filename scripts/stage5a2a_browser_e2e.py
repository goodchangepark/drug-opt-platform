#!/usr/bin/env python3
"""Focused Chromium E2E for Stage 5A-2A IVIVE inside Compound Detail > PK."""

import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


BASE_URL = os.environ.get("STAGE5A2A_BASE_URL", "http://127.0.0.1:8766")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
PROJECT_NAME = f"Stage 5A-2A Browser Acceptance {RUN_ID}"
ROOT = Path(__file__).parents[1]


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


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for argument in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1900,1800"):
        options.add_argument(argument)
    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    project_id = None
    checks = []
    result = {"stage": "5A-2A", "run_id": RUN_ID, "base_url": BASE_URL, "checks": checks}
    try:
        driver.get(BASE_URL + "/?stage5a2a-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda current: "Create New Project" in current.page_source)
        health = browser_api(driver, "GET", "/health")
        if health.get("step") != "5A-2A":
            raise AssertionError(f"Unexpected health payload: {health}")

        project = browser_api(driver, "POST", "/projects", {
            "name": PROJECT_NAME, "target": "IVIVE acceptance", "molecule_type": "Small Molecule",
            "description": "Temporary browser acceptance project; deleted after the run.",
        })
        project_id = project["id"]
        compound = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": "IVIVE-E2E-001", "name": "IVIVE-E2E-001", "smiles": "CCOc1ccccc1",
            "notes": "Stage 5A-2A browser fixture", "calculate": True,
        })
        row_id = compound["row_id"]
        version_id = compound["version"]["id"]

        browser_api(driver, "POST", f"/projects/{project_id}/admet/measurements", {
            "version_id": version_id, "endpoint": "RLM intrinsic clearance", "species": "Rat",
            "matrix": "Rat liver microsomes", "value": 10, "unit": "µL/min/mg protein",
            "method": "Microsomal depletion", "source": "E2E experimental HLM/PPB fixture",
        })
        browser_api(driver, "POST", f"/projects/{project_id}/admet/measurements", {
            "version_id": version_id, "endpoint": "Plasma Protein Binding (PPB)", "species": "Rat",
            "matrix": "Plasma", "value": 80, "unit": "% bound",
            "method": "Equilibrium dialysis", "source": "E2E experimental HLM/PPB fixture",
        })
        browser_api(driver, "POST", f"/compound-versions/{version_id}/ivive-inputs", {
            "species": "Rat", "input_endpoint": "BLOOD_PLASMA_RATIO", "input_value": 1.2,
            "unit": "ratio", "source_type": "EXPERIMENTAL", "model_source": "E2E blood partition assay",
            "confidence": "HIGH",
        })

        study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV observed systemic clearance", "species": "Rat", "route": "IV",
            "dose": 0.1, "dose_unit": "mg/kg", "matrix": "Plasma", "source": "E2E IV PK",
        })
        observations = []
        for index, time_value in enumerate((0, 1, 2, 4, 8, 12, 24)):
            observations.append({
                "subject_group_id": "Group Mean", "time_raw": time_value, "time_unit": "h",
                "concentration_raw": 100 * math.exp(-0.15 * time_value), "concentration_unit": "ng/mL",
                "blq_flag": False, "replicate": f"R{index + 1}", "notes": "",
            })
        browser_api(driver, "POST", f"/pk-studies/{study['id']}/observations", observations)
        browser_api(driver, "POST", f"/pk-studies/{study['id']}/run-nca", {})
        checks.append("Seeded exact CompoundVersion with experimental Rat microsomal Clint, PPB, B/P, and IV NCA")

        driver.refresh()
        WebDriverWait(driver, 90).until(lambda current: PROJECT_NAME in current.page_source)
        project_card = WebDriverWait(driver, 45).until(
            lambda current: current.find_element(By.XPATH, f"//article[contains(@class,'dashboard-project')][.//h3[normalize-space()={json.dumps(PROJECT_NAME)}]]")
        )
        driver.execute_script("arguments[0].click();", project_card)
        WebDriverWait(driver, 45).until(lambda current: "IVIVE-E2E-001" in current.page_source)
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'IVIVE-E2E-001')]]")
        driver.execute_script("arguments[0].click();", row.find_element(By.XPATH, ".//button[normalize-space()='Open']"))
        click_button(driver, "PK", last=True)
        WebDriverWait(driver, 60).until(lambda current: current.execute_script("return document.body.innerText.includes('E2E experimental HLM/PPB fixture') && document.body.innerText.includes('DEFAULT PHYSIOLOGY')"))
        body = driver.execute_script("return document.body.innerText")
        for text_value in (
            "Inputs · Clint", "Inputs · PPB / fu,p", "Inputs · Blood/Plasma", "Species Physiology",
            "E2E experimental HLM/PPB fixture", "RAW MICROSOMAL", "DEFAULT PHYSIOLOGY", "EXP",
        ):
            if text_value not in body:
                raise AssertionError(f"IVIVE preview missing {text_value!r}")
        checks.append("Compound Detail > PK > IVIVE preview shows experimental source precedence and Rat physiology")

        click_button(driver, "Run IVIVE")
        WebDriverWait(driver, 60).until(lambda current: current.execute_script("return document.body.innerText.includes('Hepatic Clearance Estimate (CLh)') && document.body.innerText.includes('IVIVE Provenance & Equations')"))
        body = driver.execute_script("return document.body.innerText")
        for text_value in (
            "17.2 mL/min/kg", "Hepatic Blood Flow (Qh)", "67.6", "Hepatic Extraction Ratio (Eh)",
            "Predicted Hepatic Availability (Fh)", "Experimental Comparison", "Predicted hepatic CL",
            "Estimated hepatic contribution", "Assumptions & Warnings", "Non-hepatic Clearance: Not modeled",
            "Predicted Total Clearance: Not generated", "Run confidence: HIGH",
        ):
            if text_value not in body:
                raise AssertionError(f"IVIVE result missing {text_value!r}")
        checks.append("Run IVIVE shows scaled Clint, Qh, fu,b, CLh, Eh, Fh, and no total-clearance claim")
        checks.append("Observed systemic CL comparison and estimated hepatic contribution are distinctly labeled")

        details = driver.find_element(By.XPATH, "//summary[normalize-space()='IVIVE Provenance & Equations']")
        driver.execute_script("arguments[0].click();", details)
        WebDriverWait(driver, 20).until(lambda current: current.execute_script("return document.body.innerText.includes('Input snapshot hash')"))
        body = driver.execute_script("return document.body.innerText")
        for text_value in ("CLh = (Qh * fu,b * Clint)", "compound_version_id", "selection_policy", "classification_outputs_excluded"):
            if text_value not in body:
                raise AssertionError(f"Expanded provenance missing {text_value!r}")
        checks.append("Expanded provenance exposes equation, immutable input snapshot, source policy, and hash")

        screenshot = ROOT / "validation/stage5a2a_browser_e2e.png"
        driver.save_screenshot(str(screenshot))
        result.update({"status": "PASS", "project_id": project_id, "compound_row_id": row_id,
                       "version_id": version_id, "screenshot": str(screenshot.relative_to(ROOT))})
    except Exception as exc:
        failure_screenshot = ROOT / "validation/stage5a2a_browser_e2e_failure.png"
        try:
            driver.save_screenshot(str(failure_screenshot))
            body_text = driver.execute_script("return document.body.innerText")
            console = driver.get_log("browser")
        except Exception as diagnostic_error:
            body_text = f"diagnostic unavailable: {diagnostic_error}"
            console = []
        result.update({"status": "FAIL", "error": str(exc), "current_url": driver.current_url,
                       "body_tail": body_text[-5000:], "browser_console": console,
                       "failure_screenshot": str(failure_screenshot.relative_to(ROOT))})
        raise
    finally:
        if project_id is not None:
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": PROJECT_NAME}, timeout=60)
                result["temporary_project_deleted"] = True
            except Exception as cleanup_error:
                result["temporary_project_deleted"] = False
                result["cleanup_error"] = str(cleanup_error)
        result["timestamp"] = datetime.now(timezone.utc).isoformat()
        (ROOT / "validation/stage5a2a_browser_e2e_results.json").write_text(json.dumps(result, indent=2) + "\n")
        driver.quit()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
