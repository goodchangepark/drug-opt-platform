#!/usr/bin/env python3
"""Focused Chromium E2E Acceptance Test for Stage 4C Scientific Hardening UI."""

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

BASE_URL = os.environ.get("STAGE4C_BASE_URL", "http://127.0.0.1:8766")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
ROOT = Path(__file__).parents[1]


def click_button(driver, text, timeout=60):
    def match(current):
        for item in current.find_elements(By.TAG_NAME, "button"):
            try:
                if item.is_displayed() and item.text.strip() == text:
                    return item
            except StaleElementReferenceException:
                return None
        return None

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
    checks = []
    result = {"stage": "4C", "run_id": RUN_ID, "base_url": BASE_URL, "checks": checks}

    try:
        driver.get(BASE_URL + "/?stage4c-e2e=" + RUN_ID)
        WebDriverWait(driver, 90).until(lambda d: "Create New Project" in d.page_source)

        # 1. API Endpoint verifications
        gate_res = browser_api(driver, "GET", "/evaluation/golden-gate")
        if not gate_res.get("gate_passed") or gate_res.get("passed_count") != 52:
            raise AssertionError(f"Golden gate failed: {gate_res}")
        checks.append({"name": "Golden Gate API (52 Items 100% Pass)", "status": "PASS"})

        configs = browser_api(driver, "GET", "/standardization/configurations")
        if configs.get("standardizer_name") != "CHEM_STANDARDIZER_V1":
            raise AssertionError(f"Invalid standardizer configs: {configs}")
        checks.append({"name": "Standardizer Configuration API", "status": "PASS"})

        audit = browser_api(driver, "GET", "/evaluation/lightning-audit")
        if audit.get("status") != "SECURE":
            raise AssertionError(f"PyTorch Lightning audit failed: {audit}")
        checks.append({"name": "PyTorch Lightning Security Audit API", "status": "PASS"})

        readiness = browser_api(driver, "GET", "/evaluation/rdkit-readiness")
        if readiness.get("readiness_status") != "READY_FOR_CANDIDATE_TESTING":
            raise AssertionError(f"RDKit readiness check failed: {readiness}")
        checks.append({"name": "RDKit Upgrade Readiness API", "status": "PASS"})

        # 2. Test standardization API directly
        std_res = browser_api(driver, "POST", "/standardization/standardize", {"smiles": "CC(=O)Oc1ccccc1C(=O)O.Cl"})
        if std_res.get("status") != "SUCCESS" or not std_res.get("salt_extracted"):
            raise AssertionError(f"Standardization API failed: {std_res}")
        checks.append({"name": "Standardization API Execution & Salt Removal", "status": "PASS"})

        # 3. UI Verification: Scroll to Scientific Validation
        def check_sci_val_ui(d):
            text = d.execute_script("return document.body ? document.body.innerText : '';")
            return "STAGE 4C SCIENTIFIC HARDENING" in text and "Golden Structure Set Gate" in text

        WebDriverWait(driver, 60).until(check_sci_val_ui)
        checks.append({"name": "Scientific Validation UI Rendering", "status": "PASS"})

        # 4. Click Standardize button in interactive standardizer
        click_button(driver, "Standardize")
        time.sleep(1)

        def check_std_result(d):
            text = d.execute_script("return document.body ? document.body.innerText : '';")
            return "Canonical SMILES:" in text

        WebDriverWait(driver, 60).until(check_std_result)
        checks.append({"name": "Interactive Structure Standardizer Execution", "status": "PASS"})

        # 4. Capture screenshot
        out_dir = ROOT / "validation"
        out_dir.mkdir(parents=True, exist_ok=True)
        img_path = out_dir / "stage4c_browser_e2e.png"
        driver.save_screenshot(str(img_path))
        result["screenshot"] = str(img_path)

        result["status"] = "PASS"
    except Exception as exc:
        result["status"] = "FAIL"
        result["error"] = str(exc)
        out_dir = ROOT / "validation"
        out_dir.mkdir(parents=True, exist_ok=True)
        driver.save_screenshot(str(out_dir / "stage4c_browser_e2e_failure.png"))
    finally:
        driver.quit()
        json_path = ROOT / "validation" / "stage4c_browser_e2e_results.json"
        json_path.write_text(json.dumps(result, indent=2))
        print(json.dumps(result, indent=2))
        if result["status"] != "PASS":
            sys.exit(1)


if __name__ == "__main__":
    main()
