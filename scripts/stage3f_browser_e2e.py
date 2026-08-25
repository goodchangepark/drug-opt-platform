#!/usr/bin/env python3
"""Actual Chromium Stage 3F workflow against a running project-local server."""
import json
import os
import sys
import time
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parents[1]
BASE_URL = os.environ.get("STAGE3F_BASE_URL", "http://127.0.0.1:8766")
PROJECT_NAME = "Stage 3 Acceptance — Public References"


def browser_api(driver, method, path, payload=None, timeout=180):
    driver.set_script_timeout(timeout)
    script = """
const done=arguments[arguments.length-1], method=arguments[0], path=arguments[1], payload=arguments[2];
fetch('/api'+path,{method,headers:{'Content-Type':'application/json'},body:payload===null?undefined:JSON.stringify(payload)})
 .then(async response=>{const text=await response.text();if(!response.ok)throw new Error(response.status+' '+text);return text?JSON.parse(text):null})
 .then(data=>done({ok:true,data})).catch(error=>done({ok:false,error:String(error)}));
"""
    result = driver.execute_async_script(script, method, path, payload)
    if not result["ok"]:
        raise RuntimeError(result["error"])
    return result["data"]


def click_text(driver, text, last=False):
    def match(d):
        rows = [e for e in d.find_elements(By.TAG_NAME, "button") if (e.text.strip() == text or e.text.strip().startswith(text + "\n")) and e.is_displayed()]
        return (rows[-1] if last else rows[0]) if rows else None
    element = WebDriverWait(driver, 30).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", element)
    element.click()


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for argument in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1600,1200"):
        options.add_argument(argument)
    driver = webdriver.Chrome(
        service=Service("/snap/bin/chromium.chromedriver"),
        options=options,
    )
    checks = []
    try:
        driver.get(BASE_URL + "/?stage3f-acceptance=20260826")
        WebDriverWait(driver, 30).until(lambda d: "AI Drug Optimization Platform" in d.page_source)
        projects = browser_api(driver, "GET", "/projects")
        project = next((row for row in projects if row["name"] == PROJECT_NAME), None)
        if not project:
            project = browser_api(driver, "POST", "/projects", {"name": PROJECT_NAME, "target": "Public pipeline acceptance", "description": "Public reference compounds only"})
        detail = browser_api(driver, "GET", f"/projects/{project['id']}")
        compounds = {row["compound_id"]: row for row in detail["compounds"]}
        references = [
            ("REF-HERG-001", "Dofetilide", "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"),
            ("REF-HERG-002", "Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
        ]
        for compound_id, name, smiles in references:
            if compound_id not in compounds:
                browser_api(driver, "POST", f"/projects/{project['id']}/compounds", {"compound_id": compound_id, "name": name, "smiles": smiles, "notes": "Public Stage 3 acceptance reference"})
        detail = browser_api(driver, "GET", f"/projects/{project['id']}")
        compounds = {row["compound_id"]: row for row in detail["compounds"]}
        version_id = compounds["REF-HERG-001"]["version"]["id"]
        admet = browser_api(driver, "GET", f"/projects/{project['id']}/admet")
        if not any(row["version_id"] == version_id and next((ep["name"] for ep in admet["endpoints"] if ep["id"] == row["endpoint_id"]), "") == "hERG liability" for row in admet["measurements"]):
            browser_api(driver, "POST", f"/projects/{project['id']}/admet/measurements", {"version_id": version_id, "endpoint": "hERG liability", "species": "Human", "matrix": "public reference classification", "value": 1, "unit": "class", "method": "literature annotation", "source": "PubChem/public literature"})
        browser_api(driver, "POST", f"/admet/predict/{version_id}", {})
        browser_api(driver, "POST", f"/metabolism/predict/{version_id}", {})

        driver.get(BASE_URL + "/?stage3f-acceptance=20260826")
        WebDriverWait(driver, 30).until(lambda d: PROJECT_NAME in d.page_source)
        click_text(driver, PROJECT_NAME)
        WebDriverWait(driver, 30).until(lambda d: "REF-HERG-001" in d.page_source)
        checks.append("Project → Compound")
        open_buttons = driver.find_elements(By.XPATH, "//tr[td[contains(.,'REF-HERG-001')]]//button[normalize-space()='Open']")
        open_buttons[0].click()
        WebDriverWait(driver, 30).until(lambda d: "Physicochemical properties" in d.page_source)
        checks.append("Properties")
        click_text(driver, "ADMET", last=True)
        try:
            WebDriverWait(driver, 60).until(lambda d: "Stage 3 Integrated ADMET Profile" in d.page_source)
        except Exception as exc:
            body = driver.find_element(By.TAG_NAME, "body").text
            raise AssertionError("Compound ADMET tab did not render. Visible page: " + body[-5000:]) from exc
        required = ["Absorption", "Aqueous Solubility / LogS", "Caco-2", "Distribution", "Human PPB / fu", "HLM / RLM / MLM", "Metabolism · CYP", "Metabolic Soft Spots", "Predicted Metabolites", "Transporters", "Safety", "hERG", "Ames", "DILI", "Experimental", "Predicted", "Details", "Provenance audit: PASS"]
        missing = [label for label in required if label not in driver.page_source]
        if missing:
            raise AssertionError("Missing Compound ADMET workflow labels: " + ", ".join(missing))
        checks.extend(required)
        safety_details = driver.find_elements(By.XPATH, "//tr[td[contains(.,'hERG')]]//summary[normalize-space()='Details']")
        if not safety_details:
            raise AssertionError("Safety Details control was not rendered")
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", safety_details[0])
        safety_details[0].click()
        if "License:" not in driver.page_source or "Independent validation:" not in driver.page_source:
            raise AssertionError("Safety Details provenance did not expose license/independent validation")
        checks.append("Details/provenance")
        click_text(driver, "Close")
        boxes = driver.find_elements(By.CSS_SELECTOR, "table input[type='checkbox']")[:2]
        for box in boxes:
            if not box.is_selected():
                box.click()
        click_text(driver, "Compare selected")
        WebDriverWait(driver, 30).until(lambda d: "No overall score or ranking is calculated" in d.page_source)
        checks.append("Comparison")
        workbench = driver.find_element(By.LINK_TEXT, "Open workbench")
        driver.execute_script("arguments[0].click();", workbench)
        WebDriverWait(driver, 30).until(lambda d: "Activity" in d.page_source or "Assay" in d.page_source)
        checks.append("Activity workbench")
        print(json.dumps({"status": "PASS", "base_url": BASE_URL, "project": PROJECT_NAME, "checks": checks}, indent=2))
    finally:
        driver.quit()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, indent=2), file=sys.stderr)
        raise
