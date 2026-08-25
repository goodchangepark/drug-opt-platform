#!/usr/bin/env python3
"""Actual Chromium Stage 4A strategy workflow against the production service."""

import json
import os
import sys

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


BASE_URL = os.environ.get("STAGE4A_BASE_URL", "http://127.0.0.1:8765")
PROJECT_NAME = "Stage 4A Acceptance — Public References"


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


def click_button(driver, text, last=False):
    def match(d):
        matches = [item for item in d.find_elements(By.TAG_NAME, "button") if item.is_displayed() and (item.text.strip() == text or item.text.strip().startswith(text + "\n"))]
        return (matches[-1] if last else matches[0]) if matches else None
    button = WebDriverWait(driver, 30).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    button.click()
    return button


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for argument in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1700,1400"):
        options.add_argument(argument)
    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    checks = []
    try:
        driver.get(BASE_URL + "/?stage4a-e2e=20260826")
        WebDriverWait(driver, 30).until(lambda d: "AI Drug Optimization Platform" in d.page_source)
        health = browser_api(driver, "GET", "/health")
        if health.get("step") != "4A":
            raise AssertionError(f"Unexpected health payload: {health}")

        projects = browser_api(driver, "GET", "/projects")
        project = next((row for row in projects if row["name"] == PROJECT_NAME), None)
        if not project:
            project = browser_api(driver, "POST", "/projects", {
                "name": PROJECT_NAME, "target": "Directional strategy sanity",
                "description": "Public reference structures only; no proprietary compounds",
            })
        project_detail = browser_api(driver, "GET", f"/projects/{project['id']}")
        parent = next((row for row in project_detail["compounds"] if row["compound_id"] == "OPT-LIDOCAINE"), None)
        if not parent:
            browser_api(driver, "POST", f"/projects/{project['id']}/compounds", {
                "compound_id": "OPT-LIDOCAINE", "name": "Lidocaine public direction example",
                "smiles": "CCN(CC)C(=O)c1c(C)cccc1C", "notes": "Public Stage 4A acceptance structure",
            })
            project_detail = browser_api(driver, "GET", f"/projects/{project['id']}")
            parent = next(row for row in project_detail["compounds"] if row["compound_id"] == "OPT-LIDOCAINE")
        version_id = parent["version"]["id"]
        admet = browser_api(driver, "GET", f"/projects/{project['id']}/admet")
        endpoint_by_id = {row["id"]: row["name"] for row in admet["endpoints"]}
        if not any(row["version_id"] == version_id and endpoint_by_id.get(row["endpoint_id"]) == "HLM intrinsic clearance" for row in admet["measurements"]):
            browser_api(driver, "POST", f"/projects/{project['id']}/admet/measurements", {
                "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human",
                "matrix": "HLM", "value": 2.2, "unit": "log10(mL/min/kg)",
                "method": "directional acceptance fixture", "source": "Public reference direction",
            })

        driver.get(BASE_URL + "/?stage4a-e2e=20260826")
        WebDriverWait(driver, 30).until(lambda d: PROJECT_NAME in d.page_source)
        click_button(driver, PROJECT_NAME)
        WebDriverWait(driver, 30).until(lambda d: "OPT-LIDOCAINE" in d.page_source)
        open_button = driver.find_element(By.XPATH, "//tr[td[contains(.,'OPT-LIDOCAINE')]]//button[normalize-space()='Open']")
        open_button.click()
        WebDriverWait(driver, 30).until(lambda d: "Physicochemical properties" in d.page_source)
        checks.extend(["Project", "Parent CompoundVersion", "Properties"])

        click_button(driver, "OPTIMIZATION", last=True)
        WebDriverWait(driver, 30).until(lambda d: "Optimization Run" in d.page_source and "Stage 4A deterministic strategy only" in d.page_source)
        checks.extend(["Optimization tab", "Objective selection", "Constraints"])
        click_button(driver, "Analyze strategy")
        WebDriverWait(driver, 45).until(lambda d: "Recommended transformations" in d.page_source and "HLM metabolic instability" in d.page_source)
        required = ["Current profile", "Main liabilities", "Protected regions", "Modifiable regions", "Recommended transformations", "Experimental", "METABOLIC", "no analog structures"]
        missing = [label for label in required if label.lower() not in driver.page_source.lower()]
        if missing:
            raise AssertionError("Missing Optimization result labels: " + ", ".join(missing))
        checks.extend(required)
        details = driver.find_elements(By.XPATH, "//table//summary[normalize-space()='Details']")
        if not details:
            raise AssertionError("Transformation Details/provenance was not rendered")
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", details[-1]); details[-1].click()
        if "Reaction SMARTS:" not in driver.page_source or "Source/reference:" not in driver.page_source:
            raise AssertionError("Transformation provenance did not render")
        checks.append("Transformation provenance")

        prioritize = driver.find_elements(By.XPATH, "//button[normalize-space()='Prioritize']")
        if not prioritize:
            raise AssertionError("Manual priority action was not rendered")
        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", prioritize[0]); prioritize[0].click()
        WebDriverWait(driver, 30).until(lambda d: "Manual override saved and strategy reranked" in d.page_source)
        runs = browser_api(driver, "GET", f"/projects/{project['id']}/optimization?version_id={version_id}")["runs"]
        if not runs or not runs[0]["manual_overrides"].get("prioritize_transformations"):
            raise AssertionError("Manual priority override was not persisted")
        if runs[0]["analog_generation"] != "NOT_PERFORMED":
            raise AssertionError("Stage 4A unexpectedly generated analogs")
        checks.extend(["Manual override", "Deterministic rerank", "No analog generation"])
        print(json.dumps({"status": "PASS", "base_url": BASE_URL, "project": PROJECT_NAME, "checks": checks}, indent=2))
    finally:
        driver.quit()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, indent=2), file=sys.stderr)
        raise
