#!/usr/bin/env python3
"""Chromium E2E for the compound-centred Pre-Stage 5 UI workflow."""

import json
import os
import time
from datetime import datetime, timezone

from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


BASE_URL = os.environ.get("PRE_STAGE5_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
PROJECT_NAME = f"Pre-Stage 5 UI Acceptance {RUN_ID}"


def browser_api(driver, method, path, payload=None, timeout=360):
    driver.set_script_timeout(timeout)
    result = driver.execute_async_script(
        """
const done=arguments[arguments.length-1], method=arguments[0], path=arguments[1], payload=arguments[2];
fetch('/api'+path,{method,headers:{'Content-Type':'application/json'},body:payload===null?undefined:JSON.stringify(payload)})
 .then(async response=>{const text=await response.text();if(!response.ok)throw new Error(response.status+' '+text);return text?JSON.parse(text):null})
 .then(data=>done({ok:true,data})).catch(error=>done({ok:false,error:String(error)}));
""",
        method,
        path,
        payload,
    )
    if not result["ok"]:
        raise RuntimeError(result["error"])
    return result["data"]


def click_button(driver, text, last=False, timeout=60):
    def match(d):
        matches = []
        for item in d.find_elements(By.TAG_NAME, "button"):
            try:
                item_text = item.text.strip()
                if item.is_displayed() and (item_text == text or item_text.startswith(text + "\n")):
                    matches.append(item)
            except StaleElementReferenceException:
                return None
        return (matches[-1] if last else matches[0]) if matches else None

    button = WebDriverWait(driver, timeout).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    button.click()
    return button


def set_labeled_control(driver, label, value):
    result = driver.execute_script(
        """
const wanted=arguments[0], value=String(arguments[1]);
const label=[...document.querySelectorAll('label')].find(row=>row.offsetParent!==null && row.textContent.trim()===wanted);
if(!label)return {ok:false,reason:'label missing',labels:[...document.querySelectorAll('label')].filter(row=>row.offsetParent!==null).map(row=>row.textContent.trim())};
const control=label.nextElementSibling;
if(!control || !['INPUT','SELECT','TEXTAREA'].includes(control.tagName))return {ok:false,reason:'control missing',tag:control&&control.tagName};
control.scrollIntoView({block:'center'});
const proto=control.tagName==='SELECT'?HTMLSelectElement.prototype:(control.tagName==='TEXTAREA'?HTMLTextAreaElement.prototype:HTMLInputElement.prototype);
Object.getOwnPropertyDescriptor(proto,'value').set.call(control,value);
control.dispatchEvent(new Event('input',{bubbles:true}));control.dispatchEvent(new Event('change',{bubbles:true}));
return {ok:true};
""",
        label,
        value,
    )
    if not result["ok"]:
        raise AssertionError(f"Could not set {label!r}: {result}")


def toggle_experiment(driver, label):
    element = WebDriverWait(driver, 30).until(
        lambda d: d.find_element(By.XPATH, f"//label[normalize-space(.)={json.dumps(label)}]/input[@type='checkbox']")
    )
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", element)
    element.click()


def set_endpoint_field(driver, endpoint, label, value):
    result = driver.execute_script(
        """
const endpoint=arguments[0], wanted=arguments[1], value=String(arguments[2]);
const card=[...document.querySelectorAll('.experimental-endpoint-card')].find(row=>row.querySelector('h4')?.textContent.trim()===endpoint);
if(!card)return {ok:false,reason:'endpoint missing',endpoints:[...document.querySelectorAll('.experimental-endpoint-card h4')].map(row=>row.textContent.trim())};
const label=[...card.querySelectorAll('label')].find(row=>row.textContent.trim()===wanted);
const control=label?.nextElementSibling;
if(!control)return {ok:false,reason:'field missing',labels:[...card.querySelectorAll('label')].map(row=>row.textContent.trim())};
const proto=control.tagName==='SELECT'?HTMLSelectElement.prototype:(control.tagName==='TEXTAREA'?HTMLTextAreaElement.prototype:HTMLInputElement.prototype);
Object.getOwnPropertyDescriptor(proto,'value').set.call(control,value);
control.dispatchEvent(new Event('input',{bubbles:true}));control.dispatchEvent(new Event('change',{bubbles:true}));
return {ok:true};
""",
        endpoint,
        label,
        value,
    )
    if not result["ok"]:
        raise AssertionError(f"Could not set {endpoint} / {label}: {result}")


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for argument in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1900,1800"):
        options.add_argument(argument)
    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    checks = []
    try:
        driver.get(BASE_URL + "/?pre-stage5-ui-e2e=" + RUN_ID)
        WebDriverWait(driver, 60).until(lambda d: "Create Project" in d.page_source)
        assert browser_api(driver, "GET", "/health")["step"] == "4B"
        for text in ("MAIN DASHBOARD", "Structure-first drug optimization workspace", "Default Workspace Settings", "Projects", "Compounds"):
            if text not in driver.page_source:
                raise AssertionError(f"Main Dashboard is missing {text!r}")
        driver.set_window_size(390, 844)
        WebDriverWait(driver, 10).until(lambda d: d.execute_script("return window.innerWidth") <= 500)
        columns = driver.execute_script("return getComputedStyle(document.querySelector('.dashboard-project-grid')).gridTemplateColumns")
        if " " in columns.strip():
            raise AssertionError(f"Mobile dashboard project grid is not single-column: {columns}")
        driver.set_window_size(1900, 1800)
        checks.extend(["Main Dashboard summary", "Responsive mobile dashboard"])

        set_labeled_control(driver, "Project Name *", PROJECT_NAME)
        set_labeled_control(driver, "Target *", "EGFR")
        set_labeled_control(driver, "Molecule Type", "Small Molecule")
        set_labeled_control(driver, "Description (optional)", "Public ethanol/ethylamine workflow fixture")
        click_button(driver, "Create Project")
        WebDriverWait(driver, 45).until(lambda d: PROJECT_NAME in d.page_source and "No compounds yet" in d.page_source)
        checks.append("Simplified project creation and project navigation")

        click_button(driver, "Add Compound", last=True)
        WebDriverWait(driver, 60).until(lambda d: "Draw Chemical Structure" in d.page_source and d.find_element(By.ID, "ketcher-editor"))
        set_labeled_control(driver, "Compound Name *", "UI-DRAW-001")
        WebDriverWait(driver, 120).until(
            lambda d: d.execute_script("return Boolean(document.getElementById('ketcher-editor')?.contentWindow?.ketcher)")
        )
        driver.set_script_timeout(60)
        driver.execute_async_script(
            "const done=arguments[arguments.length-1];document.getElementById('ketcher-editor').contentWindow.ketcher.setMolecule('CCO').then(()=>done(true)).catch(e=>done(String(e)));"
        )
        WebDriverWait(driver, 30).until(lambda d: d.find_element(By.XPATH, "//label[normalize-space(.)='SMILES']/following-sibling::input").get_attribute("value"))
        smiles = driver.find_element(By.XPATH, "//label[normalize-space(.)='SMILES']/following-sibling::input").get_attribute("value")
        if "CCO" not in smiles:
            raise AssertionError(f"Ketcher drawing did not synchronize to SMILES: {smiles}")
        click_button(driver, "Save Compound")
        WebDriverWait(driver, 60).until(lambda d: "UI-DRAW-001" in d.page_source and "STRUCTURE READY" in d.page_source)
        click_button(driver, "PROPERTIES", last=True)
        WebDriverWait(driver, 30).until(lambda d: "Properties have not been calculated" in d.page_source)
        click_button(driver, "Calculate Properties")
        WebDriverWait(driver, 60).until(lambda d: "Physicochemical properties" in d.page_source and "CALCULATED" in d.page_source)
        checks.extend(["Ketcher drawing and SMILES synchronization", "Save without calculation", "Later calculation"])

        click_button(driver, "Back to Compounds")
        click_button(driver, "Add Compound", last=True)
        set_labeled_control(driver, "Compound Name *", "UI-SMILES-002")
        set_labeled_control(driver, "SMILES", "CCN")
        WebDriverWait(driver, 60).until(
            lambda d: d.execute_script("return Boolean(document.getElementById('ketcher-editor')?.contentWindow?.ketcher)")
        )
        WebDriverWait(driver, 30).until(
            lambda d: "CCN" in d.execute_async_script("const done=arguments[arguments.length-1];document.getElementById('ketcher-editor').contentWindow.ketcher.getSmiles().then(done).catch(()=>done(''));" )
        )
        click_button(driver, "Save & Calculate")
        WebDriverWait(driver, 60).until(lambda d: "UI-SMILES-002" in d.page_source and "CALCULATED" in d.page_source)
        checks.append("SMILES workflow and Save & Calculate")

        click_button(driver, "Back to Compounds")
        row = WebDriverWait(driver, 30).until(lambda d: d.find_element(By.XPATH, "//tr[td[contains(.,'UI-DRAW-001')]]"))
        row.find_element(By.XPATH, ".//button[normalize-space()='Open']").click()
        click_button(driver, "ADMET", last=True)
        click_button(driver, "Add Experimental Data")
        for endpoint in ("Caco-2 Permeability", "Plasma Protein Binding (PPB)", "Human Microsomal Stability", "CYP Inhibition"):
            toggle_experiment(driver, endpoint)
        values = {
            "Caco-2 Permeability": "0.000012",
            "Plasma Protein Binding (PPB)": "85",
            "Human Microsomal Stability": "22",
            "CYP Inhibition": "3.4",
        }
        for endpoint, value in values.items():
            set_endpoint_field(driver, endpoint, "Value", value)
            set_endpoint_field(driver, endpoint, "Source", "A1-ONLY-EVIDENCE")
        click_button(driver, "Save Experimental Data")
        WebDriverWait(driver, 60).until(lambda d: d.page_source.count("A1-ONLY-EVIDENCE") >= 4 and "Prediction" in d.page_source)
        checks.append("Endpoint-selector experimental entry and Experimental/Prediction separation")

        click_button(driver, "Back to Compounds")
        row = driver.find_element(By.XPATH, "//tr[td[contains(.,'UI-SMILES-002')]]")
        row.find_element(By.XPATH, ".//button[normalize-space()='Open']").click()
        click_button(driver, "ADMET", last=True)
        WebDriverWait(driver, 30).until(lambda d: "No experimental measurement entered" in d.page_source)
        if "A1-ONLY-EVIDENCE" in driver.page_source:
            raise AssertionError("Compound A1 evidence leaked into A2 detail")
        checks.append("CompoundVersion UI isolation")

        click_button(driver, "Back to Compounds")
        for compound_name in ("UI-DRAW-001", "UI-SMILES-002"):
            checkbox = driver.find_element(By.XPATH, f"//tr[td[contains(.,'{compound_name}')]]//input[@type='checkbox']")
            checkbox.click()
        click_button(driver, "Compare Selected")
        WebDriverWait(driver, 45).until(lambda d: "Selected Compound Comparison" in d.page_source and "UI-DRAW-001" in d.page_source and "UI-SMILES-002" in d.page_source)
        checks.append("Checkbox-based Compare Selected")

        projects = browser_api(driver, "GET", "/projects")
        project = next(row for row in projects if row["name"] == PROJECT_NAME)
        project_detail = browser_api(driver, "GET", f"/projects/{project['id']}")
        first = next(row for row in project_detail["compounds"] if row["name"] == "UI-DRAW-001")
        workspace = browser_api(driver, "GET", f"/compound-versions/{first['version']['id']}/workspace")
        if workspace["scope"]["version_id"] != first["version"]["id"] or any(row["version_id"] != first["version"]["id"] for row in workspace["admet"]["measurements"]):
            raise AssertionError("Workspace API scope mismatch")

        result = {
            "status": "PASS", "base_url": BASE_URL, "project": PROJECT_NAME,
            "project_id": project["id"], "compound_version_id": first["version"]["id"],
            "checks": checks, "experimental_records": len(workspace["admet"]["measurements"]),
        }
        print(json.dumps(result, indent=2))
    finally:
        driver.quit()


if __name__ == "__main__":
    main()
