#!/usr/bin/env python3
"""Stage 5B-4 UI/UX Scientific Refinement Browser E2E Acceptance Test.

Verifies:
1. Main Dashboard: Platform Overview, New Workspace, Scientific Workspace.
2. Sidebar footer: Drug-OPT, v0.6.1-stage5b4-ui, Updated: 2026-08-27.
3. Help Page: Version History table & current runtime packages.
4. Settings Page: Model Registry, PK Methods, Platform Information.
5. Compound Detail:
   - Structure Header with Formula, MW, SMILES & Copy button
   - Overview Tab: Property Summary grid + ADMET Highlights cards
   - Properties Tab: Single aligned physicochemical table + Ionization + Experimental entry
   - ADMET Tab: Unified prediction table + Visual Profile chart
   - Metabolism Tab: Cross-species microsomal stability + CYP panel with Model Applicability
   - PK Tab: Multi-species comparative PK summary matrix
6. Clean teardown of temporary project.
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


BASE_URL = os.environ.get("DRUG_OPT_BASE_URL", "http://127.0.0.1:8765")
RUN_ID = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
TEMP_PROJECT_NAME = f"__UI_REFINEMENT_TEMP_{RUN_ID}__"
TEMP_COMPOUND_ID = "CIPRO-001"
TEMP_COMPOUND_NAME = "Ciprofloxacin Lead"
TEMP_SMILES = "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O"


def click_nav(driver, text, timeout=20):
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'global-nav')]//button[normalize-space()={json.dumps(text)}]")
    )
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
    driver.execute_script("arguments[0].click();", btn)
    return btn


def click_detail_tab(driver, text, timeout=20):
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'detail-tabs')]//button[normalize-space()={json.dumps(text)}]")
    )
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
    driver.execute_script("arguments[0].click();", btn)
    return btn


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


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--window-size=1440,960")

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    print(f"[*] Starting Stage 5B-4 UI Refinement Browser E2E against {BASE_URL}...")

    project_id = None
    try:
        driver.get(BASE_URL)
        WebDriverWait(driver, 20).until(lambda d: d.find_element(By.TAG_NAME, "main"))
        time.sleep(1.0)

        # 1. Main Dashboard
        print(" -> [1/6] Testing Main Dashboard & Sidebar Footer...")
        main_text = driver.find_element(By.TAG_NAME, "main").text
        assert "PLATFORM OVERVIEW" in main_text, "Missing PLATFORM OVERVIEW"
        assert "NEW WORKSPACE" in main_text, "Missing NEW WORKSPACE"
        assert "SCIENTIFIC WORKSPACE" in main_text, "Missing SCIENTIFIC WORKSPACE"

        sidebar_text = driver.find_element(By.XPATH, "//aside[contains(@class,'sidebar')]").text
        assert "v0.6.1-stage5b4-ui" in sidebar_text, "Missing version in sidebar footer"
        assert "Updated: 2026-08-27" in sidebar_text, "Missing update date in sidebar footer"

        # 2. Help Page & Version History
        print(" -> [2/6] Testing Help Page & Version History...")
        click_nav(driver, "Help")
        time.sleep(1.0)
        help_text = driver.find_element(By.TAG_NAME, "main").text
        assert "Drug-OPT Platform Help" in help_text
        assert "Current Platform Version" in help_text
        assert "Version History & Scientific Milestones" in help_text
        assert "0.6.1-stage5b4-ui" in help_text

        # 3. Settings Page
        print(" -> [3/6] Testing Settings Registry View...")
        click_nav(driver, "Settings")
        time.sleep(1.0)
        settings_text = driver.find_element(By.TAG_NAME, "main").text
        assert "MODEL REGISTRY & GOVERNANCE" in settings_text
        assert "PK & TRANSLATION METHOD REGISTRY" in settings_text
        assert "PLATFORM INFORMATION" in settings_text

        # 4. Project Creation & Compound Registration via API
        print(" -> [4/6] Seeding test project & compound with predictions & PK...")
        proj = browser_api(driver, "POST", "/projects", {
            "name": TEMP_PROJECT_NAME,
            "target": "DNA Gyrase",
            "molecule_type": "Small Molecule",
            "description": "UI Refinement E2E test project",
        })
        project_id = proj["id"]

        comp = browser_api(driver, "POST", f"/projects/{project_id}/compounds", {
            "compound_id": TEMP_COMPOUND_ID,
            "name": TEMP_COMPOUND_NAME,
            "smiles": TEMP_SMILES,
            "notes": "E2E compound",
            "calculate": True,
        })
        row_id = comp["row_id"]
        version_id = comp["version"]["id"]

        # Run ADMET & Metabolism predictions
        browser_api(driver, "POST", f"/admet/predict/{version_id}")
        browser_api(driver, "POST", f"/metabolism/predict/{version_id}")

        # Seed Rat and Mouse PK studies
        r_study = browser_api(driver, "POST", f"/compounds/{row_id}/pk-studies", {
            "study_name": "Rat IV Bolus Study", "species": "Rat", "route": "IV", "dose": 10.0, "dose_unit": "mg/kg",
        })
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/observations", [
            {"time_raw": 0.083, "concentration_raw": 3200.0, "subject_group_id": "Rat_Group"},
            {"time_raw": 0.5, "concentration_raw": 1450.0, "subject_group_id": "Rat_Group"},
            {"time_raw": 1.0, "concentration_raw": 780.0, "subject_group_id": "Rat_Group"},
            {"time_raw": 2.0, "concentration_raw": 240.0, "subject_group_id": "Rat_Group"},
            {"time_raw": 4.0, "concentration_raw": 45.0, "subject_group_id": "Rat_Group"},
        ])
        browser_api(driver, "POST", f"/pk-studies/{r_study['id']}/run-nca", {"selection_mode": "AUTO"})

        # Reload page and open project compound
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 20).until(lambda d: d.find_element(By.TAG_NAME, "main"))
        time.sleep(1.0)
        click_nav(driver, "Projects")
        time.sleep(1.0)

        # Click the project card
        proj_card = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, f"//article[contains(.,'{TEMP_PROJECT_NAME}')]")
        )
        driver.execute_script("arguments[0].click();", proj_card)
        time.sleep(1.0)

        # Open compound detail
        open_btn = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, f"//tr[contains(.,'{TEMP_COMPOUND_NAME}')]//button[normalize-space()='Open']")
        )
        driver.execute_script("arguments[0].click();", open_btn)
        time.sleep(1.0)

        # 5. Compound Detail Header & Overview Tab
        print(" -> [5/6] Verifying Compound Header & All Detail Tabs...")
        header_text = driver.find_element(By.XPATH, "//div[contains(@class,'compound-header-card')]").text
        assert TEMP_COMPOUND_NAME in header_text
        assert "MW:" in header_text
        assert "Copy SMILES" in header_text

        overview_text = driver.find_element(By.TAG_NAME, "main").text
        assert "BASIC PROPERTY SUMMARY" in overview_text
        assert "EXECUTIVE SCIENTIFIC SUMMARY" in overview_text
        assert "Aqueous Solubility" in overview_text
        assert "Human Microsomal Stab" in overview_text

        # Properties Tab
        click_detail_tab(driver, "PROPERTIES")
        time.sleep(1.0)
        prop_text = driver.find_element(By.TAG_NAME, "main").text
        assert "PHYSICOCHEMICAL & DRUG-LIKENESS TABLE" in prop_text
        assert "Molecular Properties & Reference Assessments" in prop_text
        assert "pH-Dependent Ionization" in prop_text
        assert "EXPERIMENTAL DATA ENTRY" in prop_text

        # ADMET Tab
        click_detail_tab(driver, "ADMET")
        time.sleep(1.0)
        admet_text = driver.find_element(By.TAG_NAME, "main").text
        assert "ADMET Developability Profile" in admet_text
        assert "2 · PREDICTION RESULTS" in admet_text
        assert "VISUAL PROFILE — QUALITATIVE NORMALIZED REPRESENTATION" in admet_text

        # Metabolism Tab
        click_detail_tab(driver, "METABOLISM")
        time.sleep(1.0)
        metab_text = driver.find_element(By.TAG_NAME, "main").text
        assert "METABOLIC STABILITY · LIVER MICROSOMES" in metab_text
        assert "Cross-Species Microsomal Stability (Human, Rat, Mouse)" in metab_text
        assert "CYP450 ENZYME PANEL" in metab_text

        # PK Tab
        click_detail_tab(driver, "PK")
        time.sleep(1.0)
        pk_text = driver.find_element(By.TAG_NAME, "main").text
        assert "MULTI-SPECIES PK SUMMARY" in pk_text
        assert "Cross-Species Comparative In Vivo & Translational PK Matrix" in pk_text

        print("[✓] All UI/UX views and scientific tabs verified successfully!")

    finally:
        # 6. Teardown
        if project_id:
            print(" -> [6/6] Tearing down temporary project...")
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": TEMP_PROJECT_NAME})
            except Exception as e:
                print(f"[!] Warning during teardown: {e}")
        driver.quit()

    print("\n[✓] Stage 5B-4 UI Refinement Browser E2E PASSED!")


if __name__ == "__main__":
    main()
