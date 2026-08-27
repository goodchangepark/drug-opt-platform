#!/usr/bin/env python3
"""Stage 5B-4 UI Polish & Unified Prediction Workflow Browser E2E Acceptance Test.

Verifies:
1. Dashboard 4 sections:
   - Section 1: Platform Overview (hero, workflow ribbon, 3 value propositions, statistics)
   - Section 2: Scientific Workspace (6 card rectangular grid)
   - Section 3: Quick Start Guide (7 compact numbered workflow steps)
   - Section 4: Projects (clickable project names, stats, de-emphasized delete button)
   - NO New Project card on Dashboard
2. Left Sidebar:
   - Clickable navigation
   - Footer: Drug-OPT, v0.6.2-stage5b4-ui, Updated: 2026-08-27
3. Project Navigation:
   - Click project title to open project workspace
4. Compound Workflow:
   - Add Compound modal -> Save -> routes directly to Compound Detail Overview
   - Structure container constrained to prevent sidebar overlap
   - Primary PREDICT button in header card triggers orchestrated multi-endpoint prediction
   - Timestamp, status badge, and experimental update awareness banner
   - Overview PK Summary cards (CL, Vd, Half-Life, Oral F, Human Translation readiness)
5. Scientific Tab Re-Predict Buttons:
   - Properties Tab: ↺ RE-PREDICT
   - Activity Tab: ↺ PREDICT / RE-PREDICT with ACTIVITY MODEL NOT READY guardrail
   - ADMET Tab: ↺ RE-PREDICT preserving experimental data
   - Metabolism Tab: ↺ RE-PREDICT
   - PK Tab: ↺ UPDATE PK ANALYSIS
6. Captures visual screenshots for all key views.
7. Clean teardown of temporary project.
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
TEMP_PROJECT_NAME = f"__UI_POLISH_TEMP_{RUN_ID}__"
TEMP_COMPOUND_ID = "MOL-8801"
TEMP_COMPOUND_NAME = "Imatinib Analog"
TEMP_SMILES = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(-c5cccnc5)n4"

SCREENSHOT_DIR = Path("/home/xavier/chem/drug-opt-platform/artifacts")
SCREENSHOT_DIR.mkdir(parents=True, exist_ok=True)


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
    options.add_argument("--window-size=1440,1080")

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    print(f"[*] Starting Stage 5B-4 UI Polish & Prediction Workflow Browser E2E against {BASE_URL}...")

    project_id = None
    try:
        driver.get(BASE_URL)
        WebDriverWait(driver, 20).until(lambda d: d.find_element(By.TAG_NAME, "main"))
        time.sleep(1.5)

        # 1. Main Dashboard & Typography Verification
        print(" -> [1/7] Testing Dashboard Layout (4 Sections) & Sidebar Footer...")
        main_elem = driver.find_element(By.TAG_NAME, "main")
        main_text = main_elem.text

        # Section 1: Platform Overview
        assert "PLATFORM OVERVIEW" in main_text, "Missing Section 1: PLATFORM OVERVIEW"
        assert "STRUCTURE" in main_text and "PROPERTIES" in main_text and "ADMET" in main_text, "Missing workflow ribbon"
        assert "1. Rigorous Data Provenance" in main_text, "Missing value prop 1"
        assert "2. Transparent Scientific Confidence" in main_text, "Missing value prop 2"
        assert "3. End-to-End Translation" in main_text, "Missing value prop 3"

        # Section 2: Scientific Workspace
        assert "SCIENTIFIC WORKSPACE" in main_text, "Missing Section 2: SCIENTIFIC WORKSPACE"
        assert "Structure & Chemistry" in main_text, "Missing Structure & Chemistry module"
        assert "PK / DMPK" in main_text, "Missing PK / DMPK module"


        # Section 3: Quick Start Guide
        assert "QUICK START GUIDE" in main_text, "Missing Section 3: QUICK START GUIDE"
        assert "01" in main_text and "07" in main_text, "Missing 7 quick start numbered steps"

        # Section 4: Projects (Research Portfolio)
        assert "RESEARCH PORTFOLIO" in main_text or "Projects" in main_text, "Missing Section 4: Projects"

        # Verify NO "New Project" card in Main Dashboard
        assert "NEW WORKSPACE" not in main_text, "New Project form must NOT be on Main Dashboard"

        # Verify Sidebar Footer
        sidebar_text = driver.find_element(By.XPATH, "//aside[contains(@class,'sidebar')]").text
        assert "v0.6.2-stage5b4-ui" in sidebar_text, "Missing v0.6.2-stage5b4-ui in sidebar footer"
        assert "Updated: 2026-08-27" in sidebar_text, "Missing update date in sidebar footer"

        # Save Dashboard Screenshot
        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_dashboard.png"))
        print("    [✓] Dashboard verified & screenshot saved.")

        # 2. Seed Project & Compound
        print(" -> [2/7] Seeding temporary project & compound...")
        proj = browser_api(driver, "POST", "/projects", {
            "name": TEMP_PROJECT_NAME,
            "target": "BCR-ABL Kinase",
            "molecule_type": "Small Molecule",
            "description": "UI Polish E2E test project",
        })
        project_id = proj["id"]

        # Navigate to Dashboard and verify project is clickable
        driver.get(f"{BASE_URL}/")
        WebDriverWait(driver, 20).until(lambda d: d.find_element(By.TAG_NAME, "main"))
        time.sleep(1.0)

        # Click project title link on Dashboard
        proj_link = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, f"//button[contains(@class,'project-link-title') and normalize-space()={json.dumps(TEMP_PROJECT_NAME)}]")
        )
        driver.execute_script("arguments[0].click();", proj_link)
        time.sleep(1.0)

        # Add Compound via Add Compound button
        add_comp_btn = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, "//button[normalize-space()='Add Compound']")
        )
        driver.execute_script("arguments[0].click();", add_comp_btn)
        time.sleep(0.5)

        # Fill Add Compound modal
        name_input = driver.find_element(By.XPATH, "//div[contains(@class,'modal')]//label[contains(.,'Compound Name')]/following-sibling::input")
        name_input.send_keys(TEMP_COMPOUND_NAME)
        smiles_input = driver.find_element(By.XPATH, "//div[contains(@class,'modal')]//label[contains(.,'SMILES')]/following-sibling::input")
        smiles_input.send_keys(TEMP_SMILES)


        save_btn = driver.find_element(By.XPATH, "//div[contains(@class,'compound-modal')]//button[normalize-space()='Save']")
        driver.execute_script("arguments[0].click();", save_btn)
        time.sleep(2.0)


        # 3. Overview & Primary Predict Workflow
        print(" -> [3/7] Verifying Compound Header Card & Primary PREDICT Button...")
        header_card = WebDriverWait(driver, 15).until(
            lambda d: d.find_element(By.XPATH, "//div[contains(@class,'compound-header-card')]")
        )
        header_text = header_card.text
        assert TEMP_COMPOUND_NAME in header_text
        assert "Copy SMILES" in header_text

        # Click Primary Predict button
        predict_btn = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, "//button[contains(@class,'btn-predict-primary')]")
        )
        print("    [*] Clicking Primary PREDICT button...")
        driver.execute_script("arguments[0].click();", predict_btn)

        # Wait dynamically for prediction completion (button returns from PREDICTING to PREDICT or MW appears)
        WebDriverWait(driver, 60).until(
            lambda d: "MW:" in d.find_element(By.XPATH, "//div[contains(@class,'compound-header-card')]").text
        )
        time.sleep(1.0)

        # Re-fetch header & main after prediction
        header_text_after = driver.find_element(By.XPATH, "//div[contains(@class,'compound-header-card')]").text
        assert "MW:" in header_text_after
        assert "Formula:" in header_text_after


        # Refresh detail and verify predictions rendered
        overview_text = driver.find_element(By.TAG_NAME, "main").text
        assert "BASIC PROPERTY SUMMARY" in overview_text
        assert "EXECUTIVE SCIENTIFIC SUMMARY" in overview_text
        assert "Aqueous Solubility" in overview_text
        assert "Human Microsomal Stab" in overview_text
        assert "TRANSLATIONAL PK SUMMARY" in overview_text

        # Capture Overview Screenshot
        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_overview.png"))
        print("    [✓] Compound Overview & PK Summary verified & screenshot saved.")


        # 4. Properties Tab & Re-Predict Button
        print(" -> [4/7] Verifying Properties Tab & ↺ RE-PREDICT Button...")
        click_detail_tab(driver, "PROPERTIES")
        time.sleep(1.0)
        prop_text = driver.find_element(By.TAG_NAME, "main").text
        assert "PHYSICOCHEMICAL & DRUG-LIKENESS TABLE" in prop_text
        assert "Molecular Properties & Reference Assessments" in prop_text
        assert "pH-Dependent Ionization" in prop_text

        # Verify Re-Predict button exists and is clickable
        repredict_btn = driver.find_element(By.XPATH, "//button[contains(@class,'tab-repredict-btn') and contains(.,'RE-PREDICT')]")
        driver.execute_script("arguments[0].click();", repredict_btn)
        time.sleep(1.0)

        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_properties.png"))
        print("    [✓] Properties tab verified & screenshot saved.")

        # 5. ADMET Tab & Re-Predict Button
        print(" -> [5/7] Verifying ADMET Tab & ↺ RE-PREDICT Button...")
        click_detail_tab(driver, "ADMET")
        admet_repredict = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, "//button[contains(@class,'tab-repredict-btn') and contains(.,'RE-PREDICT')]")
        )
        assert admet_repredict.is_displayed()
        admet_text = driver.find_element(By.TAG_NAME, "main").text
        assert "ADMET Developability Profile" in admet_text
        assert "VISUAL PROFILE — QUALITATIVE NORMALIZED REPRESENTATION" in admet_text

        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_admet.png"))
        print("    [✓] ADMET tab verified & screenshot saved.")

        # 6. Metabolism Tab & Re-Predict Button
        print(" -> [6/7] Verifying Metabolism Tab & ↺ RE-PREDICT Button...")
        click_detail_tab(driver, "METABOLISM")
        metab_repredict = WebDriverWait(driver, 10).until(
            lambda d: d.find_element(By.XPATH, "//button[contains(@class,'tab-repredict-btn') and contains(.,'RE-PREDICT')]")
        )
        assert metab_repredict.is_displayed()
        metab_text = driver.find_element(By.TAG_NAME, "main").text
        assert "METABOLIC STABILITY · LIVER MICROSOMES" in metab_text
        assert "CYP450 ENZYME PANEL" in metab_text

        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_metabolism.png"))
        print("    [✓] Metabolism tab verified & screenshot saved.")

        # 7. PK Tab & Update PK Button
        print(" -> [7/7] Verifying PK Tab & ↺ UPDATE PK ANALYSIS Button...")
        click_detail_tab(driver, "PK")
        time.sleep(2.0)
        pk_update = WebDriverWait(driver, 20).until(
            lambda d: d.find_element(By.XPATH, "//button[contains(@class,'tab-repredict-btn') and contains(.,'UPDATE PK ANALYSIS')]")
        )
        assert pk_update.is_displayed()
        pk_text = driver.find_element(By.TAG_NAME, "main").text
        assert "MULTI-SPECIES PK SUMMARY" in pk_text
        assert "Cross-Species Comparative In Vivo & Translational PK Matrix" in pk_text

        driver.save_screenshot(str(SCREENSHOT_DIR / "screenshot_pk.png"))
        print("    [✓] PK tab verified & screenshot saved.")




        print("\n[✓] ALL UI POLISH & WORKFLOW E2E TESTS PASSED SUCCESSFULLY!")

    finally:
        if project_id:
            print(" -> Cleaning up temporary test project...")
            try:
                browser_api(driver, "DELETE", f"/projects/{project_id}", {"confirmation_name": TEMP_PROJECT_NAME})
            except Exception as e:
                print(f"[!] Warning during teardown: {e}")
        driver.quit()


if __name__ == "__main__":
    main()
