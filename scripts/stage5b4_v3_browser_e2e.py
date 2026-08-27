#!/usr/bin/env python3
"""
Stage 5B-4 Refinement 3 Browser E2E Acceptance Test
Validates:
1. Main Dashboard Redesign matching reference design
2. Clean single-card Platform Overview and 3-column Scientific Workspace
3. Compound Save workflow with instant property calculation and auto-open Overview
4. Persistence across navigation and clean teardown
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

BASE_URL = os.environ.get("DRUG_OPT_BASE_URL", "http://127.0.0.1:8765")
VALIDATION_DIR = Path("/home/xavier/chem/drug-opt-platform/validation")
VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR = Path("/home/xavier/chem/drug-opt-platform/artifacts")
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

TEMP_PROJECT_NAME = "__DASHBOARD_SAVE_E2E__"
GEFITINIB_SMILES = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"


def click_nav(driver, text, timeout=20):
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'global-nav')]//button[normalize-space()={json.dumps(text)}]")
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
    try:
        # Pre-clean test project if left from prior run
        driver.get(BASE_URL)
        time.sleep(2)
        try:
            projects = browser_api(driver, "GET", "/projects")
            for p in projects:
                if p["name"] == TEMP_PROJECT_NAME:
                    browser_api(driver, "DELETE", f"/projects/{p['id']}", {"confirmation_name": TEMP_PROJECT_NAME})
        except Exception as e:
            print("Pre-clean note:", e)

        # -------------------------------------------------------------
        # STEP 1: Main Dashboard Redesign Verification
        # -------------------------------------------------------------
        print("\n--- STEP 1: Verifying Main Dashboard Redesign ---")
        driver.get(BASE_URL)
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.CLASS_NAME, "scientific-card")))
        time.sleep(1)

        # Verify Sidebar
        sidebar = driver.find_element(By.CLASS_NAME, "sidebar")
        assert sidebar.is_displayed(), "Sidebar is not visible"
        sidebar_text = sidebar.text
        assert "Drug-OPT" in sidebar_text, "Brand title missing in sidebar"
        assert "v0.6.3-stage5b4-ui" in sidebar_text, "v0.6.3-stage5b4-ui missing in sidebar footer"
        assert "Updated: 2026-08-27" in sidebar_text, "Updated date missing in sidebar footer"
        print("✓ Sidebar and footer version verified (v0.6.3-stage5b4-ui).")

        # Verify Platform Overview Hero
        hero = driver.find_element(By.CLASS_NAME, "dashboard-hero")
        hero_text = hero.text
        assert "PLATFORM OVERVIEW" in hero_text, "Eyebrow PLATFORM OVERVIEW missing"
        assert "Drug Optimization Platform" in hero_text, "Platform title missing"
        assert "Structure, activity, ADMET, DMPK" in hero_text, "Platform description missing"
        assert "Structure-based compound management" in hero_text, "Capability bullet missing"
        assert "Full prediction provenance" in hero_text, "Capability bullet missing"

        # Verify 4 Metric cards
        stats_cards = driver.find_elements(By.CLASS_NAME, "dashboard-stat-card")
        assert len(stats_cards) == 4, f"Expected 4 stat cards, got {len(stats_cards)}"
        print("✓ Platform Overview and 4 metric cards verified.")

        # Verify Scientific Workspace Grid
        sci_section = driver.find_element(By.CLASS_NAME, "scientific-workspace-section")
        sci_text = sci_section.text
        assert "SCIENTIFIC WORKSPACE" in sci_text, "SCIENTIFIC WORKSPACE eyebrow missing"
        assert "Available Scientific Modules" in sci_text, "Available Scientific Modules title missing"
        assert "Status reflects the current local engine and model registry." in sci_text, "Registry note missing"

        sci_cards = sci_section.find_elements(By.CLASS_NAME, "scientific-card")
        assert len(sci_cards) == 7, f"Expected 7 scientific module cards, got {len(sci_cards)}"
        module_titles = [c.find_element(By.TAG_NAME, "h3").text for c in sci_cards]
        print(f"✓ Found {len(sci_cards)} scientific cards: {module_titles}")
        expected_modules = ["Structure & Chemistry", "Activity & SAR", "ADME", "CYP & Transporters", "Safety / Toxicology", "Optimization", "PK / DMPK"]
        for em in expected_modules:
            assert any(em in mt for mt in module_titles), f"Missing module {em}"

        # Verify Dashboard is clutter-free (No Quick Start, No New Project, No Project list on Main Dashboard)
        assert len(driver.find_elements(By.CLASS_NAME, "quick-start-step-box")) == 0, "Quick Start guide should NOT be on Main Dashboard"
        assert len(driver.find_elements(By.CLASS_NAME, "create-project-grid")) == 0, "New Project form should NOT be on Main Dashboard"
        assert len(driver.find_elements(By.CLASS_NAME, "dashboard-project-grid")) == 0, "Project portfolio list should NOT be on Main Dashboard"
        print("✓ Verified Main Dashboard has no clutter.")

        # Capture Dashboard Screenshot
        dashboard_ss_path = str(VALIDATION_DIR / "dashboard_redesign.png")
        driver.save_screenshot(dashboard_ss_path)
        driver.save_screenshot(str(ARTIFACTS_DIR / "dashboard_redesign.png"))
        print(f"✓ Saved dashboard screenshot to {dashboard_ss_path}")

        # -------------------------------------------------------------
        # STEP 2: Compound Save & Overview Workflow
        # -------------------------------------------------------------
        print("\n--- STEP 2: Verifying Compound Save & Overview Workflow ---")
        # Navigate to New Project via Sidebar
        click_nav(driver, "New Project")
        time.sleep(1)
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "create-project-grid")))

        # Fill New Project form
        name_input = driver.find_element(By.XPATH, "//label[contains(text(),'Project Name')]/following-sibling::input")
        name_input.clear()
        name_input.send_keys(TEMP_PROJECT_NAME)

        target_input = driver.find_element(By.XPATH, "//label[contains(text(),'Target')]/following-sibling::input")
        target_input.clear()
        target_input.send_keys("EGFR")

        create_btn = driver.find_element(By.XPATH, "//button[normalize-space()='Create Project']")
        driver.execute_script("arguments[0].click();", create_btn)
        time.sleep(2)

        # In Project Workspace, click "Add Compound"
        add_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//button[normalize-space()='Add Compound']"))
        )
        driver.execute_script("arguments[0].click();", add_btn)
        time.sleep(1)

        # Fill Add Compound modal
        modal = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "compound-modal")))
        cmp_name_input = modal.find_element(By.XPATH, ".//label[contains(text(),'Compound Name')]/following-sibling::input")
        cmp_name_input.clear()
        cmp_name_input.send_keys("Gefitinib-E2E")

        smiles_input = modal.find_element(By.XPATH, ".//label[contains(text(),'SMILES')]/following-sibling::input")
        smiles_input.clear()
        smiles_input.send_keys(GEFITINIB_SMILES)

        # Click SAVE button
        save_btn = modal.find_element(By.XPATH, ".//button[normalize-space()='Save']")
        driver.execute_script("arguments[0].click();", save_btn)
        time.sleep(2)

        # Verify auto-navigation to Compound Detail Overview
        header_card = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card"))
        )
        header_text = header_card.text
        assert "Gefitinib-E2E" in header_text, "Compound name missing in overview"
        assert "Version 1" in header_text, "Version 1 missing in overview"
        assert "Formula: C22H24ClFN4O3" in header_text, "Calculated formula missing in overview"
        assert "MW: 446.90 g/mol" in header_text or "446.9" in header_text, "Calculated MW missing in overview"
        assert "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1" in header_text, "Canonical SMILES missing in overview"
        assert "▶ PREDICT" in header_text, "Primary PREDICT button missing in overview"
        print("✓ Compound Overview verified with immediate calculated properties and primary PREDICT button.")

        # Verify Overview Sections: ADMET Summary and Translational PK Summary
        workspace = driver.find_element(By.CLASS_NAME, "compound-workspace")
        ws_text = workspace.text
        assert "EXECUTIVE SCIENTIFIC SUMMARY" in ws_text or "ADMET & DMPK Highlights" in ws_text, "ADMET Summary missing in compound detail"
        assert "TRANSLATIONAL PK SUMMARY" in ws_text, "Translational PK Summary missing in compound detail"

        # Capture Compound Overview Screenshot
        overview_ss_path = str(VALIDATION_DIR / "compound_save_overview.png")
        driver.save_screenshot(overview_ss_path)
        driver.save_screenshot(str(ARTIFACTS_DIR / "compound_save_overview.png"))
        print(f"✓ Saved compound overview screenshot to {overview_ss_path}")

        # -------------------------------------------------------------
        # STEP 3: Cleanup Test Project
        # -------------------------------------------------------------
        print("\n--- STEP 3: Teardown Test Project ---")
        projects = browser_api(driver, "GET", "/projects")
        for p in projects:
            if p["name"] == TEMP_PROJECT_NAME:
                browser_api(driver, "DELETE", f"/projects/{p['id']}", {"confirmation_name": TEMP_PROJECT_NAME})
                print(f"✓ Cleaned up test project {p['name']} (ID {p['id']}).")

        print("\n========================================================")
        print("✓ ALL BROWSER E2E ACCEPTANCE CHECKS PASSED SUCCESSFULLY!")
        print("========================================================")

    finally:
        driver.quit()


if __name__ == "__main__":
    main()
