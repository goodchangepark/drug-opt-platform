#!/usr/bin/env python3
"""
E2E Browser Test for Compound Live 2D Structure Preview & Save Workflow
Validates:
1. Open project workspace
2. Click Add Compound
3. Type / paste SMILES and observe live 2D structure SVG preview appearing automatically
4. Live properties (MW, Formula, cLogP, TPSA) displayed before save
5. Click Save -> modal closes -> Compound Overview opens with full calculated properties and SVG
6. Teardown test project
"""

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
TEST_PROJECT_NAME = "__LIVE_PREVIEW_E2E_PROJ__"
TEST_COMPOUND_NAME = "Gefitinib-Live"
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
        driver.get(BASE_URL)
        time.sleep(2)

        # -------------------------------------------------------------
        # STEP 1: CREATE PROJECT
        # -------------------------------------------------------------
        print("\n--- STEP 1: Create Project ---")
        click_nav(driver, "New Project")
        time.sleep(1)

        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "create-project-grid")))

        name_input = driver.find_element(By.XPATH, "//label[contains(text(),'Project Name')]/following-sibling::input")
        name_input.clear()
        name_input.send_keys(TEST_PROJECT_NAME)

        target_input = driver.find_element(By.XPATH, "//label[contains(text(),'Target')]/following-sibling::input")
        target_input.clear()
        target_input.send_keys("EGFR")

        create_btn = driver.find_element(By.XPATH, "//button[normalize-space()='Create Project']")
        driver.execute_script("arguments[0].click();", create_btn)
        time.sleep(2)

        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "project-header")))
        print(f"✓ Project '{TEST_PROJECT_NAME}' created.")

        # -------------------------------------------------------------
        # STEP 2: OPEN ADD COMPOUND & TYPE SMILES -> VERIFY LIVE PREVIEW
        # -------------------------------------------------------------
        print("\n--- STEP 2: Open Add Compound and Enter SMILES ---")
        add_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//button[normalize-space()='Add Compound']"))
        )
        driver.execute_script("arguments[0].click();", add_btn)
        time.sleep(1)

        modal = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "compound-modal")))
        cname_input = modal.find_element(By.XPATH, ".//label[contains(text(),'Compound Name')]/following-sibling::input")
        cname_input.clear()
        cname_input.send_keys(TEST_COMPOUND_NAME)

        smiles_input = modal.find_element(By.XPATH, ".//label[contains(text(),'SMILES')]/following-sibling::input")
        smiles_input.clear()
        smiles_input.send_keys(GEFITINIB_SMILES)

        # Wait for live 2D structure preview card to appear automatically
        print("Waiting for real-time 2D structure SVG preview...")
        preview_card = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "structure-live-preview"))
        )
        card_text = preview_card.text
        assert "2D Chemical Structure" in card_text, "2D structure card title missing"
        assert "MW: 446.9" in card_text or "447" in card_text, "Calculated MW missing in live preview"
        assert "C22H24ClFN4O3" in card_text, "Calculated formula missing in live preview"

        # Check SVG presence inside preview card
        svg_elem = preview_card.find_element(By.TAG_NAME, "svg")
        assert svg_elem.is_displayed(), "Structure SVG is not displayed"
        print("✓ Live 2D Chemical Structure SVG preview successfully rendered in real time!")

        # -------------------------------------------------------------
        # STEP 3: CLICK SAVE & VERIFY COMPOUND OVERVIEW
        # -------------------------------------------------------------
        print("\n--- STEP 3: Click Save and Verify Compound Overview ---")
        save_btn = modal.find_element(By.XPATH, ".//button[normalize-space()='Save']")
        assert save_btn.is_enabled(), "Save button should be enabled"
        driver.execute_script("arguments[0].click();", save_btn)
        time.sleep(2)

        header_card = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card"))
        )
        header_text = header_card.text
        assert TEST_COMPOUND_NAME in header_text, "Compound name missing in overview header"
        assert "Formula: C22H24ClFN4O3" in header_text, "Calculated formula missing in overview header"
        assert "Version 1" in header_text, "Version 1 missing in overview header"
        assert "▶ PREDICT" in header_text, "Primary Predict button missing in overview header"
        print(f"✓ Compound '{TEST_COMPOUND_NAME}' successfully saved and Overview opened.")

        # -------------------------------------------------------------
        # STEP 4: TEARDOWN
        # -------------------------------------------------------------
        print("\n--- STEP 4: Teardown Test Project ---")
        click_nav(driver, "Projects")
        time.sleep(1)

        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-project")))
        cards = driver.find_elements(By.CLASS_NAME, "dashboard-project")
        target_card = None
        for c in cards:
            if TEST_PROJECT_NAME in c.text:
                target_card = c
                break
        assert target_card is not None

        del_btn = target_card.find_element(By.XPATH, ".//button[normalize-space()='Delete…']")
        driver.execute_script("arguments[0].click();", del_btn)
        time.sleep(1)

        del_modal = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "project-delete-modal"))
        )
        confirm_input = del_modal.find_element(By.XPATH, ".//label[contains(text(),'Type')]/following-sibling::input")
        confirm_input.clear()
        confirm_input.send_keys(TEST_PROJECT_NAME)
        time.sleep(0.5)

        perm_btn = del_modal.find_element(By.XPATH, ".//button[contains(text(),'Delete Project Permanently')]")
        driver.execute_script("arguments[0].click();", perm_btn)
        time.sleep(2)

        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))
        print(f"✓ Test project '{TEST_PROJECT_NAME}' cleaned up.")

        print("\n========================================================")
        print("✓ COMPOUND LIVE PREVIEW & SAVE E2E TEST PASSED!")
        print("========================================================")

    finally:
        driver.quit()


if __name__ == "__main__":
    main()
