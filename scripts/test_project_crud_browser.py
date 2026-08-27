#!/usr/bin/env python3
"""
E2E Browser Test for Project Creation and Deletion Workflows
Validates:
1. New Project creation via Sidebar -> Form -> Create -> Direct Workspace load
2. Compound addition within the new project
3. Project Deletion via Projects view with typed confirmation modal
4. Clean state post-deletion
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
TEST_PROJECT_NAME = "__E2E_CRUD_PROJECT__"
TEST_TARGET = "EGFR-T790M"
TEST_COMPOUND_NAME = "Osimertinib-Lead"
TEST_SMILES = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"


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
        # Pre-clean
        driver.get(BASE_URL)
        time.sleep(2)
        try:
            projs = browser_api(driver, "GET", "/projects")
            for p in projs:
                if p["name"] == TEST_PROJECT_NAME:
                    browser_api(driver, "DELETE", f"/projects/{p['id']}", {"confirmation_name": TEST_PROJECT_NAME})
        except Exception as e:
            print("Pre-clean note:", e)

        # -------------------------------------------------------------
        # STEP 1: CREATE PROJECT VIA UI
        # -------------------------------------------------------------
        print("\n--- STEP 1: Testing Project Creation via UI ---")
        driver.get(BASE_URL)
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))

        # Click New Project in Sidebar
        click_nav(driver, "New Project")
        time.sleep(1)

        create_section = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "create-project-grid"))
        )

        name_input = driver.find_element(By.XPATH, "//label[contains(text(),'Project Name')]/following-sibling::input")
        name_input.clear()
        name_input.send_keys(TEST_PROJECT_NAME)

        target_input = driver.find_element(By.XPATH, "//label[contains(text(),'Target')]/following-sibling::input")
        target_input.clear()
        target_input.send_keys(TEST_TARGET)

        create_btn = driver.find_element(By.XPATH, "//button[normalize-space()='Create Project']")
        driver.execute_script("arguments[0].click();", create_btn)
        time.sleep(2)

        # Verify direct transition to Project Workspace
        workspace_header = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.CLASS_NAME, "project-header"))
        )
        ws_text = workspace_header.text
        assert TEST_PROJECT_NAME in ws_text, f"Expected project title '{TEST_PROJECT_NAME}' in workspace, got '{ws_text}'"
        assert TEST_TARGET in ws_text, f"Expected target '{TEST_TARGET}' in workspace header"
        print(f"✓ Project '{TEST_PROJECT_NAME}' created and workspace loaded immediately.")

        # -------------------------------------------------------------
        # STEP 2: ADD COMPOUND TO NEW PROJECT
        # -------------------------------------------------------------
        print("\n--- STEP 2: Adding Compound to Project ---")
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
        smiles_input.send_keys(TEST_SMILES)

        save_btn = modal.find_element(By.XPATH, ".//button[normalize-space()='Save']")
        driver.execute_script("arguments[0].click();", save_btn)
        time.sleep(2)

        # Verify compound overview loaded
        header_card = WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card"))
        )
        assert TEST_COMPOUND_NAME in header_card.text
        print(f"✓ Compound '{TEST_COMPOUND_NAME}' added and overview displayed.")

        # -------------------------------------------------------------
        # STEP 3: DELETE PROJECT VIA PROJECTS VIEW
        # -------------------------------------------------------------
        print("\n--- STEP 3: Testing Project Deletion via Projects View ---")
        click_nav(driver, "Projects")
        time.sleep(1)

        # Wait for projects list
        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-project")))
        time.sleep(1)

        # Find the project card for our test project
        cards = driver.find_elements(By.CLASS_NAME, "dashboard-project")
        target_card = None
        for c in cards:
            if TEST_PROJECT_NAME in c.text:
                target_card = c
                break
        assert target_card is not None, f"Could not find project card for '{TEST_PROJECT_NAME}'"

        # Click Delete button in that card
        delete_btn = target_card.find_element(By.XPATH, ".//button[normalize-space()='Delete…']")
        driver.execute_script("arguments[0].click();", delete_btn)
        time.sleep(1)

        # Confirm Modal
        del_modal = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.CLASS_NAME, "project-delete-modal"))
        )
        assert "PROJECT DELETION" in del_modal.text or "Delete Project" in del_modal.text

        # Type confirmation name
        confirm_input = del_modal.find_element(By.XPATH, ".//label[contains(text(),'Type')]/following-sibling::input")
        confirm_input.clear()
        confirm_input.send_keys(TEST_PROJECT_NAME)
        time.sleep(0.5)

        # Click Delete Permanently button
        perm_btn = del_modal.find_element(By.XPATH, ".//button[contains(text(),'Delete Project Permanently')]")
        assert perm_btn.is_enabled(), "Delete permanently button should be enabled after typing name"
        driver.execute_script("arguments[0].click();", perm_btn)
        time.sleep(2)

        # Verify return to dashboard and project is deleted
        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))
        print(f"✓ Project '{TEST_PROJECT_NAME}' permanently deleted and UI returned to Dashboard.")

        # Verify project is no longer in database via API
        projs_after = browser_api(driver, "GET", "/projects")
        assert not any(p["name"] == TEST_PROJECT_NAME for p in projs_after), f"Project '{TEST_PROJECT_NAME}' still exists in database!"
        print("✓ Verified project does not exist in database.")

        print("\n========================================================")
        print("✓ PROJECT CREATION & DELETION E2E TEST PASSED!")
        print("========================================================")

    finally:
        driver.quit()


if __name__ == "__main__":
    main()
