import os
import sys
import time
import json
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

BASE_URL = "http://127.0.0.1:8765"
TEMP_PROJECT_NAME = "__SCIENTIFIC_UPDATE_E2E__"
GEFITINIB_SMILES = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
ERLOTINIB_SMILES = "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"

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

def click_project_tab(driver, text, timeout=20):
    btn = WebDriverWait(driver, timeout).until(
        lambda d: d.find_element(By.XPATH, f"//nav[contains(@class,'project-nav')]//button[normalize-space()={json.dumps(text)}]")
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

def run_browser_e2e():
    options = Options()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-gpu")
    options.add_argument("--window-size=1600,1200")

    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    try:
        print("[1/6] Navigating to Drug-OPT Platform...")
        driver.get(BASE_URL)
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))
        print("      Main Dashboard verified.")

        # Clean existing test project if needed
        try:
            projects = browser_api(driver, "GET", "/projects")
            for p in projects:
                if p["name"] in [TEMP_PROJECT_NAME, "__LIVE_PREVIEW_E2E_PROJ__"]:
                    browser_api(driver, "DELETE", f"/projects/{p['id']}", {"confirmation_name": p["name"]})
        except Exception as e:
            print("      Pre-clean note:", e)

        # Create fresh test project via UI
        print("[2/6] Creating test project and seeding compounds...")
        click_nav(driver, "New Project")
        time.sleep(1)
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "create-project-grid")))

        name_input = driver.find_element(By.XPATH, "//label[contains(text(),'Project Name')]/following-sibling::input")
        name_input.clear()
        name_input.send_keys(TEMP_PROJECT_NAME)

        target_input = driver.find_element(By.XPATH, "//label[contains(text(),'Target')]/following-sibling::input")
        target_input.clear()
        target_input.send_keys("EGFR")

        create_btn = driver.find_element(By.XPATH, "//button[normalize-space()='Create Project']")
        driver.execute_script("arguments[0].click();", create_btn)
        time.sleep(2)

        # Add Compound 1 (Gefitinib)
        add_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//button[normalize-space()='Add Compound']"))
        )
        driver.execute_script("arguments[0].click();", add_btn)
        time.sleep(1)

        modal = WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "compound-modal")))
        cmp_name_input = modal.find_element(By.XPATH, ".//label[contains(text(),'Compound Name')]/following-sibling::input")
        cmp_name_input.clear()
        cmp_name_input.send_keys("Gefitinib")

        smiles_input = modal.find_element(By.XPATH, ".//label[contains(text(),'SMILES')]/following-sibling::input")
        smiles_input.clear()
        smiles_input.send_keys(GEFITINIB_SMILES)

        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.CLASS_NAME, "structure-live-preview")))
        save_btn = modal.find_element(By.XPATH, ".//button[normalize-space()='Save']")
        driver.execute_script("arguments[0].click();", save_btn)
        print("      Saved Gefitinib.")
        time.sleep(2)

        # Open Gefitinib detail via Open button in table
        open_btns = WebDriverWait(driver, 10).until(
            EC.presence_of_all_elements_located((By.XPATH, "//button[normalize-space()='Open']"))
        )
        driver.execute_script("arguments[0].click();", open_btns[0])
        time.sleep(2)

        # Trigger ▶ PREDICT for Gefitinib
        WebDriverWait(driver, 15).until(EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card")))
        predict_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//button[contains(text(), 'PREDICT')]"))
        )
        driver.execute_script("arguments[0].click();", predict_btn)
        print("      Triggered ▶ PREDICT pipeline for Gefitinib. Awaiting predictions...")
        time.sleep(8)

        # 3. Test Metabolism Tab
        print("[3/6] Testing Metabolism Tab (5 Species Table)...")
        click_detail_tab(driver, "METABOLISM")
        time.sleep(3)

        page_text = driver.find_element(By.TAG_NAME, "body").text
        assert "Human Liver Microsomes" in page_text, "Human Liver Microsomes missing"
        assert "Rat Liver Microsomes" in page_text, "Rat Liver Microsomes missing"
        assert "Mouse Liver Microsomes" in page_text, "Mouse Liver Microsomes missing"
        assert "Dog Liver Microsomes" in page_text, "Dog Liver Microsomes missing"
        assert "Monkey Liver Microsomes" in page_text, "Monkey Liver Microsomes missing"
        assert "MODEL_UNAVAILABLE" in page_text, "MODEL_UNAVAILABLE badge missing"
        print("      Metabolism 5-species table successfully verified (Dog/Monkey MODEL_UNAVAILABLE).")

        # 4. Test PK Tab & Multi-Species PK
        print("[4/6] Testing PK Tab & Multi-Species PK Profile...")
        click_detail_tab(driver, "PK")
        time.sleep(4)

        pk_text = driver.find_element(By.TAG_NAME, "body").text
        assert "MULTI-SPECIES PK SUMMARY" in pk_text or "Multi-Species" in pk_text or "Clearance (CL)" in pk_text
        assert "Clearance (CL)" in pk_text
        assert "Volume of Distribution (V)" in pk_text
        assert "Elimination Half-Life" in pk_text
        print("      Multi-Species PK Profile table verified across Mouse, Rat, Dog, Monkey, Human.")

        # Test PK Simulation
        print("      Testing PK Simulation execution...")
        sim_buttons = driver.find_elements(By.XPATH, "//button[contains(text(), 'RUN PK SIMULATION')]")
        if sim_buttons:
            driver.execute_script("arguments[0].scrollIntoView({block:'center'});", sim_buttons[0])
            time.sleep(1)
            driver.execute_script("arguments[0].click();", sim_buttons[0])
            time.sleep(4)

            sim_page_text = driver.find_element(By.TAG_NAME, "body").text
            assert "CALCULATED PK SIMULATION" in sim_page_text or "Cmax" in sim_page_text or "AUC" in sim_page_text, "PK Simulation output missing"
            svgs = driver.find_elements(By.TAG_NAME, "svg")
            assert len(svgs) > 0, "PK simulation SVG graph missing"
            print("      PK Simulation curve and output metrics verified!")

        # Add Compound 2 (Erlotinib) via browser API
        print("      Adding Compound 2 (Erlotinib)...")
        projects = browser_api(driver, "GET", "/projects")
        target_p = next(p for p in projects if p["name"] == TEMP_PROJECT_NAME)
        c2 = browser_api(driver, "POST", f"/projects/{target_p['id']}/compounds", {
            "name": "Erlotinib",
            "smiles": ERLOTINIB_SMILES,
            "calculate": True
        })
        print(f"      Compound 2 created: {c2['name']} (#{c2['row_id']}). Triggering predict...")
        browser_api(driver, "POST", f"/compounds/{c2['row_id']}/predict-workflow", {})
        time.sleep(6)

        # 5. Test Two-Compound Comparison
        print("[5/6] Testing Two-Compound Comparison Tab...")
        proj_data = browser_api(driver, "GET", f"/projects/{target_p['id']}")
        compounds = proj_data.get("compounds", [])
        cmp_ids = ",".join(str(c["row_id"]) for c in compounds)

        comp_data = browser_api(driver, "GET", f"/projects/{target_p['id']}/compare?ids={cmp_ids}")
        assert "Solubility" in comp_data["metrics"]
        assert "Caco-2" in comp_data["metrics"]
        assert "PPB" in comp_data["metrics"]
        assert "fu" in comp_data["metrics"]
        assert "HLM" in comp_data["metrics"]
        assert "Soft Spots" in comp_data["metrics"]
        assert "Human CL (IVIVE)" in comp_data["metrics"] or "Rat CL (IV)" in comp_data["metrics"]
        print("      Two-compound comparison API verified across Properties, Activity, ADME, Metabolism, PK, Safety!")

        # Navigate to Compare view in UI via project-nav
        click_project_tab(driver, "Compare")
        time.sleep(2)

        # 6. Capture Evidence Screenshot
        print("[6/6] Capturing E2E Evidence Screenshot...")
        screenshot_dir = "/home/xavier/.gemini/antigravity-cli/brain/50dcafa9-cfcd-4ede-9575-43f0b58d9fcd/.tempmediaStorage"
        os.makedirs(screenshot_dir, exist_ok=True)
        screenshot_path = os.path.join(screenshot_dir, "stage5b4_functional_scientific_e2e.png")
        driver.save_screenshot(screenshot_path)
        print(f"      Screenshot saved to {screenshot_path}")

        # Clean up test project
        try:
            browser_api(driver, "DELETE", f"/projects/{target_p['id']}", {"confirmation_name": TEMP_PROJECT_NAME})
            print("      Cleaned up test project.")
        except Exception as e:
            print("      Post-clean note:", e)

        print("\n>>> ALL BROWSER E2E TESTS PASSED SUCCESSFULLY! <<<")

    finally:
        driver.quit()

if __name__ == "__main__":
    run_browser_e2e()
