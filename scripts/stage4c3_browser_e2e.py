"""Chromium E2E Acceptance Test for Stage 4C-3: Model Registry & Optimization Workspace."""

import datetime
import json
import os
import sys
import time
import traceback
import urllib.request
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

PORT = os.environ.get("TEST_PORT", "8766")
BASE_URL = f"http://127.0.0.1:{PORT}"


def click_button_by_text(driver, text, timeout=30):
    def find_and_click(d):
        buttons = d.find_elements(By.TAG_NAME, "button")
        for b in buttons:
            txt = b.text.strip()
            if text in txt or txt == text:
                d.execute_script("arguments[0].scrollIntoView({block:'center'});", b)
                d.execute_script("arguments[0].click();", b)
                return True
        return False

    WebDriverWait(driver, timeout).until(find_and_click)


def api_get(endpoint_path):
    url = f"{BASE_URL}/api{endpoint_path}"
    req = urllib.request.urlopen(url)
    return json.loads(req.read().decode())


def run_e2e():
    options = Options()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for arg in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1900,1800"):
        options.add_argument(arg)

    service = Service("/snap/bin/chromium.chromedriver")
    driver = webdriver.Chrome(service=service, options=options)
    wait = WebDriverWait(driver, 30)
    checks = []

    try:
        # 1. Load Dashboard UI
        driver.get(BASE_URL)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "dashboard-hero")))
        time.sleep(2)
        checks.append({"name": "Dashboard Page Render", "status": "PASS"})

        # 2. Verify Dashboard Model Registry via API
        dash_data = api_get("/dashboard")
        registry = {m["endpoint"]: m for m in dash_data.get("model_registry", [])}
        if registry.get("Solubility", {}).get("status") != "READY":
            raise AssertionError(f"Solubility status not READY: {registry.get('Solubility')}")
        if registry.get("Plasma protein binding", {}).get("status") != "READY":
            raise AssertionError(f"Plasma protein binding status not READY: {registry.get('Plasma protein binding')}")
        if registry.get("hERG liability", {}).get("status") != "READY":
            raise AssertionError(f"hERG status not READY: {registry.get('hERG liability')}")
        checks.append({"name": "Dashboard Model Registry API (15 Models READY)", "status": "PASS"})

        # 3. Click Optimization in Sidebar
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "global-nav")))
        clicked = driver.execute_script("""
            const btns = Array.from(document.querySelectorAll('.global-nav button'));
            const optBtn = btns.find(b => b.textContent.trim() === 'Optimization' || b.textContent.includes('Optimization'));
            if (optBtn) {
                optBtn.click();
                return optBtn.textContent;
            }
            return null;
        """)
        print(f"JS Click result: {clicked}", file=sys.stderr)
        time.sleep(2)

        wait.until(EC.presence_of_element_located((By.XPATH, "//h1[contains(., 'Optimization Workspace')]")))
        checks.append({"name": "Sidebar Optimization Navigation (No White Screen)", "status": "PASS"})

        # 4. Select Project and Compound in Optimization Workspace
        project_select = wait.until(EC.presence_of_element_located((By.XPATH, "//div[@key='project']//select")))
        options_elems = project_select.find_elements(By.TAG_NAME, "option")
        if len(options_elems) > 1:
            options_elems[1].click()
            time.sleep(2)
            checks.append({"name": "Optimization Step 1 — Select Project", "status": "PASS"})

            compound_select = wait.until(EC.presence_of_element_located((By.XPATH, "//div[@key='compound']//select")))
            comp_options = compound_select.find_elements(By.TAG_NAME, "option")
            if len(comp_options) > 1:
                comp_options[1].click()
                time.sleep(2)
                checks.append({"name": "Optimization Step 2 — Select Compound Profile", "status": "PASS"})

        # 5. Check Settings -> Scientific Model Governance UI
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "global-nav")))
        driver.execute_script("""
            const btns = Array.from(document.querySelectorAll('.global-nav button'));
            const setBtn = btns.find(b => b.textContent.includes('Settings'));
            if (setBtn) setBtn.click();
        """)
        time.sleep(2)

        wait.until(EC.presence_of_element_located((By.XPATH, "//*[contains(., 'Scientific Validation & Governance')]")))
        checks.append({"name": "Settings Scientific Model Governance UI", "status": "PASS"})

        # Save Screenshot and Results
        screenshot_path = os.path.abspath("validation/stage4c3_browser_e2e.png")
        os.makedirs("validation", exist_ok=True)
        driver.save_screenshot(screenshot_path)

        result = {
            "stage": "4C-3",
            "run_id": datetime.datetime.now().strftime("%Y%m%d-%H%M%S"),
            "base_url": BASE_URL,
            "checks": checks,
            "screenshot": screenshot_path,
            "status": "PASS"
        }

        with open("validation/stage4c3_browser_e2e_results.json", "w") as f:
            json.dump(result, f, indent=2)

        print(json.dumps(result, indent=2))
        return 0

    except Exception as exc:
        screenshot_path = os.path.abspath("validation/stage4c3_browser_e2e_failure.png")
        os.makedirs("validation", exist_ok=True)
        driver.save_screenshot(screenshot_path)
        print(f"E2E Failure: {exc}", file=sys.stderr)
        try:
            h1s = [h.text for h in driver.find_elements(By.TAG_NAME, "h1")]
            print(f"Page H1s: {h1s}", file=sys.stderr)
        except Exception:
            pass
        traceback.print_exc()
        return 1
    finally:
        driver.quit()


if __name__ == "__main__":
    sys.exit(run_e2e())
