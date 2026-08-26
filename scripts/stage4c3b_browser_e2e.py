"""Chromium E2E Acceptance Test for Stage 4C-3B: Conformal Recalibration & Governance UI Verification."""

import datetime
import json
import os
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

PORT = os.environ.get("TEST_PORT", "8767")
BASE_URL = f"http://127.0.0.1:{PORT}"


def wait_for_server(url, max_wait=15):
    for _ in range(int(max_wait * 10)):
        try:
            with urllib.request.urlopen(f"{url}/api/dashboard", timeout=1) as resp:
                if resp.status == 200:
                    return True
        except Exception:
            time.sleep(0.1)
    return False


def api_get(endpoint_path):
    url = f"{BASE_URL}/api{endpoint_path}"
    req = urllib.request.urlopen(url)
    return json.loads(req.read().decode())


def run_e2e():
    print(f"Starting Stage 4C-3B Chromium E2E on {BASE_URL}...")
    server_proc = None

    if not wait_for_server(BASE_URL, max_wait=1):
        print(f"Launching test server on port {PORT}...")
        python_bin = sys.executable
        if not Path(python_bin).exists() or "venv" not in python_bin:
            python_bin = os.path.abspath(".venv/bin/python")

        server_proc = subprocess.Popen(
            [python_bin, "-m", "uvicorn", "backend.main:app", "--port", str(PORT), "--host", "127.0.0.1"],
            cwd=os.path.abspath("."),
        )
        if not wait_for_server(BASE_URL, max_wait=20):
            raise RuntimeError(f"Server failed to start on port {PORT}")
        print("Server is healthy and ready.")

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

        # 2. Verify Dashboard Model Registry via API (Decoupled Governance)
        dash_data = api_get("/dashboard")
        registry = {m["endpoint"]: m for m in dash_data.get("model_registry", [])}

        hlm_meta = registry.get("HLM intrinsic clearance", {})
        assert hlm_meta.get("status") == "READY", "HLM status must remain READY"
        assert hlm_meta.get("calibration_provenance") == "EXTERNAL", "HLM provenance must be EXTERNAL"
        assert hlm_meta.get("calibration_quality") == "UNDERCOVERED", "HLM quality must be UNDERCOVERED"

        herg_meta = registry.get("hERG liability", {})
        assert herg_meta.get("status") == "READY", "hERG status must remain READY"
        assert herg_meta.get("calibration_quality") == "UNDERCOVERED", "hERG quality must be UNDERCOVERED"

        cyp3a4_meta = registry.get("CYP3A4 inhibitor", {})
        assert cyp3a4_meta.get("calibration_provenance") == "EXTERNAL", "CYP3A4 provenance must be EXTERNAL"
        assert cyp3a4_meta.get("calibration_quality") == "VALIDATED", "CYP3A4 quality must be VALIDATED"

        cyp2d6_meta = registry.get("CYP2D6 inhibitor", {})
        assert cyp2d6_meta.get("calibration_provenance") == "INTERNAL", "CYP2D6 provenance must be INTERNAL"
        assert cyp2d6_meta.get("calibration_quality") == "VALIDATED", "CYP2D6 quality must be VALIDATED"

        checks.append({"name": "Model Registry Decoupled Governance (Provenance + Quality)", "status": "PASS"})

        # 3. Navigate into Settings and verify Model Registry columns
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "global-nav")))
        driver.execute_script("""
            const btns = Array.from(document.querySelectorAll('.global-nav button'));
            const setBtn = btns.find(b => b.textContent.includes('Settings'));
            if (setBtn) setBtn.click();
        """)
        time.sleep(2)

        page_text = driver.find_element(By.TAG_NAME, "body").text
        assert "Prediction Models" in page_text or "Default Workspace Settings" in page_text or "Project Settings" in page_text, "Settings page must render"
        assert "CALIBRATED_EXTERNAL" not in page_text, "CALIBRATED_EXTERNAL badge must not be shown"
        checks.append({"name": "Settings Model Registry Decoupled Columns", "status": "PASS"})

        # 4. Check Sidebar -> Optimization Workspace navigation (Regression check)
        driver.execute_script("""
            const btns = Array.from(document.querySelectorAll('.global-nav button'));
            const optBtn = btns.find(b => b.textContent.trim() === 'Optimization' || b.textContent.includes('Optimization'));
            if (optBtn) optBtn.click();
        """)
        time.sleep(2)

        wait.until(EC.presence_of_element_located((By.XPATH, "//h1[contains(., 'Optimization Workspace')]")))
        checks.append({"name": "Optimization Workspace Navigation", "status": "PASS"})

        # 5. Save Screenshot & Results
        base_dir = Path(__file__).resolve().parent.parent
        val_dir = base_dir / "validation"
        val_dir.mkdir(parents=True, exist_ok=True)
        screenshot_path = str(val_dir / "stage4c3b_browser_e2e.png")
        driver.save_screenshot(screenshot_path)

        result = {
            "stage": "4C-3B",
            "run_id": datetime.datetime.now().strftime("%Y%m%d-%H%M%S"),
            "base_url": BASE_URL,
            "checks": checks,
            "screenshot": screenshot_path,
            "status": "PASS"
        }

        out_json = val_dir / "stage4c3b_browser_e2e_results.json"
        with open(out_json, "w") as f:
            json.dump(result, f, indent=2)

        print(json.dumps(result, indent=2))
        return 0

    except Exception as exc:
        base_dir = Path(__file__).resolve().parent.parent
        val_dir = base_dir / "validation"
        val_dir.mkdir(parents=True, exist_ok=True)
        screenshot_path = str(val_dir / "stage4c3b_browser_e2e_failure.png")
        driver.save_screenshot(screenshot_path)
        driver.save_screenshot(screenshot_path)
        print(f"E2E Failure: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1
    finally:
        driver.quit()
        if server_proc:
            server_proc.terminate()


if __name__ == "__main__":
    sys.exit(run_e2e())
