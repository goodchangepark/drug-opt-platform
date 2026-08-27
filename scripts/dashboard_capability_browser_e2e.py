#!/usr/bin/env python3
"""Focused Chromium acceptance check for backend-driven Dashboard capabilities."""

import json
import os
import sys
import threading
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


ROOT = Path(__file__).resolve().parent.parent
BASE_URL = os.environ.get("DASHBOARD_E2E_BASE_URL", "http://127.0.0.1:8765")
sys.path.insert(0, str(ROOT))


def card_text(driver, title):
    heading = driver.find_element(By.XPATH, f"//article[contains(@class,'module-card')]//h3[normalize-space()='{title}']")
    return heading.find_element(By.XPATH, "ancestor::article[contains(@class,'module-card')]").text


def start_local_server():
    import uvicorn
    from backend.main import app

    server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=8767, log_level="warning"))
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    for _ in range(180):
        try:
            with urllib.request.urlopen("http://127.0.0.1:8767/api/health", timeout=2) as response:
                if response.status == 200:
                    return server, thread, "http://127.0.0.1:8767"
        except Exception:
            time.sleep(1)
    server.should_exit = True
    thread.join(timeout=10)
    raise RuntimeError("Dashboard E2E server did not become healthy")


def main():
    server = thread = None
    base_url = BASE_URL
    if os.environ.get("DASHBOARD_E2E_SELF_HOST") == "1":
        server, thread, base_url = start_local_server()
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for arg in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1920,2400"):
        options.add_argument(arg)
    options.set_capability("goog:loggingPrefs", {"browser": "ALL"})
    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    driver.set_script_timeout(60)
    checks = []
    try:
        driver.get(f"{base_url}/?dashboard-capability-sync={datetime.now(timezone.utc).timestamp()}")
        WebDriverWait(driver, 90).until(lambda current: "Available Scientific Modules" in current.page_source)
        WebDriverWait(driver, 90).until(lambda current: len(current.find_elements(By.CSS_SELECTOR, ".module-card")) == 7)

        cyp = card_text(driver, "CYP & Transporters")
        checks.append({"name": "CYP and transporter partial availability", "passed": all(text in cyp for text in ("PARTIAL", "CYP1A2 inhibitor", "CYP3A4 substrate", "P-gp inhibitor", "READY", "BCRP inhibitor", "MODEL UNAVAILABLE"))})

        safety = card_text(driver, "Safety / Toxicology")
        checks.append({"name": "Safety partial availability", "passed": all(text in safety for text in ("PARTIAL", "hERG", "Ames", "DILI", "Structural Alerts", "READY", "MODEL UNAVAILABLE"))})

        pk = card_text(driver, "PK / DMPK")
        pk_labels = ("Experimental PK Data Management", "NCA", "IVIVE / Hepatic Clearance", "Vd / Absorption Foundation", "IV Simulation", "PO / SC / IP Simulation", "Cross-Species Scaling", "Human Translational PK", "Prospective Prediction Freeze", "Retrospective Validation")
        checks.append({"name": "PK and DMPK ready", "passed": "READY" in pk and "PLANNED" not in pk and all(label in pk for label in pk_labels)})

        adme = card_text(driver, "ADME")
        checks.append({"name": "Availability and confidence rendered separately", "passed": all(text in adme for text in ("HLM", "RLM", "MLM", "READY", "Confidence: LOW", "Conformal:"))})

        api_check = driver.execute_async_script(
            """
const done=arguments[arguments.length-1];
fetch('/api/dashboard').then(r=>r.json()).then(data=>{
 const groups=Object.fromEntries(data.capability_summary.groups.map(group=>[group.title,group]));
 const models=Object.fromEntries(data.model_registry.map(model=>[model.endpoint,model]));
 done({source:data.capability_summary.source,cyp:groups['CYP & Transporters'].status,safety:groups['Safety / Toxicology'].status,pk:groups['PK / DMPK'].status,hlm:[models['HLM intrinsic clearance'].availability,models['HLM intrinsic clearance'].confidence]});
}).catch(error=>done({error:String(error)}));
"""
        )
        checks.append({"name": "Dashboard API and model registry consistency", "passed": api_check == {"source": "BACKEND_CAPABILITY_REGISTRY", "cyp": "PARTIAL", "safety": "PARTIAL", "pk": "READY", "hlm": ["READY", "LOW"]}})

        severe = [entry for entry in driver.get_log("browser") if entry["level"] in {"SEVERE", "ERROR"}]
        checks.append({"name": "Zero JavaScript console errors", "passed": not severe, "details": severe})
        driver.save_screenshot(str(ROOT / "validation/dashboard_capability_sync.png"))
    finally:
        driver.quit()
        if server is not None:
            server.should_exit = True
            thread.join(timeout=30)

    report = {
        "suite": "Dashboard Capability Status Sync Chromium E2E",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "all_passed": all(check["passed"] for check in checks),
        "total_checks": len(checks),
        "passed_checks": sum(check["passed"] for check in checks),
        "failed_checks": sum(not check["passed"] for check in checks),
        "checks": checks,
    }
    (ROOT / "validation/dashboard_capability_sync_results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    sys.exit(0 if report["all_passed"] else 1)


if __name__ == "__main__":
    main()
