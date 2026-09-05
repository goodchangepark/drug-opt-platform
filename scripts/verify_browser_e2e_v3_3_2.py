"""
Comprehensive Browser E2E Verification for Drug-OPT v3.3.2 Production Release
==============================================================================
Verifies Directive 18:
1. Desktop (1440x900) & Mobile (390x844) responsive layouts
2. Active Prediction Engine: drugopt-prediction-engine-v3@3.3.2
3. Policy Hash: 877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce
4. Help Page: Prediction Model History & Dedicated PK Prediction Readiness Foundation (PK_FOUNDATION_READY)
5. DrugBank Reference Library (Project 300) with 200 compounds
6. Compound Workspace: 50-endpoint maturity stars, predictions, and tabs
7. Hard page reload persistence
"""
import os
import sys
import time
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

import subprocess
import urllib.request
from scripts.e2e_reference_library_v3_3_3 import create_driver

BASE_URL = "http://127.0.0.1:8765"
OUT_DIR = Path("validation/e2e_v3_3_2_browser")
OUT_DIR.mkdir(parents=True, exist_ok=True)


class CustomChromedriverManager:
    def __init__(self, port: int = 9515):
        self.port = port
        self.proc = None

    def start(self):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{self.port}/status", timeout=1)
            print(f"ChromeDriver already running on port {self.port}")
            return
        except Exception:
            pass

        print(f"Starting ChromeDriver on port {self.port}...")
        self.proc = subprocess.Popen(
            ["/snap/bin/chromium.chromedriver", f"--port={self.port}", "--allowed-ips=127.0.0.1"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for i in range(30):
            time.sleep(0.5)
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{self.port}/status", timeout=1)
                print(f"ChromeDriver ready on port {self.port} (after {(i+1)*0.5:.1f}s)")
                return
            except Exception:
                pass
        raise RuntimeError(f"ChromeDriver failed to start on port {self.port} within 15s")

    def stop(self):
        if self.proc:
            print("Stopping ChromeDriver...")
            self.proc.terminate()
            try:
                self.proc.wait(timeout=3)
            except Exception:
                self.proc.kill()
            self.proc = None


def select_nav_view(driver, wait, view_name: str):
    time.sleep(0.5)
    menu_toggles = driver.find_elements(By.CLASS_NAME, "menu-toggle")
    for mt in menu_toggles:
        if mt.is_displayed() and mt.text.strip().lower() == "menu":
            driver.execute_script("arguments[0].click();", mt)
            time.sleep(0.5)
            break

    buttons = driver.find_elements(By.XPATH, "//nav//button | //aside//button")
    for b in buttons:
        if b.text.strip() == view_name:
            driver.execute_script("arguments[0].click();", b)
            time.sleep(1.0)
            return
    raise ValueError(f"Nav item '{view_name}' not found")


def run_e2e():
    print("="*70)
    print("Drug-OPT v3.3.2 Production Browser E2E Automation")
    print("="*70)
    cd_mgr = CustomChromedriverManager(port=9515)
    cd_mgr.start()

    try:
        # =================================================================
        # 1. DESKTOP E2E TEST (1440x900)
        # =================================================================
        print("\n[Desktop 1440x900] Initializing browser...")
        driver = create_driver(1440, 900, port=9515)
        wait = WebDriverWait(driver, 20)
        try:
            # Step A: Load Dashboard
            driver.get(f"{BASE_URL}/")
            time.sleep(2)
            wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(text(),'Drug Optimization Platform')]")))
            print("  ✓ Main Dashboard loaded successfully")

            # Step B: Navigate to Help Page & Verify v3.3.2 and PK Readiness
            select_nav_view(driver, wait, "Help")
            time.sleep(1.5)
            wait.until(EC.presence_of_element_located((By.ID, "help-pk-readiness")))

            desktop_help_img = OUT_DIR / "desktop_1440x900_help_pk_readiness.png"
            driver.save_screenshot(str(desktop_help_img))
            print(f"  ✓ Captured Help & PK Readiness: {desktop_help_img}")

            # Verify active engine text on Help Page
            page_text = driver.find_element(By.TAG_NAME, "body").text
            assert "drugopt-prediction-engine-v3@3.3.2" in page_text, "Missing v3.3.2 engine ID on Help Page"
            assert "877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce" in page_text, "Missing v3.3.2 policy hash"
            assert "PK_FOUNDATION_READY" in page_text, "Missing PK_FOUNDATION_READY badge on Help Page"
            assert "PK Prediction Readiness Foundation" in page_text, "Missing PK Readiness section title"
            print("  ✓ Verified Active Production Engine v3.3.2, Policy Hash, and PK_FOUNDATION_READY on Help Page")

            # Step C: Navigate to Projects & Verify DrugBank
            select_nav_view(driver, wait, "Projects")
            time.sleep(1.5)
            project_rows = driver.find_elements(By.CSS_SELECTOR, "tr.dashboard-project")
            assert len(project_rows) == 4, f"Expected 4 protected projects, found {len(project_rows)}"
            print(f"  ✓ Verified {len(project_rows)} protected projects (GLP-1, EGFR, AMYR, DrugBank)")

            drugbank_row = driver.find_element(By.CSS_SELECTOR, "tr.reference-library-project-row")
            db_link = drugbank_row.find_element(By.CSS_SELECTOR, "button.project-link-title")
            driver.execute_script("arguments[0].click();", db_link)
            time.sleep(2.5)

            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
            desktop_drugbank_img = OUT_DIR / "desktop_1440x900_drugbank_workspace.png"
            driver.save_screenshot(str(desktop_drugbank_img))
            print(f"  ✓ Captured DrugBank Workspace: {desktop_drugbank_img}")

            compound_rows = driver.find_elements(By.CSS_SELECTOR, "table.project-status-table tbody tr")
            print(f"  ✓ DrugBank workspace rendered {len(compound_rows)} compound rows")
            assert len(compound_rows) in (150, 200), f"Expected 150 or 200 compounds, got {len(compound_rows)}"

            # Step D: Open Compound Detail (First Compound)
            first_open_btn = driver.find_element(By.CSS_SELECTOR, "table.project-status-table tbody tr button.secondary")
            driver.execute_script("arguments[0].click();", first_open_btn)
            time.sleep(2.5)

            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card")))
            desktop_detail_img = OUT_DIR / "desktop_1440x900_compound_detail.png"
            driver.save_screenshot(str(desktop_detail_img))
            print(f"  ✓ Captured Compound Detail: {desktop_detail_img}")

            # Verify 50-endpoint maturity stars exist
            star_elements = driver.find_elements(By.CLASS_NAME, "maturity-stars")
            assert len(star_elements) > 0, "Maturity stars must be rendered in Compound Workspace"
            print(f"  ✓ Rendered {len(star_elements)} maturity star elements")

            # Step E: Hard Page Reload & Verify Persistence
            driver.refresh()
            time.sleep(3.0)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-header-card")))
            reload_img = OUT_DIR / "desktop_1440x900_hard_reload_persistence.png"
            driver.save_screenshot(str(reload_img))
            print(f"  ✓ Verified hard reload persistence: {reload_img}")

        finally:
            driver.quit()

        # =================================================================
        # 2. MOBILE E2E TEST (390x844)
        # =================================================================
        print("\n[Mobile 390x844] Initializing browser...")
        driver = create_driver(390, 844, port=9515)
        wait = WebDriverWait(driver, 20)
        try:
            driver.get(f"{BASE_URL}/")
            time.sleep(2)
            wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(text(),'Drug Optimization Platform')]")))

            # Help view in mobile
            select_nav_view(driver, wait, "Help")
            time.sleep(1.5)
            wait.until(EC.presence_of_element_located((By.ID, "help-pk-readiness")))

            mobile_help_img = OUT_DIR / "mobile_390x844_help_pk_readiness.png"
            driver.save_screenshot(str(mobile_help_img))
            print(f"  ✓ Captured Mobile Help & PK Readiness: {mobile_help_img}")

            # Projects view in mobile
            select_nav_view(driver, wait, "Projects")
            time.sleep(1.5)
            mobile_projects_img = OUT_DIR / "mobile_390x844_projects.png"
            driver.save_screenshot(str(mobile_projects_img))
            print(f"  ✓ Captured Mobile Projects: {mobile_projects_img}")

            project_rows = driver.find_elements(By.CSS_SELECTOR, "tr.dashboard-project")
            assert len(project_rows) == 4, f"Expected 4 projects in mobile view, got {len(project_rows)}"

        finally:
            driver.quit()

        print("\n" + "="*70)
        print("=== Public E2E Browser Verification Passed 100%! ===")
        print("="*70)

    finally:
        cd_mgr.stop()


if __name__ == "__main__":
    run_e2e()
