"""
Browser E2E Verification for Project Recovery & Test Fixture Cleanup (v3.3.1)
Captures screenshots for Desktop (1440x900) and Mobile (390x844).
Verifies:
1. Projects view shows only 4 protected projects (DrugBank, AMYR, EGFR, GLP-1).
2. DrugBank is prominently badged as REFERENCE LIBRARY · GLOBAL MODEL DEVELOPMENT.
3. Clicking DrugBank opens 150 compound workspace with all table details.
4. Navigation and responsiveness.
"""

import os
import sys
import time
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

from scripts.e2e_reference_library_v3_3_3 import ChromedriverManager, create_driver

BASE_URL = "http://127.0.0.1:8765"
OUT_DIR = Path("validation/recovery_browser_e2e")
OUT_DIR.mkdir(parents=True, exist_ok=True)


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
            time.sleep(0.8)
            return
    raise ValueError(f"Nav item '{view_name}' not found")


def run_e2e():
    cd_mgr = ChromedriverManager(port=9515)
    cd_mgr.start()

    try:
        # 1. Desktop Test (1440x900)
        print("\n[Desktop 1440x900] Initializing browser...")
        driver = create_driver(1440, 900, port=9515)
        wait = WebDriverWait(driver, 20)
        try:
            driver.get(f"{BASE_URL}/")
            time.sleep(2)
            wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(text(),'Drug Optimization Platform')]")))

            select_nav_view(driver, wait, "Projects")
            time.sleep(1.5)

            # Capture Desktop Projects Table
            desktop_projects_img = OUT_DIR / "desktop_1440x900_projects.png"
            driver.save_screenshot(str(desktop_projects_img))
            print(f"  Captured: {desktop_projects_img}")

            # Check projects count in table
            project_rows = driver.find_elements(By.CSS_SELECTOR, "tr.dashboard-project")
            print(f"  Desktop table shows {len(project_rows)} project rows.")
            assert len(project_rows) == 4, f"Expected 4 projects in table, got {len(project_rows)}"

            # Verify DrugBank is present
            drugbank_row = driver.find_element(By.CSS_SELECTOR, "tr.reference-library-project-row")
            assert "DrugBank" in drugbank_row.text
            assert "REFERENCE LIBRARY" in drugbank_row.text
            print("  Verified DrugBank row badged and styled as REFERENCE LIBRARY.")

            # Click DrugBank to open compound workspace
            db_link = drugbank_row.find_element(By.CSS_SELECTOR, "button.project-link-title")
            driver.execute_script("arguments[0].click();", db_link)
            time.sleep(2.5)

            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
            desktop_drugbank_img = OUT_DIR / "desktop_1440x900_drugbank_workspace.png"
            driver.save_screenshot(str(desktop_drugbank_img))
            print(f"  Captured: {desktop_drugbank_img}")

            # Check compounds count in DrugBank
            compound_rows = driver.find_elements(By.CSS_SELECTOR, "table.project-status-table tbody tr")
            print(f"  DrugBank workspace rendered {len(compound_rows)} compound rows.")
            assert len(compound_rows) == 150, f"Expected 150 compounds in DrugBank, got {len(compound_rows)}"

        finally:
            driver.quit()

        # 2. Mobile Test (390x844)
        print("\n[Mobile 390x844] Initializing browser...")
        driver = create_driver(390, 844, port=9515)
        wait = WebDriverWait(driver, 20)
        try:
            driver.get(f"{BASE_URL}/")
            time.sleep(2)
            wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(text(),'Drug Optimization Platform')]")))

            select_nav_view(driver, wait, "Projects")
            time.sleep(1.5)

            mobile_projects_img = OUT_DIR / "mobile_390x844_projects.png"
            driver.save_screenshot(str(mobile_projects_img))
            print(f"  Captured: {mobile_projects_img}")

            project_rows = driver.find_elements(By.CSS_SELECTOR, "tr.dashboard-project")
            print(f"  Mobile table shows {len(project_rows)} project rows.")
            assert len(project_rows) == 4, f"Expected 4 projects in mobile view, got {len(project_rows)}"

            # Open DrugBank in mobile view
            drugbank_row = driver.find_element(By.CSS_SELECTOR, "tr.reference-library-project-row")
            db_link = drugbank_row.find_element(By.CSS_SELECTOR, "button.project-link-title")
            driver.execute_script("arguments[0].click();", db_link)
            time.sleep(2.5)

            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
            mobile_drugbank_img = OUT_DIR / "mobile_390x844_drugbank_workspace.png"
            driver.save_screenshot(str(mobile_drugbank_img))
            print(f"  Captured: {mobile_drugbank_img}")

        finally:
            driver.quit()

        print("\n=== Browser E2E Verification Succeeded 100%! ===")

    finally:
        cd_mgr.stop()


if __name__ == "__main__":
    run_e2e()
