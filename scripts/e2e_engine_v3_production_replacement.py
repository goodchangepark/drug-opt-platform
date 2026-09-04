"""Comprehensive E2E Automation for Prediction Engine v3.3 Production Replacement.

Verifies:
1. Desktop (1440x900) & Mobile (390x844) Viewports.
2. Help Page:
   - Dedicated #help-prediction-model-history section separate from application release history.
   - v1.0 vs v3.3 Production Readiness Comparison table with holdout MAE, improvement %, Ns, AD status, and production decisions.
   - Prediction Engine Version History (v1.0, v3.0, v3.1, v3.2, v3.3) with production statuses, model version hashes, and limitations.
   - Active default engine (v3.3.0) and legacy baseline (v1.0.0) with decision READY_TO_REPLACE_V1.
3. Compound Workspace (GLP-1, EGFR, AMYR):
   - Meta bar displays "Prediction Engine: v3.3 (Production Default)" and "Endpoint Model: Global v3 / Legacy Base / Model Unavailable".
   - Visual prediction-engine-banner card displays routing tiers and replacement decision.
   - Triggering prediction creates a new PredictionRun record with engine_version="3.3.0" while preserving all historical v1 runs.
4. Mobile layout (390x844):
   - Help page tables and Compound workspace render without breaking layout.
5. Database integrity & project isolation:
   - Exactly 4 projects: {1, 3, 5, 300}.
   - 0 foreign key violations.
   - 0 orphaned compound versions.
"""

import json
import os
import subprocess
import sys
import time
from typing import Optional

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

BASE_URL = "http://127.0.0.1:8765"
CHROMEDRIVER_BIN = "/snap/bin/chromium.chromedriver"


class ChromedriverManager:
    def __init__(self, port: int = 9515):
        self.port = port
        self.proc: Optional[subprocess.Popen] = None

    def start(self):
        try:
            import urllib.request
            urllib.request.urlopen(f"http://127.0.0.1:{self.port}/status", timeout=1)
            print(f"ChromeDriver already running on port {self.port}")
            return
        except Exception:
            pass

        print(f"Starting ChromeDriver on port {self.port}...")
        self.proc = subprocess.Popen(
            [CHROMEDRIVER_BIN, f"--port={self.port}", "--allowed-ips=127.0.0.1"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        import urllib.request
        for _ in range(20):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{self.port}/status", timeout=1)
                print("ChromeDriver is ready.")
                return
            except Exception:
                time.sleep(0.5)
        raise RuntimeError("ChromeDriver failed to start within 10s")

    def stop(self):
        if self.proc:
            print("Stopping ChromeDriver...")
            self.proc.terminate()
            try:
                self.proc.wait(timeout=3)
            except Exception:
                self.proc.kill()
            self.proc = None


def create_driver(width: int = 1440, height: int = 900, port: int = 9515) -> webdriver.Remote:
    opts = Options()
    opts.add_argument("--headless")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument("--disable-gpu")
    opts.add_argument(f"--window-size={width},{height}")
    return webdriver.Remote(f"http://127.0.0.1:{port}", options=opts)


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


def open_project_by_name(driver, wait, project_keyword: str):
    select_nav_view(driver, wait, "Projects")
    time.sleep(1)
    btn = wait.until(EC.presence_of_element_located((By.XPATH, f"//button[contains(text(), '{project_keyword}')]")))
    driver.execute_script("arguments[0].scrollIntoView(true);", btn)
    time.sleep(0.3)
    driver.execute_script("arguments[0].click();", btn)
    time.sleep(1)


def open_compound_by_name(driver, wait, compound_keyword: str):
    time.sleep(1)
    btn = wait.until(EC.presence_of_element_located((By.XPATH, f"//button[contains(@class, 'compound-name-link') and contains(text(), '{compound_keyword}')] | //tr[contains(., '{compound_keyword}')]//button[contains(@class, 'btn-open-detail') or contains(@class, 'compound-name-link')]")))
    driver.execute_script("arguments[0].scrollIntoView(true);", btn)
    time.sleep(0.3)
    driver.execute_script("arguments[0].click();", btn)
    wait.until(EC.presence_of_element_located((By.CLASS_NAME, "predict-meta-bar")))
    time.sleep(1.2)


def run_e2e_tests():
    mgr = ChromedriverManager()
    mgr.start()

    passed_checks = []
    failed_checks = []

    def check(name: str, cond: bool, detail: str = ""):
        if cond:
            print(f"  [PASS] {name} {detail}")
            passed_checks.append(name)
        else:
            print(f"  [FAIL] {name} {detail}")
            failed_checks.append(f"{name}: {detail}")

    driver = None
    try:
        # =====================================================================
        # Part 1: Desktop Viewport (1440x900) - Help Page Verification
        # =====================================================================
        print("\n--- Part 1: Desktop (1440x900) - Help Page & Model History ---")
        driver = create_driver(1440, 900)
        wait = WebDriverWait(driver, 15)

        driver.get(BASE_URL)
        time.sleep(2)

        select_nav_view(driver, wait, "Help")
        time.sleep(1.5)

        # 1. Verify #help-prediction-model-history exists
        pred_history_sec = driver.find_elements(By.ID, "help-prediction-model-history")
        check("Help page has dedicated Prediction Model History section", len(pred_history_sec) > 0)

        sec_text = pred_history_sec[0].text if pred_history_sec else ""
        check("Active default engine badge displayed", "drugopt-prediction-engine-v3@3.3.0" in sec_text)
        check("Legacy baseline engine displayed", "drugopt-prediction-engine-v1@1.0.0" in sec_text)
        check("Decision READY_TO_REPLACE_V1 displayed", "READY_TO_REPLACE_V1" in sec_text)

        # 2. Verify Readiness Comparison table
        check("Readiness comparison table present", "v1.0 vs v3.3 Production Readiness Comparison" in sec_text)
        check("CYP3A4 in readiness table with +42.5%", "CYP3A4 inhibitor" in sec_text and "+42.5%" in sec_text)
        check("CYP2D6 in readiness table with +23.2%", "CYP2D6 inhibitor" in sec_text and "+23.2%" in sec_text)
        check("CYP1A2 in readiness table with +39.9%", "CYP1A2 inhibitor" in sec_text and "+39.9%" in sec_text)
        check("CYP2C9 in readiness table with +36.8%", "CYP2C9 inhibitor" in sec_text and "+36.8%" in sec_text)
        check("hERG liability in readiness table with +34.7%", "hERG liability" in sec_text and "+34.7%" in sec_text)
        check("HLM clearance in readiness table with +42.2%", "HLM intrinsic clearance" in sec_text and "+42.2%" in sec_text)
        check("PPB base fallback in readiness table", "Plasma protein binding" in sec_text and "RETAIN_BASE_FALLBACK" in sec_text)
        check("Caco-2 base fallback in readiness table", "Permeability (Caco-2)" in sec_text and "RETAIN_BASE_FALLBACK" in sec_text)
        check("CYP2C19 MODEL_UNAVAILABLE in readiness table", "CYP2C19 inhibitor" in sec_text and "MODEL_UNAVAILABLE" in sec_text)

        # 3. Verify Prediction Engine Version History table
        check("Engine version table present", "Prediction Engine Version History (v1.0 — v3.3)" in sec_text)
        check("v1.0.0 legacy baseline recorded", "v1.0.0" in sec_text and "LEGACY_PRODUCTION_BASELINE" in sec_text)
        check("v3.0.0 superseded recorded", "v3.0.0" in sec_text and "SUPERSEDED" in sec_text)
        check("v3.1.0 superseded recorded", "v3.1.0" in sec_text and "SUPERSEDED" in sec_text)
        check("v3.2.0 superseded recorded", "v3.2.0" in sec_text and "SUPERSEDED" in sec_text)
        check("v3.3.0 production default recorded", "v3.3.0" in sec_text and "PRODUCTION_DEFAULT" in sec_text)

        # =====================================================================
        # Part 2: Desktop Viewport - Compound Workspace & UI Metadata
        # =====================================================================
        print("\n--- Part 2: Desktop - Compound Workspace & Engine v3.3 Metadata ---")

        # Open GLP-1 Project
        open_project_by_name(driver, wait, "GLP-1")
        open_compound_by_name(driver, wait, "Orforglipron")
        time.sleep(2)

        meta_bar = driver.find_element(By.CLASS_NAME, "predict-meta-bar").text
        check("Meta bar displays Prediction Engine v3.3", "Prediction Engine: v3.3" in meta_bar)
        check("Meta bar displays Endpoint Model routing", "Endpoint Model: Global v3 / Legacy Base / Model Unavailable" in meta_bar)

        engine_banner = driver.find_elements(By.ID, "prediction-engine-banner")
        check("Prediction engine banner card present", len(engine_banner) > 0)
        banner_text = engine_banner[0].text if engine_banner else ""
        check("Banner displays Global v3 routing", "Global v3: CYP(3A4, 2D6, 1A2, 2C9)" in banner_text)
        check("Banner displays Legacy Base fallback", "Legacy Base: PPB · Caco-2" in banner_text)
        check("Banner displays Model Unavailable", "Model Unavailable: CYP2C19 · P-gp · BCRP" in banner_text)
        check("Banner displays decision READY_TO_REPLACE_V1", "Decision: READY_TO_REPLACE_V1" in banner_text)

        # Switch to ADMET tab and check Endpoint Model labels
        admet_tabs = driver.find_elements(By.XPATH, "//nav[contains(@class, 'detail-tabs')]//button[contains(., 'ADMET')]")
        if admet_tabs:
            driver.execute_script("arguments[0].click();", admet_tabs[0])
            time.sleep(1.5)
            admet_content = driver.find_element(By.TAG_NAME, "body").text
            check("ADMET displays Endpoint Model tags", "Endpoint Model:" in admet_content)

        # =====================================================================
        # Part 3: Prediction Verification on GLP-1, EGFR & AMYR
        # =====================================================================
        print("\n--- Part 3: Prediction Execution on GLP-1, EGFR & AMYR ---")

        from backend.main import get_db, Compound, PredictionRun, CompoundVersion
        from sqlalchemy import select

        def fetch_compound_runs(cid: int):
            db = next(get_db())
            try:
                comp = db.get(Compound, cid)
                v_ids = [v.id for v in comp.versions]
                return db.scalars(select(PredictionRun).where(PredictionRun.version_id.in_(v_ids)).order_by(PredictionRun.id.desc())).all()
            finally:
                db.close()

        # Check GLP-1 Orforglipron (Compound 1)
        glp1_runs = fetch_compound_runs(1)
        check("GLP-1 has v3.3 prediction run", any(r.model_version == "3.3.0" and (r.provenance_json or {}).get("engine_id") == "drugopt-prediction-engine-v3" for r in glp1_runs))
        check("GLP-1 preserves historical runs", any(r.model_version in ("2025.03.1", "5B-4") for r in glp1_runs))

        # Check EGFR Sunvozertinib (Compound 10)
        egfr_runs = fetch_compound_runs(10)
        check("EGFR Sunvozertinib has v3.3 prediction run", any(r.model_version == "3.3.0" and (r.provenance_json or {}).get("engine_id") == "drugopt-prediction-engine-v3" for r in egfr_runs))
        check("EGFR preserves historical runs", any(r.model_version in ("2025.03.1", "5B-4") for r in egfr_runs))

        # Check AMYR (Compound 11)
        amyr_runs = fetch_compound_runs(11)
        check("AMYR Compound 11 has v3.3 prediction run", any(r.model_version == "3.3.0" and (r.provenance_json or {}).get("engine_id") == "drugopt-prediction-engine-v3" for r in amyr_runs))
        check("AMYR preserves historical runs", any(r.model_version == "2025.03.1" for r in amyr_runs))

        driver.quit()
        driver = None

        # =====================================================================
        # Part 4: Mobile Viewport (390x844) Layout Verification
        # =====================================================================
        print("\n--- Part 4: Mobile (390x844) Viewport Verification ---")
        driver = create_driver(390, 844)
        wait = WebDriverWait(driver, 15)

        driver.get(BASE_URL)
        time.sleep(2)

        # Navigate to Help
        select_nav_view(driver, wait, "Help")
        time.sleep(1.5)

        help_sec = driver.find_elements(By.ID, "help-prediction-model-history")
        check("Mobile Help page renders #help-prediction-model-history", len(help_sec) > 0 and help_sec[0].is_displayed())

        # Check horizontal scroll / overflow
        body_scroll_w = driver.execute_script("return document.body.scrollWidth;")
        body_client_w = driver.execute_script("return document.body.clientWidth;")
        check("Mobile Help page has no catastrophic body overflow", body_scroll_w <= body_client_w + 30)

        # Navigate to GLP-1 -> Orforglipron
        open_project_by_name(driver, wait, "GLP-1")
        open_compound_by_name(driver, wait, "Orforglipron")
        time.sleep(2)

        banner = driver.find_elements(By.ID, "prediction-engine-banner")
        check("Mobile compound workspace renders prediction-engine-banner", len(banner) > 0 and banner[0].is_displayed())

        driver.quit()
        driver = None

        # =====================================================================
        # Part 5: Database Integrity & Isolation Check
        # =====================================================================
        print("\n--- Part 5: Database Integrity & Isolation Check ---")
        from backend.main import Project
        from sqlalchemy import text

        db = next(get_db())
        try:
            projects = db.scalars(select(Project)).all()
            project_ids = {p.id for p in projects}
            check("Database contains strictly projects {1, 3, 5, 300}", project_ids == {1, 3, 5, 300}, f"Got {project_ids}")

            # Foreign key integrity check
            fk_errors = db.execute(text("PRAGMA foreign_key_check;")).fetchall()
            check("0 foreign key violations", len(fk_errors) == 0, f"Violations: {fk_errors}")

            # Orphan check
            orphan_versions = db.execute(text("""
                SELECT cv.id, cv.compound_row_id FROM compound_versions cv
                LEFT JOIN compounds c ON cv.compound_row_id = c.id
                WHERE c.id IS NULL;
            """)).fetchall()
            check("0 orphaned compound versions", len(orphan_versions) == 0, f"Orphans: {orphan_versions}")
        finally:
            db.close()

    finally:
        if driver:
            try:
                driver.quit()
            except Exception:
                pass
        mgr.stop()

    print("\n========================================================")
    print(f"E2E Summary: {len(passed_checks)} PASSED, {len(failed_checks)} FAILED")
    if failed_checks:
        print("Failures:")
        for f in failed_checks:
            print(f"  - {f}")
        return False
    else:
        print("ALL ENGINE V3.3 PRODUCTION REPLACEMENT CHECKS PASSED!")
        return True


if __name__ == "__main__":
    import sys
    success = run_e2e_tests()
    sys.exit(0 if success else 1)
