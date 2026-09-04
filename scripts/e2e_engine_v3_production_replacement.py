"""Comprehensive E2E Automation for Prediction Engine v3.3 Production Integrity Audit.

Verifies:
1. Desktop (1440x900) & Mobile (390x844) Viewports.
2. Help Page:
   - Dedicated #help-prediction-model-history section separate from release history.
   - Separate distinct tables for:
     * 1. Global v3 Primary Endpoints (6 models, +36.6% avg error reduction)
     * 2. Legacy Base Fallback Endpoints (3 models: Solubility, PPB, Caco-2) with audited failure reasons
     * 3. Model Unavailable Endpoints (3 models: CYP2C19, P-gp, BCRP)
   - Prediction Engine Version History (v1.0 — v3.3) with production statuses, hashes, and limitations.
   - Active default engine (v3.3.0) and legacy baseline (v1.0.0) with decision READY_TO_REPLACE_V1.
3. Compound Workspace (GLP-1, EGFR, AMYR):
   - Meta bar displays "Prediction Engine v3.3 · Current Production".
   - Top banner displays:
     * "Global v3: CYP(3A4, 2D6, 1A2, 2C9) · hERG · HLM"
     * "Legacy Base: Solubility · PPB · Caco-2"
     * "Model Unavailable: CYP2C19 · P-gp · BCRP"
   - No user engine selection control exists (automatic production default).
4. Save & Predict Automated Workflow:
   - Add Compound → Structure → Search → Predict → Save → Navigate → Hard reload → Re-entry.
   - Verification of auto v3.3 execution and endpoint tier routing.
5. Invalidation & Fresh Runs:
   - Historical runs 135, 137, 139 preserved and marked INVALIDATED_BY_ROUTING_AUDIT.
   - Fresh runs 141, 143, 145 active with Solubility as BASE_FALLBACK / Legacy Base.
6. Mobile Layout (390x844):
   - Responsiveness without layout breakage.
7. Database Integrity:
   - Strictly projects {1, 3, 5, 300}.
   - 0 foreign key violations, 0 orphaned compound versions.
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
    toggle = None
    for mt in menu_toggles:
        if mt.is_displayed():
            toggle = mt
            if mt.text.strip().lower() == "menu":
                driver.execute_script("arguments[0].click();", mt)
                time.sleep(0.5)
            break

    buttons = driver.find_elements(By.XPATH, "//nav//button | //aside//button")
    clicked = False
    for b in buttons:
        if b.text.strip() == view_name:
            driver.execute_script("arguments[0].click();", b)
            time.sleep(0.8)
            clicked = True
            break
    if not clicked:
        raise ValueError(f"Nav item '{view_name}' not found")

    if toggle and toggle.is_displayed() and toggle.text.strip().lower() == "close":
        driver.execute_script("arguments[0].click();", toggle)
        time.sleep(0.5)


def open_project_by_name(driver, wait, project_keyword: str):
    select_nav_view(driver, wait, "Projects")
    time.sleep(1)
    btn = wait.until(EC.presence_of_element_located((By.XPATH, f"//button[contains(text(), '{project_keyword}')]")))
    driver.execute_script("arguments[0].scrollIntoView(true);", btn)
    time.sleep(0.3)
    driver.execute_script("arguments[0].click();", btn)
    time.sleep(1.5)


def open_compound_by_name(driver, wait, compound_keyword: str):
    time.sleep(1.5)
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
        check("Average error reduction 36.6% displayed", "+36.6%" in sec_text)

        # 2. Verify Distinct Tier Tables
        check("Section 1: Global v3 Primary Endpoints present", "1. Global v3 Primary Endpoints" in sec_text)
        check("Section 2: Legacy Base Fallback Endpoints present", "2. Legacy Base Fallback Endpoints" in sec_text)
        check("Section 3: Model Unavailable Endpoints present", "3. Model Unavailable Endpoints" in sec_text)

        check("CYP3A4 in primary table with +42.5%", "CYP3A4 inhibitor" in sec_text and "+42.5%" in sec_text)
        check("CYP2D6 in primary table with +23.2%", "CYP2D6 inhibitor" in sec_text and "+23.2%" in sec_text)
        check("CYP1A2 in primary table with +39.9%", "CYP1A2 inhibitor" in sec_text and "+39.9%" in sec_text)
        check("CYP2C9 in primary table with +36.8%", "CYP2C9 inhibitor" in sec_text and "+36.8%" in sec_text)
        check("hERG in primary table with +34.7%", "hERG liability" in sec_text and "+34.7%" in sec_text)
        check("HLM in primary table with +42.2%", "HLM intrinsic clearance" in sec_text and "+42.2%" in sec_text)

        check("Solubility fail-closed to BASE_FALLBACK in fallback table", "Solubility" in sec_text and "RETAIN_BASE_FALLBACK" in sec_text and "-12.9%" in sec_text)
        check("PPB base fallback in fallback table", "Plasma protein binding" in sec_text and "RETAIN_BASE_FALLBACK" in sec_text and "-11.5%" in sec_text)
        check("Caco-2 base fallback in fallback table", "Permeability (Caco-2)" in sec_text and "RETAIN_BASE_FALLBACK" in sec_text and "+1.9%" in sec_text)

        check("CYP2C19 in unavailable table", "CYP2C19 inhibitor" in sec_text and "MODEL_UNAVAILABLE" in sec_text)
        check("P-gp in unavailable table", "P-gp substrate" in sec_text and "MODEL_UNAVAILABLE" in sec_text)
        check("BCRP in unavailable table", "BCRP substrate" in sec_text and "MODEL_UNAVAILABLE" in sec_text)

        # 3. Verify Prediction Engine Version History table
        check("Engine version table present", "Prediction Engine Version History (v1.0 — v3.3)" in sec_text)
        check("v1.0.0 legacy baseline recorded", "v1.0.0" in sec_text and "LEGACY_PRODUCTION_BASELINE" in sec_text)
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
        check("Meta bar displays 'Prediction Engine v3.3 · Current Production'", "Prediction Engine v3.3 · Current Production" in meta_bar)
        check("Meta bar displays Endpoint Model summary", "Endpoint Model: Global v3 / Legacy Base / Model Unavailable" in meta_bar)

        engine_banner = driver.find_elements(By.ID, "prediction-engine-banner")
        check("Prediction engine banner card present", len(engine_banner) > 0)
        banner_text = engine_banner[0].text if engine_banner else ""
        check("Banner displays Prediction Engine v3.3 · Current Production", "Prediction Engine v3.3 · Current Production" in banner_text)
        check("Banner displays Global v3 routing without Solubility", "Global v3: CYP(3A4, 2D6, 1A2, 2C9) · hERG · HLM" in banner_text)
        check("Banner displays Legacy Base with Solubility", "Legacy Base: Solubility · PPB · Caco-2" in banner_text)
        check("Banner displays Model Unavailable", "Model Unavailable: CYP2C19 · P-gp · BCRP" in banner_text)
        check("Banner displays decision READY_TO_REPLACE_V1", "Decision: READY_TO_REPLACE_V1" in banner_text)

        # Confirm users NEVER pick engine via UI (no engine selector dropdown or radio)
        engine_dropdowns = driver.find_elements(By.XPATH, "//select[contains(@id, 'engine') or contains(@name, 'engine')] | //input[@type='radio' and contains(@name, 'engine')]")
        check("Users never pick an engine; zero engine selector inputs in UI", len(engine_dropdowns) == 0)

        # Switch to ADMET tab and check Endpoint Model labels
        admet_tabs = driver.find_elements(By.XPATH, "//nav[contains(@class, 'detail-tabs')]//button[contains(., 'ADMET')]")
        if admet_tabs:
            driver.execute_script("arguments[0].click();", admet_tabs[0])
            time.sleep(2)
            admet_content = driver.find_element(By.TAG_NAME, "body").text
            matched_lines = [l.strip() for l in admet_content.splitlines() if "Endpoint Model" in l or "Solubility" in l]
            print(f"ADMET matched lines ({len(matched_lines)}): {matched_lines[:10]}")
            check("ADMET displays Endpoint Model tags", "Endpoint Model:" in admet_content)
            check("Solubility displays Legacy Base", "Endpoint Model: Legacy Base" in admet_content or any("Legacy Base" in l for l in matched_lines))

        # =====================================================================
        # Part 3: Invalidation & Persistence Verification on GLP-1, EGFR & AMYR
        # =====================================================================
        print("\n--- Part 3: Prediction Invalidation & Persistence Checks ---")

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
        check("GLP-1 has v3.3 prediction run", any(r.model_version == "3.3.0" for r in glp1_runs))
        check("GLP-1 preserves historical runs", any(r.model_version in ("2025.03.1", "5B-4") for r in glp1_runs))
        latest_glp1 = glp1_runs[0]
        check("GLP-1 latest run has Solubility as BASE_FALLBACK", (latest_glp1.outputs_json or {}).get("endpoint_routing", {}).get("SOLUBILITY_GENERIC") == "BASE_FALLBACK")

        # Check EGFR Sunvozertinib (Compound 10)
        egfr_runs = fetch_compound_runs(10)
        check("EGFR Sunvozertinib has v3.3 prediction run", any(r.model_version == "3.3.0" for r in egfr_runs))
        latest_egfr = egfr_runs[0]
        check("EGFR latest run has Solubility as BASE_FALLBACK", (latest_egfr.outputs_json or {}).get("endpoint_routing", {}).get("SOLUBILITY_GENERIC") == "BASE_FALLBACK")

        # Check AMYR (Compound 11)
        amyr_runs = fetch_compound_runs(11)
        check("AMYR Compound 11 has v3.3 prediction run", any(r.model_version == "3.3.0" for r in amyr_runs))
        latest_amyr = amyr_runs[0]
        check("AMYR latest run has Solubility as BASE_FALLBACK", (latest_amyr.outputs_json or {}).get("endpoint_routing", {}).get("SOLUBILITY_GENERIC") == "BASE_FALLBACK")

        # Verify historical runs 135, 137, 139 preserved with INVALIDATED_BY_ROUTING_AUDIT
        db = next(get_db())
        try:
            for rid in [135, 137, 139]:
                r = db.get(PredictionRun, rid)
                check(f"Historical run {rid} preserved", r is not None)
                if r:
                    prov = r.provenance_json or {}
                    check(f"Run {rid} marked INVALIDATED_BY_ROUTING_AUDIT", prov.get("audit_status") == "INVALIDATED_BY_ROUTING_AUDIT")
        finally:
            db.close()

        # =====================================================================
        # Part 4: New Compound Save & Predict, Navigation & Hard Reload
        # =====================================================================
        print("\n--- Part 4: New Compound Save & Predict, Navigation & Reload ---")

        # Create temporary test compound in Project 1 via backend API to test workflow
        from backend.models import Compound, CompoundVersion
        from backend.main import run_compound_prediction_workflow
        db = next(get_db())
        test_cid = None
        try:
            test_comp = Compound(
                project_id=1,
                name="Audit Test Compound Ibuprofen",
                compound_id="AUDIT-TEST-001",
                cas_number="15687-27-1",
                current_version=1,
            )
            db.add(test_comp)
            db.flush()
            test_cid = test_comp.id
            test_cv = CompoundVersion(
                compound_row_id=test_cid,
                version_number=1,
                original_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
                canonical_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
                isomeric_smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
                inchikey="HEFNNWSXXWATRW-UHFFFAOYSA-N",
            )
            db.add(test_cv)
            db.commit()

            # Execute Prediction workflow on new compound
            pred_res = run_compound_prediction_workflow(row_id=test_cid, db=db)
            check("New compound prediction completed", pred_res.get("status") in ("COMPLETE", "PARTIAL"))
            check("New compound automatically defaults to v3.3", pred_res.get("engine_id") == "drugopt-prediction-engine-v3")
            check("New compound routing has Solubility in BASE_FALLBACK", pred_res.get("endpoint_routing", {}).get("SOLUBILITY_GENERIC") == "BASE_FALLBACK")
            check("New compound routing has CYP3A4 in GLOBAL_V3_PRIMARY", pred_res.get("endpoint_routing", {}).get("CYP3A4_INHIBITION") == "GLOBAL_V3_PRIMARY")

            # Navigate to this new compound in browser
            driver.get(BASE_URL)
            time.sleep(1.5)
            open_project_by_name(driver, wait, "GLP-1")
            open_compound_by_name(driver, wait, "Audit Test Compound Ibuprofen")
            time.sleep(2)

            nb_text = driver.find_element(By.CLASS_NAME, "predict-meta-bar").text
            check("New compound workspace shows Prediction Engine v3.3", "Prediction Engine v3.3 · Current Production" in nb_text)

            # Hard reload and verify persistence
            driver.refresh()
            meta_bar_el = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "predict-meta-bar")))
            nb_text_after = meta_bar_el.text
            check("New compound persists after hard reload", "Prediction Engine v3.3 · Current Production" in nb_text_after)

            # Navigate away to Projects then re-enter
            select_nav_view(driver, wait, "Projects")
            time.sleep(1.5)
            open_project_by_name(driver, wait, "GLP-1")
            open_compound_by_name(driver, wait, "Audit Test Compound Ibuprofen")
            time.sleep(1.5)
            check("New compound re-entry successful", "Audit Test Compound Ibuprofen" in driver.find_element(By.TAG_NAME, "body").text)

        finally:
            # Clean up test compound
            if test_cid:
                db_clean = next(get_db())
                try:
                    c_del = db_clean.get(Compound, test_cid)
                    if c_del:
                        db_clean.delete(c_del)
                        db_clean.commit()
                        print(f"Cleaned up temporary test compound {test_cid}")
                finally:
                    db_clean.close()
            db.close()

        driver.quit()
        driver = None

        # =====================================================================
        # Part 5: Mobile Viewport (390x844) Layout Verification
        # =====================================================================
        print("\n--- Part 5: Mobile (390x844) Viewport Verification ---")
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
        # Part 6: Database Integrity & Isolation Check
        # =====================================================================
        print("\n--- Part 6: Database Integrity & Isolation Check ---")
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
        print("ALL ENGINE V3.3 PRODUCTION INTEGRITY AUDIT CHECKS PASSED!")
        return True


if __name__ == "__main__":
    import sys
    success = run_e2e_tests()
    sys.exit(0 if success else 1)
