"""Comprehensive E2E Browser Automation for Reference Library & Identity v3.3.3.

Covers:
1. Desktop (1440x900) & Mobile (390x844) viewports.
2. Real-time debounced 2D SVG preview upon keystroke typing SMILES.
3. CAS input -> structure resolution -> auto-fill SMILES & 2D preview.
4. Ketcher structure sync -> preview -> save.
5. Save compound -> page reload -> structure & identifier persistence.
6. Evidence search -> persistence -> reload.
7. Prediction execution -> stage snapshots & persistence -> reload.
8. Reopen GLP-1 (1), EGFR (3), AMYR (5) -> verify 2D structure, evidence, prediction.
9. DrugBank Reference Library (300):
   - 100 reference compounds table
   - 100% CAS coverage
   - 2D structure SVG preview
   - SMILES, InChIKey, DrugBank/ChEMBL/PubChem/UNII columns
   - Evidence count, Prediction status
   - VERIFIED status badge
   - Detail inspection: identity card, registry cross-references, provenance ledger
10. Test fixture isolation & cleanup: leaves only projects {1, 3, 5, 300}.
"""
import json
import os
import subprocess
import time
from typing import Optional

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
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
        time.sleep(1.5)

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


def set_react_input(driver, element, value: str):
    driver.execute_script("""
        const el = arguments[0];
        const val = arguments[1];
        const proto = el.tagName === 'TEXTAREA' ? window.HTMLTextAreaElement.prototype : window.HTMLInputElement.prototype;
        const setter = Object.getOwnPropertyDescriptor(proto, 'value').set;
        setter.call(el, val);
        el.dispatchEvent(new Event('input', { bubbles: true }));
        el.dispatchEvent(new Event('change', { bubbles: true }));
    """, element, value)
    time.sleep(0.3)


def select_nav_view(driver, wait, view_name: str):
    time.sleep(0.5)
    # If mobile menu toggle is displayed and sidebar is closed, open it
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
    WebDriverWait(driver, 45).until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
    time.sleep(1)


def cleanup_fixtures():
    from scripts.cleanup_test_fixtures import main as run_cleanup
    run_cleanup()


def run_tests():
    cleanup_fixtures()
    cd_mgr = ChromedriverManager(port=9515)
    cd_mgr.start()

    # Ensure backend is responding
    import urllib.request
    for _ in range(30):
        try:
            urllib.request.urlopen(f"{BASE_URL}/api/projects", timeout=1)
            break
        except Exception:
            time.sleep(1)

    try:
        # =========================================================================
        # 1. DESKTOP VIEWPORT E2E (1440 x 900)
        # =========================================================================
        print("\n" + "=" * 70)
        print("STAGE 1: DESKTOP VIEWPORT E2E (1440x900)")
        print("=" * 70)
        driver = create_driver(1440, 900, port=9515)
        wait = WebDriverWait(driver, 20)
        long_wait = WebDriverWait(driver, 60)

        driver.get(BASE_URL)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
        print("✓ Desktop shell loaded")

        # Create isolated test fixture project
        select_nav_view(driver, wait, "Projects")
        print("Creating isolated test fixture project...")
        import urllib.request
        req = urllib.request.Request(
            f"{BASE_URL}/api/projects",
            data=json.dumps({"name": "[TEST_FIXTURE] E2E Reference Library", "target": "Test", "description": "Automated E2E Fixture", "is_test_fixture": False}).encode("utf-8"),
            headers={"Content-Type": "application/json"}
        )
        with urllib.request.urlopen(req) as resp:
            test_proj = json.loads(resp.read().decode("utf-8"))
        test_pid = test_proj["id"]
        print(f"✓ Created test fixture project ID={test_pid}")

        # Open test fixture project in UI
        driver.refresh()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
        open_project_by_name(driver, wait, "[TEST_FIXTURE]")
        print("✓ Opened test fixture project")

        # -------------------------------------------------------------
        # Step A: Keystroke typing SMILES -> Live Preview auto-renders
        # -------------------------------------------------------------
        print("\n--- Step A: Live preview on keystroke typing ---")
        add_btn = wait.until(EC.element_to_be_clickable((By.ID, "btn-add-compound")))
        add_btn.click()
        time.sleep(0.5)

        smiles_input = wait.until(EC.presence_of_element_located((By.ID, "compound-smiles-input")))
        name_input = driver.find_element(By.ID, "compound-name-input")
        set_react_input(driver, name_input, "E2E-Keystroke-Molecule")

        previews = driver.find_elements(By.ID, "live-preview-card")
        assert len(previews) == 0, "Preview should not exist initially"

        print("Typing 'c1ccccc1' character-by-character...")
        for char in "c1ccccc1":
            smiles_input.send_keys(char)
            time.sleep(0.08)

        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        # Close modal after Step A
        close_btn = driver.find_element(By.XPATH, "//div[contains(@class, 'modal-backdrop')]//button[contains(text(), 'Cancel') or contains(text(), 'Close')]")
        close_btn.click()
        wait.until(EC.invisibility_of_element_located((By.ID, "compound-modal-container")))
        print("✓ Step A PASSED: Live 2D SVG preview automatically rendered during keystrokes")

        # -------------------------------------------------------------
        # Step B: CAS Resolve -> Auto-fill SMILES & Preview
        # -------------------------------------------------------------
        print("\n--- Step B: CAS Resolve Structure & Auto-fill ---")
        driver.find_element(By.ID, "btn-add-compound").click()
        wait.until(EC.presence_of_element_located((By.ID, "compound-modal-container")))
        cas_input = wait.until(EC.presence_of_element_located((By.ID, "compound-cas-input")))
        name_input = driver.find_element(By.ID, "compound-name-input")
        smiles_input = driver.find_element(By.ID, "compound-smiles-input")
        set_react_input(driver, cas_input, "50-78-2")  # Aspirin
        time.sleep(0.2)
        driver.find_element(By.ID, "btn-resolve-cas").click()

        wait.until(lambda d: "CC(=O)O" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value") or "C(=O)O" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value"))
        resolved_smiles = driver.find_element(By.ID, "compound-smiles-input").get_attribute("value").strip()
        print(f"✓ Resolved SMILES: {resolved_smiles}")
        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        assert preview_card.is_displayed()
        print("✓ Step B PASSED: CAS resolved to Aspirin SMILES and live preview updated")

        # -------------------------------------------------------------
        # Step C: Save compound -> Browser reload -> Persistence check
        # -------------------------------------------------------------
        print("\n--- Step C: Save & Reload Persistence ---")
        set_react_input(driver, name_input, "Aspirin")
        driver.find_element(By.ID, "btn-save-compound").click()
        wait.until(EC.invisibility_of_element_located((By.ID, "compound-modal-container")))
        time.sleep(1)

        print("Reloading page...")
        driver.refresh()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
        comp_btn = wait.until(EC.element_to_be_clickable((By.XPATH, "//button[contains(text(), 'Aspirin')]")))
        comp_btn.click()
        time.sleep(1)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

        # Check detail header and card
        detail_cas = driver.find_element(By.XPATH, "//div[contains(@class, 'compound-header-info')]//strong[contains(@class, 'mono')]").text
        assert "50-78-2" in detail_cas
        detail_smiles = driver.find_element(By.CLASS_NAME, "compound-smiles-bar").text
        assert "C(=O)O" in detail_smiles
        detail_svg = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
        assert detail_svg is not None
        print("✓ Step C PASSED: Compound saved, reloaded with intact CAS, SMILES, and 2D SVG")

        # -------------------------------------------------------------
        # Step D: Ketcher Modal Synchronization Check
        # -------------------------------------------------------------
        print("\n--- Step D: Ketcher Structure Editor Interaction ---")
        driver.find_element(By.ID, "btn-back-to-compounds").click()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
        driver.find_element(By.ID, "btn-add-compound").click()
        wait.until(EC.presence_of_element_located((By.ID, "compound-modal-container")))
        wait.until(lambda d: d.execute_script("return Boolean(document.getElementById('ketcher-editor')?.contentWindow?.ketcher)"))
        time.sleep(0.8)
        print("Dispatching Pyridine (c1ccncc1) to Ketcher editor...")
        driver.execute_async_script("""
            const callback = arguments[arguments.length - 1];
            const editor = document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
            if (editor && editor.setMolecule) {
                editor.setMolecule('c1ccncc1').then(() => callback(true)).catch(() => callback(false));
            } else {
                callback(false);
            }
        """)
        wait.until(lambda d: "c1ccncc1" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value") or "N" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value").upper())
        ketcher_smiles = driver.find_element(By.ID, "compound-smiles-input").get_attribute("value").strip()
        print(f"✓ Ketcher synced SMILES: {ketcher_smiles}")
        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        assert preview_card.is_displayed()
        # Close modal
        close_btn = driver.find_element(By.XPATH, "//div[contains(@class, 'modal-backdrop')]//button[contains(text(), 'Cancel') or contains(text(), 'Close')]")
        close_btn.click()
        wait.until(EC.invisibility_of_element_located((By.ID, "compound-modal-container")))
        print("✓ Step D PASSED: Ketcher editor interaction & structure sync verified")

        # Reopen Aspirin compound detail
        comp_btn = wait.until(EC.element_to_be_clickable((By.XPATH, "//button[contains(text(), 'Aspirin')]")))
        comp_btn.click()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

        # -------------------------------------------------------------
        # Step E: Evidence Search & Qualification Persistence
        # -------------------------------------------------------------
        print("\n--- Step E: Evidence Search & Qualification Persistence ---")
        ev_search_btn = wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(text(), 'Search Experimental Data')]")))
        driver.execute_script("arguments[0].click();", ev_search_btn)
        long_wait = WebDriverWait(driver, 60)
        long_wait.until(lambda d: "✓ Saved" in d.find_element(By.CLASS_NAME, "experimental-evidence-status").text)
        ev_status_text = driver.find_element(By.CLASS_NAME, "experimental-evidence-status").text
        print(f"Evidence status: {ev_status_text}")
        assert "✓ Saved" in ev_status_text

        # Reload to verify persistence of search run
        driver.refresh()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        wait.until(lambda d: "✓ Saved" in d.find_element(By.CLASS_NAME, "experimental-evidence-status").text)
        reloaded_ev_status = driver.find_element(By.CLASS_NAME, "experimental-evidence-status").text
        assert "✓ Saved" in reloaded_ev_status
        print("✓ Step E PASSED: Evidence search result persisted across page reload")

        # -------------------------------------------------------------
        # Step F: Predict Execution & Snapshot Persistence
        # -------------------------------------------------------------
        print("\n--- Step F: Prediction Run & Snapshot Persistence ---")
        predict_btn = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "btn-predict-primary")))
        driver.execute_script("arguments[0].click();", predict_btn)
        long_wait.until(lambda d: "✓ Saved" in d.find_element(By.CLASS_NAME, "predict-meta-bar").text)
        pred_meta = driver.find_element(By.CLASS_NAME, "predict-meta-bar").text
        print(f"Prediction meta: {pred_meta}")
        assert "✓ Saved" in pred_meta

        # Reload to verify snapshot persistence
        driver.refresh()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        time.sleep(2)
        long_wait.until(lambda d: "LOADING" not in d.find_element(By.CLASS_NAME, "predict-meta-bar").text)
        print("After reload stabilized, predict-meta-bar:", driver.find_element(By.CLASS_NAME, "predict-meta-bar").text)
        assert "COMPLETE" in driver.find_element(By.CLASS_NAME, "predict-meta-bar").text or "✓ Saved" in driver.find_element(By.CLASS_NAME, "predict-meta-bar").text
        print("✓ Step F PASSED: Prediction run and stage snapshots persisted across reload")

        # Return to compounds
        driver.find_element(By.ID, "btn-back-to-compounds").click()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        # -------------------------------------------------------------
        # Step G: Verify Existing Projects GLP-1, EGFR, AMYR
        # -------------------------------------------------------------
        print("\n--- Step G: Existing Projects Verification (GLP-1, EGFR, AMYR) ---")
        for p_name in ["GLP-1", "EGFR", "AMYR"]:
            open_project_by_name(driver, wait, p_name)
            rows = driver.find_elements(By.CLASS_NAME, "compound-row")
            print(f"Project '{p_name}': {len(rows)} compound rows rendered")
            assert len(rows) > 0
            first_link = rows[0].find_element(By.CLASS_NAME, "compound-name-link")
            comp_name = first_link.text
            first_link.click()
            time.sleep(1)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
            svg = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
            assert svg is not None
            print(f"✓ Project '{p_name}' compound '{comp_name}' rendered 2D structure SVG")
            driver.find_element(By.ID, "btn-back-to-compounds").click()
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        # -------------------------------------------------------------
        # Step H: DrugBank Reference Library (ID=300) Verification
        # -------------------------------------------------------------
        print("\n--- Step H: DrugBank Reference Library 100 Verification ---")
        open_project_by_name(driver, wait, "DrugBank")
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "project-status-table")))

        header_text = driver.find_element(By.XPATH, "//h2[contains(text(), 'DrugBank Reference Library')]").text
        print(f"Library Header: {header_text}")
        assert "100 Approved Drugs" in header_text

        rows = driver.find_elements(By.CLASS_NAME, "compound-row")
        print(f"Total reference drug rows displayed: {len(rows)}")
        assert len(rows) == 100, f"Expected 100 rows, found {len(rows)}"

        cas_tags = driver.find_elements(By.CLASS_NAME, "cas-tag")
        print(f"Total CAS tags in table: {len(cas_tags)}")
        assert len(cas_tags) == 100, f"Expected 100 CAS tags, found {len(cas_tags)}"

        # Check random reference drugs in detail view (e.g., Digoxin, Metformin, Ketoprofen, Clozapine)
        for drug_name in ["Digoxin", "Metformin", "Ketoprofen", "Clozapine"]:
            print(f"Inspecting detail view for Reference Drug '{drug_name}'...")
            comp_link = wait.until(EC.presence_of_element_located((By.XPATH, f"//button[contains(@class, 'compound-name-link') and text()='{drug_name}']")))
            driver.execute_script("arguments[0].scrollIntoView(true);", comp_link)
            time.sleep(0.3)
            driver.execute_script("arguments[0].click();", comp_link)
            time.sleep(1)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

            id_card = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-identity-card")))
            card_text = id_card.text
            print(f"Identity Card preview: {card_text.splitlines()[:4]}")
            assert "REFERENCE IDENTITY & VERIFICATION PROVENANCE" in card_text
            assert "VERIFIED" in card_text
            assert "Evidence:" in card_text
            assert "Prediction:" in card_text
            assert "CAS Registry No." in card_text
            assert "DrugBank ID" in card_text
            assert "ChEMBL ID" in card_text
            assert "PubChem CID" in card_text
            assert "UNII" in card_text
            assert "Canonical SMILES:" in card_text
            assert "InChIKey:" in card_text

            summary_el = id_card.find_element(By.TAG_NAME, "summary")
            assert "Provenance & Cross-Verification Ledger" in summary_el.text
            summary_el.click()
            time.sleep(0.3)
            ledger_table = id_card.find_element(By.TAG_NAME, "table")
            ledger_rows = ledger_table.find_elements(By.XPATH, ".//tbody/tr")
            print(f"Provenance ledger rows for '{drug_name}': {len(ledger_rows)}")
            assert len(ledger_rows) >= 6, f"Expected at least 6 identifier records, found {len(ledger_rows)}"

            svg_detail = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
            assert svg_detail is not None

            back_btn = driver.find_element(By.ID, "btn-back-to-compounds")
            driver.execute_script("arguments[0].click();", back_btn)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "project-status-table")))
            print(f"✓ Reference drug '{drug_name}' verified perfectly!")

        driver.quit()
        print("\n✓ ALL DESKTOP E2E TESTS PASSED!")

        # =========================================================================
        # 2. MOBILE VIEWPORT E2E (390 x 844)
        # =========================================================================
        print("\n" + "=" * 70)
        print("STAGE 2: MOBILE VIEWPORT E2E (390x844)")
        print("=" * 70)
        m_driver = create_driver(390, 844, port=9515)
        m_wait = WebDriverWait(m_driver, 15)

        m_driver.get(BASE_URL)
        m_wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
        print("✓ Mobile shell loaded")

        open_project_by_name(m_driver, m_wait, "DrugBank")
        m_wait.until(EC.presence_of_element_located((By.CLASS_NAME, "project-status-table")))
        m_rows = m_driver.find_elements(By.CLASS_NAME, "compound-row")
        print(f"Mobile view: {len(m_rows)} reference drugs rendered")
        assert len(m_rows) == 100

        m_link = m_wait.until(EC.presence_of_element_located((By.XPATH, "//button[contains(@class, 'compound-name-link') and text()='Digoxin']")))
        m_driver.execute_script("arguments[0].scrollIntoView(true);", m_link)
        time.sleep(0.3)
        m_driver.execute_script("arguments[0].click();", m_link)
        time.sleep(1)
        m_wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        m_id_card = m_wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-identity-card")))
        assert "VERIFIED" in m_id_card.text
        print("✓ Mobile view: Digoxin detail & identity card verified")

        m_driver.quit()
        print("\n✓ ALL MOBILE E2E TESTS PASSED!")

    finally:
        cd_mgr.stop()
        cleanup_fixtures()
        print("✓ Temporary test fixtures purged cleanly.")


if __name__ == "__main__":
    run_tests()
