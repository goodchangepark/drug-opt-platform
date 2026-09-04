"""E2E Browser Automation for Live Structure Preview and DrugBank Identity v3.3.2.

Executes real browser automated tests using headless Chromium & Selenium:
- Test A: Keystroke typing SMILES triggers debounced live 2D SVG preview before save
- Test B: Dynamic modification of SMILES updates 2D SVG preview in real-time
- Test C: CAS entry and resolution auto-fills SMILES and triggers live preview
- Test D: Save compound -> browser reload -> CAS, SMILES, and 2D SVG persist
- Test E: Ketcher drawing sync -> SMILES and preview auto-update -> save -> reload
- Test F: DrugBank 80 compounds -> 100% CAS coverage -> 10 random compounds detail inspection
"""
from __future__ import annotations

import json
import os
import random
import subprocess
import time
from pathlib import Path
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

def create_driver(port: int = 9515) -> webdriver.Remote:
    opts = Options()
    opts.add_argument("--headless")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument("--disable-gpu")
    opts.add_argument("--window-size=1600,1000")
    return webdriver.Remote(f"http://127.0.0.1:{port}", options=opts)


def select_nav_view(driver, wait, view_name: str):
    time.sleep(0.5)
    buttons = driver.find_elements(By.XPATH, "//nav//button | //aside//button")
    for b in buttons:
        if b.text.strip() == view_name:
            b.click()
            time.sleep(0.8)
            return
    raise ValueError(f"Nav item '{view_name}' not found")


def open_project_by_name(driver, wait, project_keyword: str):
    select_nav_view(driver, wait, "Projects")
    btn = wait.until(EC.element_to_be_clickable((By.XPATH, f"//button[contains(text(), '{project_keyword}')]")))
    btn.click()
    time.sleep(1.2)
    wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))


def cleanup_test_compounds():
    import sqlite3
    conn = sqlite3.connect("drug_opt.db")
    c = conn.cursor()
    c.execute("SELECT id FROM compounds WHERE name LIKE 'E2E-%'")
    ids = [r[0] for r in c.fetchall()]
    if ids:
        id_list = ",".join(map(str, ids))
        c.execute(f"DELETE FROM property_calculations WHERE version_id IN (SELECT id FROM compound_versions WHERE compound_row_id IN ({id_list}))")
        c.execute(f"DELETE FROM structural_alerts WHERE version_id IN (SELECT id FROM compound_versions WHERE compound_row_id IN ({id_list}))")
        c.execute(f"DELETE FROM compound_versions WHERE compound_row_id IN ({id_list})")
        c.execute(f"DELETE FROM compounds WHERE id IN ({id_list})")
        conn.commit()
        print(f"Cleaned up {len(ids)} temporary test compounds.")
    conn.close()


def run_e2e_tests():
    cleanup_test_compounds()
    cd_mgr = ChromedriverManager(port=9515)
    cd_mgr.start()

    try:
        driver = create_driver(port=9515)
        wait = WebDriverWait(driver, 15)

        print("\n" + "=" * 70)
        print("STARTING LIVE STRUCTURE PREVIEW & DRUGBANK E2E TESTS (v3.3.2)")
        print("=" * 70)

        # -------------------------------------------------------------
        # Step 0: Initial Navigation & Open Project 1 (GLP-1)
        # -------------------------------------------------------------
        print("\n[Step 0] Navigating to platform homepage...")
        driver.get(BASE_URL)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
        print("✓ Platform shell loaded successfully")

        print("Opening GLP-1 project workspace...")
        open_project_by_name(driver, wait, "GLP-1")
        print("✓ Connected to Project Workspace (GLP-1)")

        # -------------------------------------------------------------
        # Test A: Keystroke typing SMILES -> Live Preview auto-appears
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST A: Keystroke typing SMILES -> live 2D SVG preview before save")
        print("-" * 60)

        add_btn = wait.until(EC.element_to_be_clickable((By.ID, "btn-add-compound")))
        add_btn.click()
        time.sleep(0.5)

        smiles_input = wait.until(EC.presence_of_element_located((By.ID, "compound-smiles-input")))
        name_input = driver.find_element(By.ID, "compound-name-input")
        name_input.send_keys(Keys.CONTROL + "a")
        name_input.send_keys(Keys.BACKSPACE)
        name_input.send_keys("E2E-LivePreview-A")

        # Verify preview card is not present yet
        previews = driver.find_elements(By.ID, "live-preview-card")
        assert len(previews) == 0, "Preview card should not exist before typing SMILES"

        # Type SMILES character by character simulating real human typing
        test_smiles = "c1ccccc1"  # Benzene
        print(f"Typing '{test_smiles}' character-by-character into SMILES input...")
        for char in test_smiles:
            smiles_input.send_keys(char)
            time.sleep(0.08)  # Realistic typing delay (80ms per keystroke)

        # Wait for debounced preview card to appear automatically WITHOUT clicking save or validate
        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        assert preview_card.is_displayed(), "Live preview card must be displayed automatically"

        svg_el = preview_card.find_element(By.TAG_NAME, "svg")
        assert svg_el is not None, "Live preview card must contain rendered SVG"
        svg_html = svg_el.get_attribute("outerHTML")
        assert len(svg_html) > 100, f"SVG content too small: {len(svg_html)} bytes"

        print(f"✓ Test A PASSED: Live preview card auto-rendered 2D SVG ({len(svg_html)} bytes) during typing!")

        # -------------------------------------------------------------
        # Test B: Change SMILES -> Preview updates immediately
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST B: Modify SMILES -> Live preview updates dynamically")
        print("-" * 60)

        # Clear SMILES input cleanly
        set_react_input(driver, smiles_input, "")
        time.sleep(0.4)

        # Type new SMILES: Naphthalene c1ccc2ccccc2c1
        new_smiles = "c1ccc2ccccc2c1"
        print(f"Typing modified structure '{new_smiles}'...")
        for char in new_smiles:
            smiles_input.send_keys(char)
            time.sleep(0.05)

        # Wait for preview card to finish debouncing and reflect updated molecule
        def naphthalene_ready(d):
            cards = d.find_elements(By.ID, "live-preview-card")
            if not cards:
                return False
            card_text = cards[0].text
            return cards[0] if ("C10H8" in card_text or "128." in card_text) else False

        preview_card = wait.until(naphthalene_ready)
        text_content = preview_card.text
        new_svg = preview_card.find_element(By.TAG_NAME, "svg").get_attribute("outerHTML")
        assert new_svg != svg_html, "SVG depiction must change when SMILES changes"

        print("✓ Test B PASSED: Preview updated dynamically with new properties and structure!")

        # -------------------------------------------------------------
        # Test C: CAS entry -> Resolve -> Auto-fills SMILES & preview
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST C: CAS entry -> Resolve -> Auto-fills SMILES & Live Preview")
        print("-" * 60)

        # Clear SMILES completely
        set_react_input(driver, smiles_input, "")
        time.sleep(0.3)

        # Enter CAS for Ibuprofen (15687-27-1)
        cas_input = driver.find_element(By.ID, "compound-cas-input")
        set_react_input(driver, cas_input, "15687-27-1")
        time.sleep(0.2)

        # Click Resolve Structure from CAS
        resolve_btn = driver.find_element(By.ID, "btn-resolve-cas")
        resolve_btn.click()

        # Wait for SMILES input to be auto-filled with Ibuprofen
        wait.until(lambda d: "CC(C)" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value") or "C(C)C(=O)O" in d.find_element(By.ID, "compound-smiles-input").get_attribute("value"))
        resolved_smiles = driver.find_element(By.ID, "compound-smiles-input").get_attribute("value").strip()
        print(f"CAS 15687-27-1 resolved to SMILES: '{resolved_smiles}'")
        assert "CC(C)Cc1ccc(C(C)C(=O)O)cc1" in resolved_smiles or "C(C)C(=O)O" in resolved_smiles

        # Wait for live preview card to auto-render for the resolved molecule
        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        assert preview_card.is_displayed()
        cas_svg = preview_card.find_element(By.TAG_NAME, "svg").get_attribute("outerHTML")
        assert len(cas_svg) > 100

        print("✓ Test C PASSED: CAS resolved structure, auto-filled SMILES, and rendered 2D preview!")

        # -------------------------------------------------------------
        # Test D: Save compound -> Browser reload -> Persistence
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST D: Save compound -> Browser reload -> Check persistence")
        print("-" * 60)

        test_compound_name = f"E2E-Persist-Ibuprofen-{int(time.time())}"
        set_react_input(driver, name_input, test_compound_name)

        save_btn = driver.find_element(By.ID, "btn-save-compound")
        save_btn.click()

        # Wait for modal to disappear
        wait.until(EC.invisibility_of_element_located((By.ID, "compound-modal-container")))
        time.sleep(1)

        # Reload browser
        print("Reloading browser to verify full database persistence...")
        driver.refresh()
        time.sleep(1.5)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        # Locate saved compound in compound list table
        table = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
        compound_link = wait.until(EC.element_to_be_clickable((By.XPATH, f"//button[contains(text(), '{test_compound_name}')]")))
        assert compound_link is not None, f"Compound {test_compound_name} not found in table after reload"

        # Check CAS in table row
        parent_row = compound_link.find_element(By.XPATH, "./ancestor::tr")
        row_text = parent_row.text
        assert "15687-27-1" in row_text, f"CAS 15687-27-1 not displayed in table row: {row_text}"
        print(f"✓ Table row displays: {row_text.splitlines()[0]}")

        # Open compound detail view
        compound_link.click()
        time.sleep(1)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

        # Check Name, CAS, SMILES, and 2D SVG in detail view
        detail_h2 = driver.find_element(By.XPATH, "//div[contains(@class, 'compound-header-info')]//h2").text
        assert test_compound_name in detail_h2, f"Expected {test_compound_name}, got {detail_h2}"

        detail_cas = driver.find_element(By.XPATH, "//div[contains(@class, 'compound-header-info')]//strong[contains(@class, 'mono')]").text
        assert "15687-27-1" in detail_cas, f"Expected CAS 15687-27-1 in detail, got {detail_cas}"

        detail_smiles = driver.find_element(By.CLASS_NAME, "compound-smiles-bar").text
        assert "C(C)C(=O)O" in detail_smiles, f"Expected SMILES in detail bar, got {detail_smiles}"

        detail_svg = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
        assert detail_svg is not None and len(detail_svg.get_attribute("outerHTML")) > 100

        print(f"✓ Test D PASSED: Saved compound reloaded with CAS, SMILES, and 2D SVG structure fully intact!")

        # Return to compound list
        driver.find_element(By.ID, "btn-back-to-compounds").click()
        time.sleep(0.8)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        # -------------------------------------------------------------
        # Test E: Ketcher draw -> SMILES/preview sync -> Save -> Reload
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST E: Ketcher draw -> SMILES & preview sync -> Save -> Reload")
        print("-" * 60)

        # Open Add Compound modal again
        add_btn = wait.until(EC.element_to_be_clickable((By.ID, "btn-add-compound")))
        add_btn.click()
        time.sleep(0.5)

        ketcher_test_name = f"E2E-Ketcher-Pyridine-{int(time.time())}"
        name_input = driver.find_element(By.ID, "compound-name-input")
        set_react_input(driver, name_input, ketcher_test_name)

        # Wait for Ketcher iframe and editor to be ready
        print("Waiting for Ketcher structure editor...")
        wait.until(lambda d: d.execute_script("return Boolean(document.getElementById('ketcher-editor')?.contentWindow?.ketcher)"))
        time.sleep(1)

        print("Dispatching chemical structure (Pyridine c1ccncc1) to Ketcher editor...")
        driver.execute_async_script("""
            const callback = arguments[arguments.length - 1];
            const editor = document.getElementById('ketcher-editor')?.contentWindow?.ketcher;
            if (editor && editor.setMolecule) {
                editor.setMolecule('c1ccncc1').then(() => callback(true)).catch(() => callback(false));
            } else {
                callback(false);
            }
        """)

        # Wait for polling to synchronize SMILES input and preview
        wait.until(lambda d: len(d.find_element(By.ID, "compound-smiles-input").get_attribute("value").strip()) > 0)
        smiles_val = driver.find_element(By.ID, "compound-smiles-input").get_attribute("value").strip()
        print(f"Ketcher synchronized SMILES: '{smiles_val}'")
        assert "c1ccncc1" in smiles_val or "N" in smiles_val.upper(), f"Unexpected synced smiles: {smiles_val}"

        # Verify live preview rendered Pyridine
        preview_card = wait.until(EC.presence_of_element_located((By.ID, "live-preview-card")))
        assert preview_card.is_displayed()
        ketcher_svg = preview_card.find_element(By.TAG_NAME, "svg").get_attribute("outerHTML")
        assert len(ketcher_svg) > 100

        # Save Ketcher compound
        driver.find_element(By.ID, "btn-save-compound").click()
        wait.until(EC.invisibility_of_element_located((By.ID, "compound-modal-container")))
        time.sleep(1)

        # Reload and check Ketcher compound in detail view
        driver.refresh()
        time.sleep(1.5)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        ketcher_link = wait.until(EC.element_to_be_clickable((By.XPATH, f"//button[contains(text(), '{ketcher_test_name}')]")))
        ketcher_link.click()
        time.sleep(1)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

        k_smiles = driver.find_element(By.CLASS_NAME, "compound-smiles-bar").text
        assert "c1ccncc1" in k_smiles or "N" in k_smiles
        k_svg = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
        assert k_svg is not None and len(k_svg.get_attribute("outerHTML")) > 100

        print(f"✓ Test E PASSED: Ketcher drawing synchronized to SMILES & preview, saved, and persisted across reload!")

        # Return to compound list
        driver.find_element(By.ID, "btn-back-to-compounds").click()
        time.sleep(0.8)

        # Clean up temporary test compounds
        cleanup_test_compounds()

        # -------------------------------------------------------------
        # Test F: DrugBank 80 compounds -> CAS backfill -> 10 random checks
        # -------------------------------------------------------------
        print("\n" + "-" * 60)
        print("TEST F: DrugBank 80 compounds -> 100% CAS coverage & 10 random detail checks")
        print("-" * 60)

        # Open DrugBank project
        open_project_by_name(driver, wait, "DrugBank")

        # Check total compound rows in table
        table = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
        rows = table.find_elements(By.TAG_NAME, "tr")[1:]  # skip header
        print(f"DrugBank compounds rendered in table: {len(rows)}")
        assert len(rows) == 80, f"Expected exactly 80 compounds rendered in table, got {len(rows)}"

        # Verify all 80 rows display a valid CAS number
        missing_cas_rows = []
        for idx, r in enumerate(rows):
            text = r.text
            if "CAS:" not in text:
                missing_cas_rows.append((idx, text.splitlines()[0] if text else "empty"))
        assert len(missing_cas_rows) == 0, f"Found {len(missing_cas_rows)} compounds missing CAS display in table: {missing_cas_rows[:5]}"
        print("✓ All 80 compounds in table display verified CAS numbers!")

        # Sample 10 random compounds to inspect in detail view
        random.seed(42)  # Deterministic seed for reproducible verification
        sample_indices = sorted(random.sample(range(len(rows)), 10))
        print(f"\nRandomly sampling 10 compounds for in-depth Detail View inspection: indices {sample_indices}")

        for i, sample_idx in enumerate(sample_indices, 1):
            # Re-fetch rows in case DOM refreshed
            table = wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))
            cur_rows = table.find_elements(By.TAG_NAME, "tr")[1:]
            target_row = cur_rows[sample_idx]
            link = target_row.find_element(By.CLASS_NAME, "link-button")
            comp_name = link.text.strip()
            
            # Click to open detail view
            link.click()
            time.sleep(0.8)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

            # Verify detail view attributes
            h2_text = driver.find_element(By.XPATH, "//div[contains(@class, 'compound-header-info')]//h2").text
            assert comp_name == h2_text, f"Name mismatch: expected {comp_name}, got {h2_text}"

            cas_text = driver.find_element(By.XPATH, "//div[contains(@class, 'compound-header-info')]//strong[contains(@class, 'mono')]").text
            assert cas_text != "Not provided" and len(cas_text) > 4, f"CAS missing in detail for {comp_name}: {cas_text}"

            smiles_text = driver.find_element(By.CLASS_NAME, "compound-smiles-bar").text
            assert len(smiles_text) > 2, f"SMILES missing for {comp_name}"

            svg_elem = driver.find_element(By.CSS_SELECTOR, ".compound-header-structure svg, .compound-header-structure img")
            svg_content = svg_elem.get_attribute("outerHTML")
            assert len(svg_content) > 100, f"SVG structure depiction missing or invalid for {comp_name}"

            print(f"  [{i}/10] {comp_name}: CAS={cas_text} | SMILES={smiles_text[:25]}... | 2D SVG={len(svg_content)} bytes ✓")

            # Click Back to Compounds button
            driver.find_element(By.ID, "btn-back-to-compounds").click()
            time.sleep(0.6)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-list")))

        print("\n✓ Test F PASSED: 10/10 randomly sampled DrugBank compounds fully verified with Name, CAS, SMILES, and 2D SVG depiction!")

        print("\n" + "=" * 70)
        print("ALL TESTS (A, B, C, D, E, F) COMPLETED AND PASSED WITH 100% SUCCESS!")
        print("=" * 70)

        driver.quit()

    finally:
        cleanup_test_compounds()
        cd_mgr.stop()


if __name__ == "__main__":
    run_e2e_tests()
