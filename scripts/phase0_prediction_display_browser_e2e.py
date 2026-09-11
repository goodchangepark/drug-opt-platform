#!/usr/bin/env python3
"""Isolated browser regression for Predict -> current snapshot -> Compare."""

from __future__ import annotations

import json
import os
import subprocess
import time
import urllib.request
from pathlib import Path

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

ROOT = Path(__file__).resolve().parents[1]
BASE_URL = os.environ["DRUGOPT_E2E_BASE_URL"].rstrip("/")
RESULT = ROOT / "validation/phase0_prediction_display_e2e.json"
CHROMEDRIVER = "/snap/bin/chromium.chromedriver"


def api(path: str):
    with urllib.request.urlopen(BASE_URL + path, timeout=30) as response:
        return json.loads(response.read())


def click_button(driver, text: str):
    buttons = driver.find_elements(By.XPATH, f"//button[normalize-space()='{text}']")
    button = next((row for row in buttons if row.is_displayed() and row.is_enabled()), None)
    if button is None:
        raise AssertionError(f"visible enabled button {text!r} not found")
    driver.execute_script("arguments[0].click()", button)


def click_tab(driver, label: str):
    for button in driver.find_elements(By.CSS_SELECTOR, ".detail-tabs button"):
        parts = button.find_elements(By.TAG_NAME, "span")
        if parts and parts[0].text.strip() == label:
            driver.execute_script("arguments[0].click()", button)
            return
    raise AssertionError(f"tab {label!r} not found")


def main() -> None:
    service = subprocess.Popen(
        [CHROMEDRIVER, "--port=9517", "--allowed-ips=127.0.0.1"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        for _ in range(50):
            try:
                urllib.request.urlopen("http://127.0.0.1:9517/status", timeout=1)
                break
            except Exception:
                time.sleep(0.2)
        options = Options()
        for option in ("--headless", "--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu", "--window-size=1440,900"):
            options.add_argument(option)
        driver = webdriver.Remote("http://127.0.0.1:9517", options=options)
        wait = WebDriverWait(driver, 900)
        try:
            driver.get(BASE_URL)
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
            click_button(driver, "Projects")
            wait.until(lambda d: "EGFR" in d.page_source)
            click_button(driver, "EGFR")
            wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.compound-list")))
            click_button(driver, "Sunvozertinib")
            wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))

            before = api("/api/compound-versions/13/scientific-tabs/admet")
            before_rlm = next((row for row in before["rows"] if row["canonical_endpoint"] == "RLM_CLINT"), None)
            assert before_rlm is None or before_rlm["prediction"] is None

            click_button(driver, "▶ PREDICT")
            wait.until(lambda d: any(
                button.is_displayed() and button.is_enabled() and "PREDICT" in button.text and "PREDICTING" not in button.text
                for button in d.find_elements(By.CSS_SELECTOR, "button.btn-predict-primary")
            ))
            click_tab(driver, "ADMET")
            wait.until(lambda d: "RLM_CLINT" in d.page_source and "OpenADMET CheMeleon RLM" in d.page_source)
            after = api("/api/compound-versions/13/scientific-tabs/admet")
            after_rlm = next(row for row in after["rows"] if row["canonical_endpoint"] == "RLM_CLINT")
            snapshot_id = after_rlm["prediction"]["snapshot_id"]
            value = after_rlm["prediction"]["value"]

            driver.refresh()
            wait.until(lambda d: "RLM_CLINT" in d.page_source and "OpenADMET CheMeleon RLM" in d.page_source)
            reloaded = api("/api/compound-versions/13/scientific-tabs/admet")
            reloaded_rlm = next(row for row in reloaded["rows"] if row["canonical_endpoint"] == "RLM_CLINT")
            assert reloaded_rlm["prediction"]["snapshot_id"] == snapshot_id
            assert reloaded_rlm["prediction"]["value"] == value

            click_button(driver, "Back to Compounds")
            wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.compound-list")))
            rows = driver.find_elements(By.CSS_SELECTOR, "tr.compound-row")
            selected = 0
            for row in rows:
                if "Sunvozertinib" in row.text or "Zipalertinib" in row.text:
                    box = row.find_element(By.CSS_SELECTOR, "input.compound-select")
                    driver.execute_script("arguments[0].click()", box)
                    selected += 1
                    if selected == 2:
                        break
            assert selected == 2
            click_button(driver, "Compare Selected")
            wait.until(lambda d: "Selected Compound Comparison" in d.page_source and "RLM" in d.page_source)

            comparison = api("/api/projects/3/compare?ids=10,15")
            compared = next(row for row in comparison["compounds"] if row["row_id"] == 10)
            assert compared["prediction_snapshot_ids"]["RLM"] == snapshot_id
            assert compared["prediction_metadata"]["RLM"]["value"] == value
            for key in ("unit", "model_id", "model_version", "engine_version", "mode"):
                assert compared["prediction_metadata"]["RLM"][key] == after_rlm["prediction"][key]
            assert compared["prediction_metadata"]["RLM"]["species"] == after_rlm["species"]
            assert compared["prediction_metadata"]["RLM"]["context"] == after_rlm["context"]
            assert compared["prediction_metadata"]["RLM"]["maturity"] == after_rlm["maturity"]
            assert str(value) in driver.page_source
            assert compared["Solubility"] is None
            assert "Solubility" not in compared["prediction_snapshot_ids"]

            RESULT.write_text(json.dumps({
                "contract": "Phase0CanonicalPredictionDisplayE2E/v1",
                "environment": "E2E_ISOLATED_DATABASE",
                "compound": "Sunvozertinib",
                "compound_version_id": 13,
                "endpoint": "RLM_CLINT",
                "snapshot_id": snapshot_id,
                "value": value,
                "individual_display": "PASS",
                "reload_persistence": "PASS",
                "compare_snapshot_parity": "PASS",
                "legacy_only_prediction_hidden": "PASS",
            }, indent=2) + "\n", encoding="utf-8")
            print(RESULT.read_text(encoding="utf-8"))
        finally:
            driver.quit()
    finally:
        service.terminate()
        service.wait(timeout=10)


if __name__ == "__main__":
    main()
