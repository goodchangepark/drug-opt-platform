#!/usr/bin/env python3
"""Read-only desktop/mobile E2E acceptance for Stable Core v1.

Run through ``scripts/run_isolated_e2e.py`` so every browser/API request uses
the private E2E database copy.  The script never invokes prediction, evidence
search/import, edit, or delete actions.
"""

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
CHROMEDRIVER = "/snap/bin/chromium.chromedriver"
RESULT = ROOT / "validation/stable_core_v1_browser_e2e.json"


class DriverService:
    def __init__(self, port: int = 9516):
        self.port = port
        self.process: subprocess.Popen | None = None

    def __enter__(self):
        self.process = subprocess.Popen(
            [CHROMEDRIVER, f"--port={self.port}", "--allowed-ips=127.0.0.1"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for _ in range(30):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{self.port}/status", timeout=1)
                return self
            except Exception:
                time.sleep(0.2)
        raise RuntimeError("ChromeDriver did not become ready")

    def __exit__(self, *_):
        if self.process:
            self.process.terminate()
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()


def api(path: str):
    started = time.perf_counter()
    with urllib.request.urlopen(BASE_URL + path, timeout=20) as response:
        body = response.read()
    return json.loads(body), round((time.perf_counter() - started) * 1000, 3), len(body)


def new_driver(port: int, width: int, height: int):
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-gpu")
    options.add_argument(f"--window-size={width},{height}")
    return webdriver.Remote(f"http://127.0.0.1:{port}", options=options)


def click_text(driver, text: str):
    exact = driver.find_elements(By.XPATH, f"//button[normalize-space()='{text}']")
    candidates = exact or driver.find_elements(By.XPATH, f"//button[contains(normalize-space(),'{text}')]")
    visible = next((item for item in candidates if item.is_displayed() and item.is_enabled()), None)
    if not visible:
        raise AssertionError(f"Visible button containing {text!r} was not found")
    driver.execute_script("arguments[0].click()", visible)


def request_urls(driver):
    return driver.execute_script(
        "const urls=performance.getEntriesByType('resource').map(row=>row.name);"
        "performance.clearResourceTimings();return urls;"
    )


def wait_for_resource(driver, wait, fragment: str):
    wait.until(lambda d: any(
        fragment in row["name"]
        for row in d.execute_script("return performance.getEntriesByType('resource')")
    ))


def click_tab(driver, label: str):
    for button in driver.find_elements(By.CSS_SELECTOR, ".detail-tabs button"):
        parts = button.find_elements(By.TAG_NAME, "span")
        if parts and parts[0].text.strip() == label:
            driver.execute_script("arguments[0].click()", button)
            return
    raise AssertionError(f"Scientific tab {label!r} was not found")


def open_reference_project(driver, wait):
    project_buttons = driver.find_elements(By.XPATH, "//button[normalize-space()='Projects']")
    if not any(button.is_displayed() for button in project_buttons):
        menu = next((button for button in driver.find_elements(By.CSS_SELECTOR, "button.menu-toggle") if button.is_displayed()), None)
        if menu:
            driver.execute_script("arguments[0].click()", menu)
    click_text(driver, "Projects")
    wait.until(lambda d: "REFERENCE LIBRARY" in d.page_source and "DrugBank" in d.page_source)
    click_text(driver, "DrugBank")
    wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.compound-list")))


def search_reference(driver, wait, name: str):
    field = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "input.project-reference-search")))
    field.clear()
    field.send_keys(name)
    wait.until(lambda d: d.find_element(By.CSS_SELECTOR, "input.project-reference-search").get_attribute("value") == name)
    click_text(driver, "Search")
    wait.until(lambda d: len(d.find_elements(By.CSS_SELECTOR, "tr.compound-row")) == 1)


def exercise_viewport(port: int, width: int, height: int):
    driver = new_driver(port, width, height)
    wait = WebDriverWait(driver, 40)
    started = time.perf_counter()
    try:
        driver.get(BASE_URL)
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "shell")))
        open_reference_project(driver, wait)
        initial_ms = round((time.perf_counter() - started) * 1000, 3)
        initial_rows = driver.find_elements(By.CSS_SELECTOR, "tr.compound-row")
        assert len(initial_rows) == 50, f"expected bounded 50-row page, got {len(initial_rows)}"
        first_urls = request_urls(driver)
        assert any("/api/projects/300?page=1&page_size=50" in url for url in first_urls)
        assert not any("/workspace" in url for url in first_urls)

        click_text(driver, "Next")
        wait.until(lambda d: "Page 2 of 20" in d.page_source)
        page_urls = request_urls(driver)
        assert any("/api/projects/300?page=2&page_size=50" in url for url in page_urls)
        assert len(driver.find_elements(By.CSS_SELECTOR, "tr.compound-row")) == 50
        click_text(driver, "Previous")
        wait.until(lambda d: "Page 1 of 20" in d.page_source)
        request_urls(driver)

        search_reference(driver, wait, "Warfarin")
        assert "Warfarin" in driver.page_source
        click_text(driver, "Warfarin")
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        click_text(driver, "Back to Compounds")
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "input.project-reference-search")))
        assert driver.find_element(By.CSS_SELECTOR, "input.project-reference-search").get_attribute("value") == "Warfarin"

        project_buttons = driver.find_elements(By.XPATH, "//button[normalize-space()='Projects']")
        if not any(button.is_displayed() for button in project_buttons):
            menu = next((button for button in driver.find_elements(By.CSS_SELECTOR, "button.menu-toggle") if button.is_displayed()), None)
            if menu:
                driver.execute_script("arguments[0].click()", menu)
        click_text(driver, "Projects")
        wait.until(lambda d: "GLP-1 (small molecule)" in d.page_source)
        click_text(driver, "GLP-1 (small molecule)")
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.compound-list")))
        click_text(driver, "Orforglipron")
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        assert "drugopt-prediction-engine-v3@3.3.3" in driver.page_source
        assert "4 accepted observations" in driver.page_source
        detail_urls = request_urls(driver)
        assert any("/api/compounds/1/summary" in url for url in detail_urls)
        assert not any("/workspace" in url for url in detail_urls)
        assert not any("experimental-harvest" in url for url in detail_urls)

        tab_requests = {}
        for label, route in (
            ("PROPERTIES", "/scientific-tabs/properties"),
            ("ACTIVITY", "/scientific-tabs/activity"),
            ("ADMET", "/scientific-tabs/admet"),
            ("METABOLISM", "/scientific-tabs/metabolism"),
        ):
            click_tab(driver, label)
            wait_for_resource(driver, wait, route)
            urls = request_urls(driver)
            matching = [url for url in urls if route in url]
            assert len(matching) == 1, matching
            assert not any("/workspace" in url for url in urls)
            tab_requests[label] = len(matching)

        click_tab(driver, "PK")
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.stable-core-scientific-table")))
        pk_text = driver.find_element(By.CSS_SELECTOR, "table.stable-core-scientific-table").text
        for expected in ("HUMAN_PK_F_ORAL", "HUMAN_PK_VD_IV", "HUMAN_PK_CL_UNSPECIFIED", "HUMAN_PK_CMAX_UNSPECIFIED"):
            assert expected in pk_text
        assert "77.0 %" in pk_text and "285.0 L" in pk_text and "7.15 L/h" in pk_text and "149.0 ng/mL" in pk_text
        assert "dose: 36" in pk_text and "route: UNSPECIFIED" in pk_text
        pk_urls = request_urls(driver)
        canonical_pk = [url for url in pk_urls if "/scientific-tabs/pk" in url]
        assert len(canonical_pk) == 1, canonical_pk
        assert not any("/pk-studies" in url or "/ivive" in url or "/pk-simulation/run" in url for url in pk_urls)

        driver.refresh()
        wait.until(EC.presence_of_element_located((By.CLASS_NAME, "compound-workspace")))
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "table.stable-core-scientific-table")))
        reloaded_pk = driver.find_element(By.CSS_SELECTOR, "table.stable-core-scientific-table").text
        assert "77.0 %" in reloaded_pk and "149.0 ng/mL" in reloaded_pk
        reload_urls = request_urls(driver)
        assert not any("/workspace" in url or "experimental-harvest" in url for url in reload_urls)

        click_tab(driver, "EVIDENCE")
        wait_for_resource(driver, wait, "/scientific-tabs/evidence")
        wait.until(lambda d: "CANONICAL EVIDENCE STORE" in d.page_source and "FDA" in d.page_source)
        evidence_urls = request_urls(driver)
        canonical_evidence = [url for url in evidence_urls if "/scientific-tabs/evidence" in url]
        assert len(canonical_evidence) == 1, canonical_evidence

        click_tab(driver, "HISTORY")
        wait_for_resource(driver, wait, "/scientific-tabs/history")
        wait.until(lambda d: "IMMUTABLE SCIENTIFIC HISTORY" in d.page_source)
        history_urls = request_urls(driver)
        canonical_history = [url for url in history_urls if "/scientific-tabs/history" in url]
        assert len(canonical_history) == 1, canonical_history

        return {
            "viewport": f"{width}x{height}",
            "initial_project_ui_ms": initial_ms,
            "initial_rows": len(initial_rows),
            "pagination_page2_rows": 50,
            "pk_contract_rows": 4,
            "canonical_tab_requests": tab_requests,
            "workspace_requests": len([url for url in first_urls + detail_urls + pk_urls if "/workspace" in url]),
            "external_search_requests_on_open": len([url for url in detail_urls if "experimental-harvest" in url]),
            "reload_persistence": "PASS",
            "browser_console_severe": "NOT_AVAILABLE_FROM_REMOTE_DRIVER",
        }
    finally:
        driver.quit()


def main():
    projects, projects_ms, _ = api("/api/projects")
    assert {1, 3, 5, 300}.issubset({row["id"] for row in projects})
    current, current_ms, _ = api("/api/prediction-engine/current")
    assert current["current_production_engine"]["engine_id"] == "drugopt-prediction-engine-v3@3.3.3"
    pk, pk_ms, pk_bytes = api("/api/compound-versions/11/scientific-tabs/pk")
    accepted = {row["canonical_endpoint"] for row in pk["rows"] if row.get("experimental")}
    assert {"HUMAN_PK_F_ORAL", "HUMAN_PK_VD_IV", "HUMAN_PK_CL_UNSPECIFIED", "HUMAN_PK_CMAX_UNSPECIFIED"} <= accepted
    with DriverService() as service:
        viewports = [exercise_viewport(service.port, 1440, 900), exercise_viewport(service.port, 390, 844)]
    result = {
        "contract": "StableCoreBrowserE2E/v1",
        "environment": "E2E_ISOLATED_DATABASE",
        "base_url": BASE_URL,
        "api_ms": {"projects": projects_ms, "current_engine": current_ms, "orforglipron_pk": pk_ms},
        "orforglipron_pk_payload_bytes": pk_bytes,
        "viewports": viewports,
        "status": "PASS",
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
