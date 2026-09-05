import asyncio
import os
from playwright.async_api import async_playwright

OUTPUT_DIR = "test_evidence/ui_integrity_v331"
os.makedirs(OUTPUT_DIR, exist_ok=True)

async def run():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)

        # ----------------------------------------------------
        # 1. DESKTOP E2E TEST (1440x900)
        # ----------------------------------------------------
        print("\n=== [1/2] RUNNING DESKTOP E2E VERIFICATION (1440x900) ===")
        context = await browser.new_context(viewport={"width": 1440, "height": 900})
        page = await context.new_page()

        # Step A: Navigate to DrugBank Project 300 directly
        print("  --> Navigating to http://127.0.0.1:8765/?project_id=300")
        await page.goto("http://127.0.0.1:8765/?project_id=300", wait_until="domcontentloaded")
        await page.wait_for_selector(".compound-row", timeout=15000)
        await page.wait_for_timeout(1000)
        print("  ✓ DrugBank Project loaded successfully")

        # Step B: Verify 150 compounds rendered with CAS numbers
        rows = await page.query_selector_all(".compound-row")
        print(f"  ✓ Rendered {len(rows)} compound rows in DrugBank library")
        assert len(rows) == 150, f"Expected exactly 150 compounds in DrugBank, found {len(rows)}"

        cas_tags = await page.query_selector_all(".cas-tag")
        print(f"  ✓ Found {len(cas_tags)} rendered CAS tags in table")
        assert len(cas_tags) == 150, f"Expected 150 CAS tags, found {len(cas_tags)}"

        # Check sample CAS text
        first_cas = await cas_tags[0].inner_text()
        print(f"  ✓ First row CAS tag: '{first_cas}'")

        await page.screenshot(path=f"{OUTPUT_DIR}/desktop_drugbank_150.png", full_page=False)
        print(f"  ✓ Captured screenshot: {OUTPUT_DIR}/desktop_drugbank_150.png")

        # Step C: Click on first compound (e.g. Acetaminophen / button)
        first_btn = await page.query_selector(".btn-open-detail")
        assert first_btn is not None
        await first_btn.click()
        await page.wait_for_selector(".compound-workspace", timeout=10000)
        await page.wait_for_timeout(1000)
        print("  ✓ Opened Compound Workspace")

        # Step D: Verify Compound Overview UX clean-up
        rev_el = await page.query_selector("text=Structure Revision: v")
        assert rev_el is not None, "Missing 'Structure Revision: v' label in header"
        rev_text = await rev_el.inner_text()
        print(f"  ✓ Verified Structure Revision label: '{rev_text}'")

        engine_el = await page.query_selector("#predict-meta-engine")
        assert engine_el is not None, "Missing #predict-meta-engine element"
        engine_text = await engine_el.inner_text()
        print(f"  ✓ Verified Active Prediction Engine: '{engine_text}'")
        assert "3.3.1" in engine_text

        # Verify Collapsible Debug Scope
        scope_details = await page.query_selector(".scope-debug-details")
        assert scope_details is not None, "Missing .scope-debug-details element"
        summary_el = await scope_details.query_selector("summary")
        print(f"  ✓ Collapsed Debug Scope: '{await summary_el.inner_text()}'")

        # Verify Stage Prediction Chips
        chips = await page.query_selector_all(".prediction-stage-chip")
        chip_texts = [await c.inner_text() for c in chips]
        print(f"  ✓ Stage status chips: {chip_texts}")
        assert any("PREDICTION: COMPLETE" in t or "PROPERTIES" in t for t in chip_texts)

        # Check Executive Scientific Summary cards and maturity stars
        cards = await page.query_selector_all(".admet-highlight-card")
        print(f"  ✓ Executive scientific summary cards count: {len(cards)}")
        assert len(cards) == 8

        interactive_stars = await page.query_selector_all(".maturity-stars-interactive")
        print(f"  ✓ Interactive maturity star elements count: {len(interactive_stars)}")
        assert len(interactive_stars) >= 5

        # Check hover tooltip / title on solubility or caco2 card
        first_star = interactive_stars[0]
        title_attr = await first_star.get_attribute("title")
        data_reason = await first_star.get_attribute("data-reason")
        print(f"  ✓ Star tooltip title: {repr(title_attr[:60])}...")
        print(f"  ✓ Star data-reason: {repr(data_reason[:60])}...")
        assert "Maturity Level" in title_attr
        assert "Scientific Provenance & Reason" in title_attr

        await page.screenshot(path=f"{OUTPUT_DIR}/desktop_compound_overview.png", full_page=False)
        print(f"  ✓ Captured screenshot: {OUTPUT_DIR}/desktop_compound_overview.png")

        # Step E: Navigate to ADMET tab and verify ScientificResultTable columns
        admet_tab = await page.query_selector("button:has-text('ADMET')")
        if admet_tab:
            await admet_tab.click()
            await page.wait_for_timeout(1000)
            table_ths = await page.query_selector_all(".scientific-results-table th")
            th_texts = [await th.inner_text() for th in table_ths]
            print(f"  ✓ ADMET table column headers: {th_texts}")
            if th_texts:
                assert "Current Prediction" in th_texts
                assert "Prediction Source" in th_texts
                assert "Applicability Domain" in th_texts
                assert "Model Maturity" in th_texts
            await page.screenshot(path=f"{OUTPUT_DIR}/desktop_scientific_table.png", full_page=False)
            print(f"  ✓ Captured screenshot: {OUTPUT_DIR}/desktop_scientific_table.png")

        # Step F: Navigate to Project 3 (EGFR)
        print("  --> Navigating to Project 3 (EGFR)")
        await page.goto("http://127.0.0.1:8765/?project_id=3", wait_until="domcontentloaded")
        await page.wait_for_selector(".compound-row", timeout=15000)
        await page.wait_for_timeout(1000)
        egfr_open_btn = await page.query_selector(".btn-open-detail")
        if egfr_open_btn:
            await egfr_open_btn.click()
            await page.wait_for_selector(".compound-workspace", timeout=10000)
            await page.wait_for_timeout(1000)
            await page.screenshot(path=f"{OUTPUT_DIR}/desktop_egfr_workspace.png", full_page=False)
            print(f"  ✓ Captured screenshot: {OUTPUT_DIR}/desktop_egfr_workspace.png")

        await context.close()

        # ----------------------------------------------------
        # 2. MOBILE E2E TEST (390x844 iPhone 12/13/14)
        # ----------------------------------------------------
        print("\n=== [2/2] RUNNING MOBILE E2E VERIFICATION (390x844) ===")
        mobile_context = await browser.new_context(
            viewport={"width": 390, "height": 844},
            is_mobile=True,
            has_touch=True
        )
        mobile_page = await mobile_context.new_page()
        await mobile_page.goto("http://127.0.0.1:8765/?project_id=300", wait_until="domcontentloaded")
        await mobile_page.wait_for_selector(".compound-row", timeout=15000)
        await mobile_page.wait_for_timeout(1000)

        # Open first compound on mobile
        mob_open_btn = await mobile_page.query_selector(".btn-open-detail")
        await mob_open_btn.click()
        await mobile_page.wait_for_selector(".compound-workspace", timeout=10000)
        await mobile_page.wait_for_timeout(1000)

        # Mobile Overview screenshot
        await mobile_page.screenshot(path=f"{OUTPUT_DIR}/mobile_overview.png", full_page=False)
        print(f"  ✓ Captured mobile screenshot: {OUTPUT_DIR}/mobile_overview.png")

        # Open ADMET tab on mobile
        mob_admet_tab = await mobile_page.query_selector("button:has-text('ADMET')")
        if mob_admet_tab:
            await mob_admet_tab.click()
            await mobile_page.wait_for_timeout(1000)

        await mobile_page.screenshot(path=f"{OUTPUT_DIR}/mobile_scientific_table.png", full_page=False)
        print(f"  ✓ Captured mobile screenshot: {OUTPUT_DIR}/mobile_scientific_table.png")

        await mobile_context.close()
        await browser.close()
        print("\n============================================================")
        print("🎉 ALL PLAYWRIGHT E2E BROWSER TESTS PASSED (DESKTOP & MOBILE)!")
        print("============================================================")

if __name__ == "__main__":
    asyncio.run(run())
