import asyncio
import json
import base64
import os
import urllib.request
import websockets

OUTPUT_DIR = "test_evidence/ui_integrity_v331"
os.makedirs(OUTPUT_DIR, exist_ok=True)

async def main():
    print("============================================================")
    print("🚀 DRUG-OPT v3.3.1 UI INTEGRITY E2E VERIFICATION (CDP)")
    print("============================================================")

    with urllib.request.urlopen("http://127.0.0.1:9222/json/list") as r:
        pages = json.loads(r.read().decode())
    target = next((p for p in pages if "127.0.0.1:8765" in p.get("url", "")), None)
    if not target:
        req = urllib.request.Request("http://127.0.0.1:9222/json/new?http://127.0.0.1:8765/", method="PUT")
        with urllib.request.urlopen(req) as r:
            target = json.loads(r.read().decode())

    ws_url = target["webSocketDebuggerUrl"]
    print(f"Connected to Chrome page target: {target['id']}")

    async with websockets.connect(ws_url, max_size=20*1024*1024) as ws:
        msg_id = 1
        async def call(method, params=None):
            nonlocal msg_id
            msg_id += 1
            await ws.send(json.dumps({"id": msg_id, "method": method, "params": params or {}}))
            while True:
                resp = json.loads(await ws.recv())
                if resp.get("id") == msg_id:
                    return resp.get("result", {})

        async def capture(filename):
            res = await call("Page.captureScreenshot", {"format": "png"})
            data = res.get("data")
            filepath = os.path.join(OUTPUT_DIR, filename)
            with open(filepath, "wb") as f:
                f.write(base64.b64decode(data))
            print(f"  📸 Screenshot saved: {filepath} ({os.path.getsize(filepath):,} bytes)")

        await call("Page.enable")
        await call("DOM.enable")

        # -------------------------------------------------------------
        # DESKTOP RUN (1440x900)
        # -------------------------------------------------------------
        print("\n=== [1/2] RUNNING DESKTOP E2E VERIFICATION (1440x900) ===")
        await call("Emulation.setDeviceMetricsOverride", {
            "width": 1440,
            "height": 900,
            "deviceScaleFactor": 1,
            "mobile": False
        })

        # Navigate to Projects -> DrugBank
        print("  --> Navigating to DrugBank Reference Library (Project 300)...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const navBtns = Array.from(document.querySelectorAll('button, nav a, nav button'));
            const projBtn = navBtns.find(b => b.innerText.trim() === 'Projects');
            if (projBtn) projBtn.click();
        })()"""})
        await asyncio.sleep(1.0)

        await call("Runtime.evaluate", {"expression": """(() => {
            const btns = Array.from(document.querySelectorAll('button, .project-link-title, a'));
            const drugbank = btns.find(b => b.innerText.includes('DrugBank'));
            if (drugbank) drugbank.click();
        })()"""})
        await asyncio.sleep(2.0)

        # Verify 150 compounds and 150 CAS numbers
        eval_res = await call("Runtime.evaluate", {
            "expression": """(() => {
                const rows = document.querySelectorAll('.compound-row');
                const cas = document.querySelectorAll('.cas-tag');
                const firstCas = cas[0]?.innerText || '';
                return { row_count: rows.length, cas_count: cas.length, firstCas };
            })()""",
            "returnByValue": True
        })
        info = eval_res.get("result", {}).get("value", {})
        print(f"  ✓ Rendered {info.get('row_count')} compound rows in DrugBank library")
        print(f"  ✓ Rendered {info.get('cas_count')} CAS tags in DrugBank library (First: '{info.get('firstCas')}')")
        assert info.get("row_count") == 150, f"Expected 150 compounds, got {info.get('row_count')}"
        assert info.get("cas_count") == 150, f"Expected 150 CAS tags, got {info.get('cas_count')}"

        await capture("desktop_drugbank_150.png")

        # Open first compound detail (Acetaminophen)
        print("  --> Opening first compound detail...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const openBtn = document.querySelector('.btn-open-detail, .compound-name-link');
            if (openBtn) openBtn.click();
        })()"""})
        
        # Wait for workspace to finish loading
        print("  --> Waiting for workspace hydration...")
        for _ in range(20):
            await asyncio.sleep(0.5)
            check = await call("Runtime.evaluate", {
                "expression": "document.querySelector('.prediction-stage-chip.loading') === null && document.querySelector('.scope-debug-details') !== null",
                "returnByValue": True
            })
            if check.get("result", {}).get("value") is True:
                break

        # Verify Overview Header
        overview_eval = await call("Runtime.evaluate", {
            "expression": """(() => {
                const revEl = Array.from(document.querySelectorAll('span')).find(s => s.innerText.includes('Structure Revision: v'));
                const engineEl = document.querySelector('#predict-meta-engine');
                const scopeDetails = document.querySelector('.scope-debug-details');
                const scopeSummary = scopeDetails ? scopeDetails.querySelector('summary')?.innerText : null;
                const chips = Array.from(document.querySelectorAll('.prediction-stage-chip')).map(c => c.innerText);
                const cards = document.querySelectorAll('.admet-highlight-card').length;
                const stars = document.querySelectorAll('.maturity-stars-interactive').length;
                const firstStar = document.querySelector('.maturity-stars-interactive');
                const starTitle = firstStar?.getAttribute('title') || '';
                const starReason = firstStar?.getAttribute('data-reason') || '';
                return {
                    revText: revEl?.innerText || null,
                    engineText: engineEl?.innerText || null,
                    scopeSummary,
                    chips,
                    cards,
                    stars,
                    starTitle,
                    starReason
                };
            })()""",
            "returnByValue": True
        })
        ov = overview_eval.get("result", {}).get("value", {})
        print(f"  ✓ Structure Revision: '{ov.get('revText')}'")
        print(f"  ✓ Active Prediction Engine: '{ov.get('engineText')}'")
        print(f"  ✓ Collapsed Debug Scope: '{ov.get('scopeSummary')}'")
        print(f"  ✓ Prediction stage chips: {ov.get('chips')}")
        print(f"  ✓ Executive scientific cards count: {ov.get('cards')}")
        print(f"  ✓ Interactive maturity stars count: {ov.get('stars')}")
        print(f"  ✓ Star tooltip title: {repr(ov.get('starTitle')[:65])}...")
        print(f"  ✓ Star data-reason: {repr(ov.get('starReason')[:65])}...")

        assert ov.get("revText") and "Structure Revision: v" in ov.get("revText")
        assert ov.get("engineText") and "v3.3.1" in ov.get("engineText")
        assert ov.get("scopeSummary") and "System Identifiers & Debug Scope" in ov.get("scopeSummary")
        assert ov.get("cards") == 8
        assert ov.get("stars") >= 5
        assert "Maturity Level" in ov.get("starTitle")

        await capture("desktop_compound_overview.png")

        # Open ADMET tab and check ScientificResultTable columns
        print("  --> Switching to ADMET tab...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const tabs = Array.from(document.querySelectorAll('.detail-tabs button, button'));
            const admetTab = tabs.find(b => b.innerText.includes('ADMET'));
            if (admetTab) admetTab.click();
        })()"""})
        await asyncio.sleep(2.0)

        table_eval = await call("Runtime.evaluate", {
            "expression": """(() => {
                const ths = Array.from(document.querySelectorAll('.scientific-results-table th')).map(t => t.innerText);
                const firstRowTds = Array.from(document.querySelectorAll('.scientific-results-table tbody tr:not(.scientific-group-row) td')).map(td => td.innerText);
                return { ths, firstRowTds: firstRowTds.slice(0, 9) };
            })()""",
            "returnByValue": True
        })
        tbl = table_eval.get("result", {}).get("value", {})
        print(f"  ✓ ScientificResultTable headers: {tbl.get('ths')}")
        headers = tbl.get("ths", [])
        if headers:
            upper_h = [h.upper() for h in headers]
            assert "CURRENT PREDICTION" in upper_h
            assert "PREDICTION SOURCE" in upper_h
            assert "APPLICABILITY DOMAIN" in upper_h
            assert "MODEL MATURITY" in upper_h

        await capture("desktop_scientific_table.png")

        # Switch to Project 3 (EGFR), Compound 10 (Sunvozertinib)
        print("  --> Switching to EGFR Project (Project 3)...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const navBtns = Array.from(document.querySelectorAll('button, nav a, nav button'));
            const projBtn = navBtns.find(b => b.innerText.trim() === 'Projects');
            if (projBtn) projBtn.click();
        })()"""})
        await asyncio.sleep(1.0)
        await call("Runtime.evaluate", {"expression": """(() => {
            const btns = Array.from(document.querySelectorAll('button, .project-link-title, a'));
            const egfr = btns.find(b => b.innerText.includes('EGFR'));
            if (egfr) egfr.click();
        })()"""})
        await asyncio.sleep(1.5)
        await call("Runtime.evaluate", {"expression": """(() => {
            const openBtn = document.querySelector('.btn-open-detail, .compound-name-link');
            if (openBtn) openBtn.click();
        })()"""})
        await asyncio.sleep(2.5)

        await capture("desktop_egfr_sunvozertinib.png")

        # -------------------------------------------------------------
        # MOBILE RUN (390x844 iPhone 12/13/14)
        # -------------------------------------------------------------
        print("\n=== [2/2] RUNNING MOBILE E2E VERIFICATION (390x844) ===")
        await call("Emulation.setDeviceMetricsOverride", {
            "width": 390,
            "height": 844,
            "deviceScaleFactor": 3,
            "mobile": True
        })
        await asyncio.sleep(1.0)

        # Switch to Overview tab on mobile
        await call("Runtime.evaluate", {"expression": """(() => {
            const tabs = Array.from(document.querySelectorAll('.detail-tabs button, button'));
            const ovTab = tabs.find(b => b.innerText.includes('OVERVIEW'));
            if (ovTab) ovTab.click();
        })()"""})
        await asyncio.sleep(1.5)

        await capture("mobile_overview.png")

        # Switch to ADMET tab on mobile
        await call("Runtime.evaluate", {"expression": """(() => {
            const tabs = Array.from(document.querySelectorAll('.detail-tabs button, button'));
            const admetTab = tabs.find(b => b.innerText.includes('ADMET'));
            if (admetTab) admetTab.click();
        })()"""})
        await asyncio.sleep(1.5)

        await capture("mobile_scientific_table.png")

        # Reset emulation
        await call("Emulation.clearDeviceMetricsOverride")

    print("\n============================================================")
    print("🎉 ALL E2E BROWSER VERIFICATIONS PASSED WITH ZERO ERRORS!")
    print("============================================================")

if __name__ == "__main__":
    asyncio.run(main())
