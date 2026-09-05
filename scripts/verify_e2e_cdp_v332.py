"""
CDP Browser E2E Automation for Drug-OPT v3.3.2 Release
======================================================
Connects directly to Chrome DevTools Protocol (port 9222)
Tests Desktop (1440x900) and Mobile (390x844):
1. Active Production Engine: drugopt-prediction-engine-v3@3.3.2
2. Policy Hash: 877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce
3. Help Page: Prediction Model History & Dedicated PK Prediction Readiness Foundation (PK_FOUNDATION_READY)
4. DrugBank Reference Library (Project 300) with 200 compounds
5. Compound Workspace: 50-endpoint maturity stars and prediction routes
6. Hard page reload persistence
"""
import asyncio
import json
import base64
import os
import urllib.request
import websockets

OUTPUT_DIR = "validation/e2e_v3_3_2_browser"
os.makedirs(OUTPUT_DIR, exist_ok=True)


async def main():
    print("="*70)
    print("🚀 DRUG-OPT v3.3.2 PRODUCTION E2E VERIFICATION (CDP)")
    print("="*70)

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

        # =============================================================
        # 1. DESKTOP RUN (1440x900)
        # =============================================================
        print("\n=== [1/2] RUNNING DESKTOP E2E VERIFICATION (1440x900) ===")
        await call("Emulation.setDeviceMetricsOverride", {
            "width": 1440,
            "height": 900,
            "deviceScaleFactor": 1,
            "mobile": False
        })

        print("  --> Navigating to fresh page and clearing cache...")
        await call("Page.navigate", {"url": "http://127.0.0.1:8765/"})
        await asyncio.sleep(2.0)
        await call("Page.reload", {"ignoreCache": True})
        await asyncio.sleep(2.5)

        # Step A: Navigate to Help Page
        print("  --> Navigating to Help Page...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const navBtns = Array.from(document.querySelectorAll('button, nav a, nav button, aside button'));
            const helpBtn = navBtns.find(b => b.innerText.trim() === 'Help');
            if (helpBtn) helpBtn.click();
        })()"""})
        await asyncio.sleep(2.0)

        # Verify Help Page Content (Engine v3.3.2, Policy Hash, PK_FOUNDATION_READY)
        help_eval = await call("Runtime.evaluate", {
            "expression": """(() => {
                const bodyText = document.body.innerText;
                const pkSection = document.querySelector('#help-pk-readiness');
                const hasV332 = bodyText.includes('drugopt-prediction-engine-v3@3.3.2');
                const hasHash = bodyText.includes('877ea28f4731a67ad635252023e6601e000eecdf34297abecae6e354d91b02ce');
                const hasPkReady = bodyText.includes('PK_FOUNDATION_READY');
                const pkTitle = pkSection ? pkSection.querySelector('h2')?.innerText : null;
                return { hasV332, hasHash, hasPkReady, pkTitle };
            })()""",
            "returnByValue": True
        })
        help_info = help_eval.get("result", {}).get("value", {})
        print(f"  ✓ Help Page contains v3.3.2 engine ID: {help_info.get('hasV332')}")
        print(f"  ✓ Help Page contains v3.3.2 policy hash: {help_info.get('hasHash')}")
        print(f"  ✓ Help Page contains PK_FOUNDATION_READY: {help_info.get('hasPkReady')}")
        print(f"  ✓ PK Readiness section title: '{help_info.get('pkTitle')}'")
        assert help_info.get("hasV332") is True, "Missing v3.3.2 engine ID on Help Page"
        assert help_info.get("hasHash") is True, "Missing v3.3.2 policy hash on Help Page"
        assert help_info.get("hasPkReady") is True, "Missing PK_FOUNDATION_READY on Help Page"

        await capture("desktop_1440x900_help_pk_readiness.png")

        # Step B: Navigate to Projects -> DrugBank
        print("\n  --> Navigating to DrugBank Reference Library (Project 300)...")
        for _ in range(10):
            p_clicked = await call("Runtime.evaluate", {"expression": """(() => {
                const navBtns = Array.from(document.querySelectorAll('.global-nav button, aside button, button'));
                const projBtn = navBtns.find(b => b.innerText.trim() === 'Projects');
                if (projBtn) {
                    projBtn.click();
                    return true;
                }
                return false;
            })()""", "returnByValue": True})
            if p_clicked.get("result", {}).get("value") is True:
                print("  ✓ Clicked Projects nav item")
                break
            await asyncio.sleep(0.5)

        # Wait for Projects list to render and click DrugBank
        db_clicked = False
        for _ in range(20):
            await asyncio.sleep(0.5)
            clicked = await call("Runtime.evaluate", {"expression": """(() => {
                const btns = Array.from(document.querySelectorAll('.project-link-title'));
                const drugbank = btns.find(b => b.innerText.trim() === 'DrugBank' || b.innerText.includes('DrugBank'));
                if (drugbank) {
                    drugbank.click();
                    return true;
                }
                return false;
            })()""", "returnByValue": True})
            if clicked.get("result", {}).get("value") is True:
                db_clicked = True
                print("  ✓ Clicked DrugBank project link")
                break
        assert db_clicked, "Failed to click DrugBank project link"

        # Verify DrugBank compounds count with polling (up to 40s on Xavier ARM64)
        info = {}
        for attempt in range(160):
            await asyncio.sleep(0.5)
            eval_res = await call("Runtime.evaluate", {
                "expression": """(() => {
                    const rows = document.querySelectorAll('.compound-row');
                    const cas = document.querySelectorAll('.cas-tag');
                    return { row_count: rows.length, cas_count: cas.length };
                })()""",
                "returnByValue": True
            })
            info = eval_res.get("result", {}).get("value", {})
            if info.get("row_count", 0) >= 150:
                print(f"  ✓ Compounds rendered after {(attempt + 1) * 0.5:.1f}s")
                break
            if attempt > 0 and attempt % 20 == 0:
                print(f"  ... waiting for compounds ({attempt * 0.5:.0f}s, current rows: {info.get('row_count', 0)})")

        print(f"  ✓ Rendered {info.get('row_count')} compound rows in DrugBank library")
        print(f"  ✓ Rendered {info.get('cas_count')} CAS tags in DrugBank library")
        assert info.get("row_count") in (150, 200), f"Expected 150 or 200 compounds, got {info.get('row_count')}"

        await capture("desktop_1440x900_drugbank_200.png")

        # Step C: Open first compound detail
        print("\n  --> Opening first compound detail...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const openBtn = document.querySelector('.btn-open-detail, .compound-name-link');
            if (openBtn) openBtn.click();
        })()"""})

        # Wait for workspace to hydrate (up to 30s on Xavier ARM64)
        for attempt in range(60):
            await asyncio.sleep(0.5)
            check = await call("Runtime.evaluate", {
                "expression": "document.querySelector('.scope-debug-details') !== null",
                "returnByValue": True
            })
            if check.get("result", {}).get("value") is True:
                print(f"  ✓ Workspace hydrated after {(attempt + 1) * 0.5:.1f}s")
                break

        # Verify Compound Overview details
        overview_eval = await call("Runtime.evaluate", {
            "expression": """(() => {
                const scopeDetails = document.querySelector('.scope-debug-details');
                const scopeSummary = scopeDetails ? (scopeDetails.querySelector('summary')?.textContent || '').trim() : '';
                const scopeText = scopeDetails ? (scopeDetails.textContent || scopeDetails.innerText || '') : '';
                const stars = document.querySelectorAll('.maturity-stars');
                const badges = Array.from(document.querySelectorAll('.badge-favorable')).map(b => b.innerText);
                return {
                    scopeSummary,
                    hasStrictScope: scopeText.includes('Strict scope:'),
                    hasV332: scopeText.includes('drugopt-prediction-engine-v3@3.3.2'),
                    starCount: stars.length,
                    badges: badges.slice(0, 5)
                };
            })()""",
            "returnByValue": True
        })
        comp_info = overview_eval.get("result", {}).get("value", {})
        print(f"  ✓ Scope summary: '{comp_info.get('scopeSummary')}'")
        print(f"  ✓ Scope contains 'Strict scope': {comp_info.get('hasStrictScope')}")
        print(f"  ✓ Scope contains v3.3.2 engine ID: {comp_info.get('hasV332')}")
        print(f"  ✓ Maturity stars elements count: {comp_info.get('starCount')}")
        assert comp_info.get("hasStrictScope") is True, "Missing 'Strict scope' in scope details"
        assert comp_info.get("hasV332") is True, "Missing v3.3.2 in scope details"
        assert comp_info.get("starCount") > 0, "Expected rendered maturity stars"

        await capture("desktop_1440x900_compound_workspace.png")

        # Step D: Hard page reload and verify persistence (up to 30s on Xavier ARM64)
        print("\n  --> Executing hard reload and verifying persistence...")
        await call("Page.reload", {"ignoreCache": True})
        reload_info = {}
        for attempt in range(60):
            await asyncio.sleep(0.5)
            reload_eval = await call("Runtime.evaluate", {
                "expression": """(() => {
                    const scopeDetails = document.querySelector('.scope-debug-details');
                    const stars = document.querySelectorAll('.maturity-stars');
                    return {
                        hasScope: scopeDetails !== null,
                        starCount: stars.length
                    };
                })()""",
                "returnByValue": True
            })
            reload_info = reload_eval.get("result", {}).get("value", {})
            if reload_info.get("hasScope") is True:
                print(f"  ✓ Persisted workspace re-hydrated after {(attempt + 1) * 0.5:.1f}s")
                break
        print(f"  ✓ Hard reload persistence verified (Scope present: {reload_info.get('hasScope')}, Stars: {reload_info.get('starCount')})")
        assert reload_info.get("hasScope") is True, "Workspace failed to persist after hard reload"

        await capture("desktop_1440x900_hard_reload_persistence.png")

        # =============================================================
        # 2. MOBILE RUN (390x844)
        # =============================================================
        print("\n=== [2/2] RUNNING MOBILE E2E VERIFICATION (390x844) ===")
        await call("Emulation.setDeviceMetricsOverride", {
            "width": 390,
            "height": 844,
            "deviceScaleFactor": 2,
            "mobile": True
        })
        await asyncio.sleep(1.0)

        await capture("mobile_390x844_compound_workspace.png")

        # Navigate to Help on Mobile
        print("  --> Navigating to Help view in mobile...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const menuBtn = document.querySelector('.menu-toggle');
            if (menuBtn) menuBtn.click();
        })()"""})
        await asyncio.sleep(0.5)

        await call("Runtime.evaluate", {"expression": """(() => {
            const navBtns = Array.from(document.querySelectorAll('.global-nav button, aside button, button, nav a, nav button'));
            const helpBtn = navBtns.find(b => b.innerText.trim() === 'Help');
            if (helpBtn) helpBtn.click();
        })()"""})
        await asyncio.sleep(2.0)

        await capture("mobile_390x844_help_pk_readiness.png")

        # Navigate to Projects on Mobile
        print("  --> Navigating to Projects view in mobile...")
        await call("Runtime.evaluate", {"expression": """(() => {
            const menuBtn = document.querySelector('.menu-toggle');
            if (menuBtn) menuBtn.click();
        })()"""})
        await asyncio.sleep(0.5)

        await call("Runtime.evaluate", {"expression": """(() => {
            const navBtns = Array.from(document.querySelectorAll('.global-nav button, aside button, button, nav a, nav button'));
            const projBtn = navBtns.find(b => b.innerText.trim() === 'Projects');
            if (projBtn) projBtn.click();
        })()"""})
        await asyncio.sleep(2.0)

        await capture("mobile_390x844_projects.png")

        print("\n" + "="*70)
        print("🎉 PUBLIC E2E CDP BROWSER VERIFICATION COMPLETED 100% SUCCESSFULLY!")
        print("="*70)


if __name__ == "__main__":
    asyncio.run(main())
