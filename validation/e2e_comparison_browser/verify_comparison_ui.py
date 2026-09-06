import os
import time
from playwright.sync_api import sync_playwright

def main():
    os.makedirs("validation/e2e_comparison_browser", exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        
        # 1. Desktop Viewport (1440x900)
        context_desktop = browser.new_context(viewport={"width": 1440, "height": 900})
        page_desktop = context_desktop.new_page()
        print("Navigating to Desktop URL...")
        page_desktop.goto("http://127.0.0.1:8765/#/projects/3/compounds/10", timeout=30000)
        time.sleep(3)
        page_desktop.screenshot(path="validation/e2e_comparison_browser/sunvozertinib_desktop.png", full_page=True)
        print("Desktop screenshot captured successfully!")

        # 2. Mobile Viewport (390x844)
        context_mobile = browser.new_context(viewport={"width": 390, "height": 844})
        page_mobile = context_mobile.new_page()
        print("Navigating to Mobile URL...")
        page_mobile.goto("http://127.0.0.1:8765/#/projects/3/compounds/10", timeout=30000)
        time.sleep(3)
        page_mobile.screenshot(path="validation/e2e_comparison_browser/sunvozertinib_mobile.png", full_page=True)
        print("Mobile screenshot captured successfully!")

        # 3. Help Page Desktop
        print("Navigating to Help Page URL...")
        page_desktop.goto("http://127.0.0.1:8765/#/help", timeout=30000)
        time.sleep(2)
        page_desktop.screenshot(path="validation/e2e_comparison_browser/help_comparison_section.png", full_page=True)
        print("Help section screenshot captured successfully!")

        browser.close()
        print("ALL_SCREENSHOTS_CAPTURED")

if __name__ == "__main__":
    main()
