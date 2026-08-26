#!/usr/bin/env python3
"""Actual Chromium Stage 4B proposal workflow against the production service."""

import json
import os
import sys
import time

from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


BASE_URL = os.environ.get("STAGE4B_BASE_URL", "http://127.0.0.1:8765")
PROJECT_NAME = "Stage 4B Acceptance — Public References"


def browser_api(driver, method, path, payload=None, timeout=300):
    driver.set_script_timeout(timeout)
    script = """
const done=arguments[arguments.length-1], method=arguments[0], path=arguments[1], payload=arguments[2];
fetch('/api'+path,{method,headers:{'Content-Type':'application/json'},body:payload===null?undefined:JSON.stringify(payload)})
 .then(async response=>{const text=await response.text();if(!response.ok)throw new Error(response.status+' '+text);return text?JSON.parse(text):null})
 .then(data=>done({ok:true,data})).catch(error=>done({ok:false,error:String(error)}));
"""
    result = driver.execute_async_script(script, method, path, payload)
    if not result["ok"]:
        raise RuntimeError(result["error"])
    return result["data"]


def click_button(driver, text, last=False, timeout=60):
    def match(d):
        matches = []
        for item in d.find_elements(By.TAG_NAME, "button"):
            try:
                item_text = item.text.strip()
                if item.is_displayed() and (item_text == text or item_text.startswith(text + "\n")):
                    matches.append(item)
            except StaleElementReferenceException:
                return None
        return (matches[-1] if last else matches[0]) if matches else None
    button = WebDriverWait(driver, timeout).until(match)
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", button)
    button.click()
    return button


def set_labeled_input(driver, label, value):
    result = driver.execute_script("""
const wanted=arguments[0], value=arguments[1];
const label=[...document.querySelectorAll('label')].find(row=>row.textContent.trim()===wanted);
if(!label)return {ok:false,labels:[...document.querySelectorAll('label')].filter(row=>row.offsetParent!==null).map(row=>row.textContent.trim())};
const input=label.nextElementSibling;
if(!input || input.tagName!=='INPUT')return {ok:false,labels:['label found but no sibling input']};
input.scrollIntoView({block:'center'});
Object.getOwnPropertyDescriptor(HTMLInputElement.prototype,'value').set.call(input,String(value));
input.dispatchEvent(new Event('input',{bubbles:true}));
input.dispatchEvent(new Event('change',{bubbles:true}));
return {ok:true};
""", label, value)
    if not result["ok"]:
        raise AssertionError(f"Input label {label!r} not found; available labels: {result['labels']}")


def wait_for_proposal(driver, proposal_id, timeout=360):
    deadline = time.time() + timeout
    statuses = []
    while time.time() < deadline:
        proposal = browser_api(driver, "GET", f"/proposals/{proposal_id}?view=all")
        if not statuses or statuses[-1] != proposal["status"]:
            statuses.append(proposal["status"])
        if proposal["status"] == "COMPLETED":
            return proposal, statuses
        if proposal["status"] == "FAILED":
            raise AssertionError("Proposal job failed: " + proposal.get("stage_message", ""))
        time.sleep(1.5)
    raise TimeoutError(f"Proposal {proposal_id} did not complete; observed {statuses}")


def main():
    options = webdriver.ChromeOptions()
    options.binary_location = "/snap/chromium/current/usr/lib/chromium-browser/chrome"
    for argument in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage", "--window-size=1800,1600"):
        options.add_argument(argument)
    driver = webdriver.Chrome(service=Service("/snap/bin/chromium.chromedriver"), options=options)
    checks = []
    try:
        driver.get(BASE_URL + "/?stage4b-e2e=20260826")
        WebDriverWait(driver, 45).until(lambda d: "AI Drug Optimization Platform" in d.page_source)
        health = browser_api(driver, "GET", "/health")
        if health.get("step") != "4B":
            raise AssertionError(f"Unexpected health payload: {health}")

        projects = browser_api(driver, "GET", "/projects")
        project = next((row for row in projects if row["name"] == PROJECT_NAME), None)
        if not project:
            project = browser_api(driver, "POST", "/projects", {
                "name": PROJECT_NAME, "target": "Deterministic proposal workflow",
                "description": "Public lidocaine structure; no proprietary compounds",
            })
        detail = browser_api(driver, "GET", f"/projects/{project['id']}")
        parent = next((row for row in detail["compounds"] if row["compound_id"] == "OPT-LIDOCAINE"), None)
        if not parent:
            browser_api(driver, "POST", f"/projects/{project['id']}/compounds", {
                "compound_id": "OPT-LIDOCAINE", "name": "Lidocaine public direction example",
                "smiles": "CCN(CC)C(=O)c1c(C)cccc1C", "notes": "Public Stage 4B browser fixture",
            })
            detail = browser_api(driver, "GET", f"/projects/{project['id']}")
            parent = next(row for row in detail["compounds"] if row["compound_id"] == "OPT-LIDOCAINE")
        version_id = parent["version"]["id"]
        admet = browser_api(driver, "GET", f"/projects/{project['id']}/admet")
        endpoint_by_id = {row["id"]: row["name"] for row in admet["endpoints"]}
        if not any(row["version_id"] == version_id and endpoint_by_id.get(row["endpoint_id"]) == "HLM intrinsic clearance" for row in admet["measurements"]):
            browser_api(driver, "POST", f"/projects/{project['id']}/admet/measurements", {
                "version_id": version_id, "endpoint": "HLM intrinsic clearance", "species": "Human",
                "matrix": "HLM", "value": 2.2, "unit": "log10(mL/min/kg)",
                "method": "directional browser fixture", "source": "Public reference direction",
            })

        driver.get(BASE_URL + "/?stage4b-e2e=20260826")
        WebDriverWait(driver, 45).until(lambda d: PROJECT_NAME in d.page_source)
        click_button(driver, PROJECT_NAME)
        WebDriverWait(driver, 30).until(lambda d: "OPT-LIDOCAINE" in d.page_source)
        driver.find_element(By.XPATH, "//tr[td[contains(.,'OPT-LIDOCAINE')]]//button[normalize-space()='Open']").click()
        WebDriverWait(driver, 30).until(lambda d: "Physicochemical properties" in d.page_source)
        click_button(driver, "OPTIMIZATION", last=True)
        WebDriverWait(driver, 45).until(lambda d: "Optimization Run" in d.page_source)
        set_labeled_input(driver, "TPSA minimum Å²", "")
        click_button(driver, "Analyze strategy")
        WebDriverWait(driver, 60).until(lambda d: "Recommended transformations" in d.page_source and "HLM metabolic instability" in d.page_source)
        set_labeled_input(driver, "Maximum raw candidates (1–200)", 12)
        checks.extend(["Parent selection", "Objective/constraints", "Stage 4A strategy"])

        click_button(driver, "Generate analogs")
        WebDriverWait(driver, 30).until(lambda d: "Analog proposal job queued" in d.page_source)
        optimization_data = browser_api(driver, "GET", f"/projects/{project['id']}/optimization?version_id={version_id}")
        optimization_id = optimization_data["runs"][0]["id"]
        runs = browser_api(driver, "GET", f"/optimization/runs/{optimization_id}/proposals")
        proposal_id = runs["proposal_runs"][0]["id"]
        proposal, statuses = wait_for_proposal(driver, proposal_id)
        if proposal["raw_candidate_count"] < 1 or proposal["accepted_count"] < 1 or proposal["top_count"] < 1:
            reasons = [
                (row["candidate_number"], [reason["code"] for reason in row["rejection_reasons"]])
                for row in proposal["candidates"]
            ]
            raise AssertionError(
                f"Unexpected completed counts: raw={proposal['raw_candidate_count']} "
                f"accepted={proposal['accepted_count']} top={proposal['top_count']}; reasons={reasons}"
            )
        if statuses[0] not in {"PENDING", "GENERATING", "FILTERING", "PREDICTING", "RANKING"}:
            raise AssertionError(f"Unexpected job lifecycle: {statuses}")
        checks.extend(["Background job", "Chemical filtering", "Stage 1/2/3 rescoring", "Pareto/ranking", "Top 10"])

        driver.refresh()
        WebDriverWait(driver, 45).until(lambda d: PROJECT_NAME in d.page_source)
        click_button(driver, PROJECT_NAME)
        WebDriverWait(driver, 30).until(lambda d: "OPT-LIDOCAINE" in d.page_source)
        driver.find_element(By.XPATH, "//tr[td[contains(.,'OPT-LIDOCAINE')]]//button[normalize-space()='Open']").click()
        click_button(driver, "OPTIMIZATION", last=True)
        WebDriverWait(driver, 60).until(lambda d: "Candidate Filtering" in d.page_source and "Show Top 10" in d.page_source)
        for label in ("Show all generated", "Show accepted", "Show rejected", "Show Pareto front", "Show Top 10"):
            click_button(driver, label)
        click_button(driver, "Details")
        WebDriverWait(driver, 30).until(lambda d: "Parent vs Candidate" in d.page_source and "Parent difference" in d.page_source and "Candidate difference" in d.page_source)
        required = ["Activity prediction", "Stage 1 property changes", "Soft spot changes", "Safety flags", "Synthetic complexity", "Information Value", "Ranking formula and prediction provenance"]
        missing = [label for label in required if label not in driver.page_source]
        if missing:
            raise AssertionError("Candidate details missing: " + ", ".join(missing))
        checks.extend(["Candidate filters", "Candidate detail", "Parent comparison", "Structure difference"])

        click_button(driver, "Promote")
        WebDriverWait(driver, 30).until(lambda d: "Candidate decision saved" in d.page_source)
        selected = browser_api(driver, "GET", f"/proposals/{proposal_id}?view=top10")["candidates"]
        promoted = next((row for row in selected if row["user_decision"] == "PROMOTED"), None)
        if not promoted:
            raise AssertionError("Manual promotion was not persisted")
        browser_api(driver, "PATCH", f"/proposal-candidates/{promoted['id']}/decision", {"decision": "REJECTED", "reason": "Browser E2E manual rejection fixture"})
        rejected = browser_api(driver, "GET", f"/proposals/{proposal_id}?view=rejected")["candidates"]
        if not any(row["id"] == promoted["id"] and any(reason["code"] == "USER_REJECTED" for reason in row["rejection_reasons"]) for row in rejected):
            raise AssertionError("Manual rejection and reason were not persisted")
        checks.extend(["Manual promote", "Manual reject reason"])

        user_candidate = browser_api(driver, "POST", f"/proposals/{proposal_id}/candidates", {
            "smiles": "CCN(CC)C(=O)c1c(C)cc(Cl)cc1C", "reason": "Browser E2E user-authored analog",
        }, timeout=300)
        if not user_candidate["user_added"] or not user_candidate["stage1"] or not user_candidate["admet"]:
            raise AssertionError("User-added analog did not receive the common rescoring workflow")
        checks.append("User-added analog common workflow")

        print(json.dumps({
            "status": "PASS", "base_url": BASE_URL, "project": PROJECT_NAME,
            "proposal_id": proposal_id, "job_states_observed": statuses,
            "raw": proposal["raw_candidate_count"], "accepted": proposal["accepted_count"],
            "rejected": proposal["rejected_count"], "top": proposal["top_count"], "checks": checks,
        }, indent=2))
    finally:
        driver.quit()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(json.dumps({"status": "FAIL", "error": str(exc)}, indent=2), file=sys.stderr)
        raise
