#!/usr/bin/env python3
"""Restartable, sequential completion campaign for Prediction Engine v3.3.3.

The runner intentionally executes one bounded task at a time.  It persists
state before and after every command, refuses concurrent invocations, and
stops on a failed scientific or validation gate.  It is suitable for a single
user-systemd service and never modifies the runtime database directly.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = ROOT / "validation" / "v333_completion_campaign.json"
LOCK_PATH = Path("/tmp/drugopt-v333-completion-campaign.lock")


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_state() -> dict:
    return json.loads(STATE_PATH.read_text(encoding="utf-8"))


def save_state(state: dict) -> None:
    state["updated_at"] = now()
    temporary = STATE_PATH.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, STATE_PATH)


def next_task(state: dict) -> dict | None:
    return next((task for task in state["tasks"] if task["status"] in {"PENDING", "RETRY"}), None)


def run_task(state: dict, task: dict) -> bool:
    task["status"] = "RUNNING"
    task["started_at"] = now()
    state["phase"] = task["phase"]
    state["current_task"] = task["id"]
    state["exact_next_action"] = f"complete {task['id']}"
    save_state(state)
    completed = subprocess.run(task["command"], cwd=ROOT, text=True, capture_output=True)
    task["finished_at"] = now()
    task["returncode"] = completed.returncode
    task["stdout_tail"] = completed.stdout[-4000:]
    task["stderr_tail"] = completed.stderr[-4000:]
    task["status"] = "COMPLETE" if completed.returncode == 0 else "FAILED"
    if completed.returncode == 0:
        state.setdefault("completed_tasks", []).append(task["id"])
        pending = next_task(state)
        state["exact_next_action"] = f"run {pending['id']}" if pending else "final release acceptance audit"
    else:
        state["exact_next_action"] = f"triage {task['id']} failure; then set status RETRY"
    save_state(state)
    return completed.returncode == 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--max-tasks", type=int, default=1)
    args = parser.parse_args()
    state = load_state()
    if args.retry_failed:
        failed = next((task for task in state["tasks"] if task["status"] == "FAILED"), None)
        if failed:
            failed["status"] = "RETRY"
            failed["retry_requested_at"] = now()
            state["exact_next_action"] = f"retry {failed['id']}"
            save_state(state)
    if args.status or not args.run:
        print(json.dumps({
            "campaign": state["campaign"],
            "phase": state["phase"],
            "current_task": state.get("current_task"),
            "completed_tasks": state.get("completed_tasks", []),
            "next_task": (next_task(state) or {}).get("id"),
            "exact_next_action": state["exact_next_action"],
        }, indent=2))
        return 0
    LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("another v3.3.3 campaign runner is active")
            return 2
        completed = 0
        while completed < max(1, args.max_tasks):
            task = next_task(state)
            if not task:
                break
            if not run_task(state, task):
                return 1
            completed += 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
