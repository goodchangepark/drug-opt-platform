#!/usr/bin/env bash
set -euo pipefail

PORT=8766
export STAGE5A2B_BASE_URL="http://127.0.0.1:${PORT}"

echo "=== Starting background server on port ${PORT} for Stage 5A-2B E2E ==="
.venv/bin/python3 -m uvicorn backend.main:app --host 127.0.0.1 --port "${PORT}" > /dev/null 2>&1 &
SERVER_PID=$!

cleanup() {
  echo "=== Cleaning up server PID ${SERVER_PID} ==="
  kill "${SERVER_PID}" || true
}
trap cleanup EXIT

sleep 5

echo "=== Running Stage 5A-2B Chromium E2E ==="
.venv/bin/python3 scripts/stage5a2b_browser_e2e.py
