#!/usr/bin/env bash
set -euo pipefail

PORT=8766
export STAGE4C_BASE_URL="http://127.0.0.1:${PORT}"

echo "Starting isolated test backend on port ${PORT}..."
fuser -k ${PORT}/tcp || true
sleep 1

.venv/bin/python3 -m uvicorn backend.main:app --host 127.0.0.1 --port ${PORT} > /tmp/stage4c_e2e_uvicorn.log 2>&1 &
SERVER_PID=$!

cleanup() {
    echo "Stopping test backend server PID ${SERVER_PID}..."
    kill -9 ${SERVER_PID} 2>/dev/null || true
    fuser -k ${PORT}/tcp || true
}
trap cleanup EXIT

echo "Waiting for server to respond on port ${PORT}..."
for i in {1..30}; do
    if curl -s "${STAGE4C_BASE_URL}/api/health" | grep -q '"status":"ok"'; then
        echo "Server is ready!"
        break
    fi
    sleep 0.5
done

echo "Running Stage 4C Chromium E2E script..."
.venv/bin/python3 scripts/stage4c_browser_e2e.py

echo "Stage 4C Chromium E2E Completed Successfully!"
