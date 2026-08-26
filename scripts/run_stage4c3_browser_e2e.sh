#!/usr/bin/env bash
set -e

PORT=8766
echo "Starting isolated test backend on port ${PORT}..."
.venv/bin/python3 -m uvicorn backend.main:app --host 127.0.0.1 --port ${PORT} > /tmp/stage4c3_e2e_uvicorn.log 2>&1 &
SERVER_PID=$!

cleanup() {
    echo "Stopping test backend server PID ${SERVER_PID}..."
    kill ${SERVER_PID} || true
}
trap cleanup EXIT

echo "Waiting for server to respond on port ${PORT}..."
for i in {1..30}; do
    if curl -s http://127.0.0.1:${PORT}/api/health > /dev/null; then
        echo "Server is ready!"
        break
    fi
    sleep 1
done

echo "Running Stage 4C-3 Chromium E2E script..."
TEST_PORT=${PORT} .venv/bin/python3 scripts/stage4c3_browser_e2e.py

echo "Stage 4C-3 Chromium E2E Completed Successfully!"
