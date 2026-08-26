#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
e2e_port="${PRE_STAGE5_E2E_PORT:-8766}"
server_log="${TMPDIR:-/tmp}/pre-stage5-ui-server.log"

cd "${project_root}"
project_openmp="${project_root}/.venv/lib/python3.11/site-packages/torch.libs/libgomp-947d5fa1.so.1.0.0:${project_root}/.venv/lib/python3.11/site-packages/scikit_learn.libs/libgomp-d22c30c5.so.1.0.0"
env LD_PRELOAD="${project_openmp}" \
  "${project_root}/.venv/bin/uvicorn" backend.main:app --host 127.0.0.1 --port "${e2e_port}" >"${server_log}" 2>&1 &
server_pid=$!
cleanup(){ kill "${server_pid}" 2>/dev/null || true; wait "${server_pid}" 2>/dev/null || true; }
trap cleanup EXIT

for _ in $(seq 1 90); do
  if curl -fsS "http://127.0.0.1:${e2e_port}/api/health" >/dev/null; then
    PRE_STAGE5_BASE_URL="http://127.0.0.1:${e2e_port}" "${project_root}/.venv/bin/python" scripts/pre_stage5_ui_browser_e2e.py
    exit $?
  fi
  sleep 1
done

echo "Pre-Stage 5 UI test server did not become healthy" >&2
tail -80 "${server_log}" >&2
exit 1
