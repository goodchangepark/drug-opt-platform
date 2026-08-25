#!/bin/sh
set -eu

project_venv="/home/xavier/chem/drug-opt-platform/.venv"
project_openmp="${project_venv}/lib/python3.11/site-packages/torch.libs/libgomp-947d5fa1.so.1.0.0:${project_venv}/lib/python3.11/site-packages/scikit_learn.libs/libgomp-d22c30c5.so.1.0.0"

# Both ARM64 wheels carry different SONAMEs. Loading both at process start avoids
# the older host glibc failing later with "cannot allocate memory in static TLS block".
exec env LD_PRELOAD="${project_openmp}" "${project_venv}/bin/uvicorn" backend.main:app --host 127.0.0.1 --port 8765
