# AI Drug Optimization Platform — Stage 1

Stage 1 implements project management, compound versioning, SMILES validation, RDKit molecular properties, drug-likeness, structural alerts, rule-based medicinal-chemistry assessment, comparison, visualization and immutable calculation provenance.

## Run
```bash
cd /home/xavier/chem/drug-opt-platform
.venv/bin/uvicorn backend.main:app --host 127.0.0.1 --port 8765
```
Open http://127.0.0.1:8765

## Test
```bash
cd /home/xavier/chem/drug-opt-platform
.venv/bin/pytest -q
node --check frontend/static/app.js
```

## Architecture
- FastAPI + Uvicorn REST API.
- SQLAlchemy 2 / SQLite file database (`drug_opt.db`).
- RDKit 2025.03.1 chemistry engine.
- React 18 UMD dashboard served by FastAPI.
- Versioned `Project → Compound → CompoundVersion` data model.
- `PropertyCalculation`, `StructuralAlert`, and `PredictionRun` preserve provenance and audit history.

Stage 2–5 entities can attach experimental/prediction records to `CompoundVersion.id` without schema redesign.
