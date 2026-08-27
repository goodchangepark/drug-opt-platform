# Stage 5B-4 Platform Stabilization

This release stabilizes the implemented platform before Stage 5C or Stage 6. It does not add a scientific stage, retrain a model, introduce a new endpoint, or modernize dependencies.

## Stabilization changes

- Unified Dashboard capability state with the backend model and feature registries.
- Added `/api/help/registry`, a runtime-backed researcher Help inventory.
- Rebuilt the sidebar Help view around installed modules, models, versions, workflows, terminology, and limitations.
- Added repeatable model-runtime, project-cleanup, database-integrity, API/workflow, and Chromium acceptance audits.
- Conservatively classified development projects and used the production cascade-delete service only for confirmed test fixtures.
- Repaired one legacy orphaned Stage 5B-3 browser-fixture snapshot after positive fixture identification.
- Set the stabilization application version to `0.6.0-stage5b4-stable`; service stage remains `5B-4`.

## Evidence artifacts

- `validation/stage5b4_stabilization_model_audit.json`
- `validation/test_projects_cleanup.json`
- `validation/stage5b4_stabilization_database_audit.json`
- `validation/stage5b4_stabilization_browser_e2e_results.json`

The model audit requires registry presence, packaged assets, successful loader execution, finite CPU inference, and endpoint/version/unit agreement. The cleanup artifact records every pre-cleanup project classification and preserves ambiguous projects.

## Validation policy

The release must pass targeted stabilization tests, the complete regression suite, Chromium navigation and workflow acceptance, zero uncaught JavaScript errors, SQLite integrity and foreign-key checks, and `/api/health`. The controlled project `__STABILIZATION_E2E_TEMP__` must be removed through the production project-delete workflow after acceptance.

Warnings from optional performance/logging suggestions are classified as known nonblocking. Scientific availability is never inferred from confidence or conformal status.

