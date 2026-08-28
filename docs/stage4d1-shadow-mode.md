# Stage 4D-1: Shadow Mode & Production Isolation Contract

## 1. Production UI Design Freeze

Under the mandatory Stage 4D-1 **DESIGN FREEZE**, no visual modifications, new UI cards, consensus tables, layout alterations, or CSS changes are made to the Drug-OPT frontend.

The system guarantees that:
* Existing views (Dashboard, Projects, Compounds, Overview, Properties, ADMET, Metabolism, PK, Optimization, Comparison, Help) render with zero visual divergence.
* Primary visible prediction values displayed to users remain 100% identical to their historical reference values.

---

## 2. Consensus Mode Lifecycle

The ensemble consensus engine operates across three discrete operational modes:

| Mode | Database Storage | Consensus Calculation | User-Facing Production Value | Audit / Provenance API |
| :--- | :--- | :--- | :--- | :--- |
| `OFF` | Individual models only | Disabled | Primary Model Prediction | Individual model records |
| **`SHADOW` (Default)** | **All individual models** | **Computed & Persisted** | **Primary Model Prediction** | **Full consensus + provenance** |
| `ACTIVE` (Future 4D-2) | All individual models | Computed & Persisted | Multi-Model Consensus | Full consensus + provenance |

### Why Shadow Mode is Mandatory in Stage 4D-1
1. Allows thorough stress testing of multi-model execution, adapter concurrency, and database persistence across all chemical series without risk of destabilizing drug discovery programs.
2. Enables validation of model disagreement metrics ($\sigma_w$) against prospective experimental measurements before promoting consensus predictions to active production status.
3. Provides an audit trail linking individual model predictions, contract versions, and consensus weights.

---

## 3. Backward-Compatible API Architecture

Existing API endpoints continue to return identical structures:
* `POST /api/compounds/{compound_id}/predict` -> returns `{status: "COMPLETE", predictions: [...], consensus_predictions: [...]}`
* `GET /api/compound-versions/{version_id}/admet` -> returns standard ADMET payload with historical properties.

New multi-model introspection is exposed through dedicated audit routes:
* `GET /api/compound-versions/{version_id}/multimodel-provenance`:
  Returns full multi-model provenance, individual model outputs, weights, dispersion metrics, vote patterns, and shadow consensus records.
