# Stage 4D-1: Multi-Model Execution & Adapter Architecture

## 1. Overview and Design Principles

Stage 4D-1 establishes the production-grade foundation that allows Drug-OPT to run, store, compare, and later ensemble **multiple qualified models for the same scientifically compatible endpoint**.

The system is governed by five non-negotiable principles:
1. **Preserve Every Individual Model Prediction**: Every model result is stored in its own database record with model ID, model version, execution status, runtime latency, and raw outputs. The system never saves only the consensus.
2. **Deterministic Endpoint Contracts**: Models can only be combined if they satisfy strict scientific compatibility rules (species, unit, physical assay definition, and directionality) defined in `backend/endpoint_contracts.py`.
3. **Resource-Aware Scheduling on Xavier ARM64**: Heavy PyTorch/Chemprop deep learning models are executed sequentially or under strict worker controls to prevent memory spikes and out-of-memory (OOM) faults on the 32GB Jetson AGX Xavier platform.
4. **Fault-Tolerant Execution**: A runtime failure, missing checkpoint, or out-of-domain flag in Model $M_k$ must never crash the endpoint execution or the compound Save & Predict workflow.
5. **Zero UI Redesign / Production Design Freeze**: All multi-model infrastructure executes in **SHADOW MODE** by default, ensuring existing user-facing predictions, cards, tables, and workflows remain 100% stable and visually unchanged.

---

## 2. Standardized BaseModelAdapter Framework

Every installed model is encapsulated inside a subclass of `BaseModelAdapter` (`backend/endpoint_contracts.py` & `backend/multimodel.py`), exposing a uniform lifecycle:

```python
class BaseModelAdapter(ABC):
    model_id: str
    model_name: str
    model_family: str
    model_version: str
    supported_endpoints: Set[str]
    execution_tier: ExecutionTier
    arm64_status: ARM64Status
    standardizer_version: str

    @abstractmethod
    def is_available(self) -> Tuple[bool, str]: ...

    @abstractmethod
    def execute(self, canonical_smiles: str, contract: EndpointContract) -> ModelExecutionPayload: ...
```

### Registered Adapters
* `AdmeticaChempropAdapter`: 13 ADMET endpoints (Solubility, Caco-2, PPB, CYP Panel: 1A2/2C9/2C19/2D6/3A4 inhibitors, 2C9/2D6/3A4 substrates, P-gp inhibitor, hERG liability).
* `OpenADMETClearanceAdapter`: 3 microsomal clearance endpoints (HLM, RLM, MLM continuous clearance in scaled $\log_{10}(\text{mL/min/kg})$).
* `ADMETAISafetyAdapter`: 2 toxicological safety endpoints (Ames mutagenicity, DILI clinical liability via 5-model Chemprop ensemble).
* `SyGMaMetabolismAdapter`: Phase I and Phase II SMARTS rule engine for atom-level metabolic soft spot prediction.

---

## 3. Standardized Execution Status Codes

Every model execution produces a deterministic `ExecutionStatus`:

| Status Code | Meaning | Consensus Inclusion |
| :--- | :--- | :--- |
| `SUCCESS` | Model ran cleanly, produced numerical/categorical prediction. | Included with full weight. |
| `MODEL_UNAVAILABLE` | Checkpoint missing, corrupt, or uninstalled on local system. | Excluded; weights renormalized. |
| `INCOMPATIBLE_ENDPOINT` | Model physics/units do not match requested endpoint contract. | Excluded; failure logged. |
| `OUT_OF_DOMAIN` | Compound exceeds applicability domain distance threshold. | Downweighted ($0.10\times$) or excluded. |
| `RUNTIME_ERROR` | Python exception, PyTorch inference failure, or timeout. | Excluded; isolated gracefully. |
| `INVALID_INPUT` | SMILES unparseable by RDKit sanitization engine. | Aborted with diagnostic warning. |
| `SKIPPED` | Model omitted due to execution tier or user constraint. | Excluded cleanly. |

---

## 4. Jetson Xavier ARM64 Resource Isolation & Benchmark Findings

Real measured hardware benchmarks on the 8-core ARM64 Jetson Xavier CPU (32GB Unified LPDDR4x):

| Model Family | Cold Load / Graph Init | Warm Single Compound | Batch (10 Mol) | Cache Hit | Execution Tier |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Admetica Chemprop (Single D-MPNN)** | $732.3\text{ ms}$ | $63.6\text{ ms}$ | $142.5\text{ ms}$ | $< 0.1\text{ ms}$ | `TIER_1_LOCAL_FAST` |
| **OpenADMET CheMeleon (Multi-Task MPNN)** | $610.1\text{ ms}$ | $106.7\text{ ms}$ | $225.1\text{ ms}$ | $< 0.1\text{ ms}$ | `TIER_1_LOCAL_FAST` |
| **ADMET-AI (5-Ensemble Chemprop)** | $418.3\text{ ms}$ | $195.4\text{ ms}$ | $540.2\text{ ms}$ | $< 0.1\text{ ms}$ | `TIER_2_LOCAL_HEAVY` |
| **SyGMa SMARTS Rule Engine** | $70.2\text{ ms}$ | $54.9\text{ ms}$ | $180.0\text{ ms}$ | $< 0.1\text{ ms}$ | `TIER_1_LOCAL_FAST` |
| **Full 18-Endpoint ADMET Panel** | $2,954\text{ ms}$ (Cold) | $1,570\text{ ms}$ (Warm) | $4,820\text{ ms}$ | $0.85\text{ ms}$ | Sequential Scheduling |

---

## 5. Granular Collision-Resistant Cache Key

To ensure predictions are never improperly cross-pollinated across model versions, standardizers, or chemical structures, cache keys are computed as:

$$\text{CacheKey} = \text{SHA256}(\text{version\_id} \mathbin{\Vert} \text{canonical\_smiles} \mathbin{\Vert} \text{endpoint\_id} \mathbin{\Vert} \text{model\_id} \mathbin{\Vert} \text{model\_version} \mathbin{\Vert} \text{standardizer\_version})$$
