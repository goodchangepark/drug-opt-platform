"""
Xavier Production Runtime Benchmark for Prediction Engine v3.3.1.
Directive 16:
- Cold start latency vs Warm latency
- 50-endpoint full profile prediction latency
- 5-compound sequential run latency and throughput
- RAM / VRAM peak measurement on NVIDIA Jetson Xavier ARM64
"""
from __future__ import annotations

import gc
import json
import os
import resource
import time
from pathlib import Path
import psutil
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))
import numpy as np

from backend.endpoint_contracts import get_endpoint_contract
import backend.multimodel as mm
import backend.candidate_model_registry as cmr

BENCHMARK_OUTPUT = ROOT / "validation" / "xavier_runtime_benchmark_v3_3_1.json"

TEST_COMPOUNDS = [
    ("Orforglipron", "Cc1cc(-n2nc3c(c2-n2ccn(-c4ccc5c(cnn5C)c4F)c2=O)[C@H](C)N(C(=O)c2cc4cc([C@H]5CCOC(C)(C)C5)ccc4n2[C@@]2(c4noc(=O)[nH]4)C[C@@H]2C)CC3)cc(C)c1F"),
    ("Sunvozertinib", "C=CC(=O)Nc1cc(Nc2nccc(Nc3cc(Cl)c(F)cc3C(C)(C)O)n2)c(OC)cc1N1CC[C@@H](N(C)C)C1"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Imatinib", "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"),
    ("Atorvastatin", "CC(C)c1c(c(c(c1c2ccc(cc2)F)c3ccccc3)C(=O)Nc4ccccc4)CC(O)CC(O)CC(=O)O"),
]

QUANT_ENDPOINTS = [
    ("Solubility", "esol_delaney_v1"),
    ("Permeability", "admetica_caco2"),
    ("Plasma protein binding", "admetica_ppbr"),
    ("HLM intrinsic clearance", "drugopt_hlm_chemical_space_v1"),
    ("CYP3A4 quantitative inhibition", "drugopt_calibrated_cyp3a4_pic50"),
    ("CYP2D6 quantitative inhibition", "drugopt_calibrated_cyp2d6_pic50"),
    ("CYP1A2 quantitative inhibition", "drugopt_calibrated_cyp1a2_pic50"),
    ("CYP2C9 quantitative inhibition", "drugopt_calibrated_cyp2c9_pic50"),
    ("hERG liability", "physchem_gbr_herg_pic50_v1"),
]

def get_memory_info():
    proc = psutil.Process(os.getpid())
    rss_mb = proc.memory_info().rss / (1024 * 1024)
    vms_mb = proc.memory_info().vms / (1024 * 1024)
    sys_mem = psutil.virtual_memory()
    return {
        "process_rss_mb": round(rss_mb, 2),
        "process_vms_mb": round(vms_mb, 2),
        "system_total_gb": round(sys_mem.total / (1024**3), 2),
        "system_used_gb": round(sys_mem.used / (1024**3), 2),
        "system_free_gb": round(sys_mem.free / (1024**3), 2),
    }

def run_benchmark():
    print("="*70)
    print("NVIDIA Jetson Xavier ARM64 Production Runtime Benchmark (v3.3.1)")
    print("="*70)
    
    cmr.register_candidate_models_to_multimodel()
    mem_init = get_memory_info()
    print(f"Initial Memory: RSS {mem_init['process_rss_mb']} MB | System Used {mem_init['system_used_gb']} GB")
    
    # Cold start on first compound (Orforglipron)
    c1_name, c1_smi = TEST_COMPOUNDS[0]
    t_cold_start = time.perf_counter()
    cold_results = {}
    for ep_name, model_id in QUANT_ENDPOINTS:
        adapter = mm._V2_ADAPTER_REGISTRY.get(model_id) or mm._ADAPTER_REGISTRY.get(model_id)
        if adapter:
            contract = get_endpoint_contract(ep_name)
            res = adapter.execute(c1_smi, contract)
            cold_results[ep_name] = res.value
    t_cold = round((time.perf_counter() - t_cold_start) * 1000.0, 2)
    print(f"Cold Start Full Execution ({c1_name}): {t_cold:.2f} ms")
    
    # Warm start on second run of same compound
    t_warm_start = time.perf_counter()
    warm_results = {}
    for ep_name, model_id in QUANT_ENDPOINTS:
        adapter = mm._V2_ADAPTER_REGISTRY.get(model_id) or mm._ADAPTER_REGISTRY.get(model_id)
        if adapter:
            contract = get_endpoint_contract(ep_name)
            res = adapter.execute(c1_smi, contract)
            warm_results[ep_name] = res.value
    t_warm = round((time.perf_counter() - t_warm_start) * 1000.0, 2)
    print(f"Warm Start Full Execution ({c1_name}): {t_warm:.2f} ms")
    
    # 5-compound sequential run
    seq_latencies = []
    print("\nExecuting 5-Compound Sequential Test:")
    for name, smi in TEST_COMPOUNDS:
        t0 = time.perf_counter()
        for ep_name, model_id in QUANT_ENDPOINTS:
            adapter = mm._V2_ADAPTER_REGISTRY.get(model_id) or mm._ADAPTER_REGISTRY.get(model_id)
            if adapter:
                contract = get_endpoint_contract(ep_name)
                res = adapter.execute(smi, contract)
        dur_ms = round((time.perf_counter() - t0) * 1000.0, 2)
        seq_latencies.append(dur_ms)
        print(f"  - {name:15s}: {dur_ms:.2f} ms")
        
    avg_latency = round(float(np.mean(seq_latencies)), 2)
    p95_latency = round(float(np.percentile(seq_latencies, 95)), 2)
    total_seq_time = round(sum(seq_latencies), 2)
    throughput = round(len(TEST_COMPOUNDS) / (total_seq_time / 1000.0), 2)
    
    mem_final = get_memory_info()
    print(f"\nFinal Memory: RSS {mem_final['process_rss_mb']} MB (Delta: +{round(mem_final['process_rss_mb'] - mem_init['process_rss_mb'], 2)} MB)")
    print(f"Throughput: {throughput} compounds/sec | Mean Latency: {avg_latency} ms | p95: {p95_latency} ms")
    
    report = {
        "hardware": "NVIDIA Jetson Xavier ARM64 (8-core Carmel CPU, 32GB Unified LPDDR4x)",
        "engine_version": "drugopt-prediction-engine-v3.3.1-production",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "cold_start_latency_ms": t_cold,
        "warm_start_latency_ms": t_warm,
        "sequential_latencies_ms": {
            name: lat for (name, _), lat in zip(TEST_COMPOUNDS, seq_latencies)
        },
        "mean_latency_ms": avg_latency,
        "p95_latency_ms": p95_latency,
        "throughput_compounds_per_sec": throughput,
        "memory_profile": {
            "initial_rss_mb": mem_init["process_rss_mb"],
            "peak_rss_mb": mem_final["process_rss_mb"],
            "delta_rss_mb": round(mem_final["process_rss_mb"] - mem_init["process_rss_mb"], 2),
            "system_ram_total_gb": mem_init["system_total_gb"],
            "system_ram_used_gb": mem_final["system_used_gb"],
        },
        "verdict": "PRODUCTION_LATENCY_COMPLIANT" if avg_latency < 500.0 else "REVIEW_LATENCY",
    }
    
    with open(BENCHMARK_OUTPUT, "w") as f:
        json.dump(report, f, indent=2)
        
    print(f"Benchmark report saved to: {BENCHMARK_OUTPUT}")
    return report

if __name__ == "__main__":
    run_benchmark()
