import base64
import io
import math
import pickle
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors, rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold, cross_val_predict

from .chemistry import ENGINE_VERSION, parse_smiles

FINGERPRINT_CONFIG = {"radius": 2, "bit_size": 2048}
DESCRIPTOR_NAMES = ["MW", "ALogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "AromaticRings"]
SEED = 42


def normalize_concentration(value: float, unit: str):
    factors = {"nM": 1.0, "µM": 1000.0, "uM": 1000.0, "mM": 1_000_000.0}
    if unit not in factors:
        raise ValueError(f"Unsupported concentration unit: {unit}")
    return float(value) * factors[unit], {"raw_value": value, "raw_unit": unit,
                                            "normalized_unit": "nM", "factor_to_nM": factors[unit]}


def pactivity(normalized_nm: float):
    if normalized_nm <= 0:
        raise ValueError("Activity value must be positive")
    molar = normalized_nm * 1e-9
    return -math.log10(molar)


def value_from_pactivity(p_value: float):
    return (10 ** (-p_value)) / 1e-9


def fingerprint_and_descriptors(smiles: str):
    mol = parse_smiles(smiles)
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=FINGERPRINT_CONFIG["radius"], fpSize=FINGERPRINT_CONFIG["bit_size"]
    )
    fingerprint = list(generator.GetFingerprint(mol).GetOnBits())
    descriptors = {
        "MW": Descriptors.MolWt(mol), "ALogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol), "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol), "RotB": Lipinski.NumRotatableBonds(mol),
        "Fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "AromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }
    murcko = Chem.MurckoDecompose(mol)
    scaffold = Chem.MolToSmiles(murcko) if murcko is not None else ""
    return mol, fingerprint, descriptors, scaffold


def feature_vector(fingerprint, descriptors):
    bits = np.zeros(FINGERPRINT_CONFIG["bit_size"], dtype=float)
    for bit in fingerprint: bits[bit] = 1.0
    descriptor_values = np.array([descriptors[name] for name in DESCRIPTOR_NAMES], dtype=float)
    descriptor_values[0] /= 500; descriptor_values[2] /= 150
    return np.concatenate([bits, descriptor_values])


def tanimoto_similarity(fp_a, fp_b):
    a, b = set(fp_a), set(fp_b)
    union = len(a | b)
    return len(a & b) / union if union else 1.0


def _metrics(y_true, y_pred):
    n = len(y_true)
    rho = None
    if n >= 3 and np.std(y_true) > 0 and np.std(y_pred) > 0:
        rho = float(spearmanr(y_true, y_pred).statistic)
    return {
        "N": int(n), "MAE": round(float(mean_absolute_error(y_true, y_pred)), 4),
        "RMSE": round(float(math.sqrt(mean_squared_error(y_true, y_pred))), 4),
        "R2": round(float(r2_score(y_true, y_pred)), 4), "Spearman": round(rho, 4) if rho is not None else None,
    }


def select_model(n_samples: int):
    return n_samples


def train_model(dataset):
    x = np.vstack(dataset["features"]); y=np.array(dataset["targets"], dtype=float); groups=np.array(dataset["scaffolds"])
    n=len(y)
    if n < 15:
        raise ValueError("Insufficient data for formal QSAR training")
    if n >= 30: candidates=[("Ridge baseline",Ridge(alpha=10.0)),("RandomForest",RandomForestRegressor(n_estimators=250,random_state=SEED,n_jobs=-1)),("ExtraTrees",ExtraTreesRegressor(n_estimators=250,random_state=SEED,n_jobs=-1))]; kfold_splits=5
    else: candidates=[("Ridge baseline",Ridge(alpha=5.0)),("RandomForest",RandomForestRegressor(n_estimators=120,random_state=SEED,n_jobs=-1))]; kfold_splits=4
    kf=KFold(n_splits=kfold_splits,shuffle=True,random_state=SEED)
    comparisons={}
    for name,model in candidates:
        pred=cross_val_predict(model,x,y,cv=kf)
        comparisons[name]=_metrics(y,pred)
    best_name=max(comparisons,key=lambda key:(comparisons[key]["R2"] if comparisons[key]["R2"] is not None else -999))
    best_model=next(model for name,model in candidates if name==best_name)
    unique_groups=len(set(groups))
    scaffold_metrics=None
    if unique_groups>=4:
        gkf=GroupKFold(n_splits=min(4,unique_groups)); scaffold_pred=cross_val_predict(best_model,x,y,cv=gkf,groups=groups)
        scaffold_metrics=_metrics(y,scaffold_pred)
    best_model.fit(x,y)
    buffer=io.BytesIO(); pickle.dump({"model":best_model,"name":best_name},buffer)
    encoded=base64.b64encode(buffer.getvalue()).decode()
    metrics={"random_cv":comparisons[best_name],"all_models_random_cv":comparisons,
             "scaffold_cv":scaffold_metrics,"selection":"lowest RMSE among validated models"}
    reason=f"Selected {best_name} by validation performance; N={n}; policy threshold={'N>=30' if n>=30 else '15<=N<30'}"
    return encoded,best_name,metrics,reason,n


def nearest_neighbors(target_fp, dataset, top_k=5):
    rows=[]
    for i,row in enumerate(dataset["rows"]):
        sim=tanimoto_similarity(target_fp,dataset["fingerprints"][i]); rows.append({**row,"similarity":round(sim,3)})
    return sorted(rows,key=lambda row:row["similarity"],reverse=True)[:top_k]


def applicability(similarities, descriptors, dataset):
    max_sim=max([row["similarity"] for row in similarities],default=0)
    known=np.array(dataset["descriptors"])
    target=np.array([descriptors[name] for name in DESCRIPTOR_NAMES])
    if known.shape[0]==0:
        return "OUT OF DOMAIN","INSUFFICIENT DATA",max_sim,True
    ranges=np.ptp(known,axis=0)
    outside=bool(np.any((target < known.min(axis=0)-.05*ranges) | (target > known.max(axis=0)+.05*ranges)))
    if max_sim>=.55 and not outside: domain="IN DOMAIN"
    elif max_sim>=.35: domain="BORDERLINE"
    else: domain="OUT OF DOMAIN"
    confidence="INSUFFICIENT DATA" if not similarities else ("HIGH" if max_sim>=.7 and domain=="IN DOMAIN" else ("MEDIUM" if max_sim>=.55 and domain!="OUT OF DOMAIN" else "LOW"))
    return domain,confidence,max_sim,outside
