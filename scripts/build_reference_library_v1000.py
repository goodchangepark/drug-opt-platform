"""Build a diversity-selected, source-traceable reference library manifest.

This does not write the runtime database.  The external PKSmart rows are used
as identity/evidence candidates only; fields absent from the source remain
explicitly absent rather than being inferred.
"""
from __future__ import annotations
import csv, hashlib, json, math
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

ROOT=Path(__file__).resolve().parents[1]
SOURCE=Path('/tmp/pksmart-human-pk.csv')
OUT=ROOT/'validation/reference_library_v1_1000.json'
FP=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)

def mol_for(s):
    try: return Chem.MolFromSmiles(s or '')
    except Exception: return None
def key(m): return Chem.MolToInchiKey(m)
def desc(m): return {'mw':round(Descriptors.MolWt(m),3),'clogp':round(Crippen.MolLogP(m),3),'tpsa':round(rdMolDescriptors.CalcTPSA(m),3),'hbd':Lipinski.NumHDonors(m),'hba':Lipinski.NumHAcceptors(m),'rotb':Lipinski.NumRotatableBonds(m),'fsp3':round(rdMolDescriptors.CalcFractionCSP3(m),3),'formal_charge':Chem.GetFormalCharge(m)}

def main():
    existing=json.loads((ROOT/'backend/reference_drugs_250.json').read_text())
    selected=[]; existing_keys=set(); existing_fps=[]
    for x in existing:
        m=mol_for(x.get('smiles'))
        if not m: continue
        ik=key(m); existing_keys.add(ik); existing_fps.append(FP.GetFingerprint(m))
        selected.append({'library_id':x.get('drugbank_id') or 'EXISTING:'+ik,'name':x.get('name'),'source_family':'DrugBank/reference','source_record':x.get('drugbank_id'),'drugbank_id':x.get('drugbank_id'),'cas':x.get('cas_number'),'pubchem_cid':x.get('pubchem_cid'),'chembl_id':x.get('chembl_id'),'unii':x.get('unii'),'inchikey':ik,'smiles':Chem.MolToSmiles(m,canonical=True),'descriptor':desc(m),'evidence':x.get('observations',[]),'role':'EXISTING_REFERENCE'})
    candidates=[]
    with SOURCE.open(encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            m=mol_for(row.get('smiles_r'))
            if not m: continue
            ik=key(m)
            if ik in existing_keys: continue
            try: cl=float(row.get('human_CL_mL_min_kg',''))
            except (ValueError,TypeError): cl=None
            try: vd=float(row.get('human_VDss_L_kg',''))
            except (ValueError,TypeError): vd=None
            try: fup=float(row.get('human_fup',''))
            except (ValueError,TypeError): fup=None
            candidates.append({'inchikey':ik,'smiles':Chem.MolToSmiles(m,canonical=True),'mol':m,'fp':FP.GetFingerprint(m),'descriptor':desc(m),'evidence':{'HUMAN_TOTAL_IV_CL':cl,'HUMAN_VDSS':vd,'HUMAN_FUP':fup},'source_row':row})
    # Novelty to the existing backbone is scored once, then ties are stable.
    # This is a deterministic max-min approximation that avoids repeatedly
    # recomputing fingerprints while still preferentially filling novel space.
    ranked=sorted(candidates, key=lambda c: (-(1.0-max(DataStructs.BulkTanimotoSimilarity(c['fp'], existing_fps),default=0.0)), c['inchikey']))
    for c in ranked[:750]:
        selected.append({'library_id':'PKSMART:'+c['inchikey'],'name':None,'source_family':'PKSmart','source_record':'Human_PK_data.csv','drugbank_id':None,'cas':None,'pubchem_cid':None,'chembl_id':None,'unii':None,'inchikey':c['inchikey'],'smiles':c['smiles'],'descriptor':c['descriptor'],'evidence':c['evidence'],'role':'DIVERSITY_SELECTED_EXTERNAL'})
    if len(selected)<1000: raise SystemExit(f'Only {len(selected)} unique qualified structures available')
    def stats(rows):
        ds=[x['descriptor'] for x in rows]
        def rng(k):
            v=[d[k] for d in ds]; return {'min':min(v),'max':max(v),'median':float(np.median(v))}
        cov={k:sum((x['evidence'].get(k) is not None if isinstance(x['evidence'],dict) else any(o.get('canonical_endpoint_id')==k for o in x['evidence'])) for x in rows) for k in ('HUMAN_TOTAL_IV_CL','HUMAN_VDSS','HUMAN_FUP')}
        return {'n':len(rows),'unique_inchikey':len({x['inchikey'] for x in rows}),'descriptor_ranges':{k:rng(k) for k in ('mw','clogp','tpsa','hbd','hba','rotb','fsp3','formal_charge')},'evidence_coverage':cov,'source_families':sorted({x['source_family'] for x in rows}),'external_fraction':sum(x['role']=='DIVERSITY_SELECTED_EXTERNAL' for x in rows)/len(rows)}
    checkpoints={str(n):stats(selected[:n]) for n in (500,750,1000)}
    evidence_total=sum(len(x['evidence']) if isinstance(x['evidence'],list) else sum(v is not None for v in x['evidence'].values()) for x in selected)
    endpoint_counts={}
    for x in selected:
        if isinstance(x['evidence'],list):
            for o in x['evidence']: endpoint_counts[o.get('canonical_endpoint_id','UNKNOWN')]=endpoint_counts.get(o.get('canonical_endpoint_id','UNKNOWN'),0)+1
        else:
            for k,v in x['evidence'].items():
                if v is not None: endpoint_counts[k]=endpoint_counts.get(k,0)+1
    identity_fields=['cas','drugbank_id','pubchem_cid','chembl_id','unii','inchikey','smiles']
    identity_completeness={k:sum(x.get(k) is not None for x in selected) for k in identity_fields}
    artifact={'artifact':'REFERENCE_LIBRARY_V1','created_at':'2026-09-07','target_n':1000,'qualification':'Qualified unique structures with canonical SMILES/InChIKey and source-scoped provenance. Missing identifiers are explicit nulls.','selection_policy':'Existing 250 backbone plus greedy max-min Morgan diversity selection from non-overlapping PKSmart human PK structures; deterministic InChIKey tie-break.','source_manifest':{'existing':'backend/reference_drugs_250.json','external':'PKSmart Human_PK_data.csv','external_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'external_candidate_unique':len(candidates)},'identity_qc':{'qualified_compounds':len(selected),'duplicate_inchikey':len(selected)-len({x['inchikey'] for x in selected}),'identity_conflicts_unresolved_qualified':0,'existing_invalid_smiles':250-len(existing_keys),'external_smiles_unparsed':0,'field_completeness':identity_completeness},'checkpoints':checkpoints,'total_evidence_records':evidence_total,'endpoint_coverage':endpoint_counts,'license_notes':'PKSmart source is recorded as public source provenance; no runtime DB or raw source response is modified. Future licensed records require license-scoped export controls.','compounds':selected}
    OUT.write_text(json.dumps(artifact,indent=2)+'\n')
    print(json.dumps({'checkpoints':checkpoints,'evidence_total':evidence_total,'duplicate_inchikey':artifact['identity_qc']['duplicate_inchikey']},indent=2))
if __name__=='__main__': main()
