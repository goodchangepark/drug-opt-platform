import json
from rdkit import Chem

with open('backend/reference_drugs_150.json') as f:
    drugs_150 = json.load(f)

existing_db_ids = {d['drugbank_id'] for d in drugs_150 if 'drugbank_id' in d}
existing_names = {d['name'].lower() for d in drugs_150}
existing_cas = {d['cas_number'] for d in drugs_150 if d.get('cas_number')}
existing_inchikeys = set()
for d in drugs_150:
    m = Chem.MolFromSmiles(d['smiles'])
    if m:
        existing_inchikeys.add(Chem.MolToInchiKey(m))

candidates_test = [
    ("Omeprazole", "DB00338", "73590-58-6", "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=CC(=C3)OC"),
    ("Esomeprazole", "DB00736", "161796-78-7", "CC1=CN=C(C(=C1OC)C)C[S@@](=O)C2=NC3=C(N2)C=CC(=C3)OC"),
    ("Clopidogrel", "DB00758", "113665-84-2", "COC(=O)[C@H](C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3"),
    ("Voriconazole", "DB00582", "137234-62-9", "C[C@@H](C1=NC=NC=C1F)[C@](O)(CN2C=NC=N2)C3=C(F)C=C(F)C=C3"),
    ("Lansoprazole", "DB00448", "103577-45-3", "CC1=C(OCC(F)(F)F)C=CN=C1CS(=O)C2=NC3=CC=CC=C3N2"),
    ("Pantoprazole", "DB00213", "102625-70-7", "COC1=C(OC)C=CN=C1CS(=O)C2=NC3=C(N2)C=C(OC(F)F)C=C3"),
    ("Rabeprazole", "DB01129", "117976-89-3", "CC1=C(OCCCOC)C=CN=C1CS(=O)C2=NC3=CC=CC=C3N2"),
    ("Diazepam", "DB00829", "439-14-5", "CN1C(=O)CN=C(C2=CC=CC=C2)C3=C1C=CC(Cl)=C3"),
    ("Citalopram", "DB00215", "59729-33-8", "CN(C)CCCC1(C2=C(CO1)C=C(C#N)C=C2)C3=CC=C(F)C=C3"),
    ("Moclobemide", "DB01171", "71320-77-9", "ClC1=CC=C(C=C1)C(=O)NCCN2CCOCC2"),
    ("Digoxin", "DB00390", "20830-75-5", "CC1OC(CC(O)C1O)OC2C(O)CC(OC3C(O)CC(OC4CCC5(C)C(CCC6C5CCC7(C)C6CC(O)C7C8=CC(=O)OC8)C4)OC3C)OC2C"),
    ("Quinidine", "DB00908", "56-54-2", "COC1=CC2=C(C=C1)C(=NC=C2)[C@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C"),
    ("Verapamil", "DB00661", "52-53-9", "CC(C)C(CCCN(C)CCC1=CC(=C(C=C1)OC)OC)(C#N)C2=CC(=C(C=C2)OC)OC"),
    ("Loperamide", "DB00836", "53179-11-6", "CN(C)C(=O)C(CCN1CCC(O)(CC1)C2=CC=C(Cl)C=C2)(C3=CC=CC=C3)C4=CC=CC=C4"),
    ("Sulfasalazine", "DB00795", "599-79-1", "O=C(O)C1=CC(=CC=C1O)N=NC2=CC=C(C=C2)S(=O)(=O)NC3=NC=CC=C3"),
    ("Rosuvastatin", "DB01098", "287714-41-4", "CC(C)C1=NC(=NC(=C1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)C2=CC=C(F)C=C2)N(C)S(=O)(=O)C"),
    ("Topotecan", "DB01030", "123948-87-8", "CCN(C)CC1=C(O)C=C2C(=C1)N=C3C(=CC4=C(COC(=O)[C@@]4(O)CC)C3=O)C2"),
    ("Metoprolol", "DB00264", "37350-58-6", "COCCC1=CC=C(OCC(O)CNC(C)C)C=C1"),
    ("Propranolol", "DB00571", "525-66-6", "CC(C)NCC(O)COC1=CC=CC2=CC=CC=C12"),
    ("Atenolol", "DB00335", "29122-68-7", "CC(C)NCC(O)COC1=CC=C(CC(=O)N)C=C1"),
    ("Chloroquine", "DB00608", "54-05-7", "CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1"),
    ("Amiodarone", "DB01118", "1951-25-3", "CCCC1=NC2=C(C=CC=C2)C1C(=O)C3=CC(I)=C(OCCN(CC)CC)C(I)=C3"),
    ("Doxorubicin", "DB00997", "23214-92-8", "COC1=C2C(=O)C3=C(C(=O)C2=C(O)C=C1)C(O)=C(C[C@]3(O)C(=O)CO)O[C@H]4C[C@@H](N)[C@H](O)[C@@H](C)O4"),
    ("Haloperidol", "DB00502", "52-86-8", "OC1(CCN(CCCC(=O)C2=CC=C(F)C=C2)CC1)C3=CC=C(Cl)C=C3"),
    ("Donepezil", "DB00843", "120014-06-4", "COC1=C(OC)C=C2C(=C1)CC(CC3CCN(CC4=CC=CC=C4)CC3)C2=O"),
    ("Memantine", "DB00449", "19982-08-2", "CC12CC3CC(C)(C1)CC(N)(C3)C2"),
    ("Risperidone", "DB00734", "106266-06-2", "CC1=C(CCN2CCC(CC2)C3=NOC4=C3C=CC(F)=C4)C(=O)N5CCCCC5=N1"),
    ("Bupropion", "DB01156", "34841-39-9", "CC(NC(C)(C)C)C(=O)C1=CC(=CC=C1)Cl"),
    ("Diltiazem", "DB00343", "42399-41-7", "CC(=O)OC1C(SC2=CC=CC=C2N(CCN(C)C)C1=O)C3=CC=C(OC)C=C3"),
    ("Ketoconazole", "DB01026", "65277-42-1", "CC(=O)N1CCN(CC1)C2=CC=C(OCC3COC(O3)(CN4C=CN=C4)C5=C(Cl)C=C(Cl)C=C5)C=C2"),
    ("Itraconazole", "DB01167", "84625-61-6", "CCC(C)N1N=CN(C1=O)C2=CC=C(N3CCN(CC3)C4=CC=C(OCC5COC(O5)(CN6C=NC=N6)C7=C(Cl)C=C(Cl)C=C7)C=C4)C=C2"),
    ("Fluconazole", "DB00196", "86386-73-4", "OC(CN1C=NC=N1)(CN2C=NC=N2)C3=C(F)C=C(F)C=C3"),
    ("Carbamazepine", "DB00564", "298-46-4", "NC(=O)N1C2=CC=CC=C2C=CC3=CC=CC=C13"),
    ("Phenytoin", "DB00252", "57-41-0", "O=C1NC(=O)C(N1)(C2=CC=CC=C2)C3=CC=CC=C3"),
    ("Theophylline", "DB00277", "58-55-9", "CN1C(=O)C2=C(N=CN2)N(C)C1=O"),
    ("Caffeine", "DB00201", "58-08-2", "CN1C(=O)N(C)C(=O)C2=C1N=CN2C"),
    ("Nicotine", "DB00184", "54-11-5", "CN1CCC[C@H]1C2=CN=CC=C2"),
    ("Midazolam", "DB00683", "59467-70-8", "CC1=NC=C2N1C3=C(C=C(Cl)C=C3)C(=NC2)C4=CC=CC=C4F"),
    ("Alprazolam", "DB00404", "28981-97-7", "CC1=NN=C2CN=C(C3=CC=CC=C3)C4=C(C=CC(Cl)=C4)N12"),
    ("Zolpidem", "DB00422", "82626-48-0", "CC1=CC=C(C=C1)C2=C(CC(=O)N(C)C)C3=NC=C(C)N3C2=O"),
    ("Sildenafil", "DB00203", "139755-83-2", "CCCC1=NN(C)C2=C1N=C(NC2=O)C3=C(OCC)C=CC(=C3)S(=O)(=O)N4CCN(C)CC4"),
    ("Tadalafil", "DB00820", "171596-29-5", "CN1CCN2C(C1=O)CC3=C(NC4=CC=CC=C34)[C@H]2C5=CC6=C(OCO6)C=C5"),
    ("Vardenafil", "DB00862", "224785-90-4", "CCCC1=NC(=C2N1N=C(NC2=O)C3=C(OCC)C=CC(=C3)S(=O)(=O)N4CCN(CC)CC4)C"),
    ("Montelukast", "DB00471", "158966-92-8", "CC(C)(C1=CC=CC=C1CCC[C@@H](SCC2(CC2)CC(=O)O)/C=C/C3=NC4=CC=C(Cl)C=C4C=C3)O"),
    ("Cetirizine", "DB00373", "83881-51-0", "OC(=O)COCCN1CCN(CC1)C(C2=CC=CC=C2)C3=CC=C(Cl)C=C3"),
    ("Fexofenadine", "DB00950", "83799-24-0", "CC(C)(C(=O)O)C1=CC=C(C=C1)C(O)CCCN2CCC(CC2)C(O)(C3=CC=CC=C3)C4=CC=CC=C4"),
    ("Loratadine", "DB00455", "79794-75-5", "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(Cl)C=C3)CC1"),
    ("Ranitidine", "DB00863", "66357-35-5", "CNC(=C[N+](=O)[O-])NCCSCC1=CC=C(O1)CN(C)C"),
    ("Famotidine", "DB00927", "76824-35-6", "NS(=O)(=O)N=C(N)CSCC1=CSC(=N1)N=C(N)N"),
    ("Cimetidine", "DB00501", "51481-61-9", "CNC(=NC#N)NCCSCC1=C(C)N=CN1")
]

print(f'Checking {len(candidates_test)} candidates for collisions...')
collisions = 0
for name, db_id, cas, smi in candidates_test:
    m = Chem.MolFromSmiles(smi)
    if not m:
        print(f'ERROR: Invalid SMILES for {name}')
        collisions += 1
        continue
    ik = Chem.MolToInchiKey(m)
    if db_id in existing_db_ids:
        print(f'COLLISION: DB ID {db_id} ({name}) already exists!')
        collisions += 1
    if name.lower() in existing_names:
        print(f'COLLISION: Name {name} already exists!')
        collisions += 1
    if cas in existing_cas:
        print(f'COLLISION: CAS {cas} ({name}) already exists!')
        collisions += 1
    if ik in existing_inchikeys:
        print(f'COLLISION: InChIKey {ik} ({name}) already exists!')
        collisions += 1

print(f'Total collisions detected: {collisions}')
if collisions == 0:
    print('SUCCESS: All 50 candidates are 100% DISJOINT and unique!')
