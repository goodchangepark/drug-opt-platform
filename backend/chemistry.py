import hashlib
from collections import OrderedDict

from rdkit import Chem, RDLogger
from io import BytesIO

from PIL import Image
from rdkit.Chem import Crippen, Descriptors, Draw, FilterCatalog, Lipinski, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ENGINE = "RDKit"
try:
    from rdkit import __version__ as ENGINE_VERSION
except ImportError:
    ENGINE_VERSION = "unknown"


class ChemistryError(ValueError):
    pass


_PAINS_FACTORY = FilterCatalog.FilterCatalogParams()
_PAINS_FACTORY.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_BRENS_FACTORY = FilterCatalog.FilterCatalogParams()
_BRENS_FACTORY.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)

_ALERT_SETS = [
    ("PAINS", FilterCatalog.FilterCatalog(_PAINS_FACTORY), "Frequent-hitter/reactive PAINS motif; assay interference risk."),
    ("Brenk", FilterCatalog.FilterCatalog(_BRENS_FACTORY), "Undesirable/reactive medicinal-chemistry motif from Brenk et al."),
]

_SMARTS_ALERTS = [
    ("Reactive acyl halide", "[C,Cl;X4](=O)Cl", "Acyl halides are highly electrophilic and generally unsuitable leads."),
    ("Aldehyde", "[CX3H1](=O)[#6]", "Aldehydes can form covalent adducts; selectivity/toxicity concern."),
    ("Michael acceptor", "[C]=[C][C,S]=O", "Electrophilic Michael acceptor may react with biological nucleophiles."),
    ("Epoxide", "[OX2r3]1[#6][#6]1", "Strained epoxides are potentially genotoxic/electrophilic."),
    ("Anhydride", "[CX3](=[OX1])-[OX2]-[CX3]=[OX1]", "Reactive hydrolysis-prone electrophile."),
    ("Nitro aromatic", "[c][$([N+](=O)[O-])]", "Aromatic nitro groups can be associated with mutagenicity/reduction toxicity."),
]


def parse_smiles(smiles: str):
    if not smiles or not smiles.strip():
        raise ChemistryError("SMILES is empty")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ChemistryError("Invalid SMILES: RDKit sanitization failed")
    return mol


def _round(value, digits=3):
    if value is None or isinstance(value, str):
        return value
    return round(float(value), digits)


def calculate_properties(mol) -> dict:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    descriptors = {
        "molecular_formula": formula,
        "molecular_weight": _round(Descriptors.MolWt(mol), 2),
        "exact_molecular_weight": _round(Descriptors.ExactMolWt(mol), 4),
        "clogp": _round(Crippen.MolLogP(mol)),
        "tpsa": _round(rdMolDescriptors.CalcTPSA(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
        "heteroatom_count": int(rdMolDescriptors.CalcNumHeteroatoms(mol)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(mol)),
        "aromatic_ring_count": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fraction_csp3": _round(rdMolDescriptors.CalcFractionCSP3(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
        "molar_refractivity": _round(Crippen.MolMR(mol), 2),
        "aromatic_proportion": _round(sum(atom.GetIsAromatic() for atom in mol.GetAtoms()) / max(mol.GetNumHeavyAtoms(), 1), 3),
        "molecular_flexibility": _round(int(Lipinski.NumRotatableBonds(mol)) / max(mol.GetNumHeavyAtoms(), 1), 3),
    }
    qed_props = QED.properties(mol)._asdict()
    descriptors["qed"] = _round(QED.qed(mol))
    descriptors["qed_components"] = {name: _round(value) for name, value in qed_props.items()}
    return descriptors


def drug_likeness(props: dict):
    failures = {
        "Lipinski Rule of Five": [],
        "Veber": [],
        "Ghose": [],
        "Egan": [],
        "Rule of Three": [],
    }
    mw, logp, tpsa = props["molecular_weight"], props["clogp"], props["tpsa"]
    hbd, hba, rotb = props["hbd"], props["hba"], props["rotatable_bonds"]

    if mw > 500: failures["Lipinski Rule of Five"].append("MW > 500")
    if logp > 5: failures["Lipinski Rule of Five"].append("cLogP > 5")
    if hbd > 5: failures["Lipinski Rule of Five"].append("HBD > 5")
    if hba > 10: failures["Lipinski Rule of Five"].append("HBA > 10")

    if tpsa > 140: failures["Veber"].append("TPSA > 140 Å²")
    if rotb > 10: failures["Veber"].append("rotatable bonds > 10")

    if not (160 <= mw <= 480): failures["Ghose"].append("MW outside 160–480")
    if not (-0.4 <= logp <= 5.6): failures["Ghose"].append("cLogP outside -0.4 to 5.6")
    if not (40 <= props["molar_refractivity"] <= 130): failures["Ghose"].append("molar refractivity outside 40–130")
    if not (20 <= props["heavy_atom_count"] <= 70): failures["Ghose"].append("heavy atoms outside 20–70")

    if logp > 5.88: failures["Egan"].append("cLogP > 5.88")
    if tpsa > 131.6: failures["Egan"].append("TPSA > 131.6 Å²")

    if mw > 300: failures["Rule of Three"].append("MW > 300")
    if logp > 3: failures["Rule of Three"].append("cLogP > 3")
    if hbd > 3: failures["Rule of Three"].append("HBD > 3")
    if hba > 3: failures["Rule of Three"].append("HBA > 3")
    if rotb > 3: failures["Rule of Three"].append("rotatable bonds > 3")

    rules = {}
    for name, reasons in failures.items():
        rules[name] = {"result": "PASS" if not reasons else "FAIL", "reasons": reasons}
    return rules


def structural_alerts(mol):
    alerts = []
    for set_name, catalog, default_reason in _ALERT_SETS:
        entries = catalog.GetMatches(mol)
        seen = set()
        for entry in entries:
            name = entry.GetDescription()
            key = (set_name, name)
            if key in seen:
                continue
            seen.add(key)
            atoms = []
            for filter_match in entry.GetFilterMatches(mol):
                for pair in filter_match.atomPairs:
                    try:
                        atoms.append(int(pair.target))
                    except (TypeError, AttributeError):
                        pass
            alerts.append({
                "alert_set": set_name,
                "alert_name": name,
                "matched_atoms": sorted(set(atoms)),
                "matched_smiles": "",
                "reason": default_reason,
            })
    for name, pattern, reason in _SMARTS_ALERTS:
        patt = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(patt)
        flat = sorted({atom_id for match in matches for atom_id in match})
        if flat:
            submol = Chem.PathToSubmol(mol, [bond.GetIdx() for bond in mol.GetBonds()
                                             for pair in matches if bond.GetBeginAtomIdx() in pair and bond.GetEndAtomIdx() in pair])
            alerts.append({
                "alert_set": "Functional group",
                "alert_name": name,
                "matched_atoms": flat,
                "matched_smiles": Chem.MolToSmiles(submol),
                "reason": reason,
            })
    return alerts


def assessment(props: dict, alerts: list, rules: dict):
    strengths, concerns = [], []
    mw, logp, tpsa, fsp3 = props["molecular_weight"], props["clogp"], props["tpsa"], props["fraction_csp3"]
    qed = props["qed"]
    if 250 <= mw <= 500: strengths.append("acceptable molecular weight range")
    elif mw < 200: concerns.append("very low molecular weight may limit target engagement")
    else: concerns.append("high molecular weight")
    if 0 <= logp <= 4: strengths.append("reasonable lipophilicity")
    elif logp > 5: concerns.append("high lipophilicity; solubility and clearance risk")
    else: concerns.append("low lipophilicity; permeability risk possible")
    if tpsa <= 120: strengths.append("acceptable TPSA for passive permeability")
    else: concerns.append("excessive polar surface area")
    if fsp3 >= 0.35: strengths.append("good Fsp3 relative to common oral-drug space")
    else: concerns.append("low Fsp3/high aromatic proportion; poor-solubility risk")
    if props["aromatic_proportion"] >= 0.65: concerns.append("excessive aromaticity")
    if props["rotatable_bonds"] >= 10: concerns.append("excessive molecular flexibility")
    if alerts: concerns.append(f"{len(alerts)} structural alert(s)")
    if qed >= 0.55: strengths.append(f"good QED ({qed:.2f})")
    elif qed <= 0.30: concerns.append(f"low QED ({qed:.2f})")
    failed = [name for name, result in rules.items() if result["result"] == "FAIL"]
    if failed: concerns.append("drug-likeness rule failure: " + ", ".join(failed))
    return {"strengths": sorted(set(strengths)), "concerns": sorted(set(concerns))}


def structure_images(mol, highlight_atoms=None):
    highlight_atoms = highlight_atoms or []
    Draw.rdDepictor.Compute2DCoords(mol)
    normal = Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(420, 320), useSVG=True)
    png = Draw.MolToImage(mol, size=(420, 320), highlightAtoms=highlight_atoms)
    buffer = BytesIO(); png.save(buffer, format="PNG"); buffer.seek(0)
    image = Image.open(buffer); svg_buffer = BytesIO()
    image.save(svg_buffer, format="PNG")
    # Return a data URI for highlighted output so the UI always renders true highlighting
    # across RDKit builds whose Python SVG API omits atom-highlight controls.
    import base64
    data_uri = "data:image/png;base64," + base64.b64encode(svg_buffer.getvalue()).decode("ascii")
    return str(normal), data_uri


def analyze_smiles(smiles: str) -> dict:
    mol = parse_smiles(smiles)
    canonical = Chem.MolToSmiles(mol)
    props = calculate_properties(mol)
    rules = drug_likeness(props)
    alerts = structural_alerts(mol)
    highlight = sorted({atom_id for alert in alerts for atom_id in alert.get("matched_atoms", [])})
    svg, highlighted_svg = structure_images(mol, highlight)
    identity = {
        "canonical_smiles": canonical,
        "isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
    }
    try:
        identity["inchi"] = Chem.inchi.MolToInchi(mol)
        identity["inchikey"] = Chem.inchi.InchiToInchiKey(identity["inchi"])
    except Exception as exc:
        identity["inchi"], identity["inchikey"] = "", ""
        identity["inchi_error"] = str(exc)
    inputs_hash = hashlib.sha256(canonical.encode()).hexdigest()
    provenance = {"type": "Calculated", "engine": ENGINE, "engine_version": ENGINE_VERSION,
                  "methods": ["Crippen cLogP/MR", "Ertl TPSA", "RDKit QED", "RDKit FilterCatalog"],
                  "calculated_at_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat()}
    return {
        "identity": identity, "properties": props, "rules": rules, "alerts": alerts,
        "assessment": assessment(props, alerts, rules), "svg": svg, "highlighted_svg": highlighted_svg,
        "provenance": provenance, "inputs_hash": inputs_hash,
    }
