"""Canonical Structure Standardization & Reproducibility Engine (Stage 4C-1).

Standardizer Pipeline: CHEM_STANDARDIZER_V1
1. Parse & Sanitize
2. Normalize Functional Groups (rdMolStandardize.Normalizer)
3. Fragment & Salt Handling (rdMolStandardize.MetalDisconnector & Salt/Parent Extractor)
4. Charge-Parent Handling (rdMolStandardize.Uncharger for neutralizable acids/bases)
5. Isotope Policy (Preserve isotopic labels [2H], [13C], etc.)
6. Tautomer Canonicalization (rdMolStandardize.TautomerEnumerator)
7. Stereochemistry Preservation (AssignStereochemistry with defined E/Z & tetrahedral preservation)
8. Export Canonical Isomeric SMILES, InChI, and InChIKey
"""

from __future__ import annotations

import logging
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")
logger = logging.getLogger(__name__)

STANDARDIZER_NAME = "CHEM_STANDARDIZER_V1"
STANDARDIZER_VERSION = "1.0"
RDKIT_VERSION = str(getattr(Chem, "__version__", getattr(Chem, "rdkitVersion", "2025.03.1")))

# Global Versioned Fingerprint Configuration
GLOBAL_FINGERPRINT_CONFIG = {
    "name": "Morgan",
    "radius": 2,
    "nBits": 2048,
    "useChirality": True,
    "useFeatures": False,
    "rdkit_version": RDKIT_VERSION,
    "disclaimer": "Global canonical fingerprint setting. Model-specific overrides apply if trained representation is fixed.",
}

# Global Versioned Descriptor Configuration
GLOBAL_DESCRIPTOR_CONFIG = {
    "version": "1.0",
    "mw_definition": "Descriptors.MolWt (Average Molecular Weight)",
    "exact_mw_definition": "Descriptors.ExactMolWt (Monoisotopic Mass)",
    "clogp_definition": "Crippen.MolLogP (Crippen SlogP implementation)",
    "tpsa_definition": "rdMolDescriptors.CalcTPSA (Ertl et al. 2000 polar surface area)",
    "hbd_definition": "Lipinski.NumHDonors (NH and OH count)",
    "hba_definition": "Lipinski.NumHAcceptors (N and O count)",
    "rotatable_bonds_definition": "Lipinski.NumRotatableBonds",
    "aromaticity_model": "RDKit Default Aromaticity Model",
    "fraction_csp3_definition": "rdMolDescriptors.CalcFractionCSP3",
    "disclaimer": "Calculated descriptors are in silico molecular properties, not experimental laboratory measurements.",
}

# Common inorganic and organic counterions / salts / solvents
RECOGNIZED_SALTS = {
    "[Cl-]", "[Br-]", "[I-]", "[F-]", "[Na+]", "[K+]", "[Ca+2]", "[Mg+2]", "[Li+]",
    "[NH4+]", "[Zn+2]", "O", "CC(=O)O", "OS(=O)(=O)O", "O=C(O)C(=O)O", "C(=O)(C(=O)O)O",
    "O=P(O)(O)O", "CS(=O)(=O)O", "Cc1ccc(S(=O)(=O)O)cc1", "O=C(O)c1ccccc1O",
    "OC(C(=O)O)C(O)C(=O)O", "OC(C(=O)O)CC(=O)O", "O=C(O)C=CC(=O)O", "Cl", "Br",
}

# Standardizer Components Initialization
_normalizer = rdMolStandardize.Normalizer()
_metal_disconnector = rdMolStandardize.MetalDisconnector()
_uncharger = rdMolStandardize.Uncharger()
_tautomer_enumerator = rdMolStandardize.TautomerEnumerator()


class StandardizationError(ValueError):
    pass


def standardize_molecule(smiles_input: str) -> dict[str, Any]:
    """Execute canonical CHEM_STANDARDIZER_V1 pipeline on input SMILES.

    Returns structured record containing canonical SMILES, isomeric SMILES, InChIKey,
    salt/fragment metadata, stereochemistry flags, and provenance.
    """
    if not smiles_input or not isinstance(smiles_input, str) or not smiles_input.strip():
        return {
            "status": "INVALID_INPUT",
            "original_smiles": smiles_input,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "inchikey": None,
            "warnings": ["Input SMILES is empty or invalid."],
            "provenance": _get_provenance(),
        }

    raw_smiles = smiles_input.strip()
    warnings = []

    # 1. Parse & Initial Sanitize
    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
        return {
            "status": "PARSING_FAILED",
            "original_smiles": raw_smiles,
            "canonical_smiles": None,
            "isomeric_smiles": None,
            "inchikey": None,
            "warnings": ["RDKit MolFromSmiles failed to parse structure."],
            "provenance": _get_provenance(),
        }

    # Record initial stereochemistry features
    has_defined_chiral = any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms())
    has_defined_ez = any(bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds())
    has_isotopes = any(atom.GetIsotope() != 0 for atom in mol.GetAtoms())
    if has_isotopes:
        warnings.append("ISOTOPIC_LABEL_PRESERVED: Structure contains isotopic labels (e.g. 2H, 13C).")

    # 2. Normalize Functional Groups
    try:
        mol = _normalizer.normalize(mol)
    except Exception as err:
        warnings.append(f"Normalizer warning: {err}")

    # 3. Disconnect Metals
    try:
        mol = _metal_disconnector.Disconnect(mol)
    except Exception as err:
        warnings.append(f"Metal disconnector warning: {err}")

    # 4. Fragment & Salt Handling
    frags = Chem.GetMolFrags(mol, asMols=True)
    status = "SUCCESS"
    salt_extracted = False
    parent_mol = mol

    if len(frags) > 1:
        frag_smiles_list = [Chem.MolToSmiles(f) for f in frags]
        organic_frags = []
        salt_frags = []

        for f, f_smi in zip(frags, frag_smiles_list):
            if f_smi in RECOGNIZED_SALTS or f.GetNumHeavyAtoms() <= 2:
                salt_frags.append(f_smi)
            else:
                organic_frags.append((f.GetNumHeavyAtoms(), f, f_smi))

        if salt_frags and organic_frags:
            # Recognized salt / solvent stripped
            organic_frags.sort(key=lambda x: x[0], reverse=True)
            parent_mol = organic_frags[0][1]
            salt_extracted = True
            warnings.append(f"SALT_REMOVED: Extracted main organic parent. Removed counterion/solvent: {', '.join(salt_frags)}")
        elif len(organic_frags) > 1:
            # Ambiguous multi-component organic mixture
            status = "MULTICOMPONENT_REVIEW_REQUIRED"
            warnings.append(f"MULTICOMPONENT_REVIEW_REQUIRED: Input contains multiple organic fragments ({', '.join(f[2] for f in organic_frags)}). No arbitrary parent selection executed.")
            # Select largest fragment for preview only, but flag status
            organic_frags.sort(key=lambda x: x[0], reverse=True)
            parent_mol = organic_frags[0][1]

    # 5. Charge-Parent Handling (Neutralize acids/bases, preserve quaternary ammonium)
    try:
        uncharged_mol = _uncharger.uncharge(parent_mol)
        parent_mol = uncharged_mol
    except Exception as err:
        warnings.append(f"Uncharger warning: {err}")

    # 6. Tautomer Canonicalization
    try:
        tautomer_mol = _tautomer_enumerator.Canonicalize(parent_mol)
        parent_mol = tautomer_mol
        warnings.append("CANONICAL_TAUTOMER_REPRESENTATION: Canonical tautomer generated for 2D registration reproducibility, not thermodynamic equilibrium prediction.")
    except Exception as err:
        warnings.append(f"Tautomer canonicalization warning: {err}")

    # 7. Stereochemistry Preservation & Assignment
    try:
        Chem.AssignStereochemistry(parent_mol, cleanIt=True, force=True)
    except Exception as err:
        warnings.append(f"Stereochemistry assignment warning: {err}")

    # 8. Export Representations
    canonical_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=False)
    isomeric_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=True)
    inchi = Chem.MolToInchi(parent_mol) if parent_mol.GetNumHeavyAtoms() > 0 else ""
    inchikey = Chem.MolToInchiKey(parent_mol) if parent_mol.GetNumHeavyAtoms() > 0 else ""

    return {
        "status": status,
        "original_smiles": raw_smiles,
        "canonical_smiles": canonical_smiles,
        "isomeric_smiles": isomeric_smiles,
        "inchi": inchi,
        "inchikey": inchikey,
        "num_heavy_atoms": parent_mol.GetNumHeavyAtoms(),
        "has_chiral_centers": has_defined_chiral,
        "has_ez_bonds": has_defined_ez,
        "has_isotopes": has_isotopes,
        "salt_extracted": salt_extracted,
        "warnings": warnings,
        "provenance": _get_provenance(),
    }


def _get_provenance() -> dict[str, str]:
    return {
        "standardizer_name": STANDARDIZER_NAME,
        "standardizer_version": STANDARDIZER_VERSION,
        "rdkit_version": RDKIT_VERSION,
    }
