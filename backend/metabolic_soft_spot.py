"""SyGMa-backed, atom-mapped metabolic hypotheses with strict scientific labels."""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
SYGMA_VERSION = "1.1.0"
RULE_ROOT = Path(__file__).resolve().parents[1] / "models" / "sygma"
ENGINE_NAME = "SyGMa empirical rules + RDKit chemical validation"
ENGINE_VERSION = f"stage3d-sygma-{SYGMA_VERSION}-rdkit-2025.03.1-v1"
ENGINE_SOURCE = "https://github.com/3D-e-Chem/sygma"
ENGINE_LICENSE = "GPL-3.0"
PREDICTED_LABEL = "PREDICTED METABOLITE HYPOTHESIS"
NO_CYP_ASSIGNMENT = "CYP isoform not assigned"

MODEL_STATUS = {
    "status": "MODEL_UNAVAILABLE",
    "reason": (
        "No atom-level ML checkpoint qualified for this installation: FAME3/FAME3R weights are "
        "non-commercial-research restricted; ATTNSOM publishes MIT code/data but no pretrained "
        "checkpoint; XenoSite has no reproducible local checkpoint. No model probability is generated."
    ),
    "evaluated_models": ["FAME3/FAME3R", "ATTNSOM", "XenoSite", "SMARTCyp 2.4.2"],
}

PUBLISHER_VALIDATION = {
    "SyGMa_full_phase1_phase2": "68% metabolite recall on 175 parent compounds; 45% of known metabolites within top 10",
    "SyGMa_single_step_CYP": "84% metabolite recall on 127 reactions; 66% of known metabolites within top 3",
    "scope": "Publisher-reported metabolite-level validation, not atom-probability calibration",
}

SUPPORTED_TRANSFORMATIONS = (
    "Aromatic hydroxylation", "Benzylic oxidation", "Aliphatic oxidation",
    "N-dealkylation", "O-dealkylation", "N-oxidation", "S-oxidation",
    "Ester hydrolysis", "Amide hydrolysis", "Glucuronidation", "Sulfation",
)

MITIGATION_STRATEGIES = {
    "Aromatic hydroxylation": ["block the reactive aryl position", "reduce local electron density", "reduce local lipophilicity"],
    "Benzylic oxidation": ["steric shielding", "fluorination", "remove the benzylic CH", "heteroatom replacement", "reduce local lipophilicity"],
    "Aliphatic oxidation": ["steric shielding", "fluorination where appropriate", "reduce exposed lipophilicity", "conformational constraint"],
    "N-dealkylation": ["steric shielding", "reduce N-basicity where appropriate", "constrain the substituent", "replace the metabolically labile alkyl group"],
    "O-dealkylation": ["steric shielding", "replace the labile ether", "constrain the O-alkyl substituent", "reduce local lipophilicity"],
    "N-oxidation": ["reduce N-basicity where appropriate", "steric shielding", "replace or constrain the tertiary amine"],
    "S-oxidation": ["replace the thioether", "steric shielding", "evaluate an oxidized bioisostere"],
    "Ester hydrolysis": ["steric shielding", "replace the ester with a hydrolysis-resistant bioisostere", "reduce esterase accessibility"],
    "Amide hydrolysis": ["steric shielding", "constrain amide geometry", "evaluate an amide bioisostere"],
    "Glucuronidation": ["sterically shield the conjugation handle", "modulate phenol/alcohol acidity", "evaluate a non-conjugating bioisostere"],
    "Sulfation": ["sterically shield the conjugation handle", "modulate phenol/alcohol electronics", "evaluate a non-conjugating bioisostere"],
}


@dataclass(frozen=True)
class RuleRecord:
    phase: str
    transformation: str
    name: str
    prior: float
    smirks: str
    count_evidence: str
    reaction: object


def _classify_rule(name: str, phase: str) -> str | None:
    lower = name.lower()
    if phase == "Phase II":
        if "glucuronidation" in lower:
            return "Glucuronidation"
        if "sulfation" in lower:
            return "Sulfation"
        return None
    if "benzylic_hydroxylation" in lower:
        return "Benzylic oxidation"
    if "aromatic_hydroxylation" in lower:
        return "Aromatic hydroxylation"
    if "aliphatic_hydroxylation" in lower:
        return "Aliphatic oxidation"
    if any(term in lower for term in ("n-demethylation", "n-dealkylation", "n-depropylation")):
        return "N-dealkylation"
    if any(term in lower for term in ("o-demethylation", "o-dealkylation")):
        return "O-dealkylation"
    if lower.startswith("n-oxidation"):
        return "N-oxidation"
    if "sulfide_oxidation" in lower or "sulfoxide_oxidation" in lower:
        return "S-oxidation"
    if lower in {"hydrolysis_(methoxyester)", "hydrolysis_(ester)"}:
        return "Ester hydrolysis"
    if lower.startswith("hydrolysis_") and "amide" in lower:
        return "Amide hydrolysis"
    return None


@lru_cache(maxsize=1)
def selected_rules() -> tuple[RuleRecord, ...]:
    selected = []
    for key, phase in (("phase1", "Phase I"), ("phase2", "Phase II")):
        with (RULE_ROOT / f"{key}.txt").open(encoding="utf-8") as stream:
            for line in stream:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                parts = line.rstrip().split("\t")
                if len(parts) < 3:
                    continue
                smirks, prior_text, name = parts[:3]
                transformation = _classify_rule(name, phase)
                if not transformation:
                    continue
                reaction = AllChem.ReactionFromSmarts(smirks)
                if reaction is None:
                    continue
                selected.append(RuleRecord(
                    phase=phase, transformation=transformation, name=name,
                    prior=float(prior_text), smirks=smirks,
                    count_evidence=parts[3].strip() if len(parts) > 3 else "",
                    reaction=reaction,
                ))
    return tuple(selected)


def _atom_environment(mol, atom_index: int) -> str:
    bonds = set(Chem.FindAtomEnvironmentOfRadiusN(mol, 1, atom_index))
    atoms = {atom_index}
    for bond_index in bonds:
        bond = mol.GetBondWithIdx(bond_index)
        atoms.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    return Chem.MolFragmentToSmiles(
        mol, atomsToUse=sorted(atoms), bondsToUse=sorted(bonds),
        canonical=True, isomericSmiles=True,
    ) or mol.GetAtomWithIdx(atom_index).GetSymbol()


def _mapped_parent_index(atom) -> int | None:
    return atom.GetIntProp("react_atom_idx") if atom.HasProp("react_atom_idx") else None


def _changed_parent_atoms(parent, product) -> set[int]:
    changed = set()
    mapped = {}
    for atom in product.GetAtoms():
        parent_index = _mapped_parent_index(atom)
        if parent_index is None:
            for neighbor in atom.GetNeighbors():
                neighbor_parent = _mapped_parent_index(neighbor)
                if neighbor_parent is not None:
                    changed.add(neighbor_parent)
            continue
        mapped[parent_index] = atom.GetIdx()
        original = parent.GetAtomWithIdx(parent_index)
        if atom.GetAtomicNum() != original.GetAtomicNum() or atom.GetFormalCharge() != original.GetFormalCharge():
            changed.add(parent_index)
    for bond in parent.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if first not in mapped or second not in mapped:
            continue
        product_bond = product.GetBondBetweenAtoms(mapped[first], mapped[second])
        if product_bond is None or product_bond.GetBondType() != bond.GetBondType():
            changed.update((first, second))
    return changed


def _source_atom(parent, product, transformation: str) -> int | None:
    changed = _changed_parent_atoms(parent, product)
    mapped_indices = {
        index for atom in product.GetAtoms()
        if (index := _mapped_parent_index(atom)) is not None
    }
    candidates = changed or mapped_indices
    if not candidates:
        return None

    def first_where(predicate):
        preferred = next((index for index in sorted(changed) if predicate(parent.GetAtomWithIdx(index))), None)
        return preferred if preferred is not None else next(
            (index for index in sorted(candidates) if predicate(parent.GetAtomWithIdx(index))), None,
        )

    if transformation == "Aromatic hydroxylation":
        return first_where(lambda atom: atom.GetAtomicNum() == 6 and atom.GetIsAromatic())
    if transformation == "Benzylic oxidation":
        return first_where(lambda atom: atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and any(n.GetIsAromatic() for n in atom.GetNeighbors()))
    if transformation == "Aliphatic oxidation":
        return first_where(lambda atom: atom.GetAtomicNum() == 6 and not atom.GetIsAromatic())
    if transformation in {"N-dealkylation", "N-oxidation"}:
        return first_where(lambda atom: atom.GetAtomicNum() == 7)
    if transformation == "O-dealkylation":
        return first_where(lambda atom: atom.GetAtomicNum() == 8)
    if transformation == "S-oxidation":
        return first_where(lambda atom: atom.GetAtomicNum() == 16)
    if transformation in {"Ester hydrolysis", "Amide hydrolysis"}:
        carbonyl = first_where(lambda atom: atom.GetAtomicNum() == 6 and any(
            bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8
            for bond in atom.GetBonds()
        ))
        return carbonyl if carbonyl is not None else min(changed)
    if transformation in {"Glucuronidation", "Sulfation"}:
        return first_where(lambda atom: atom.GetAtomicNum() in {7, 8})
    return min(candidates)


def _sanitized_fragments(product, parent_heavy_atoms: int, source_atom: int) -> list[dict]:
    structures = []
    for fragment in Chem.GetMolFrags(product, asMols=True, sanitizeFrags=False):
        candidate = copy.copy(fragment)
        try:
            Chem.SanitizeMol(candidate)
        except Exception:
            continue
        if candidate.GetNumHeavyAtoms() < max(2, math.ceil(parent_heavy_atoms * 0.40)):
            continue
        if not any(_mapped_parent_index(atom) == source_atom for atom in candidate.GetAtoms()):
            continue
        if any(atom.GetAtomicNum() == 0 for atom in candidate.GetAtoms()):
            continue
        canonical = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=False)
        isomeric = Chem.MolToSmiles(candidate, canonical=True, isomericSmiles=True)
        reparsed = Chem.MolFromSmiles(isomeric)
        if not canonical or reparsed is None or Chem.DetectChemistryProblems(reparsed):
            continue
        structures.append({"canonical_smiles": canonical, "isomeric_smiles": isomeric})
    return structures


def _reaction_products(parent, rule: RuleRecord) -> list[dict]:
    products = []
    parent_canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=False)
    for outcome in rule.reaction.RunReactants((parent,)):
        combined = None
        for product in outcome:
            combined = copy.copy(product) if combined is None else Chem.CombineMols(combined, product)
        if combined is None:
            continue
        try:
            Chem.SanitizeMol(combined)
        except Exception:
            # Disconnected products may sanitize only after fragmentation.
            pass
        source_atom = _source_atom(parent, combined, rule.transformation)
        if source_atom is None:
            continue
        for structure in _sanitized_fragments(combined, parent.GetNumHeavyAtoms(), source_atom):
            if structure["canonical_smiles"] == parent_canonical:
                continue
            products.append({**structure, "source_atom": source_atom})
    return products


def render_soft_spot_svg(smiles: str, spots: list[dict]) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")
    mol = copy.copy(mol)
    rdDepictor.Compute2DCoords(mol)
    colors = {1: (0.95, 0.20, 0.18), 2: (1.0, 0.60, 0.10), 3: (0.16, 0.48, 0.90)}
    atom_colors, atom_radii, highlight_atoms = {}, {}, []
    for spot in spots[:3]:
        atom_index = int(spot["atom_index"])
        if atom_index < 0 or atom_index >= mol.GetNumAtoms() or atom_index in highlight_atoms:
            continue
        rank = int(spot["rank"])
        highlight_atoms.append(atom_index)
        atom_colors[atom_index] = colors.get(rank, (0.55, 0.55, 0.55))
        atom_radii[atom_index] = 0.42
        mol.GetAtomWithIdx(atom_index).SetProp("atomNote", f"Rank {rank}")
    drawer = rdMolDraw2D.MolDraw2DSVG(560, 380)
    drawer.DrawMolecule(
        mol, highlightAtoms=highlight_atoms, highlightAtomColors=atom_colors,
        highlightAtomRadii=atom_radii,
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    rank_metadata = " · ".join(
        f"Rank {spot['rank']}: atom {spot['atom_index']} {spot['transformation']}"
        for spot in spots[:3]
    )
    return svg.replace("</svg>", f"<metadata id='soft-spot-ranks'>{rank_metadata}</metadata>\n</svg>")


def predict_soft_spots(smiles: str, context: dict | None = None, max_spots: int = 24) -> dict:
    """Generate atom-ranked reaction hypotheses; the score is a rule prior, not ML probability."""
    parent = Chem.MolFromSmiles(smiles)
    if parent is None:
        raise ValueError("Invalid SMILES")
    context = context or {}
    candidates: dict[tuple[int, str, str], dict] = {}
    for rule in selected_rules():
        for product in _reaction_products(parent, rule):
            key = (product["source_atom"], rule.transformation, rule.phase)
            row = candidates.setdefault(key, {
                "atom_index": product["source_atom"],
                "atom_environment": _atom_environment(parent, product["source_atom"]),
                "transformation": rule.transformation,
                "phase": rule.phase,
                "cyp_isoform": NO_CYP_ASSIGNMENT,
                "model_evidence": {**MODEL_STATUS, "linked_cyp_context": context.get("cyp", [])},
                "rule_evidence": {
                    "tool": "SyGMa", "tool_version": SYGMA_VERSION, "source": ENGINE_SOURCE,
                    "license": ENGINE_LICENSE, "rules": [],
                    "ranking_basis": "Maximum published SyGMa empirical occurrence prior among matching rules",
                    "not_atom_probability": True,
                    "mitigation_strategies": MITIGATION_STRATEGIES.get(rule.transformation, []),
                },
                "score": rule.prior,
                "score_type": "SyGMa empirical reaction-rule occurrence prior; not an atom probability",
                "confidence": "LOW",
                "products": {},
            })
            row["score"] = max(row["score"], rule.prior)
            evidence = {"rule_name": rule.name, "empirical_prior": rule.prior, "smirks": rule.smirks, "count_evidence": rule.count_evidence}
            if evidence not in row["rule_evidence"]["rules"]:
                row["rule_evidence"]["rules"].append(evidence)
            row["products"][product["canonical_smiles"]] = product

    ranked = sorted(candidates.values(), key=lambda row: (-row["score"], row["phase"] != "Phase I", row["atom_index"], row["transformation"]))[:max_spots]
    metabolites, seen_metabolites = [], set()
    for rank, spot in enumerate(ranked, 1):
        spot["rank"] = rank
        for structure in spot.pop("products").values():
            canonical = structure["canonical_smiles"]
            if canonical in seen_metabolites:
                continue
            seen_metabolites.add(canonical)
            metabolites.append({
                "canonical_smiles": canonical,
                "isomeric_smiles": structure["isomeric_smiles"],
                "transformation": spot["transformation"],
                "source_atom": spot["atom_index"],
                "phase": spot["phase"],
                "rank": rank,
                "confidence": spot["confidence"],
                "label": PREDICTED_LABEL,
                "evidence": {
                    "rule_evidence": spot["rule_evidence"],
                    "chemical_validation": "RDKit sanitization, valence/aromaticity/charge parse, parent and duplicate exclusion passed",
                },
            })

    summary = None
    if ranked:
        primary = ranked[0]
        summary = {
            "primary_predicted_liability": f"{primary['transformation']} at atom {primary['atom_index']}",
            "primary": {"rank": 1, "atom_index": primary["atom_index"], "transformation": primary["transformation"], "confidence": primary["confidence"]},
            "supporting_evidence": "SyGMa empirical reaction rule plus RDKit structural match; no qualified atom-level ML probability",
            "microsomal_evidence": context.get("microsomal", []),
            "cyp_evidence": context.get("cyp", []),
            "cyp_attribution_limit": "Compound-level CYP substrate evidence does not assign this atom or transformation to an isoform.",
            "mitigation_strategies": primary["rule_evidence"]["mitigation_strategies"],
            "strategy_scope": "General medicinal-chemistry strategies only; no analog structure or compound proposal generated.",
        }
    return {
        "engine": ENGINE_NAME,
        "engine_version": ENGINE_VERSION,
        "model_status": MODEL_STATUS,
        "publisher_validation": PUBLISHER_VALIDATION,
        "supported_transformations": list(SUPPORTED_TRANSFORMATIONS),
        "spots": ranked,
        "metabolites": sorted(metabolites, key=lambda row: (row["rank"], row["canonical_smiles"])),
        "highlighted_svg": render_soft_spot_svg(smiles, ranked),
        "liability_summary": summary or {
            "primary_predicted_liability": "No supported SyGMa transformation matched",
            "microsomal_evidence": context.get("microsomal", []),
            "cyp_evidence": context.get("cyp", []),
            "strategy_scope": "No analog structure or compound proposal generated.",
        },
    }
