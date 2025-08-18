from dimorphite_dl.protonate.run import protonate_smiles
from rdkit import Chem


def enum_ionization_states_dimorphite(smi: str, pH: float, dPH: float, max_variants: int = 256):
    lo, hi = round(pH - dPH, 1), round(pH + dPH, 1)
    protomers = protonate_smiles(smiles_input=smi, ph_min=min(lo, hi, pH), ph_max=max(lo, hi, pH),
                                 max_variants=max_variants)
    return protomers


def ionization_key(mol):
    chg = Chem.GetFormalCharge(mol)
    charges = tuple(sorted((a.GetIdx(), a.GetFormalCharge()) for a in mol.GetAtoms()))
    return f"q{chg}|" + ",".join(f"{i}:{c}" for i, c in charges)


def ti_enables_new_tautomer_rules(mol, cfg) -> bool:
    patts = ["[O-]-C=O", "c[N-]", "C=[O-]"]  # extend with your charge-gated rule triggers
    for p in patts:
        q = Chem.MolFromSmarts(p)
        if mol.HasSubstructMatch(q):
            return True
    return False
