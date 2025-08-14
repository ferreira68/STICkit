from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def enum_tautomers(mol, cfg, charge_aware=False):
    te = rdMolStandardize.TautomerEnumerator()
    tset = te.Enumerate(mol)
    out, seen = [], set()
    for t in tset:
        smi = Chem.MolToSmiles(Chem.AddHs(t), isomericSmiles=False)
        if smi not in seen:
            seen.add(smi)
            out.append(t)
    return out

def tautomer_key(mol):
    return Chem.MolToSmiles(Chem.AddHs(mol), isomericSmiles=False)

