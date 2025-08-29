from typing import List
from rdkit import Chem
from rdkit.Chem import rdMolAlign


def rmsd_between_confs(mol: Chem.Mol, i: int, j: int) -> float:
    return rdMolAlign.GetBestRMS(mol, mol, prbId=i, refId=j)


def unique_confs_by_rmsd(mol: Chem.Mol, ids: List[int], thresh: float = 0.5) -> List[int]:
    kept: List[int] = []
    for i in ids:
        if all(rmsd_between_confs(mol, i, k) >= thresh for k in kept):
            kept.append(i)
    return kept
