from rdkit import Chem
from stickit.conformers import make_conformers
from tests.helpers import unique_confs_by_rmsd


def test_rdkit_conformers_embed_and_dedup(cfg):
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
    mol3d, ids = make_conformers(m, cfg)
    assert len(ids) > 0
    uniq = unique_confs_by_rmsd(mol3d, ids, 0.3)
    assert len(uniq) <= len(ids)
