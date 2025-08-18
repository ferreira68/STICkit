from rdkit import Chem
from stickit.conformers import make_conformers
from stickit.openmm_min import minimize_openmm

def test_openmm_min_returns_finite(cfg):
    m = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    m = Chem.AddHs(m)
    mol3d, ids = make_conformers(m, cfg)
    mol3d, e = minimize_openmm(mol3d, cfg, ids)
    assert isinstance(mol3d, Chem.Mol)
    assert all(k in e.keys() for k in ids)
    assert all(isinstance(t, tuple) and len(t) == 2 for t in e.values())
    for (a, b) in e.values():
        assert isinstance(a, float)
        assert isinstance(b, float)