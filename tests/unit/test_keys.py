from rdkit import Chem
from stickit.utils import canonical_parent_key
from stickit.tautomer import tautomer_key
from stickit.stereo import stereo_key_for

def test_parent_key_strips_stereo():
    m = Chem.MolFromSmiles("C[C@H](O)C")
    # Compare to RDKit's canonical non-isomeric SMILES to avoid brittle literals
    expected = Chem.MolToSmiles(Chem.MolFromSmiles("CC(O)C"), isomericSmiles=False)
    assert canonical_parent_key(m) == expected

def test_tautomer_key_includes_explicit_H():
    m = Chem.MolFromSmiles("OC=CN")   # imine/enamine toy
    k1 = tautomer_key(m)
    k2 = tautomer_key(Chem.AddHs(m))
    assert k1 == k2  # our implementation adds Hs internally

def test_stereo_key_is_isomeric_smiles():
    m = Chem.MolFromSmiles("F/C=C/F")
    s = stereo_key_for(m)
    assert "/" in s or "\\" in s
