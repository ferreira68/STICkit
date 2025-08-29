import pytest
from rdkit import Chem
from stickit.tautomer import enum_tautomers


@pytest.mark.parametrize(
    "smi, min_count",
    [
        ("O=C1NC(=O)CC1", 1),
        ("ClC(C[C@@](Cl)(C)F)(F)C(C(C)C)=O", 4),
        ("O=C1NC=NC(=O)N1", 5),  # uracil
    ],
)
def test_enumerate_tautomers_returns_unique_variants(smi, min_count, test_cfg_base):
    m = Chem.MolFromSmiles(smi)
    ts = enum_tautomers(m, test_cfg_base)
    # unique by our key inside enum; ensure at least min_count
    assert len(ts) >= min_count
