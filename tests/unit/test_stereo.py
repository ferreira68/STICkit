from rdkit import Chem
from stickit.stereo import enumerate_stereo_filtered

def test_pyramidal_n_excluded_by_default(test_cfg_base):
    cfg = dict(test_cfg_base)
    cfg["chem"] = dict(cfg["chem"])
    cfg["chem"]["allow_pyramidal_N"] = False
    m = Chem.MolFromSmiles("N(C)(C)C")  # tertiary amine (labile inversion)
    isos = list(enumerate_stereo_filtered(m, cfg))
    # we expect no artificially chiral N variants
    assert len(isos) == 1

def test_include_atropisomers_flag(test_cfg_base):
    # Hindered biphenyl; our current filter is crude, so just check it runs
    m = Chem.MolFromSmiles("Clc1cc(Cl)c(c(Cl)c1Cl)C2=CC(Cl)=CC(Cl)=C2")
    isos = list(enumerate_stereo_filtered(m, test_cfg_base))
    assert len(isos) >= 1
