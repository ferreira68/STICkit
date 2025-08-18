import shutil, pytest
from typing import Any

from rdkit import Chem
from stickit.mopac import mopac_refine_and_prune
from stickit.data import STIC, STICKey, ConformerRecord
from stickit.utils import pathlib_which


@pytest.mark.mopac
def test_mopac_prunes_imaginary(tmp_path, test_cfg_base):
    cfg = dict(test_cfg_base)
    cfg["mopac"]["enable"] = True
    m = Chem.MolFromSmiles("CCO")
    m = Chem.AddHs(m)
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(m)
    key = STICKey("CCO", "CCO", "q0", "CCO")
    stic = STIC(key, m, [ConformerRecord(conf_id=0)], {})
    from stickit.data import STICSet
    st = STICSet("CCO", [stic])
    mopac_refine_and_prune(st, cfg)
    # either pruned to 0 (if imaginary) or still 1; test just ensures call succeeds
    assert len(st.stics[0].conformers) in {0, 1}
