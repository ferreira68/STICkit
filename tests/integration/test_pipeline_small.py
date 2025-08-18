# tests/integration/test_pipeline_small.py
import pytest
from stickit.pipeline import stic_generation

@pytest.mark.parametrize("dataset_fixture", ["small_mols"])  # add "medium_mols" if you like
def test_pipeline_param(cfg, request, dataset_fixture):
    mols = request.getfixturevalue(dataset_fixture)
    sticsets = stic_generation(molecules=mols, config=cfg)
    assert sticsets
    n_inputs = len(mols)
    n_stics = sum(len(s.stics) for s in sticsets)
    n_confs = sum(len(st.conformers) for s in sticsets for st in s.stics)
    assert n_stics >= n_inputs
    assert n_confs >= n_inputs
