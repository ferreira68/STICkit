import pytest
from stickit.pipeline import stic_generation


@pytest.mark.medium
def test_pipeline_medium(cfg, medium_mols):
    # Increase conformers a bit for medium set
    cfg["conformers"]["num_confs"] = 8
    sticsets = stic_generation(molecules=medium_mols, config=cfg)
    assert sticsets
    # Ensure caps are respected and runtime stays bounded
    assert sum(len(st.stics) for st in sticsets) < 200
