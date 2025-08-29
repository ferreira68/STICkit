import pytest
from stickit.pipeline import stic_generation


@pytest.mark.parametrize(
    "smi,name",
    [
        ("N[C@@H](CC(=O)[O-])C(=O)O", "glycine_zwitterion"),
        ("N(C)(C)C", "tert_amine"),
    ],
)
def test_edge_zwitterion_and_amine(cfg, smi, name):
    outs = stic_generation(molecules=[(smi, name)], config=cfg)
    assert outs and outs[0].stics
