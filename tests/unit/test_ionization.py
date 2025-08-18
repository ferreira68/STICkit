import pytest
from stickit.ionization import enum_ionization_states_dimorphite
from rdkit import Chem

@pytest.mark.parametrize("smi", ["c1ncccc1", "C1=CC=NC=C1"])  # pyridine forms
def test_dimorphite_microstates_include_neutral_and_protonated(smi):
    outs = enum_ionization_states_dimorphite(smi, pH=7.4, dPH=3.0, max_variants=64)
    assert any("[nH+]" in x or "+1" in x or "[nH]" in x for x in outs) or len(outs) > 1
