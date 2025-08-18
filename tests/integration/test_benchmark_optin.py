import os
import pytest
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from stickit.pipeline import stic_generation
from tests.helpers import rmsd_between_confs

@pytest.mark.bench
def test_against_external_conformers(cfg, bench_data_dir):
    # Opt-in via STICKIT_BENCH_DATA=/path/to/sdf_dir (each SDF has 1+ conformers)
    if not bench_data_dir or not bench_data_dir.exists():
        pytest.skip("STICKIT_BENCH_DATA not set")
    sdf_files = list(bench_data_dir.glob("*.sdf"))[:5]
    assert sdf_files, "No SDFs in benchmark dir"

    for sdf in sdf_files:
        suppl = SDMolSupplier(str(sdf), removeHs=False)
        ref = [m for m in suppl if m]
        if not ref:
            continue
        ref_mol = ref[0]
        smi = Chem.MolToSmiles(Chem.RemoveHs(ref_mol))
        outs = stic_generation(molecules=[(smi, sdf.stem)], config=cfg)
        assert outs
        stic = outs[0].stics[0]
        # crude: compare best RMSD of any generated conf to the first ref conf
        best = min(rmsd_between_confs(stic.mol, c.conf_id, stic.mol.GetConformer(0).GetId())
                   for c in stic.conformers)
        assert best < 2.0  # loose threshold
