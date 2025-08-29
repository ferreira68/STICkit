import os
import json
import pytest
from pathlib import Path
from rdkit import RDLogger


@pytest.fixture(autouse=True, scope="session")
def silence_rdkit_warnings():
    # Hide lines like "Explicit valence for atom # ... is greater than permitted"
    RDLogger.DisableLog("rdApp.warning")
    # If you also want to hide RDKit errors (not recommended), uncomment:
    # RDLogger.DisableLog("rdApp.error")


# ---------- Config fixtures tuned for tests ----------


@pytest.fixture(scope="session")
def test_cfg_base():
    # Minimal, fast defaults; MOPAC disabled in CI
    return {
        "io": {"output_dir": "./_test_out", "write_json": False, "write_sdf": False},
        "chem": {
            "ph": 7.4,
            "ph_delta": 2.0,
            "max_ti_iters": 1,
            "max_stereoisomers": 16,
            "include_atropisomers": True,
            "allow_pyramidal_N": False,
        },
        "conformers": {"engine": "rdkit", "num_confs": 5, "rmsd_thresh": 0.5},
        "minimization": {
            "engine": "openmm",
            "forcefield_spec": "openff-2.2.0.offxml",
            "implicit_solvent": None,
            "max_steps": 200,
            "gradient": 0.01,
        },
        "mopac": {"method": "PM6-D3H4", "enable": False},
        "parallel": {"backend": "none"},
    }


@pytest.fixture(scope="session")
def tmp_outdir(tmp_path_factory):
    p = tmp_path_factory.mktemp("_stickit")
    yield p
    # cleanup handled by pytest


@pytest.fixture
def cfg(test_cfg_base, tmp_outdir):
    cfg = json.loads(json.dumps(test_cfg_base))  # deep copy
    cfg["io"]["output_dir"] = str(tmp_outdir)
    return cfg


# ---------- Helpers to load SMILES ----------


def _load_smi(path):
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            smi, *rest = line.strip().split()
            name = rest[0] if rest else smi
            yield smi, name


@pytest.fixture(scope="session")
def small_mols():
    return list(_load_smi(Path(__file__).parent / "data" / "small.smi"))


@pytest.fixture(scope="session")
def medium_mols():
    return list(_load_smi(Path(__file__).parent / "data" / "medium.smi"))


@pytest.fixture(scope="session")
def large_mols():
    return list(_load_smi(Path(__file__).parent / "data" / "large.smi"))


@pytest.fixture(scope="session")
def edge_mols():
    return list(_load_smi(Path(__file__).parent / "data" / "edge.smi"))


# Optional external benchmark (opt-in)
@pytest.fixture(scope="session")
def bench_data_dir():
    p = os.environ.get("STICKIT_BENCH_DATA")
    return Path(p) if p else None
