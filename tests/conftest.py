# - Unifies GPU options/markers and existing fixtures
# - Sets conservative default parallelism for tests to 4
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Generator

import pytest
from rdkit import RDLogger

# ---------------------- Runtime stability knobs ----------------------
# Prefer 'spawn' so workers don't inherit a bunch of open FDs from the parent.
try:
    import multiprocessing as mp

    mp.set_start_method("spawn", force=True)
except Exception:
    pass

# Bump file descriptor soft limit if possible (Linux/macOS).
try:
    import resource  # type: ignore[attr-defined]

    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)  # type: ignore[attr-defined]
    new_soft = min(max(soft, 8192), hard)  # raise to 8192 if allowed
    if new_soft != soft:
        resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))  # type: ignore[attr-defined]
except Exception:
    pass

# Clamp hidden thread pools to avoid oversubscription during tests.
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

# Let our code (and wrappers like Gypsum) know to keep process parallelism modest.
os.environ.setdefault("STICKIT_PARALLEL_PROCS", "4")


def has_cuda() -> bool:
    if shutil.which("nvidia-smi"):
        try:
            smi = shutil.which("nvidia-smi")
            if not smi:
                return False
            out = subprocess.check_output(
                [smi], stderr=subprocess.STDOUT, timeout=5
            )  # nosec B603
            return b"CUDA" in out or b"NVIDIA-SMI" in out
        except Exception:
            return False
    try:
        import torch  # type: ignore

        return torch.cuda.is_available()
    except Exception:
        return False


# ----------------------------- PyTest CLI integration -----------------------------
def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--require-gpu",
        action="store_true",
        default=False,
        help="Fail tests marked 'gpu' if CUDA is unavailable (otherwise skip).",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    require_gpu = config.getoption("--require-gpu")
    for item in items:
        if "gpu" in item.keywords and not has_cuda():
            if require_gpu:
                item.add_marker(
                    pytest.mark.xfail(reason="CUDA not available", run=False)
                )
            else:
                item.add_marker(pytest.mark.skip(reason="CUDA not available"))


# ------------------------------ Global test fixtures ------------------------------
@pytest.fixture(autouse=True, scope="session")
def silence_rdkit_warnings():
    """Reduce RDKit noise in test output."""
    RDLogger.DisableLog("rdApp.warning")
    # To also hide RDKit errors (not recommended), uncomment:
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
        # Tests default to single-process or modest procs; your code can read STICKIT_PARALLEL_PROCS
        "parallel": {"backend": "none"},
    }


@pytest.fixture(scope="session")
def tmp_outdir(tmp_path_factory: pytest.TempPathFactory) -> Generator[Path, None, None]:
    p = tmp_path_factory.mktemp("_stickit")
    yield p
    # pytest cleans up


@pytest.fixture
def cfg(test_cfg_base, tmp_outdir: Path):
    cfg = json.loads(json.dumps(test_cfg_base))  # deep copy
    cfg["io"]["output_dir"] = str(tmp_outdir)
    return cfg


# ---------- Helpers to load SMILES ----------
def _load_smi(path: Path):
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


# Optional external benchmark (opt-in via env)
@pytest.fixture(scope="session")
def bench_data_dir():
    p = os.environ.get("STICKIT_BENCH_DATA")
    return Path(p) if p else None
