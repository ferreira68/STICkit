# Optional pytest add-on (won't conflict if you have your own conftest.py).
import pytest
from stickit.utils.gpu import has_cuda

def pytest_addoption(parser):
    parser.addoption("--require-gpu", action="store_true", default=False,
                     help="Fail tests marked gpu if CUDA is unavailable")

def pytest_collection_modifyitems(config, items):
    require_gpu = config.getoption("--require-gpu")
    for item in items:
        if "gpu" in item.keywords and not has_cuda():
            if require_gpu:
                item.add_marker(pytest.mark.xfail(reason="CUDA not available", run=False))
            else:
                item.add_marker(pytest.mark.skip(reason="CUDA not available"))
