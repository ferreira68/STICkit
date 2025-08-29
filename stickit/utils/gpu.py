from __future__ import annotations
import shutil
import subprocess


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
