from __future__ import annotations
import shutil, subprocess

def has_cuda() -> bool:
    if shutil.which("nvidia-smi"):
        try:
            out = subprocess.check_output(["nvidia-smi"], stderr=subprocess.STDOUT)
            return b"CUDA" in out or b"NVIDIA-SMI" in out
        except Exception:
            return False
    try:
        import torch  # type: ignore
        return torch.cuda.is_available()
    except Exception:
        return False
