from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import yaml

try:
    import boto3  # optional
except Exception:
    boto3 = None  # type: ignore

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "configs"


@dataclass
class DBConfig:
    url: str = "postgresql+psycopg://user:pass@localhost:5432/stickit"


@dataclass
class SecretsManagerConfig:
    enabled: bool = False
    region_name: str = "us-east-1"
    secret_name: str = "stickit/default"


@dataclass
class Settings:
    db: DBConfig = DBConfig()
    secrets_manager: SecretsManagerConfig = SecretsManagerConfig()
    mlflow_tracking_uri: str = "file:./mlruns"

    @classmethod
    def load(cls) -> "Settings":
        data: Dict[str, Any] = {}

        plain = CONFIG_DIR / "settings.yaml"
        if plain.exists():
            data.update(yaml.safe_load(plain.read_text()) or {})

        sm = data.get("secrets_manager", {})
        if sm.get("enabled") and boto3 is not None:
            try:
                client = boto3.client(
                    "secretsmanager", region_name=sm.get("region_name", "us-east-1")
                )
                resp = client.get_secret_value(SecretId=sm["secret_name"])
                secret_str = resp.get("SecretString") or "{}"
                secret_data = json.loads(secret_str)
                data = _deep_merge(data, secret_data)
            except Exception as e:
                print(f"{e}")

        enc = CONFIG_DIR / "settings.enc.yaml"
        if enc.exists():
            try:
                import shutil
                import subprocess

                sops_path = shutil.which("sops")
                if not sops_path:
                    raise FileNotFoundError("sops not found on PATH")
                out = subprocess.check_output(
                    [sops_path, "-d", str(enc)], text=True, timeout=30
                )  # nosec B603
                enc_data = yaml.safe_load(out) or {}
                data = _deep_merge(data, enc_data)
            except Exception as e:
                print(f"{e}")

        db_conf = DBConfig(**(data.get("db") or {}))
        sm_conf = SecretsManagerConfig(**(data.get("secrets_manager") or {}))
        return Settings(
            db=db_conf,
            secrets_manager=sm_conf,
            mlflow_tracking_uri=data.get("mlflow_tracking_uri", "file:./mlruns"),
        )


def _deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)  # type: ignore
        else:
            out[k] = v
    return out
