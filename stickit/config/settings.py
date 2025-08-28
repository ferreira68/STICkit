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

        # 1) Non-secret defaults
        plain = CONFIG_DIR / "settings.yaml"
        if plain.exists():
            data.update(yaml.safe_load(plain.read_text()) or {})

        # 2) AWS Secrets Manager overlay (preferred secret store)
        sm = data.get("secrets_manager", {})
        if sm.get("enabled") and boto3 is not None:
            try:
                client = boto3.client("secretsmanager", region_name=sm.get("region_name", "us-east-1"))
                resp = client.get_secret_value(SecretId=sm["secret_name"])
                secret_str = resp.get("SecretString") or "{}"
                secret_data = json.loads(secret_str)
                data = _deep_merge(data, secret_data)
            except Exception:
                pass  # soft-fail

        # 3) SOPS-encrypted YAML overlay (optional)
        enc = CONFIG_DIR / "settings.enc.yaml"
        if enc.exists():
            try:
                import subprocess
                out = subprocess.check_output(["sops", "-d", str(enc)], text=True)
                enc_data = yaml.safe_load(out) or {}
                data = _deep_merge(data, enc_data)
            except Exception:
                pass

        # dataclass mapping
        db_conf = DBConfig(**(data.get("db") or {}))
        sm_conf = SecretsManagerConfig(**(data.get("secrets_manager") or {}))
        return Settings(db=db_conf, secrets_manager=sm_conf, mlflow_tracking_uri=data.get("mlflow_tracking_uri", "file:./mlruns"))

def _deep_merge(a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)  # type: ignore
        else:
            out[k] = v
    return out
