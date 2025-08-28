from __future__ import annotations
import mlflow
from typing import Any
from ..config.settings import Settings

_settings = Settings.load()
mlflow.set_tracking_uri(_settings.mlflow_tracking_uri)

def register_model(model: Any, artifact_path: str, params: dict | None = None) -> str:
    with mlflow.start_run():
        if params:
            mlflow.log_params(params)
        mlflow.sklearn.log_model(model, artifact_path=artifact_path)  # change flavor as needed
        run = mlflow.active_run()
        assert run is not None
        return run.info.run_id

def load_model(model_uri: str) -> Any:
    return mlflow.pyfunc.load_model(model_uri)
