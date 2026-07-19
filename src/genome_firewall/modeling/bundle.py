from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib

from genome_firewall.modeling.baseline import feature_schema_sha256


class ModelBundle:
    """Validated, eagerly loaded collection of per-antibiotic model artifacts."""

    def __init__(self, directory: Path) -> None:
        self.directory = directory.resolve()
        manifest_path = self.directory / "bundle-manifest.json"
        if not manifest_path.is_file():
            raise FileNotFoundError(f"Model bundle manifest is missing: {manifest_path}")
        self.manifest: dict[str, Any] = json.loads(manifest_path.read_text(encoding="utf-8"))
        if self.manifest.get("schema_version") != "genome-firewall-bundle-v1":
            raise ValueError("Unsupported model bundle schema")

        columns = list(self.manifest.get("feature_columns", []))
        expected_hash = self.manifest.get("feature_schema_sha256")
        if not columns or feature_schema_sha256(columns) != expected_hash:
            raise ValueError("Model bundle feature schema hash is invalid")

        self.artifacts: dict[str, dict[str, Any]] = {}
        for antibiotic in self.manifest.get("antibiotics", []):
            entry = self.manifest.get("models", {}).get(antibiotic)
            if not entry:
                raise ValueError(f"Manifest entry is missing for {antibiotic}")
            model_path = self.directory / entry["path"]
            artifact = joblib.load(model_path)
            if artifact.get("schema_version") != "genome-firewall-model-v1":
                raise ValueError(f"Unsupported artifact schema for {antibiotic}")
            if artifact.get("antibiotic") != antibiotic:
                raise ValueError(f"Artifact antibiotic mismatch for {antibiotic}")
            if artifact.get("feature_columns") != columns:
                raise ValueError(f"Feature order mismatch for {antibiotic}")
            if artifact.get("feature_schema_sha256") != expected_hash:
                raise ValueError(f"Feature hash mismatch for {antibiotic}")
            self.artifacts[antibiotic] = artifact

        if not self.artifacts:
            raise ValueError("Model bundle contains no trained models")

    @property
    def antibiotics(self) -> list[str]:
        return list(self.artifacts)

    @property
    def feature_columns(self) -> list[str]:
        return list(self.manifest["feature_columns"])

    def artifact(self, antibiotic: str) -> dict[str, Any]:
        try:
            return self.artifacts[antibiotic]
        except KeyError as error:
            raise KeyError(f"Unsupported antibiotic: {antibiotic}") from error

    def public_metadata(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in self.manifest.items()
            if key not in {"feature_columns"}
        }
