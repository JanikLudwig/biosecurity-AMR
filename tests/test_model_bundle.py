import hashlib
import json
from pathlib import Path

import joblib
import pytest

from genome_firewall.modeling.bundle import ModelBundle


def _bundle(tmp_path: Path, *, artifact_hash: str | None = None) -> Path:
    columns = ["gene::mecA", "mutation::gyrA::S84L"]
    schema_hash = hashlib.sha256("\n".join(columns).encode()).hexdigest()
    joblib.dump(
        {
            "schema_version": "genome-firewall-model-v1",
            "antibiotic": "cefoxitin",
            "feature_columns": columns,
            "feature_schema_sha256": artifact_hash or schema_hash,
        },
        tmp_path / "cefoxitin.joblib",
    )
    (tmp_path / "bundle-manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "genome-firewall-bundle-v1",
                "antibiotics": ["cefoxitin"],
                "feature_count": 2,
                "feature_columns": columns,
                "feature_schema_sha256": schema_hash,
                "models": {"cefoxitin": {"path": "cefoxitin.joblib"}},
            }
        ),
        encoding="utf-8",
    )
    return tmp_path


def test_model_bundle_loads_validated_artifacts_once(tmp_path: Path) -> None:
    bundle = ModelBundle(_bundle(tmp_path))
    assert bundle.antibiotics == ["cefoxitin"]
    assert bundle.feature_columns == ["gene::mecA", "mutation::gyrA::S84L"]


def test_model_bundle_rejects_feature_hash_mismatch(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Feature hash mismatch"):
        ModelBundle(_bundle(tmp_path, artifact_hash="wrong"))
