import importlib.util
from pathlib import Path

import pandas as pd


def _explorer_module():
    path = Path("scripts/explore_amrfinder.py")
    spec = importlib.util.spec_from_file_location("explore_amrfinder", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_feature_vector_uses_feature_union_and_binary_presence() -> None:
    evidence = pd.DataFrame(
        {
            "genome_id": ["g1", "g1", "g1"],
            "feature_key": ["gene::mecA", "mutation::gyrA_S84L", "gene::mecA"],
            "feature_value": [1, 1, 1],
        }
    )
    vector = _explorer_module().feature_vector(evidence, genome_id="g1")
    assert vector.loc[0, "gene::mecA"] == 1
    assert vector.loc[0, "mutation::gyrA_S84L"] == 1
    assert len(vector.columns) == 3
