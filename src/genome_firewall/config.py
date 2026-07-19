from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

DEFAULT_CONFIG = Path("configs/staphylococcus_aureus.toml")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load the checked-in experiment configuration."""
    with path.open("rb") as handle:
        return tomllib.load(handle)

