import yaml
from pathlib import Path
from typing import Any, Dict

def get_repo_root() -> Path:
    """Returns the absolute path to the repository root."""
    # Assuming config.py is in src/genome_firewall/
    return Path(__file__).resolve().parent.parent.parent

def load_config(config_path: str = "config/pipeline.yaml") -> Dict[str, Any]:
    """
    Loads and validates the pipeline YAML configuration.
    Paths are resolved relative to the repository root.
    """
    repo_root = get_repo_root()
    full_path = repo_root / config_path

    if not full_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {full_path}")

    try:
        with open(full_path, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML schema in {full_path}: {e}")

    if not isinstance(config, dict):
        raise ValueError(f"Configuration in {full_path} must be a YAML dictionary.")

    # Validate required sections
    if "amrfinder" not in config:
        raise ValueError("Missing 'amrfinder' section in configuration.")
    if "paths" not in config:
        raise ValueError("Missing 'paths' section in configuration.")

    # Ensure critical defaults are not silently missing
    amr_config = config["amrfinder"]
    if "docker_image" not in amr_config or not amr_config["docker_image"]:
        raise ValueError("Configuration must explicitly define amrfinder.docker_image.")
    if "organism" not in amr_config or not amr_config["organism"]:
        raise ValueError("Configuration must explicitly define amrfinder.organism.")

    return config
