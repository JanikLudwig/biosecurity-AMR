from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class HealthResponse(BaseModel):
    status: str
    species: str
    antibiotics: list[str]
    amrfinder_version: str
    amrfinder_database: str
    model_bundle_schema: str


class ModelsResponse(BaseModel):
    species: str
    antibiotics: list[str]
    feature_count: int
    feature_schema_sha256: str
    bundle: dict[str, Any]
    warning: str


class AnalysisResponse(BaseModel):
    analysis_id: str
    schema_version: str
    generated_at: str
    genome_id: str
    species_scope: str
    warning: str
    defensive_use_only: bool
    qc: dict[str, Any]
    provenance: dict[str, Any]
    evidence_sources: dict[str, Any]
    workflows: dict[str, Any]
    lineage: dict[str, Any]
    decisions: list[dict[str, Any]] = Field(default_factory=list)
