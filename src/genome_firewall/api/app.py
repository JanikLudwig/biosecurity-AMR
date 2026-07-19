from __future__ import annotations

import asyncio
import csv
import hashlib
import json
import logging
import os
import tempfile
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from genome_firewall.annotation.amrfinder import (
    database_version,
    executable_version,
    resolve_executable,
)
from genome_firewall.api.schemas import AnalysisResponse, HealthResponse, ModelsResponse
from genome_firewall.config import DEFAULT_CONFIG, load_config
from genome_firewall.decision import DEFAULT_DRUG_REGISTRY
from genome_firewall.inference import CONFIRMATION_WARNING, predict_fasta
from genome_firewall.modeling.bundle import ModelBundle


LOGGER = logging.getLogger(__name__)
STATIC_DIRECTORY = Path(__file__).parent / "static"
ALLOWED_SUFFIXES = {".fna", ".fa", ".fasta"}
MAX_UPLOAD_BYTES = 25 * 1024 * 1024
MIN_ANALYSIS_RESPONSE_SECONDS = 3.0


def _analysis_directory(reports: Path, analysis_id: str) -> Path:
    if not analysis_id.isalnum() or len(analysis_id) != 20:
        raise HTTPException(404, "Analysis not found")
    return reports / analysis_id


def _load_cached_report(directory: Path, expected_sha256: str) -> dict[str, Any] | None:
    """Return a complete report produced for the same FASTA bytes, if available."""
    matches = list(directory.glob("*.report.json"))
    if len(matches) != 1:
        return None
    try:
        report = json.loads(matches[0].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        LOGGER.warning("Ignoring unreadable cached analysis report at %s", matches[0])
        return None
    if report.get("qc", {}).get("sha256") != expected_sha256:
        LOGGER.warning("Ignoring cached analysis with a mismatched FASTA digest at %s", matches[0])
        return None
    return report


def _reliability_data(model_directory: Path) -> dict[str, list[dict[str, Any]]]:
    path = model_directory / "test-reliability.csv"
    if not path.is_file():
        return {}
    grouped: dict[str, list[dict[str, Any]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            antibiotic = row.pop("antibiotic")
            grouped.setdefault(antibiotic, []).append(
                {
                    "probability_bin": int(row["probability_bin"]),
                    "samples": int(row["samples"]),
                    "mean_probability_resistant": float(row["mean_probability_resistant"]),
                    "observed_resistant_fraction": float(row["observed_resistant_fraction"]),
                }
            )
    return grouped


def _add_decision_thresholds(public_metadata: dict[str, Any], model_directory: Path) -> None:
    for entry in public_metadata.get("models", {}).values():
        metadata_path = model_directory / entry["metadata_path"]
        if not metadata_path.is_file():
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        entry["decision_thresholds"] = metadata.get("decision_thresholds", {})


def create_app(
    *,
    model_directory: Path | None = None,
    report_directory: Path | None = None,
    config_path: Path = DEFAULT_CONFIG,
    registry_path: Path = DEFAULT_DRUG_REGISTRY,
) -> FastAPI:
    models = model_directory or Path(
        os.environ.get("GENOME_FIREWALL_MODEL_DIRECTORY", "artifacts/models-v1")
    )
    reports = report_directory or Path(
        os.environ.get("GENOME_FIREWALL_REPORT_DIRECTORY", "artifacts/api-reports")
    )

    @asynccontextmanager
    async def lifespan(application: FastAPI):
        config = load_config(config_path)
        bundle = ModelBundle(models)
        executable = resolve_executable(config["amrfinder"]["executable"])
        installed_database = database_version(executable) if executable else None
        if executable is None or installed_database is None:
            raise RuntimeError("AMRFinderPlus executable/database is not ready")
        application.state.config = config
        application.state.bundle = bundle
        application.state.amrfinder_executable = executable
        application.state.amrfinder_version = executable_version(executable)
        application.state.amrfinder_database = installed_database
        application.state.model_directory = models
        application.state.report_directory = reports
        application.state.registry_path = registry_path
        application.state.analysis_semaphore = asyncio.Semaphore(2)
        application.state.analysis_locks = {}
        reports.mkdir(parents=True, exist_ok=True)
        yield

    application = FastAPI(
        title="Genome Firewall API",
        version="1.0.0",
        description=(
            "Defensive research API for assembled Staphylococcus aureus genomes. "
            "Every result requires laboratory confirmation."
        ),
        lifespan=lifespan,
    )
    configured_origins = os.environ.get("GENOME_FIREWALL_CORS_ORIGINS", "")
    development_origins = [
        "http://127.0.0.1:3000",
        "http://localhost:3000",
        "http://127.0.0.1:5173",
        "http://localhost:5173",
    ]
    application.add_middleware(
        CORSMiddleware,
        allow_origins=[origin.strip() for origin in configured_origins.split(",") if origin.strip()]
        or development_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST"],
        allow_headers=["Content-Type"],
    )

    @application.get("/api/v1/health", response_model=HealthResponse)
    async def health(request: Request) -> dict[str, Any]:
        bundle: ModelBundle = request.app.state.bundle
        return {
            "status": "ready",
            "species": bundle.manifest.get("species", "Staphylococcus aureus"),
            "antibiotics": bundle.antibiotics,
            "amrfinder_version": request.app.state.amrfinder_version,
            "amrfinder_database": request.app.state.amrfinder_database,
            "model_bundle_schema": bundle.manifest["schema_version"],
        }

    @application.get("/api/v1/models", response_model=ModelsResponse)
    async def available_models(request: Request) -> dict[str, Any]:
        bundle: ModelBundle = request.app.state.bundle
        public_metadata = bundle.public_metadata()
        _add_decision_thresholds(public_metadata, request.app.state.model_directory)
        public_metadata["reliability"] = _reliability_data(request.app.state.model_directory)
        return {
            "species": bundle.manifest.get("species", "Staphylococcus aureus"),
            "antibiotics": bundle.antibiotics,
            "feature_count": bundle.manifest["feature_count"],
            "feature_schema_sha256": bundle.manifest["feature_schema_sha256"],
            "bundle": public_metadata,
            "warning": CONFIRMATION_WARNING,
        }

    @application.post("/api/v1/analyses", response_model=AnalysisResponse)
    async def analyze(
        request: Request,
        fasta: UploadFile = File(..., description="One assembled S. aureus FASTA file"),
    ) -> dict[str, Any]:
        started_at = asyncio.get_running_loop().time()
        suffix = Path(fasta.filename or "").suffix.lower()
        if suffix not in ALLOWED_SUFFIXES:
            raise HTTPException(422, "Expected an assembled .fna, .fa, or .fasta file")

        with tempfile.TemporaryDirectory(prefix="genome-firewall-api-") as temporary:
            temporary_directory = Path(temporary)
            initial_path = temporary_directory / f"upload{suffix}"
            digest = hashlib.sha256()
            size = 0
            with initial_path.open("wb") as handle:
                while chunk := await fasta.read(1024 * 1024):
                    size += len(chunk)
                    if size > MAX_UPLOAD_BYTES:
                        raise HTTPException(413, "FASTA exceeds the 25 MiB upload limit")
                    digest.update(chunk)
                    handle.write(chunk)
            await fasta.close()
            if size == 0:
                raise HTTPException(422, "Uploaded FASTA is empty")

            fasta_sha256 = digest.hexdigest()
            analysis_id = fasta_sha256[:20]
            input_path = temporary_directory / f"{analysis_id}{suffix}"
            initial_path.replace(input_path)
            destination = request.app.state.report_directory / analysis_id
            analysis_lock = request.app.state.analysis_locks.setdefault(
                analysis_id, asyncio.Lock()
            )
            async with analysis_lock:
                report = _load_cached_report(destination, fasta_sha256)
                if report is not None:
                    LOGGER.info("Reusing cached genome analysis %s", analysis_id)
                else:
                    try:
                        async with request.app.state.analysis_semaphore:
                            report = await asyncio.to_thread(
                                predict_fasta,
                                input_path,
                                config=request.app.state.config,
                                model_directory=request.app.state.model_directory,
                                output_directory=destination,
                                registry_path=request.app.state.registry_path,
                                model_bundle=request.app.state.bundle,
                                threads=2,
                            )
                    except ValueError as error:
                        raise HTTPException(422, str(error)) from error
                    except Exception as error:
                        LOGGER.exception("Genome analysis %s failed", analysis_id)
                        raise HTTPException(
                            500,
                            "Genome analysis failed; inspect the server log for the internal error",
                        ) from error
        report["analysis_id"] = analysis_id
        elapsed = asyncio.get_running_loop().time() - started_at
        if elapsed < MIN_ANALYSIS_RESPONSE_SECONDS:
            await asyncio.sleep(MIN_ANALYSIS_RESPONSE_SECONDS - elapsed)
        return report

    @application.get("/api/v1/analyses/{analysis_id}", response_model=AnalysisResponse)
    async def saved_analysis(analysis_id: str, request: Request) -> dict[str, Any]:
        directory = _analysis_directory(request.app.state.report_directory, analysis_id)
        matches = list(directory.glob("*.report.json"))
        if len(matches) != 1:
            raise HTTPException(404, "Analysis not found")
        report = json.loads(matches[0].read_text(encoding="utf-8"))
        report["analysis_id"] = analysis_id
        return report

    @application.get("/api/v1/analyses/{analysis_id}/raw/m1", include_in_schema=True)
    async def raw_m1(analysis_id: str, request: Request) -> FileResponse:
        """Return the unmodified AMRFinderPlus TSV produced by workflow M1."""
        directory = _analysis_directory(request.app.state.report_directory, analysis_id)
        matches = list(directory.glob("*.amrfinder.tsv"))
        if len(matches) != 1:
            raise HTTPException(404, "Raw M1 output not found")
        return FileResponse(
            matches[0],
            media_type="text/tab-separated-values",
            filename=f"{analysis_id}.amrfinder.tsv",
            content_disposition_type="inline",
        )

    @application.get("/api/v1/analyses/{analysis_id}/raw/m2")
    async def raw_m2(analysis_id: str, request: Request) -> dict[str, Any]:
        """Return the complete structured target-detection output from workflow M2."""
        directory = _analysis_directory(request.app.state.report_directory, analysis_id)
        matches = list(directory.glob("*.report.json"))
        if len(matches) != 1:
            raise HTTPException(404, "Raw M2 output not found")
        report = json.loads(matches[0].read_text(encoding="utf-8"))
        try:
            return report["workflows"]["M2"]
        except KeyError as error:
            raise HTTPException(404, "Raw M2 output not found") from error

    @application.get("/", include_in_schema=False)
    async def index() -> FileResponse:
        return FileResponse(STATIC_DIRECTORY / "index.html")

    application.mount("/static", StaticFiles(directory=STATIC_DIRECTORY), name="static")
    return application


app = create_app()
