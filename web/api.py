"""Genome Firewall web backend (FastAPI).

Serves the single-page dashboard and a small JSON API:

    GET /                      -> the SPA
    GET /api/panel             -> the data-driven drug panel (tiers, targets)
    GET /api/metrics           -> per-drug hidden-test metrics (+ synthetic flag)
    GET /api/genomes           -> index of precomputed genome reports
    GET /api/report/{gid}      -> a genome report (precomputed if present, else live)
    GET /api/predict/{gid}     -> force a live end-to-end prediction
    GET /reports/{file}        -> static plots (reliability.png, performance.png)

Run:  uvicorn web.api:app --host 0.0.0.0 --port 8000   (from repo root)
"""

from __future__ import annotations

import json
import os
import sys

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw import SAFETY_NOTICE, SPECIES, __version__
from gfw.config import PANEL_JSON, REPORTS_DIR
from gfw.engine import Engine
from gfw.panel import build_panel, load_panel

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(HERE, "static")
GENOME_REPORTS = os.path.join(REPORTS_DIR, "genomes")

app = FastAPI(title="Genome Firewall", version=__version__)
_engine = None


def engine() -> Engine:
    global _engine
    if _engine is None:
        _engine = Engine()
    return _engine


@app.get("/api/meta")
def meta():
    features_synthetic = getattr(engine(), "features_synthetic", False)
    return {"version": __version__, "species": SPECIES,
            "safety_notice": SAFETY_NOTICE,
            "features_synthetic": bool(features_synthetic)}


@app.get("/api/panel")
def panel():
    try:
        p = load_panel()
    except Exception:
        p = build_panel()
    return [e.as_dict() for e in p]


@app.get("/api/metrics")
def metrics():
    path = os.path.join(REPORTS_DIR, "metrics.json")
    if not os.path.exists(path):
        return JSONResponse({"metrics": [], "synthetic_features": False})
    with open(path) as fh:
        return json.load(fh)


@app.get("/api/genomes")
def genomes():
    path = os.path.join(REPORTS_DIR, "index.json")
    if not os.path.exists(path):
        return {"genomes": [], "count": 0}
    with open(path) as fh:
        return json.load(fh)


@app.get("/api/report/{gid}")
def report(gid: str):
    cached = os.path.join(GENOME_REPORTS, f"{gid}.json")
    if os.path.exists(cached):
        with open(cached) as fh:
            return json.load(fh)
    return predict(gid)


@app.get("/api/predict/{gid}")
def predict(gid: str):
    try:
        rep = engine().predict_genome(gid)
    except FileNotFoundError:
        raise HTTPException(404, f"genome '{gid}' not found in genomes/")
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(500, str(exc))
    return rep.as_dict()


@app.get("/reports/{name}")
def report_asset(name: str):
    if name not in ("reliability.png", "performance.png"):
        raise HTTPException(404)
    path = os.path.join(REPORTS_DIR, name)
    if not os.path.exists(path):
        raise HTTPException(404)
    return FileResponse(path)


@app.get("/")
def index():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))


if os.path.isdir(STATIC_DIR):
    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
