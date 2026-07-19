from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable


def _fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def _database_directory(amrfinder_executable: Path) -> Path:
    return amrfinder_executable.parent.parent / "share/amrfinderplus/data/latest"


def _protein_probes(
    database: Path, accessions: list[str], symbols: list[str]
) -> list[tuple[str, str]]:
    wanted = dict(zip(accessions, symbols, strict=True))
    probes: dict[str, tuple[str, str]] = {}
    for header, sequence in _fasta_records(database / "AMRProt.fa"):
        accession = header.split("|", 1)[0]
        if accession in wanted and accession not in probes:
            probes[accession] = (wanted[accession], sequence)
    missing = set(accessions).difference(probes)
    if missing:
        raise RuntimeError(f"AMRFinder target references missing: {sorted(missing)}")
    return [probes[accession] for accession in accessions]


def _rrna_probe(database: Path, kind: str) -> list[tuple[str, str]]:
    if kind == "staphylococcus_23s_rrna":
        path = database / "AMR_DNA-Staphylococcus_aureus.fa"
        marker = "@23S"
        symbol = "23S rRNA"
    elif kind == "conserved_16s_rrna":
        # AMRFinder's curated 16S reference is used only as a conserved target-presence
        # probe; resistance calls still come from the species-specific AMRFinder run.
        path = database / "AMR_DNA-Neisseria_gonorrhoeae.fa"
        marker = "@16S"
        symbol = "16S rRNA"
    else:
        raise ValueError(f"Unsupported rRNA target probe kind: {kind}")
    for header, sequence in _fasta_records(path):
        if marker in header:
            return [(symbol, sequence)]
    raise RuntimeError(f"Target reference {marker} missing from {path}")


def _run_target_search(
    fasta: Path,
    probes: list[tuple[str, str]],
    *,
    tool: Path,
    protein: bool,
) -> set[str]:
    with tempfile.TemporaryDirectory(prefix="genome-firewall-targets-") as temporary:
        query = Path(temporary) / "targets.fasta"
        query_symbols: dict[str, str] = {}
        with query.open("w", encoding="ascii") as handle:
            for index, (symbol, sequence) in enumerate(probes):
                query_id = f"probe_{index}"
                query_symbols[query_id] = symbol
                handle.write(f">{query_id}\n{sequence}\n")
        command = [
            str(tool),
            "-query",
            str(query),
            "-subject",
            str(fasta.resolve()),
            "-evalue",
            "1e-10",
            "-max_target_seqs",
            "20",
            "-outfmt",
            "6 qseqid qlen length pident evalue",
        ]
        if protein:
            command.extend(["-seg", "no"])
        else:
            command.extend(["-task", "blastn", "-dust", "no"])
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, "BLAST_USAGE_REPORT": "false"},
        )
    detected: set[str] = set()
    for line in result.stdout.splitlines():
        query_id, query_length, alignment_length, identity, _ = line.split("\t")
        coverage = int(alignment_length) / int(query_length)
        minimum_identity = 70.0 if protein else 70.0
        minimum_coverage = 0.80 if not protein else 0.85
        if coverage >= minimum_coverage and float(identity) >= minimum_identity:
            detected.add(query_symbols[query_id])
    return detected


def verify_drug_targets(
    fasta: Path,
    *,
    amrfinder_executable: Path,
    drugs: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    database = _database_directory(amrfinder_executable)
    search_cache: dict[tuple[Any, ...], set[str]] = {}
    statuses: dict[str, dict[str, Any]] = {}
    for antibiotic, drug in drugs.items():
        kind = drug["target_probe_kind"]
        if kind == "protein":
            key = (kind, *drug["target_accessions"])
            if key not in search_cache:
                probes = _protein_probes(
                    database, drug["target_accessions"], drug["target_symbols"]
                )
                search_cache[key] = _run_target_search(
                    fasta,
                    probes,
                    tool=amrfinder_executable.parent / "tblastn",
                    protein=True,
                )
        else:
            key = (kind,)
            if key not in search_cache:
                probes = _rrna_probe(database, kind)
                search_cache[key] = _run_target_search(
                    fasta,
                    probes,
                    tool=amrfinder_executable.parent / "blastn",
                    protein=False,
                )
        detected = sorted(search_cache[key])
        required = int(drug["minimum_targets_detected"])
        statuses[antibiotic] = {
            "status": "present" if len(detected) >= required else "not_verified",
            "detected": detected,
            "required_count": required,
        }
    return statuses
