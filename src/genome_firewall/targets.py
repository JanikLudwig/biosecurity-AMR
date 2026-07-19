from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pyrodigal
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import Background, Pipeline


TARGET_REFERENCE_DIRECTORY = Path(__file__).parent / "resources" / "targets"
MIN_PROTEIN_IDENTITY = 0.80
MIN_PROTEIN_REFERENCE_COVERAGE = 0.60
MAX_PROTEIN_EVALUE = 1e-10
_AMINO = Alphabet.amino()


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


def _call_proteins(fasta: Path) -> list[tuple[str, str, str]]:
    """Call ORFs once and return (ORF id, contig, translated sequence)."""
    records = list(_fasta_records(fasta))
    sequences = [sequence for _, sequence in records if sequence]
    total = sum(map(len, sequences))
    if not sequences:
        return []
    finder = pyrodigal.GeneFinder(meta=total < 20_000)
    if total >= 20_000:
        finder.train(*[sequence.encode("ascii") for sequence in sequences])
    proteins: list[tuple[str, str, str]] = []
    for header, sequence in records:
        contig = header.split()[0]
        genes = finder.find_genes(sequence.encode("ascii"))
        for index, gene in enumerate(genes, start=1):
            translation = str(gene.translate()).rstrip("*")
            if translation:
                proteins.append((f"{contig}_{index}", contig, translation))
    return proteins


def _load_pyhmmer_references(
    directory: Path = TARGET_REFERENCE_DIRECTORY,
) -> dict[str, tuple[str, str]]:
    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"M2 target manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    references: dict[str, tuple[str, str]] = {}
    for path in sorted(directory.glob("*.fasta")):
        expected = manifest.get("references", {}).get(path.stem, {}).get("sha256")
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected != observed:
            raise RuntimeError(f"M2 target reference checksum mismatch: {path}")
        records = list(_fasta_records(path))
        if len(records) != 1:
            raise RuntimeError(f"Expected one target reference in {path}")
        header, sequence = records[0]
        references[path.stem] = (header, sequence)
    if not references:
        raise RuntimeError(f"No PyHMMER target references found in {directory}")
    return references


def _identity_and_coverage(hit: Any, reference_length: int) -> tuple[float, float]:
    alignment = hit.best_domain.alignment
    identity_sequence = alignment.identity_sequence or ""
    aligned = len(identity_sequence)
    identical = sum(1 for character in identity_sequence if character.isalpha())
    identity = identical / aligned if aligned else 0.0
    coverage = (
        (alignment.hmm_to - alignment.hmm_from + 1) / reference_length
        if reference_length
        else 0.0
    )
    return identity, coverage


def _search_pyhmmer_targets(
    proteins: list[tuple[str, str, str]],
    references: dict[str, tuple[str, str]],
) -> dict[str, dict[str, Any]]:
    contigs = {orf_id: contig for orf_id, contig, _ in proteins}
    database = DigitalSequenceBlock(
        _AMINO,
        [
            TextSequence(name=orf_id.encode(), sequence=sequence).digitize(_AMINO)
            for orf_id, _, sequence in proteins
        ],
    )
    background = Background(_AMINO)
    results: dict[str, dict[str, Any]] = {}
    for symbol, (header, sequence) in references.items():
        query = TextSequence(name=symbol.encode(), sequence=sequence).digitize(_AMINO)
        hits = Pipeline(_AMINO, background=background).search_seq(query, database)
        if not hits:
            results[symbol] = {
                "symbol": symbol,
                "reference": header,
                "present": False,
                "identity": 0.0,
                "reference_coverage": 0.0,
                "evalue": None,
                "bitscore": 0.0,
                "orf_id": None,
                "contig": None,
            }
            continue
        hit = hits[0]
        identity, coverage = _identity_and_coverage(hit, len(sequence))
        raw_name = hit.name
        orf_id = raw_name.decode() if isinstance(raw_name, bytes) else str(raw_name)
        present = (
            float(hit.evalue) <= MAX_PROTEIN_EVALUE
            and identity >= MIN_PROTEIN_IDENTITY
            and coverage >= MIN_PROTEIN_REFERENCE_COVERAGE
        )
        results[symbol] = {
            "symbol": symbol,
            "reference": header,
            "present": present,
            "identity": round(identity, 4),
            "reference_coverage": round(coverage, 4),
            "evalue": float(hit.evalue),
            "bitscore": round(float(hit.score), 2),
            "orf_id": orf_id,
            "contig": contigs.get(orf_id),
        }
    return results


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
    return analyze_drug_targets(
        fasta, amrfinder_executable=amrfinder_executable, drugs=drugs
    )["drugs"]


def analyze_drug_targets(
    fasta: Path,
    *,
    amrfinder_executable: Path,
    drugs: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Run M2 once and return both drug gates and auditable target evidence."""
    database = _database_directory(amrfinder_executable)
    protein_symbols = sorted(
        {
            symbol
            for drug in drugs.values()
            if drug["target_probe_kind"] in {"protein", "pyhmmer_protein"}
            for symbol in drug["target_symbols"]
        }
    )
    proteins: list[tuple[str, str, str]] = []
    protein_hits: dict[str, dict[str, Any]] = {}
    if protein_symbols:
        references = _load_pyhmmer_references()
        missing = sorted(set(protein_symbols).difference(references))
        if missing:
            raise RuntimeError(f"M2 target references missing: {missing}")
        proteins = _call_proteins(fasta)
        if proteins:
            protein_hits = _search_pyhmmer_targets(
                proteins, {symbol: references[symbol] for symbol in protein_symbols}
            )
        else:
            protein_hits = {
                symbol: {"symbol": symbol, "present": False, "reason": "no_orfs_predicted"}
                for symbol in protein_symbols
            }

    nucleotide_cache: dict[str, set[str]] = {}
    statuses: dict[str, dict[str, Any]] = {}
    for antibiotic, drug in drugs.items():
        kind = drug["target_probe_kind"]
        if kind in {"protein", "pyhmmer_protein"}:
            evidence = [protein_hits[symbol] for symbol in drug["target_symbols"]]
            detected = sorted(hit["symbol"] for hit in evidence if hit.get("present"))
        else:
            if kind not in nucleotide_cache:
                probes = _rrna_probe(database, kind)
                nucleotide_cache[kind] = _run_target_search(
                    fasta,
                    probes,
                    tool=amrfinder_executable.parent / "blastn",
                    protein=False,
                )
            detected = sorted(nucleotide_cache[kind])
            evidence = [
                {"symbol": symbol, "present": symbol in detected, "method": "blastn"}
                for symbol in drug["target_symbols"]
            ]
        required = int(drug["minimum_targets_detected"])
        statuses[antibiotic] = {
            "status": "present" if len(detected) >= required else "not_verified",
            "detected": detected,
            "required_count": required,
            "probe_kind": kind,
            "evidence": evidence,
        }
    return {
        "workflow": "M2",
        "method": "Pyrodigal ORF calling + PyHMMER protein homology; BLASTN for RNA targets",
        "predicted_proteins": len(proteins),
        "protein_thresholds": {
            "minimum_identity": MIN_PROTEIN_IDENTITY,
            "minimum_reference_coverage": MIN_PROTEIN_REFERENCE_COVERAGE,
            "maximum_evalue": MAX_PROTEIN_EVALUE,
        },
        "protein_hits": protein_hits,
        "drugs": statuses,
    }
