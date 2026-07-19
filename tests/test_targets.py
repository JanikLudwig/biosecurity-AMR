from pathlib import Path

from genome_firewall.targets import (
    TARGET_REFERENCE_DIRECTORY,
    _call_proteins,
    _load_pyhmmer_references,
    _search_pyhmmer_targets,
)


def test_pyhmmer_target_workflow_detects_reference_orf(tmp_path: Path) -> None:
    references = _load_pyhmmer_references(TARGET_REFERENCE_DIRECTORY)
    protein = references["gyrA"][1]
    # Reverse translate with a deterministic codon per amino acid so Pyrodigal can call it.
    codons = {
        "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
        "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
        "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
        "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
    }
    sequence = "ATG" + "".join(codons[amino] for amino in protein[1:]) + "TAA"
    fasta = tmp_path / "target.fna"
    fasta.write_text(f">target\n{sequence}\n", encoding="ascii")
    proteins = _call_proteins(fasta)
    hits = _search_pyhmmer_targets(proteins, {"gyrA": references["gyrA"]})
    assert hits["gyrA"]["present"] is True
