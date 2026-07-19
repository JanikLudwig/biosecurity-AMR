from pathlib import Path

from genome_firewall.data.qc import inspect_fasta


def test_fasta_metrics(tmp_path: Path) -> None:
    fasta = tmp_path / "tiny.fna"
    fasta.write_text(">long\nAAAAANNNNN\n>short\nACGT\n", encoding="ascii")
    metrics = inspect_fasta(fasta)
    assert metrics.genome_length == 14
    assert metrics.contigs == 2
    assert metrics.contig_n50 == 10
    assert metrics.ambiguous_bases == 5
