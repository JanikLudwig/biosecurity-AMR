from pathlib import Path

from genome_firewall.annotation.amrfinder import parse_output


def test_parser_separates_gene_and_mutation_evidence(tmp_path: Path) -> None:
    output = tmp_path / "sample.tsv"
    output.write_text(
        "Element symbol\tElement name\tType\tSubtype\tClass\tSubclass\tMethod\t"
        "% Coverage of reference\t% Identity to reference\n"
        "mecA\tPBP2a\tAMR\tAMR\tBETA-LACTAM\tMETHICILLIN\tEXACTX\t100\t100\n"
        "gyrA_S84L\tGyrA\tAMR\tPOINT\tQUINOLONE\tQUINOLONE\tPOINTX\t100\t99.8\n",
        encoding="utf-8",
    )
    evidence = parse_output(output, genome_id="1280.test").set_index("element_symbol")
    assert evidence.loc["mecA", "feature_key"] == "gene::mecA"
    assert evidence.loc["mecA", "evidence_category"] == "known_resistance_gene"
    assert evidence.loc["gyrA_S84L", "feature_key"] == "mutation::gyrA::S84L"
    assert evidence.loc["gyrA_S84L", "evidence_category"] == "known_resistance_mutation"
