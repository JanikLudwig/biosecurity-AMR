#!/usr/bin/env python3
"""Download explicit BV-BRC genome assemblies, reusing the local cache when possible.

Example:
    uv run python scripts/download_bvbrc_fastas.py \
        1280.51926 1280.51872 1280.51741 \
        --output-directory data/demo-data
"""

from __future__ import annotations

import argparse
import asyncio
import re
import shutil
from pathlib import Path

from genome_firewall.config import DEFAULT_CONFIG, load_config
from genome_firewall.data.bvbrc import BvbrcClient
from genome_firewall.data.qc import inspect_fasta


GENOME_ID_PATTERN = re.compile(r"^[0-9]+\.[0-9]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download selected assembled BV-BRC genomes as nucleotide FASTA files."
    )
    parser.add_argument("genome_ids", nargs="+", help="BV-BRC genome IDs, e.g. 1280.51926")
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=Path("data/demo-data"),
        help="Destination for <genome_id>.fna files.",
    )
    parser.add_argument(
        "--cache-directory",
        type=Path,
        default=Path("data/raw/genomes"),
        help="Existing assembly cache checked before making a network request.",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Ignore cached and destination files and request fresh copies from BV-BRC.",
    )
    return parser.parse_args()


async def obtain_fastas(args: argparse.Namespace) -> None:
    invalid = [
        genome_id for genome_id in args.genome_ids if not GENOME_ID_PATTERN.fullmatch(genome_id)
    ]
    if invalid:
        raise ValueError(f"Invalid BV-BRC genome ID(s): {', '.join(invalid)}")

    config = load_config(args.config)
    dataset = config["dataset"]
    bvbrc = config["bvbrc"]
    args.output_directory.mkdir(parents=True, exist_ok=True)

    async with BvbrcClient(
        bvbrc["api_base_url"],
        timeout_seconds=bvbrc["timeout_seconds"],
        sequence_result_limit=bvbrc["sequence_result_limit"],
    ) as client:
        for genome_id in dict.fromkeys(args.genome_ids):
            destination = args.output_directory / f"{genome_id}.fna"
            cached = args.cache_directory / f"{genome_id}.fna"
            source = "existing destination"

            if args.force_download:
                metadata = await client.metadata(genome_id)
                if (
                    metadata.get("species") != dataset["species"]
                    or int(metadata.get("taxon_id", -1)) != dataset["taxon_id"]
                ):
                    raise ValueError(f"{genome_id} is outside the configured species scope")
                await client.fasta(genome_id, destination)
                source = "BV-BRC API"
            elif destination.is_file():
                pass
            elif cached.is_file():
                shutil.copy2(cached, destination)
                source = f"local cache ({cached})"
            else:
                metadata = await client.metadata(genome_id)
                if (
                    metadata.get("species") != dataset["species"]
                    or int(metadata.get("taxon_id", -1)) != dataset["taxon_id"]
                ):
                    raise ValueError(f"{genome_id} is outside the configured species scope")
                await client.fasta(genome_id, destination)
                source = "BV-BRC API"

            metrics = await asyncio.to_thread(inspect_fasta, destination)
            print(
                f"{genome_id}: {destination} | {source} | "
                f"{metrics.genome_length:,} bp | {metrics.contigs:,} contigs | "
                f"N50 {metrics.contig_n50:,}"
            )


def main() -> None:
    asyncio.run(obtain_fastas(parse_args()))


if __name__ == "__main__":
    main()
