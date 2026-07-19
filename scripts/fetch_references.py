#!/usr/bin/env python3
"""Fetch the M2 target reference proteins from UniProt into data/references/targets/.

For each gene in ``gfw.targets.specs.TARGET_PROTEINS`` we pull the best
*S. aureus* entry (reviewed/SwissProt preferred), writing one ``<gene>.fasta``
plus a combined ``targets.fasta`` and a ``manifest.json`` recording provenance.

Stdlib-only HTTP (urllib) so there is no extra dependency. Network required.
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gfw.config import REFERENCES_DIR, ensure_dirs
from gfw.targets.specs import TARGET_PROTEINS

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def _fetch(query: str, reviewed: bool) -> str:
    q = query + (" AND reviewed:true" if reviewed else "")
    params = {"query": q, "format": "fasta", "size": "1"}
    url = f"{UNIPROT}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": "genome-firewall/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8")


def fetch_gene(gene: str, query: str) -> str | None:
    for reviewed in (True, False):
        try:
            fasta = _fetch(query, reviewed=reviewed)
        except Exception as exc:  # noqa: BLE001
            print(f"  ! {gene}: request failed ({exc})")
            continue
        if fasta.strip().startswith(">"):
            tag = "reviewed" if reviewed else "unreviewed"
            print(f"  ✓ {gene}: {tag}  ({len(fasta.split(chr(10))[0])} char header)")
            return fasta
    print(f"  ✗ {gene}: no sequence found")
    return None


def main() -> int:
    ensure_dirs()
    manifest = {}
    combined = []
    for gene, spec in TARGET_PROTEINS.items():
        fasta = fetch_gene(gene, spec["query"])
        time.sleep(0.3)  # be polite to UniProt
        if not fasta:
            continue
        # Re-tag the header so downstream code knows the gene symbol.
        header, *body = fasta.strip().split("\n")
        acc = header.split("|")[1] if "|" in header else header[1:].split()[0]
        new_header = f">{gene}|{acc} {spec['desc']}"
        record = "\n".join([new_header] + body) + "\n"
        with open(os.path.join(REFERENCES_DIR, f"{gene}.fasta"), "w") as fh:
            fh.write(record)
        combined.append(record)
        manifest[gene] = {"accession": acc, "desc": spec["desc"],
                          "length": sum(len(b) for b in body)}
    with open(os.path.join(REFERENCES_DIR, "targets.fasta"), "w") as fh:
        fh.write("".join(combined))
    with open(os.path.join(REFERENCES_DIR, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\nFetched {len(manifest)}/{len(TARGET_PROTEINS)} target references "
          f"-> {REFERENCES_DIR}")
    return 0 if manifest else 1


if __name__ == "__main__":
    raise SystemExit(main())
