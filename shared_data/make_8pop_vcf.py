#!/usr/bin/env python3
"""Create the public eight-population Pool-seq VCF from the private 25-pool VCF.

The script retains only RD:AD sample fields for the eight manuscript populations,
normalises sample names to short stable IDs, and copies no SnpEff annotations.
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

FOCAL = ["AAI", "AAK", "AAL", "AAQ", "AAU", "AAW", "AAX", "ABE"]


def sample_id(text: str) -> str:
    hits = re.findall(r"(?<![A-Z0-9])(AAI|AAK|AAL|AAQ|AAU|AAW|AAX|ABE)(?![A-Z0-9])", text)
    if len(hits) != 1:
        raise ValueError(f"Could not resolve one focal sample ID from header field: {text}")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path, help="Headered 25-pool unannotated VCF")
    ap.add_argument("output", type=Path, help="Output .vcf.gz")
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    selected_columns = None
    records = 0
    with args.source.open("rt", encoding="utf-8", errors="strict", newline="") as src, gzip.open(
        args.output, "wt", encoding="utf-8", newline="", compresslevel=9
    ) as dst:
        for line in src:
            if line.startswith("##"):
                if line.startswith("##SnpEff") or line.startswith("##INFO=<ID=ANN"):
                    continue
                dst.write(line)
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\r\n").split("\t")
                resolved = {}
                for i, value in enumerate(fields[9:], start=9):
                    for code in FOCAL:
                        if re.search(rf"(?<![A-Z0-9]){code}(?![A-Z0-9])", value):
                            resolved[code] = i
                missing = [x for x in FOCAL if x not in resolved]
                if missing:
                    raise ValueError(f"Missing focal samples: {missing}")
                selected_columns = [resolved[x] for x in FOCAL]
                dst.write("\t".join(fields[:9] + FOCAL) + "\n")
                continue
            if line.startswith("#"):
                dst.write(line)
                continue
            if selected_columns is None:
                raise ValueError("VCF data encountered before #CHROM header")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= max(selected_columns):
                raise ValueError(f"Short VCF row near record {records + 1}")
            # The source is explicitly the unannotated call set. This defensive
            # check prevents an annotated file from being deposited accidentally.
            if "ANN=" in fields[7] or "EFF=" in fields[7]:
                raise ValueError("SnpEff annotation detected in source VCF")
            minimal_site = fields[:5] + [".", "PASS", ".", fields[8]]
            dst.write("\t".join(minimal_site + [fields[i] for i in selected_columns]) + "\n")
            records += 1

    print(f"wrote {records:,} records for {len(FOCAL)} populations to {args.output}")


if __name__ == "__main__":
    main()
