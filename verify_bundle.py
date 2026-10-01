#!/usr/bin/env python3
"""Lightweight structural and result checks for the public bundle."""
from __future__ import annotations
import argparse, csv, gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXPECTED_SAMPLES = ["AAI", "AAK", "AAL", "AAQ", "AAU", "AAW", "AAX", "ABE"]


def check_vcf(full: bool) -> None:
    path = ROOT / "shared_data/8pop_poolseq.vcf.gz"
    samples = None; records = 0; annotated = False
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
            elif not line.startswith("#"):
                records += 1
                annotated |= "ANN=" in line or "EFF=" in line
                if not full and records >= 1000: break
    assert samples == EXPECTED_SAMPLES, samples
    assert not annotated, "SnpEff annotation found"
    if full: assert records == 6_444_023, records
    print(f"VCF OK: samples={len(samples)}, records={'6,444,023' if full else 'first 1,000 checked'}")


def tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--full-vcf", action="store_true")
    args = ap.parse_args(); check_vcf(args.full_vcf)
    counts = {r["selection_class"]: int(r["n_loci"])
              for r in tsv(ROOT / "01_selective_sweep_signatures/results/selection_locus_counts.tsv")}
    assert counts == {"farmed_only": 672, "farmed_and_wild": 774, "wild_only": 213}, counts
    intro = tsv(ROOT / "02_adaptive_introgression/results/overlap_permutation_test.tsv")[0]
    assert int(intro["observed_introgression_regions_overlapping"]) == 68
    assert int(intro["observed_selection_loci_overlapping"]) == 96
    qc = tsv(ROOT / "03_qtl_vs_genome_signatures/results/run_QC.tsv")
    assert qc and all(r["pass"].lower() == "true" for r in qc)
    print("Selection counts, adaptive-introgression result, and QTL QC: OK")


if __name__ == "__main__": main()
