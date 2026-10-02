#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SAMPLES = ["AAI", "AAK", "AAL", "AAQ", "AAU", "AAW", "AAX", "ABE"]
ZENODO = "https://doi.org/10.5281/zenodo.23099934"


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_vcf(path, full):
    samples = None
    n_records = 0
    annotated = False
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
            elif not line.startswith("#"):
                n_records += 1
                annotated |= "ANN=" in line or "EFF=" in line
            if not full and n_records >= 1000:
                break
    assert samples == SAMPLES and not annotated
    assert not full or n_records == 6_444_023
    print(f"VCF check: OK ({path})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=Path, help="downloaded Zenodo VCF")
    parser.add_argument("--full-vcf", action="store_true")
    args = parser.parse_args()

    vcf = args.vcf
    if vcf is not None and not vcf.is_absolute():
        vcf = ROOT / vcf
    default_vcf = ROOT / "shared_data/8pop_poolseq.vcf.gz"
    if vcf is None and default_vcf.exists():
        vcf = default_vcf
    if vcf is not None:
        if not vcf.exists():
            parser.error(f"VCF not found: {vcf}")
        check_vcf(vcf, args.full_vcf)
    elif args.full_vcf:
        parser.error(f"download the VCF from {ZENODO} and provide it with --vcf")
    else:
        print(f"VCF check skipped; data are archived at {ZENODO}")

    counts = {
        row["selection_class"]: int(row["n_loci"])
        for row in rows(
            ROOT / "01_selective_sweep_signatures/results/selection_locus_counts.tsv"
        )
    }
    assert counts == {"farmed_only": 672, "farmed_and_wild": 774, "wild_only": 213}
    intro = rows(
        ROOT / "02_adaptive_introgression/results/overlap_permutation_test.tsv"
    )[0]
    assert int(intro["observed_introgression_regions_overlapping"]) == 68
    assert int(intro["observed_selection_loci_overlapping"]) == 96
    qc = rows(ROOT / "03_qtl_vs_genome_signatures/results/run_QC.tsv")
    assert qc and all(row["pass"].lower() == "true" for row in qc)
    print("Selection, introgression, and QTL checks: OK")


if __name__ == "__main__":
    main()
