#!/usr/bin/env python3
"""Test overlap between FdM-5% regions and farm-only low-D5% loci."""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def overlaps(query: pd.DataFrame, target: pd.DataFrame) -> tuple[int, pd.DataFrame]:
    rows, hit_ids = [], set()
    by_chr = {str(c): g for c, g in target.groupby("chromosome", sort=False)}
    for q in query.itertuples(index=False):
        g = by_chr.get(str(q.chromosome))
        if g is None: continue
        hits = g[(g.start <= q.end) & (g.end >= q.start)]
        if len(hits): hit_ids.add(q.introgression_region_id)
        for t in hits.itertuples(index=False):
            s, e = max(q.start, t.start), min(q.end, t.end)
            rows.append({"introgression_region_id": q.introgression_region_id,
                         "selection_region_id": t.selection_region_id,
                         "chromosome": q.chromosome, "overlap_start": s,
                         "overlap_end": e, "overlap_bp": e - s + 1})
    return len(hit_ids), pd.DataFrame(rows)


def main() -> None:
    here, root = Path(__file__).resolve().parent, Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--permutations", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=20260907)
    ap.add_argument("--outdir", type=Path, default=here / "results")
    args = ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    selection_file = root / "01_selective_sweep_signatures/results/overlapaware_selection_loci_5pct.tsv"
    if not selection_file.exists():
        raise SystemExit("Run ../01_selective_sweep_signatures/run_analysis.py first")
    selection = pd.read_csv(selection_file, sep="\t", dtype={"chromosome": str})
    selection = selection[selection.selection_class.eq("farmed_only")].copy()
    intro = pd.read_csv(root / "shared_data/introgression/fdm_merged_regions.tsv", sep="\t", dtype={"chromosome": str})
    intro = intro[np.isclose(intro.fdm_tail, 0.05)].copy()
    lengths = pd.read_csv(root / "shared_data/genome/chromosome_lengths.tsv", sep="\t", dtype={"chromosome": str})
    length_map = dict(zip(lengths.chromosome.str.replace("chr", "", regex=False), lengths.length.astype(int)))
    if "X" in length_map:
        length_map["10"] = length_map["X"]
    observed, pairs = overlaps(intro, selection)
    pairs.to_csv(args.outdir / "introgression_selection_overlap_pairs.tsv", sep="\t", index=False)

    rng = np.random.default_rng(args.seed)
    null = np.zeros(args.permutations, dtype=np.int32)
    target_by_chr = {str(c): g[["start", "end"]].to_numpy(int) for c, g in selection.groupby("chromosome")}
    for q in intro.itertuples(index=False):
        chrom, width = str(q.chromosome), int(q.end - q.start + 1)
        max_start = length_map[chrom] - width + 1
        starts = rng.integers(1, max_start + 1, size=args.permutations)
        ends = starts + width - 1
        hit = np.zeros(args.permutations, dtype=bool)
        for s, e in target_by_chr.get(chrom, np.empty((0, 2), int)):
            hit |= (starts <= e) & (ends >= s)
        null += hit
    expected = float(null.mean())
    p_enrichment = (1 + int((null >= observed).sum())) / (args.permutations + 1)
    result = pd.DataFrame([{
        "n_introgression_regions": len(intro), "n_farm_only_selection_loci": len(selection),
        "observed_introgression_regions_overlapping": observed,
        "observed_selection_loci_overlapping": pairs.selection_region_id.nunique(),
        "null_mean_introgression_regions_overlapping": expected,
        "fold_enrichment": observed / expected if expected else np.nan,
        "empirical_p_enrichment": p_enrichment, "n_permutations": args.permutations,
        "seed": args.seed,
    }])
    result.to_csv(args.outdir / "overlap_permutation_test.tsv", sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
