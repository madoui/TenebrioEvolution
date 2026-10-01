#!/usr/bin/env python3
"""Trait-wise QTL overlap sensitivity analysis for D, pi, and fDM.

Lower-tail Tajima's D and nucleotide-diversity windows and upper-tail positive
fDM windows are merged into loci at 10%, 5%, 1%, and 0.1% cutoffs.  Pairwise
QTL overlaps and exact shared-base multipartite overlaps are compared with a
chromosome- and QTL-width-preserving random-repositioning null.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOCAL = HERE / ".python_packages"
if LOCAL.exists():
    sys.path.insert(0, str(LOCAL))
CACHE = HERE / ".matplotlib_cache"
CACHE.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(CACHE))

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

import run_qtl_overlap_sensitivity as core
from cluster_low_diversity import read_metric


ROOT = HERE.parents[0]
DATA = ROOT / "shared_data"
DEFAULT_OUT = HERE / "results_traitwise"
FARM_SPECIFIC_OUT = HERE / "results_farm_specific"
MARKER50_FARM_SPECIFIC_OUT = (
    HERE / "results_marker50_farm_specific"
)
OUT = DEFAULT_OUT
FARM_SPECIFIC_MODE = False
MARKER50_MODE = False
PI_FILE = DATA / "selection" / "pi_windows.tsv.gz"
LIFTED_FDM_FILE = DATA / "introgression" / "lifted_fdm_windows.csv.gz"
SIGNIFICANT_MARKER_FILE = (
    DATA / "qtl" / "significant_markers.csv"
)
TAILS = [0.10, 0.05, 0.01, 0.001]
TAIL_LABEL = {0.10: "10%", 0.05: "5%", 0.01: "1%", 0.001: "0.1%"}
PERMUTATIONS = int(os.environ.get("N_PERMUTATIONS", "100000"))
SEED = 20260910
DISPLAY_TRAITS = [*core.TRAIT_ORDER, "ALL"]
TRAIT_LABEL = {"ALL": "All QTLs", **core.TRAIT_LABEL}
PAIRWISE_SIGNALS = ["D", "pi", "fDM"]
MULTIPARTITE_TYPES = ["D_pi", "D_fDM", "pi_fDM", "D_pi_fDM"]
MULTIPARTITE_LABEL = {
    "D_pi": "QTL-D-pi",
    "D_fDM": "QTL-D-fDM",
    "pi_fDM": "QTL-pi-fDM",
    "D_pi_fDM": "QTL-D-pi-fDM",
}
FARMED_POPULATIONS = [
    "Dole7",
    "Starfood Holland",
    "Kingsect Belgium",
    "Pronutrix+ Belgium",
    "Papek Czechia",
    "Mystik Canada",
]
WILD_POPULATION = "Wild Serbian"
RAW_TRAIT_MAP = {
    "EHR_i": "EHR",
    "TELC_i": "TELC",
    "DELR2_i": "DELR2",
    "DELR4_i": "DELR4",
    "DELR5_i": "DELR5",
    "PW_M2_4quantiles": "PW",
    "DTEP_cat": "DTEP categorical",
}


def safe(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")


def load_qtls_and_marker_windows(marker_50kb: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load annotated QTLs and, when requested, independently verify their marker build."""
    qtl = pd.read_csv(core.QTL_FILE).rename(columns={"CHR": "chromosome"})
    qtl["chromosome"] = qtl.chromosome.map(core.chrom)
    qtl[["start", "end"]] = qtl[["start", "end"]].astype(int)
    if qtl.QTL_ID.nunique() != 305 or qtl.QTL_ID.duplicated().any():
        raise ValueError("Expected 305 unique QTL intervals")
    if set(qtl.Trait) != set(core.TRAIT_ORDER):
        raise ValueError("Unexpected QTL trait set")
    if not marker_50kb:
        return qtl, pd.DataFrame()

    markers = pd.read_csv(SIGNIFICANT_MARKER_FILE).copy()
    markers["Trait_raw"] = markers["Trait"]
    markers["Trait"] = markers.Trait.replace(RAW_TRAIT_MAP)
    markers["chromosome"] = markers.CHR.map(core.chrom)
    markers["POS"] = pd.to_numeric(markers.POS, errors="raise").astype(int)
    markers["window_start"] = (markers.POS - 50_000).clip(lower=1).astype(int)
    markers["window_end"] = (markers.POS + 50_000).astype(int)

    # Independently rebuild the merged coordinates from the significant-marker table.
    built_rows: list[dict[str, object]] = []
    for (trait, chromosome), group in markers.groupby(["Trait", "chromosome"], sort=True):
        current_start = current_end = None
        n_markers = 0
        for row in group.sort_values(["window_start", "window_end", "SNP"]).itertuples(index=False):
            start, end = int(row.window_start), int(row.window_end)
            if current_end is None or start > current_end + 1:
                if current_end is not None:
                    built_rows.append({
                        "Trait": trait, "chromosome": chromosome,
                        "start": current_start, "end": current_end,
                        "n_sig_SNPs": n_markers,
                    })
                current_start, current_end, n_markers = start, end, 1
            else:
                current_end = max(current_end, end)
                n_markers += 1
        if current_end is not None:
            built_rows.append({
                "Trait": trait, "chromosome": chromosome,
                "start": current_start, "end": current_end,
                "n_sig_SNPs": n_markers,
            })
    built = pd.DataFrame(built_rows)
    key = ["Trait", "chromosome", "start", "end", "n_sig_SNPs"]
    if len(built) != len(qtl) or set(map(tuple, built[key].to_numpy())) != set(
        map(tuple, qtl[key].to_numpy())
    ):
        raise ValueError("Fresh ±50-kb marker-window build does not match the annotated QTL table")

    # Assign every marker window to its one verified merged, trait-specific QTL.
    membership = markers.merge(
        qtl[["QTL_ID", "Trait", "chromosome", "start", "end"]],
        on=["Trait", "chromosome"], how="left",
    )
    membership = membership[
        (membership.POS >= membership.start) & (membership.POS <= membership.end)
    ].copy()
    counts = membership.groupby(["Trait", "SNP"], dropna=False).QTL_ID.nunique()
    if len(counts) != len(markers) or not counts.eq(1).all():
        raise ValueError("A significant trait–SNP record did not map to exactly one merged QTL")
    marker_columns = [
        "QTL_ID", "SNP", "Trait_raw", "Trait", "family", "chromosome", "POS",
        "window_start", "window_end", "n_models_significant", "methods_found",
        "Genes", "Consequence", "Amino acid substitution", "Function (DMEL)",
        "Function (TCAS)", "InterPro",
    ]
    return qtl, membership[marker_columns].sort_values(
        ["Trait", "chromosome", "POS", "SNP"]
    ).reset_index(drop=True)


def q_text(value: float) -> str:
    if value < 0.001:
        return "<.001"
    if value < 0.01:
        return f"{value:.3f}"
    return f"{value:.2f}"


def p_text(value: float) -> str:
    """Compactly show the empirical P value at its useful precision."""
    if value < 0.001:
        return f"{value:.5f}"
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.3f}"


def stars(value: float) -> str:
    if value <= 0.001:
        return "***"
    if value <= 0.01:
        return "**"
    if value <= 0.05:
        return "*"
    return ""


def make_pi_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = read_metric(PI_FILE, "pi").rename(
        columns={"chrom": "chromosome", "window_start": "start", "window_end": "end"}
    )
    raw["chromosome"] = raw.chromosome.map(core.chrom)
    raw = raw.dropna(subset=["pi", "pop", "chromosome", "start", "end"]).copy()
    raw[["start", "end"]] = raw[["start", "end"]].astype(int)
    threshold_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    for tail in TAILS:
        for population, group in raw.groupby("pop", sort=True):
            threshold = float(group.pi.quantile(tail))
            selected = group[group.pi <= threshold].copy()
            threshold_rows.append({
                "pi_tail": tail,
                "pi_tail_label": TAIL_LABEL[tail],
                "population": population,
                "pi_threshold": threshold,
                "n_callable_windows": len(group),
                "n_outlier_windows": len(selected),
                "actual_outlier_fraction": len(selected) / len(group),
            })
            for chromosome, cg in selected.groupby("chromosome", sort=False):
                block = None
                for row in cg.sort_values(["start", "end"]).itertuples(index=False):
                    start, end = int(row.start), int(row.end)
                    if block is not None and start <= block["end"] + 1:
                        block["end"] = max(block["end"], end)
                        block["n_10kb_windows"] += 1
                        block["min_pi"] = min(block["min_pi"], float(row.pi))
                    else:
                        if block is not None:
                            region_rows.append(block)
                        block = {
                            "pi_tail": tail,
                            "pi_tail_label": TAIL_LABEL[tail],
                            "population": str(population),
                            "chromosome": str(chromosome),
                            "start": start,
                            "end": end,
                            "n_10kb_windows": 1,
                            "min_pi": float(row.pi),
                        }
                if block is not None:
                    region_rows.append(block)
    regions = pd.DataFrame(region_rows)
    regions["pi_region_id"] = [
        f"PI{TAIL_LABEL[t].replace('%', 'pct')}_{safe(p)}_{c}_R{i + 1:05d}"
        for i, (t, p, c) in enumerate(zip(
            regions.pi_tail, regions.population, regions.chromosome
        ))
    ]
    return pd.DataFrame(threshold_rows), regions


def load_metric_windows(metric: str) -> pd.DataFrame:
    """Load D or pi as aligned 10-kb population windows."""
    if metric == "D":
        raw = pd.read_csv(core.D_FILE, sep="\t", na_values=["NA", "NaN"]).reset_index(drop=True)
        raw["chromosome"] = raw["chr"].map(core.chrom)
        raw["start"] = pd.to_numeric(raw["pos"], errors="coerce")
        raw["end"] = pd.to_numeric(raw["pos2"], errors="coerce")
        raw["value"] = pd.to_numeric(raw["D"], errors="coerce")
    elif metric == "pi":
        raw = read_metric(PI_FILE, "pi").rename(
            columns={"chrom": "chromosome", "window_start": "start", "window_end": "end"}
        )
        raw["chromosome"] = raw.chromosome.map(core.chrom)
        raw["value"] = pd.to_numeric(raw["pi"], errors="coerce")
    else:
        raise ValueError(metric)
    raw = raw.rename(columns={"pop": "population"})
    raw = raw.dropna(subset=["population", "chromosome", "start", "end", "value"]).copy()
    raw[["start", "end"]] = raw[["start", "end"]].astype(int)
    return raw[["population", "chromosome", "start", "end", "value"]]


def intersect_interval_sets(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """Return exact shared bases between two chromosome-indexed interval sets."""
    rows = []
    for chromosome in sorted(set(left.chromosome.astype(str)) & set(right.chromosome.astype(str))):
        a = core.merge_simple(left[left.chromosome.astype(str) == chromosome])
        b = core.merge_simple(right[right.chromosome.astype(str) == chromosome])
        av = list(a[["start", "end"]].itertuples(index=False, name=None))
        bv = list(b[["start", "end"]].itertuples(index=False, name=None))
        i = j = 0
        while i < len(av) and j < len(bv):
            a_start, a_end = map(int, av[i])
            b_start, b_end = map(int, bv[j])
            start, end = max(a_start, b_start), min(a_end, b_end)
            if start <= end:
                rows.append({"chromosome": chromosome, "start": start, "end": end})
            if a_end < b_end:
                i += 1
            elif b_end < a_end:
                j += 1
            else:
                i += 1
                j += 1
    return pd.DataFrame(rows, columns=["chromosome", "start", "end"])


def subtract_interval_sets(base: pd.DataFrame, exclusion: pd.DataFrame) -> pd.DataFrame:
    """Subtract exclusion bases from base intervals, chromosome by chromosome."""
    rows = []
    exclusion_by_chrom = {
        str(chromosome): core.merge_simple(group).sort_values(["start", "end"])
        for chromosome, group in exclusion.groupby("chromosome", sort=False)
    }
    for chromosome, group in core.merge_simple(base).groupby("chromosome", sort=False):
        excluded = exclusion_by_chrom.get(str(chromosome), pd.DataFrame(columns=["start", "end"]))
        ex = list(excluded[["start", "end"]].itertuples(index=False, name=None))
        for start, end in group[["start", "end"]].itertuples(index=False, name=None):
            segments = [(int(start), int(end))]
            for ex_start, ex_end in ex:
                ex_start, ex_end = int(ex_start), int(ex_end)
                if ex_end < start:
                    continue
                if ex_start > end:
                    break
                next_segments = []
                for seg_start, seg_end in segments:
                    if ex_end < seg_start or ex_start > seg_end:
                        next_segments.append((seg_start, seg_end))
                    else:
                        if seg_start < ex_start:
                            next_segments.append((seg_start, ex_start - 1))
                        if ex_end < seg_end:
                            next_segments.append((ex_end + 1, seg_end))
                segments = next_segments
                if not segments:
                    break
            rows.extend({"chromosome": str(chromosome), "start": s, "end": e}
                        for s, e in segments if s <= e)
    if not rows:
        return pd.DataFrame(columns=["chromosome", "start", "end"])
    return core.merge_simple(pd.DataFrame(rows))


def make_farm_specific_metric_regions(metric: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select low-tail farm windows with callable, non-low wild counterparts."""
    raw = load_metric_windows(metric)
    available_farms = [x for x in FARMED_POPULATIONS if x in set(raw.population)]
    required = set(available_farms + [WILD_POPULATION])
    raw = raw[raw.population.isin(required)].copy()
    if WILD_POPULATION not in set(raw.population):
        raise ValueError(f"{metric}: missing {WILD_POPULATION}")
    threshold_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    for tail in TAILS:
        thresholds = raw.groupby("population").value.quantile(tail).to_dict()
        for population, group in raw.groupby("population", sort=True):
            threshold = float(thresholds[population])
            threshold_rows.append({
                "metric": metric,
                "cutoff": tail,
                "cutoff_label": TAIL_LABEL[tail],
                "population": population,
                "population_type": "wild" if population == WILD_POPULATION else "farmed",
                "lower_tail_threshold": threshold,
                "n_callable_windows": len(group),
                "n_low_tail_windows": int((group.value <= threshold).sum()),
                "included_as_farm": population in available_farms,
            })
        wild = raw[raw.population == WILD_POPULATION].copy()
        wild_callable = core.merge_simple(wild[["chromosome", "start", "end"]])
        wild_low = core.merge_simple(
            wild[wild.value <= float(thresholds[WILD_POPULATION])][["chromosome", "start", "end"]]
        )
        for farm in available_farms:
            farm_low = raw[(raw.population == farm) &
                           (raw.value <= float(thresholds[farm]))].copy()
            farm_low_union = core.merge_simple(farm_low[["chromosome", "start", "end"]])
            callable_farm_low = intersect_interval_sets(farm_low_union, wild_callable)
            selected = subtract_interval_sets(callable_farm_low, wild_low)
            for row in selected.itertuples(index=False):
                block = {
                    "population": farm,
                    "chromosome": str(row.chromosome),
                    "start": int(row.start),
                    "end": int(row.end),
                    "width_bp": int(row.end) - int(row.start) + 1,
                    "wild_threshold": float(thresholds[WILD_POPULATION]),
                    "farm_specific_definition": "farm_low_bases_with_wild_callable_low_bases_subtracted",
                }
                if metric == "D":
                    block["selection_tail"] = tail
                    block["selection_tail_label"] = TAIL_LABEL[tail]
                else:
                    block["pi_tail"] = tail
                    block["pi_tail_label"] = TAIL_LABEL[tail]
                region_rows.append(block)
    regions = pd.DataFrame(region_rows)
    tail_col = "selection_tail" if metric == "D" else "pi_tail"
    id_col = "selection_region_id" if metric == "D" else "pi_region_id"
    regions[id_col] = [
        f"FARM_{metric.upper()}{TAIL_LABEL[t].replace('%', 'pct')}_{safe(p)}_{c}_R{i + 1:05d}"
        for i, (t, p, c) in enumerate(zip(
            regions[tail_col], regions.population, regions.chromosome
        ))
    ]
    return pd.DataFrame(threshold_rows), regions


def load_lifted_fdm_windows() -> pd.DataFrame:
    """Parse the alternating source/liftover records in the full fDM file."""
    rows: list[dict[str, object]] = []
    comparison_counts: dict[str, int] = {}
    opener = gzip.open if LIFTED_FDM_FILE.suffix == ".gz" else open
    with opener(LIFTED_FDM_FILE, "rt", encoding="utf-8", errors="replace") as handle:
        names_source = next(csv.reader([next(handle).rstrip("\n")]))
        names_lifted = next(csv.reader([next(handle).rstrip("\n")]))[1:]
        while True:
            source_line = handle.readline()
            lifted_line = handle.readline()
            if not source_line or not lifted_line:
                break
            source_values = next(csv.reader([source_line.rstrip("\n")]))
            lifted_values = next(csv.reader([lifted_line.rstrip("\n")]))[1:]
            if len(source_values) != len(names_source) or len(lifted_values) != len(names_lifted):
                continue
            row = dict(zip(names_source, source_values))
            row.update(dict(zip(names_lifted, lifted_values)))
            if row.get("window_bp") not in {"1e+05", "100000"}:
                continue
            comparison = str(row["comparison"])
            comparison_counts[comparison] = comparison_counts.get(comparison, 0) + 1
            lifted = row.get("liftover_status") == "LIFTED"
            midpoint = int(float(row["LiftedMidPos"])) if lifted else None
            rows.append({
                "window_id": f"{comparison}_100kb_{comparison_counts[comparison]:08d}",
                "comparison": comparison,
                "P2": str(row["P2"]),
                "recipient_population": core.RECIPIENT_MAP.get(str(row["P2"]), str(row["P2"])),
                "chromosome": core.chrom(row["LiftedChr"]) if lifted else None,
                "start": max(1, midpoint - 50_000) if lifted else None,
                "end": midpoint + 50_000 if lifted else None,
                "nsnp": int(float(row["nsnp"])),
                "D": float(row["D"]),
                "FdM": float(row["FdM"]),
                "valid_window": str(row["valid_window"]).upper() == "TRUE",
                "liftover_status": str(row.get("liftover_status", "")),
            })
    return pd.DataFrame(rows)


def make_fdm_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = load_lifted_fdm_windows()
    raw = raw[(raw.valid_window == True) & (raw.nsnp >= 100)].copy()  # noqa: E712
    threshold_rows: list[dict[str, object]] = []
    region_rows: list[dict[str, object]] = []
    for tail in TAILS:
        for comparison, group in raw.groupby("comparison", sort=True):
            positive = group.loc[group.FdM > 0, "FdM"].dropna()
            threshold = float(positive.quantile(1 - tail))
            # The established candidate definition uses a strict upper-tail
            # comparison (FdM > quantile), which also resolves tied thresholds.
            selected_all = group[(group.D < 0) & (group.FdM > threshold)].copy()
            selected = selected_all[selected_all.liftover_status == "LIFTED"].copy()
            threshold_rows.append({
                "fdm_tail": tail,
                "fdm_tail_label": TAIL_LABEL[tail],
                "comparison": comparison,
                "fdm_threshold": threshold,
                "n_valid_windows": len(group),
                "n_positive_fdm_windows": len(positive),
                "n_selected_windows_after_D_filter": len(selected_all),
                "n_selected_windows_with_lifted_coordinates": len(selected),
            })
            for chromosome, cg in selected.groupby("chromosome", sort=False):
                block = None
                block_ids: list[str] = []
                for row in cg.sort_values(["start", "end"]).itertuples(index=False):
                    start, end = int(row.start), int(row.end)
                    if block is not None and start <= block["end"] + 1:
                        block["end"] = max(block["end"], end)
                        block["n_100kb_windows"] += 1
                        block["max_fdm"] = max(block["max_fdm"], float(row.FdM))
                        block_ids.append(str(row.window_id))
                    else:
                        if block is not None:
                            block["window_ids"] = ",".join(block_ids)
                            region_rows.append(block)
                        block = {
                            "fdm_tail": tail,
                            "fdm_tail_label": TAIL_LABEL[tail],
                            "comparison": comparison,
                            "P2": str(row.P2),
                            "recipient_population": str(row.recipient_population),
                            "chromosome": str(chromosome),
                            "start": start,
                            "end": end,
                            "n_100kb_windows": 1,
                            "max_fdm": float(row.FdM),
                        }
                        block_ids = [str(row.window_id)]
                if block is not None:
                    block["window_ids"] = ",".join(block_ids)
                    region_rows.append(block)
    regions = pd.DataFrame(region_rows)
    regions["introgression_region_id"] = [
        f"FDM{TAIL_LABEL[t].replace('%', 'pct')}_{safe(comp)}_{c}_R{i + 1:05d}"
        for i, (t, comp, c) in enumerate(zip(
            regions.fdm_tail, regions.comparison, regions.chromosome
        ))
    ]
    return pd.DataFrame(threshold_rows), regions


def matched_intersection(
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_population: str,
    right_population: str,
) -> pd.DataFrame:
    """Return exact shared bases, matching population and chromosome."""
    rows: list[dict[str, object]] = []
    left_pops = set(left[left_population].astype(str)) if not left.empty else set()
    right_pops = set(right[right_population].astype(str)) if not right.empty else set()
    for population in sorted(left_pops & right_pops):
        a_pop = left[left[left_population].astype(str) == population]
        b_pop = right[right[right_population].astype(str) == population]
        for chromosome in sorted(set(a_pop.chromosome.astype(str)) & set(b_pop.chromosome.astype(str))):
            a = a_pop[a_pop.chromosome.astype(str) == chromosome].sort_values(["start", "end"])
            b = b_pop[b_pop.chromosome.astype(str) == chromosome].sort_values(["start", "end"])
            av = list(a[["start", "end"]].itertuples(index=False, name=None))
            bv = list(b[["start", "end"]].itertuples(index=False, name=None))
            i = j = 0
            while i < len(av) and j < len(bv):
                a_start, a_end = map(int, av[i])
                b_start, b_end = map(int, bv[j])
                start, end = max(a_start, b_start), min(a_end, b_end)
                if start <= end:
                    rows.append({
                        "population": population,
                        "chromosome": chromosome,
                        "start": start,
                        "end": end,
                    })
                if a_end < b_end:
                    i += 1
                elif b_end < a_end:
                    j += 1
                else:
                    i += 1
                    j += 1
    if not rows:
        return pd.DataFrame(columns=["population", "chromosome", "start", "end"])
    raw = pd.DataFrame(rows)
    merged_rows = []
    for (population, chromosome), group in raw.groupby(["population", "chromosome"], sort=False):
        merged = core.merge_simple(group)
        merged["population"] = population
        merged_rows.append(merged[["population", "chromosome", "start", "end"]])
    return pd.concat(merged_rows, ignore_index=True)


def permutation_results_by_trait(
    qtl: pd.DataFrame,
    target: pd.DataFrame,
    lengths: dict[str, int],
    key: str,
) -> list[dict[str, float]]:
    """Test all trait strata in one set of random QTL placements."""
    target = core.merge_simple(target)
    index = core.target_index(target)
    n_rows = len(DISPLAY_TRAITS)
    trait_index = {trait: i for i, trait in enumerate(DISPLAY_TRAITS)}
    observed = np.zeros(n_rows, dtype=np.int32)
    null = np.zeros((n_rows, PERMUTATIONS), dtype=np.uint16)
    rng = np.random.default_rng(SEED + zlib.crc32(key.encode("utf-8")))
    for row in qtl.itertuples(index=False):
        chromosome = str(row.chromosome)
        start, end = int(row.start), int(row.end)
        width = end - start + 1
        max_start = int(lengths[chromosome]) - width + 1
        if max_start < 1:
            raise ValueError(f"QTL wider than chromosome: {chromosome}:{start}-{end}")
        target_starts, target_ends = index.get(
            (chromosome,),
            (np.array([], dtype=np.int64), np.array([], dtype=np.int64)),
        )
        if len(target_starts):
            observed_j = np.searchsorted(target_starts, end, side="right") - 1
            is_observed = int(observed_j >= 0 and target_ends[observed_j] >= start)
            random_starts = rng.integers(1, max_start + 1, size=PERMUTATIONS)
            random_ends = random_starts + width - 1
            j = np.searchsorted(target_starts, random_ends, side="right") - 1
            valid = j >= 0
            jj = np.maximum(j, 0)
            hits = (valid & (target_ends[jj] >= random_starts)).astype(np.uint16)
        else:
            is_observed = 0
            hits = np.zeros(PERMUTATIONS, dtype=np.uint16)
        ti = trait_index[str(row.Trait)]
        observed[0] += is_observed
        observed[ti] += is_observed
        null[0] += hits
        null[ti] += hits

    rows: list[dict[str, float]] = []
    qtl_counts = qtl.Trait.value_counts().to_dict()
    for i, trait in enumerate(DISPLAY_TRAITS):
        obs = int(observed[i])
        values = null[i]
        p_upper = (1 + int((values >= obs).sum())) / (PERMUTATIONS + 1)
        p_lower = (1 + int((values <= obs).sum())) / (PERMUTATIONS + 1)
        mean = float(values.mean())
        rows.append({
            "Trait": trait,
            "n_query_loci": len(qtl) if trait == "ALL" else int(qtl_counts.get(trait, 0)),
            "observed_query_loci_overlapping": obs,
            "percent_query_loci_overlapping": 100 * obs / (len(qtl) if trait == "ALL" else qtl_counts[trait]),
            "null_mean_query_loci_overlapping": mean,
            "null_sd_query_loci_overlapping": float(values.std(ddof=1)),
            "fold_enrichment": obs / mean if mean > 0 else np.nan,
            "empirical_p_enrichment": p_upper,
            "empirical_p_depletion": p_lower,
            "empirical_p_two_sided": min(1.0, 2 * min(p_upper, p_lower)),
            "n_permutations": PERMUTATIONS,
        })
    return rows


def add_statistics(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["fdr_q_enrichment"] = core.bh(frame.empirical_p_enrichment)
    frame["fdr_q_depletion"] = core.bh(frame.empirical_p_depletion)
    frame["fdr_q_two_sided"] = core.bh(frame.empirical_p_two_sided)
    frame["log2_observed_expected"] = np.log2(
        (frame.observed_query_loci_overlapping + 0.5) /
        (frame.null_mean_query_loci_overlapping + 0.5)
    )
    frame["direction"] = np.where(
        frame.observed_query_loci_overlapping >= frame.null_mean_query_loci_overlapping,
        "enrichment", "depletion",
    )
    return frame


def count_intervals_not_contained(strict: pd.DataFrame, relaxed: pd.DataFrame) -> int:
    """Count strict union intervals not wholly contained in a relaxed union interval."""
    index = core.target_index(core.merge_simple(relaxed))
    failures = 0
    for row in core.merge_simple(strict).itertuples(index=False):
        starts, ends = index.get(
            (str(row.chromosome),),
            (np.array([], dtype=np.int64), np.array([], dtype=np.int64)),
        )
        j = np.searchsorted(starts, int(row.start), side="right") - 1
        failures += int(j < 0 or ends[j] < int(row.end))
    return failures


def heatmap_panel(ax: plt.Axes, frame: pd.DataFrame, columns: list[tuple[str, float]],
                  column_labels: list[str], title: str, norm: TwoSlopeNorm):
    matrix = np.full((len(DISPLAY_TRAITS), len(columns)), np.nan)
    labels = [["" for _ in columns] for _ in DISPLAY_TRAITS]
    significance_labels = [["" for _ in columns] for _ in DISPLAY_TRAITS]
    for i, trait in enumerate(DISPLAY_TRAITS):
        for j, (kind, tail) in enumerate(columns):
            row = frame[(frame.analysis_type == kind) & np.isclose(frame.cutoff, tail)].loc[
                lambda x: x.Trait == trait
            ].iloc[0]
            direction = 1.0 if row.observed_query_loci_overlapping >= row.null_mean_query_loci_overlapping else -1.0
            matrix[i, j] = direction * -np.log10(max(float(row.empirical_p_two_sided), 1 / (PERMUTATIONS + 1)))
            labels[i][j] = f"{int(row.observed_query_loci_overlapping)}"
            significance_labels[i][j] = stars(float(row.fdr_q_two_sided))
    image = ax.imshow(matrix, cmap="RdBu_r", norm=norm, aspect="auto")
    column_fontsize = 5.2 if len(columns) > 15 else 6.2
    cell_fontsize = 5.5 if len(columns) > 15 else 6.5
    ax.set_xticks(range(len(columns)), column_labels, fontsize=column_fontsize)
    ax.xaxis.tick_top()
    y_fontsize = 5.5 if len(columns) > 15 else 7
    ax.set_yticks(range(len(DISPLAY_TRAITS)), [TRAIT_LABEL.get(x, x) for x in DISPLAY_TRAITS], fontsize=y_fontsize)
    ax.tick_params(length=0, pad=4)
    for i in range(len(DISPLAY_TRAITS)):
        for j in range(len(columns)):
            colour = "white" if abs(matrix[i, j]) >= 2.0 else "#202020"
            x_count = j - (0.11 if significance_labels[i][j] else 0.0)
            ax.text(x_count, i, labels[i][j], ha="center", va="center", fontsize=cell_fontsize,
                    color=colour, linespacing=1.0)
            if significance_labels[i][j]:
                ax.text(j + 0.16, i - 0.03, significance_labels[i][j], ha="center", va="center",
                        fontsize=7.5, fontweight="bold", color=colour, linespacing=1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=9, fontweight="bold", pad=20)
    return image


def make_heatmap(pairwise: pd.DataFrame, multipartite: pd.DataFrame) -> None:
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "pdf.fonttype": 42})
    norm = TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5)
    pair_columns = [(signal, tail) for signal in PAIRWISE_SIGNALS for tail in TAILS]
    pair_labels = [f"{signal}\n{TAIL_LABEL[tail]}" for signal, tail in pair_columns]
    multi_columns = [(kind, tail) for kind in MULTIPARTITE_TYPES for tail in TAILS]
    multi_labels = [f"{MULTIPARTITE_LABEL[kind].replace('QTL-', '')}\n{TAIL_LABEL[tail]}"
                    for kind, tail in multi_columns]

    # A single continuous heatmap at a fixed journal-column width of 18 cm.
    all_columns = [*pair_columns, *multi_columns]
    all_labels = [TAIL_LABEL[tail] for _, tail in all_columns]
    all_results = pd.concat([pairwise, multipartite], ignore_index=True)
    fig = plt.figure(figsize=(18.0 / 2.54, 7.0 / 2.54))
    gs = fig.add_gridspec(1, 1, left=0.10, right=0.985, bottom=0.20, top=0.80)
    ax = fig.add_subplot(gs[0, 0])
    image = heatmap_panel(
        ax, all_results, all_columns, all_labels,
        "", norm,
    )
    # True 1-mm gutters separate statistic-by-trait-family blocks.
    axes_width_mm = 180.0 * (0.985 - 0.10)
    axes_height_mm = 70.0 * (0.80 - 0.20)
    x_gap = len(all_columns) / axes_width_mm
    y_gap = len(DISPLAY_TRAITS) / axes_height_mm
    column_bounds = [-0.5, *[x - 0.5 for x in range(len(TAILS), len(all_columns), len(TAILS))],
                     len(all_columns) - 0.5]
    row_bounds = [-0.5, 4.5, 11.5, 12.5, len(DISPLAY_TRAITS) - 0.5]
    for boundary in column_bounds[1:-1]:
        ax.axvspan(boundary - x_gap / 2, boundary + x_gap / 2,
                   color="white", linewidth=0, zorder=4)
    for boundary in row_bounds[1:-1]:
        ax.axhspan(boundary - y_gap / 2, boundary + y_gap / 2,
                   color="white", linewidth=0, zorder=4)
    for column_index, (left, right) in enumerate(zip(column_bounds[:-1], column_bounds[1:])):
        x0 = left + (x_gap / 2 if column_index else 0)
        x1 = right - (x_gap / 2 if column_index < len(column_bounds) - 2 else 0)
        for row_index, (top, bottom) in enumerate(zip(row_bounds[:-1], row_bounds[1:])):
            y0 = top + (y_gap / 2 if row_index else 0)
            y1 = bottom - (y_gap / 2 if row_index < len(row_bounds) - 2 else 0)
            ax.add_patch(mpl.patches.Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor="#242424", linewidth=0.55,
                clip_on=False, zorder=6,
            ))
    ax.text(4.0, 1.22, "PAIRWISE", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#444444")
    ax.text(14.5, 1.22, "MULTIPARTITE", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#444444")
    for centre, label in zip(
        (1, 4, 7, 10, 13, 16, 19),
        ("D", "pi", "fDM", "D-pi", "D-fDM", "pi-fDM", "D-pi-fDM"),
    ):
        ax.text(centre, 1.08, label, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=5.8, fontweight="bold")
    legend_width = 2.0 / 18.0
    cax = fig.add_axes([(1.0 - legend_width) / 2.0, 0.12, legend_width, 0.025])
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    cbar.set_ticks([-5, 0, 5])
    cbar.set_label("signed -log10(P)", fontsize=4.8)
    cbar.ax.tick_params(labelsize=4.8, length=2)
    fig.text(0.10, 0.018,
             "Blue = depletion; red = enrichment. * q<=0.05, ** q<=0.01, *** q<=0.001. pi denotes nucleotide diversity.",
             fontsize=4.3, ha="left")
    for suffix, kwargs in [(".png", {"dpi": 600}), (".pdf", {}), (".svg", {})]:
        fig.savefig((OUT / "QTL_D_pi_fDM_overlap_sensitivity_heatmap").with_suffix(suffix),
                    facecolor="white", **kwargs)
    plt.close(fig)


def write_report(pairwise: pd.DataFrame, multipartite: pd.DataFrame,
                 region_counts: pd.DataFrame, qtl_count: int,
                 marker_count: int | None = None) -> None:
    if MARKER50_MODE and FARM_SPECIFIC_MODE:
        scope_title = "50-kb marker-window QTL overlap with farm-specific low D and pi"
    elif FARM_SPECIFIC_MODE:
        scope_title = "Trait-wise QTL overlap with farm-specific low D and pi"
    else:
        scope_title = "Trait-wise QTL overlap sensitivity across D, pi, and fDM thresholds"
    cutoff_text = ", ".join(TAIL_LABEL[tail] for tail in TAILS)
    qtl_method = (
        f"This analysis rebuilt {qtl_count} QTL intervals from {marker_count} significant "
        "trait-SNP records by expanding every SNP by 50 kb on each side and merging "
        "overlapping or adjacent windows only within the same trait and chromosome. "
        if MARKER50_MODE else
        f"This analysis used {qtl_count} LD-defined QTLs. "
    )
    lines = [
        f"# {scope_title}",
        "",
        qtl_method + f"Each test used {PERMUTATIONS:,} chromosome- and QTL-width-preserving "
        "random placements per target set.",
        "",
        "Tajima's D and nucleotide diversity (pi) were tested at their population-specific "
        f"lower {cutoff_text} tails. fDM was tested at the corresponding comparison-specific upper tails "
        "of positive fDM among valid 100-kb windows with at least 100 SNPs and D < 0. "
        "Adjacent or overlapping outlier windows were merged before analysis. fDM thresholds "
        "were estimated from all valid windows; selected windows without lifted coordinates "
        "were recorded in the threshold table but could not enter physical-overlap tests.",
        "",
        "Multipartite tests require exact shared genomic bases. D-pi intersections match "
        "population identity. Intersections involving fDM match the D or pi population to the "
        "fDM recipient population. The four-way test requires a QTL, D, pi, and fDM to share "
        "at least one base. Cutoff levels are matched within each displayed multipartite test.",
        "",
        "Benjamini-Hochberg correction was applied separately across all cells in the pairwise "
        "panel and all cells in the multipartite panel, including the all-QTL summary rows.",
        "",
        "## Region counts",
        "",
    ]
    if FARM_SPECIFIC_MODE:
        lines[8:8] = [
            "Farm-specific D and pi windows were required to fall in the lower tail in at least "
            "one focal farm while the coordinate-matched Wild Serbian window was callable and "
            "above its lower-tail threshold. USA and non-focal populations neither define nor "
            "disqualify farm specificity. Pronutrix+ Belgium was included for D but could not be "
            "included for pi because it is absent from the pi input.",
            "",
        ]
    if MARKER50_MODE:
        lines[8:8] = [
            "The rebuilt QTL coordinates were independently checked against the annotated "
            "50-kb QTL table: all marker records mapped exactly once and all merged coordinates "
            "and marker counts agreed.",
            "",
        ]
    for row in region_counts.itertuples(index=False):
        lines.append(
            f"- {row.analysis_type}, {row.cutoff_label}: {int(row.n_regions):,} merged regions; "
            f"{int(row.n_union_regions):,} genomic-union regions spanning {row.union_bp / 1e6:.2f} Mb."
        )
    lines += ["", "## FDR-significant departures from chance", ""]
    significant = pd.concat([
        pairwise.assign(panel="pairwise"),
        multipartite.assign(panel="multipartite"),
    ], ignore_index=True)
    significant = significant[significant.fdr_q_two_sided <= 0.05]
    if significant.empty:
        lines.append("No displayed test differed from chance at 5% FDR.")
    else:
        for row in significant.sort_values(["panel", "analysis_type", "cutoff", "Trait"]).itertuples(index=False):
            lines.append(
                f"- {row.panel}, {row.analysis_type}, {row.cutoff_label}, "
                f"{TRAIT_LABEL.get(row.Trait, row.Trait)}: observed "
                f"{int(row.observed_query_loci_overlapping)}, expected "
                f"{row.null_mean_query_loci_overlapping:.2f}, {row.direction} fold "
                f"{row.fold_enrichment:.2f}, P={row.empirical_p_two_sided:.4g}, "
                f"BH q={row.fdr_q_two_sided:.4g}."
            )
    lines += [
        "",
        "Complete observed counts, null distributions, P values, and q values are provided in "
        "the pairwise and multipartite result tables. Heatmap colour is signed -log10 of the "
        "two-sided empirical P value (blue depletion, red enrichment). Cells display that P "
        "value, while stars are assigned from BH q (* q<=0.05, ** q<=0.01, *** q<=0.001).",
        "",
    ]
    (OUT / "ANALYSIS_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main(farm_specific: bool = False, marker_50kb: bool = False) -> None:
    global OUT, FARM_SPECIFIC_MODE, MARKER50_MODE, TAILS, TAIL_LABEL
    FARM_SPECIFIC_MODE = farm_specific
    MARKER50_MODE = marker_50kb
    if marker_50kb:
        TAILS = [0.10, 0.05, 0.01]
        TAIL_LABEL = {0.10: "10%", 0.05: "5%", 0.01: "1%"}
    else:
        TAILS = [0.10, 0.05, 0.01, 0.001]
        TAIL_LABEL = {0.10: "10%", 0.05: "5%", 0.01: "1%", 0.001: "0.1%"}
    OUT = (MARKER50_FARM_SPECIFIC_OUT if marker_50kb and farm_specific else
           FARM_SPECIFIC_OUT if farm_specific else DEFAULT_OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    core.TAILS = TAILS
    core.TAIL_LABEL = TAIL_LABEL

    qtl, marker_windows = load_qtls_and_marker_windows(marker_50kb)

    if farm_specific:
        d_thresholds, d_regions = make_farm_specific_metric_regions("D")
        pi_thresholds, pi_regions = make_farm_specific_metric_regions("pi")
    else:
        d_thresholds, d_regions = core.make_selection_regions()
        pi_thresholds, pi_regions = make_pi_regions()
    fdm_thresholds, fdm_regions = make_fdm_regions()
    d_thresholds.to_csv(OUT / "TajimaD_population_thresholds.tsv", sep="\t", index=False)
    pi_thresholds.to_csv(OUT / "pi_population_thresholds.tsv", sep="\t", index=False)
    fdm_thresholds.to_csv(OUT / "fDM_comparison_thresholds.tsv", sep="\t", index=False)
    d_regions.to_csv(OUT / "TajimaD_merged_regions.tsv", sep="\t", index=False)
    pi_regions.to_csv(OUT / "pi_merged_regions.tsv", sep="\t", index=False)
    fdm_regions.to_csv(OUT / "fDM_merged_regions.tsv", sep="\t", index=False)
    qtl.to_csv(OUT / "QTL_intervals_used.tsv", sep="\t", index=False)
    if marker_50kb:
        marker_windows.to_csv(OUT / "significant_SNP_50kb_windows.tsv", sep="\t", index=False)

    lengths = core.chromosome_lengths(core.GFF)
    for frame in (qtl, d_regions, pi_regions, fdm_regions):
        for chromosome, end in frame.groupby("chromosome").end.max().items():
            lengths[str(chromosome)] = max(lengths.get(str(chromosome), 0), int(end))

    pairwise_rows: list[dict[str, object]] = []
    multipartite_rows: list[dict[str, object]] = []
    shared_parts: list[pd.DataFrame] = []
    pair_parts: list[pd.DataFrame] = []
    count_rows: list[dict[str, object]] = []
    target_cache: dict[tuple[str, float], pd.DataFrame] = {}

    for tail in TAILS:
        signal_regions = {
            "D": d_regions[np.isclose(d_regions.selection_tail, tail)].copy(),
            "pi": pi_regions[np.isclose(pi_regions.pi_tail, tail)].copy(),
            "fDM": fdm_regions[np.isclose(fdm_regions.fdm_tail, tail)].copy(),
        }
        for signal, regions in signal_regions.items():
            union = core.merge_simple(regions)
            target_cache[(signal, tail)] = union
            count_rows.append({
                "analysis_type": signal,
                "cutoff": tail,
                "cutoff_label": TAIL_LABEL[tail],
                "n_regions": len(regions),
                "n_union_regions": len(union),
                "union_bp": int((union.end - union.start + 1).sum()) if len(union) else 0,
            })
            for result in permutation_results_by_trait(qtl, union, lengths, f"pair|{signal}|{tail}"):
                pairwise_rows.append({
                    "analysis_type": signal,
                    "cutoff": tail,
                    "cutoff_label": TAIL_LABEL[tail],
                    **result,
                })
            pairs = core.interval_pairs(qtl, union)
            if len(pairs):
                pairs.insert(0, "analysis_type", signal)
                pairs.insert(1, "cutoff", tail)
                pairs.insert(2, "cutoff_label", TAIL_LABEL[tail])
                pair_parts.append(pairs)

        d = signal_regions["D"]
        pi = signal_regions["pi"]
        fdm = signal_regions["fDM"]
        d_pi = matched_intersection(d, pi, "population", "population")
        d_fdm = matched_intersection(fdm, d, "recipient_population", "population")
        pi_fdm = matched_intersection(fdm, pi, "recipient_population", "population")
        d_pi_fdm = matched_intersection(fdm, d_pi, "recipient_population", "population")
        shared_sets = {
            "D_pi": d_pi,
            "D_fDM": d_fdm,
            "pi_fDM": pi_fdm,
            "D_pi_fDM": d_pi_fdm,
        }
        for kind, shared in shared_sets.items():
            if len(shared):
                z = shared.copy()
                z.insert(0, "analysis_type", kind)
                z.insert(1, "cutoff", tail)
                z.insert(2, "cutoff_label", TAIL_LABEL[tail])
                z["shared_region_id"] = [
                    f"{kind}_{TAIL_LABEL[tail].replace('%', 'pct')}_{i + 1:06d}"
                    for i in range(len(z))
                ]
                shared_parts.append(z)
            union = core.merge_simple(shared)
            target_cache[(kind, tail)] = union
            count_rows.append({
                "analysis_type": kind,
                "cutoff": tail,
                "cutoff_label": TAIL_LABEL[tail],
                "n_regions": len(shared),
                "n_union_regions": len(union),
                "union_bp": int((union.end - union.start + 1).sum()) if len(union) else 0,
            })
            for result in permutation_results_by_trait(qtl, union, lengths, f"multi|{kind}|{tail}"):
                multipartite_rows.append({
                    "analysis_type": kind,
                    "cutoff": tail,
                    "cutoff_label": TAIL_LABEL[tail],
                    **result,
                })
            pairs = core.interval_pairs(qtl, union)
            if len(pairs):
                pairs.insert(0, "analysis_type", kind)
                pairs.insert(1, "cutoff", tail)
                pairs.insert(2, "cutoff_label", TAIL_LABEL[tail])
                pair_parts.append(pairs)

    pairwise = add_statistics(pd.DataFrame(pairwise_rows))
    multipartite = add_statistics(pd.DataFrame(multipartite_rows))
    pairwise.to_csv(OUT / "pairwise_traitwise_permutation_tests.tsv", sep="\t", index=False)
    multipartite.to_csv(OUT / "multipartite_traitwise_permutation_tests.tsv", sep="\t", index=False)
    pd.concat([pairwise.assign(panel="pairwise"), multipartite.assign(panel="multipartite")],
              ignore_index=True).to_csv(OUT / "all_traitwise_permutation_tests.tsv", sep="\t", index=False)
    region_counts = pd.DataFrame(count_rows)
    region_counts.to_csv(OUT / "merged_region_counts.tsv", sep="\t", index=False)
    if shared_parts:
        pd.concat(shared_parts, ignore_index=True).to_csv(
            OUT / "exact_shared_regions.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(OUT / "exact_shared_regions.tsv", sep="\t", index=False)
    if pair_parts:
        all_pairs = pd.concat(pair_parts, ignore_index=True)
        all_pairs.to_csv(OUT / "QTL_overlap_pairs_all_targets.tsv", sep="\t", index=False)
    else:
        all_pairs = pd.DataFrame()
        all_pairs.to_csv(OUT / "QTL_overlap_pairs_all_targets.tsv", sep="\t", index=False)

    # Reconcile every observed count against the detailed interval-pair table.
    qc_rows = []
    all_results = pd.concat([pairwise, multipartite], ignore_index=True)
    for row in all_results.itertuples(index=False):
        if all_pairs.empty:
            actual = 0
        else:
            z = all_pairs[(all_pairs.analysis_type == row.analysis_type) &
                          np.isclose(all_pairs.cutoff, row.cutoff)]
            if row.Trait != "ALL":
                z = z[z.left_Trait == row.Trait]
            actual = int(z.left_QTL_ID.nunique())
        expected = int(row.observed_query_loci_overlapping)
        qc_rows.append({
            "check": "observed overlap reconciliation",
            "analysis_type": row.analysis_type,
            "cutoff": row.cutoff,
            "Trait": row.Trait,
            "permutation_observed": expected,
            "pair_table_distinct_QTLs": actual,
            "pass": expected == actual,
        })
    for name, frame in (("QTL", qtl), ("D", d_regions), ("pi", pi_regions), ("fDM", fdm_regions)):
        normalized = bool(frame.chromosome.astype(str).str.match(
            r"^(?:[0-9]+|Tmol_scf_[0-9]+)$"
        ).all())
        coordinates_valid = bool(((frame.start >= 1) & (frame.end >= frame.start)).all())
        qc_rows.append({
            "check": "normalized chromosome labels and valid coordinates",
            "analysis_type": name,
            "cutoff": np.nan,
            "Trait": "ALL",
            "permutation_observed": np.nan,
            "pair_table_distinct_QTLs": np.nan,
            "pass": normalized and coordinates_valid,
        })
    # A farm-vs-wild contrast need not be nested across percentiles: relaxing
    # the cutoff can add newly wild-low bases to the exclusion set. Only fDM
    # retains a mathematically required nesting relationship in this mode.
    nesting_kinds = ["fDM"] if FARM_SPECIFIC_MODE else [
        *PAIRWISE_SIGNALS, *MULTIPARTITE_TYPES
    ]
    for kind in nesting_kinds:
        for strict, relaxed in zip(TAILS[1:], TAILS[:-1]):
            failures = count_intervals_not_contained(
                target_cache[(kind, strict)], target_cache[(kind, relaxed)]
            )
            qc_rows.append({
                "check": "strict cutoff nested within relaxed cutoff",
                "analysis_type": kind,
                "cutoff": strict,
                "Trait": "ALL",
                "permutation_observed": failures,
                "pair_table_distinct_QTLs": 0,
                "pass": failures == 0,
            })
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "run_QC.tsv", sep="\t", index=False)
    if not qc["pass"].all():
        raise ValueError("Observed overlap reconciliation failed")

    make_heatmap(pairwise, multipartite)
    write_report(
        pairwise, multipartite, region_counts, qtl.QTL_ID.nunique(),
        len(marker_windows) if marker_50kb else None,
    )
    significant = pd.concat([
        pairwise.assign(panel="pairwise"),
        multipartite.assign(panel="multipartite"),
    ], ignore_index=True)
    significant = significant[significant.fdr_q_two_sided <= 0.05]
    print(f"pairwise_tests={len(pairwise)} multipartite_tests={len(multipartite)} "
          f"significant={len(significant)} qc_failures={(~qc['pass']).sum()}")
    if len(significant):
        print(significant[["panel", "analysis_type", "cutoff_label", "Trait",
                           "observed_query_loci_overlapping",
                           "null_mean_query_loci_overlapping", "fold_enrichment",
                           "empirical_p_two_sided", "fdr_q_two_sided"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--farm-specific",
        action="store_true",
        help="restrict low-D and low-pi loci to farm-low, callable non-low wild windows",
    )
    parser.add_argument(
        "--marker-50kb",
        action="store_true",
        help="rebuild QTLs as merged ±50-kb windows around every significant trait-SNP record",
    )
    args = parser.parse_args()
    main(farm_specific=args.farm_specific, marker_50kb=args.marker_50kb)
