#!/usr/bin/env python3
"""Dole7-specific QTL overlap tests conditional on callability."""
from __future__ import annotations

import os
import sys
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
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

import run_qtl_D_pi_fdm_traitwise_sensitivity as source
import run_qtl_overlap_sensitivity as core

OUT = HERE / "results"
GWAS_MARKERS = ROOT / "shared_data" / "qtl" / "gwas_marker_positions.csv.gz"
TAILS = [0.10, 0.05, 0.01]
TAIL_LABEL = {0.10: "10%", 0.05: "5%", 0.01: "1%"}
PERMUTATIONS = int(os.environ.get("N_PERMUTATIONS", "100000"))
SEED = 20260917
TRAITS = list(core.TRAIT_ORDER)
DISPLAY_TRAITS = [*TRAITS, "ALL"]
TRAIT_LABEL = {"ALL": "All QTLs", **core.TRAIT_LABEL}
ANALYSES = ["D", "pi", "fDM", "D_pi", "D_fDM", "pi_fDM", "D_pi_fDM"]
PAIRWISE = {"D", "pi", "fDM"}
LABEL = {"D": "D", "pi": "pi", "fDM": "fDM", "D_pi": "D-pi",
         "D_fDM": "D-fDM", "pi_fDM": "pi-fDM", "D_pi_fDM": "D-pi-fDM"}
DOLE7, WILD = "Dole7", "Wild Serbian"


def merge(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=["chromosome", "start", "end"])
    return core.merge_simple(frame[["chromosome", "start", "end"]].copy())


def intersect_many(*frames: pd.DataFrame) -> pd.DataFrame:
    if not frames:
        return pd.DataFrame(columns=["chromosome", "start", "end"])
    answer = merge(frames[0])
    for frame in frames[1:]:
        answer = source.intersect_interval_sets(answer, merge(frame))
        if answer.empty:
            break
    return merge(answer)


def interval_bp(frame: pd.DataFrame) -> int:
    z = merge(frame)
    return int((z.end - z.start + 1).sum()) if len(z) else 0


def metric_signal_and_mask(metric: str):
    raw = source.load_metric_windows(metric)
    farm = raw[raw.population == DOLE7].copy()
    wild = raw[raw.population == WILD].copy()
    callable_mask = intersect_many(merge(farm), merge(wild))
    thresholds, selected = [], {}
    for tail in TAILS:
        farm_threshold = float(farm.value.quantile(tail))
        wild_threshold = float(wild.value.quantile(tail))
        farm_low = merge(farm[farm.value <= farm_threshold])
        wild_low = merge(wild[wild.value <= wild_threshold])
        selected[tail] = intersect_many(
            source.subtract_interval_sets(intersect_many(farm_low, callable_mask), wild_low),
            callable_mask,
        )
        for population, threshold, group in (
            (DOLE7, farm_threshold, farm), (WILD, wild_threshold, wild)
        ):
            thresholds.append({
                "analysis_type": metric, "cutoff": tail,
                "cutoff_label": TAIL_LABEL[tail], "population": population,
                "lower_tail_threshold": threshold,
                "n_callable_windows": len(group),
                "n_low_tail_windows": int((group.value <= threshold).sum()),
            })
    return callable_mask, selected, pd.DataFrame(thresholds)


def fdm_signal_and_mask():
    raw = source.load_lifted_fdm_windows()
    raw = raw[raw.recipient_population == DOLE7].copy()
    comparisons = sorted(raw.comparison.unique())
    callable_by_comparison, selected_raw, threshold_rows = {}, {}, []
    for comparison in comparisons:
        group = raw[raw.comparison == comparison].copy()
        valid_all = group[(group.valid_window == True) & (group.nsnp >= 100)].copy()  # noqa
        valid_lifted = valid_all[
            (valid_all.liftover_status == "LIFTED") & valid_all.chromosome.notna()
            & valid_all.start.notna() & valid_all.end.notna()
        ].copy()
        callable_by_comparison[comparison] = merge(valid_lifted)
        for tail in TAILS:
            positive = valid_all.loc[valid_all.FdM > 0, "FdM"].dropna()
            threshold = float(positive.quantile(1 - tail))
            chosen = valid_lifted[(valid_lifted.D < 0) &
                                  (valid_lifted.FdM > threshold)].copy()
            selected_raw[(tail, comparison)] = merge(chosen)
            threshold_rows.append({
                "analysis_type": "fDM", "cutoff": tail,
                "cutoff_label": TAIL_LABEL[tail], "comparison": comparison,
                "recipient_population": DOLE7, "upper_tail_threshold": threshold,
                "n_valid_windows_all_coordinates": len(valid_all),
                "n_valid_lifted_windows": len(valid_lifted),
                "n_selected_lifted_windows": len(chosen),
            })
    common_mask = intersect_many(*[callable_by_comparison[x] for x in comparisons])
    selected, comparison_parts = {}, []
    for tail in TAILS:
        pieces = []
        for comparison in comparisons:
            z = intersect_many(selected_raw[(tail, comparison)], common_mask)
            pieces.append(z)
            if len(z):
                q = z.copy()
                q.insert(0, "comparison", comparison)
                q.insert(0, "cutoff", tail)
                comparison_parts.append(q)
        selected[tail] = intersect_many(pd.concat(pieces, ignore_index=True), common_mask)
    comparison_regions = pd.concat(comparison_parts, ignore_index=True)
    return (common_mask, selected, pd.DataFrame(threshold_rows),
            comparison_regions, callable_by_comparison)


def interval_index(frame: pd.DataFrame):
    answer = {}
    for chromosome, group in merge(frame).groupby("chromosome", sort=False):
        z = group.sort_values("start")
        starts, ends = z.start.to_numpy(np.int64), z.end.to_numpy(np.int64)
        answer[str(chromosome)] = (starts, ends, np.cumsum(ends - starts + 1))
    return answer


def covered_to(x, starts, ends, prefix):
    idx = np.searchsorted(starts, x, side="right") - 1
    safe = np.maximum(idx, 0)
    before = np.where(idx > 0, prefix[safe - 1], 0)
    partial = np.where(idx >= 0,
                       np.clip(x - starts[safe] + 1, 0, ends[safe] - starts[safe] + 1), 0)
    return before + partial


def covered_bp(index, chromosome, starts, ends):
    values = index.get(str(chromosome))
    if values is None:
        return np.zeros(len(starts), dtype=np.int64)
    ms, me, prefix = values
    return covered_to(ends, ms, me, prefix) - covered_to(starts - 1, ms, me, prefix)


def interval_hit(index, chromosome, starts, ends):
    values = index.get((str(chromosome),))
    if values is None:
        return np.zeros(len(starts), dtype=bool)
    ts, te = values
    j = np.searchsorted(ts, ends, side="right") - 1
    valid, safe = j >= 0, np.maximum(j, 0)
    return valid & (te[safe] >= starts)


def gwas_marker_universe():
    frame = pd.read_csv(GWAS_MARKERS, usecols=["CHR", "POS"])
    frame["chromosome"] = frame.CHR.map(core.chrom)
    frame["POS"] = pd.to_numeric(frame.POS, errors="raise").astype(int)
    return {str(c): np.sort(g.POS.drop_duplicates().to_numpy(np.int64))
            for c, g in frame.groupby("chromosome", sort=False)}


def build_catalog(qtl, callable_mask, markers, lengths, analysis_type):
    mask_index = interval_index(callable_mask)
    entries, rows = [], []
    for row in qtl.itertuples(index=False):
        chromosome = str(row.chromosome)
        start, end = int(row.start), int(row.end)
        width = end - start + 1
        callable_obs = int(covered_bp(mask_index, chromosome,
                                      np.array([start]), np.array([end]))[0])
        fraction_obs = callable_obs / width
        record = {
            "analysis_type": analysis_type, "QTL_ID": row.QTL_ID,
            "Trait": row.Trait, "chromosome": chromosome,
            "start": start, "end": end, "width_bp": width,
            "callable_bp": callable_obs, "callable_fraction": fraction_obs,
            "eligible": callable_obs > 0,
        }
        if callable_obs <= 0:
            record.update({"candidate_pool_size": 0, "matching_tolerance": np.nan,
                           "candidate_mean_callable_fraction": np.nan})
            rows.append(record)
            continue
        positions = markers.get(chromosome, np.array([], dtype=np.int64))
        candidate_starts = positions - ((width - 1) // 2)
        candidate_ends = candidate_starts + width - 1
        valid = (candidate_starts >= 1) & (candidate_ends <= int(lengths[chromosome]))
        candidate_starts, candidate_ends = candidate_starts[valid], candidate_ends[valid]
        fractions = covered_bp(mask_index, chromosome, candidate_starts, candidate_ends) / width
        differences = np.abs(fractions - fraction_obs)
        chosen, tolerance_used = None, np.nan
        for tolerance in (0.025, 0.05, 0.10, 0.20, 0.35):
            test = np.flatnonzero(differences <= tolerance)
            if len(test) >= 200:
                chosen, tolerance_used = test, tolerance
                break
        if chosen is None:
            chosen = np.argsort(differences)[:min(200, len(candidate_starts))]
            tolerance_used = float(differences[chosen].max()) if len(chosen) else np.nan
        pool = candidate_starts[chosen].astype(np.int64)
        if not len(pool):
            raise ValueError(f"No callable-matched marker placement for {row.QTL_ID}")
        entries.append({"row": row, "pool": pool, "width": width,
                        "observed_start": start, "observed_end": end})
        record.update({
            "candidate_pool_size": len(pool), "matching_tolerance": tolerance_used,
            "candidate_mean_callable_fraction": float(fractions[chosen].mean()),
        })
        rows.append(record)
    return entries, pd.DataFrame(rows)


def permutation_test(entries, target, analysis_type, cutoff):
    target_index = core.target_index(merge(target))
    n, all_index = len(DISPLAY_TRAITS), len(DISPLAY_TRAITS) - 1
    trait_index = {trait: i for i, trait in enumerate(TRAITS)}
    observed = np.zeros(n, dtype=np.int32)
    null = np.zeros((n, PERMUTATIONS), dtype=np.uint16)
    query_counts = np.zeros(n, dtype=np.int32)
    for entry in entries:
        row = entry["row"]
        chromosome, ti = str(row.chromosome), trait_index[str(row.Trait)]
        query_counts[[ti, all_index]] += 1
        observed_hit = int(interval_hit(
            target_index, chromosome, np.array([entry["observed_start"]]),
            np.array([entry["observed_end"]]))[0])
        observed[[ti, all_index]] += observed_hit
        seed_key = f"{analysis_type}|{row.QTL_ID}|callable_marker_placement"
        rng = np.random.default_rng(SEED + zlib.crc32(seed_key.encode("utf-8")))
        sampled_starts = entry["pool"][rng.integers(
            0, len(entry["pool"]), size=PERMUTATIONS
        )]
        hits = interval_hit(target_index, chromosome, sampled_starts,
                            sampled_starts + entry["width"] - 1).astype(np.uint16)
        null[ti] += hits
        null[all_index] += hits
    rows = []
    for i, trait in enumerate(DISPLAY_TRAITS):
        obs, values, count = int(observed[i]), null[i], int(query_counts[i])
        p_upper = (1 + int((values >= obs).sum())) / (PERMUTATIONS + 1)
        p_lower = (1 + int((values <= obs).sum())) / (PERMUTATIONS + 1)
        mean = float(values.mean())
        rows.append({
            "analysis_type": analysis_type, "cutoff": cutoff,
            "cutoff_label": TAIL_LABEL[cutoff], "Trait": trait,
            "n_query_loci_callable": count,
            "observed_query_loci_overlapping": obs,
            "percent_callable_query_loci_overlapping": 100 * obs / count if count else np.nan,
            "null_mean_query_loci_overlapping": mean,
            "null_sd_query_loci_overlapping": float(values.std(ddof=1)),
            "fold_enrichment": obs / mean if mean else np.nan,
            "empirical_p_enrichment": p_upper, "empirical_p_depletion": p_lower,
            "empirical_p_two_sided": min(1.0, 2 * min(p_upper, p_lower)),
            "n_permutations": PERMUTATIONS,
            "null_model": "same chromosome and width; GWAS-marker centred; callable-fraction matched",
        })
    return pd.DataFrame(rows)


def star(p):
    return "***" if p <= .001 else "**" if p <= .01 else "*" if p <= .05 else ""


def make_heatmap(results):
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "pdf.fonttype": 42})
    columns = [(kind, tail) for kind in ANALYSES for tail in TAILS]
    matrix = np.full((len(DISPLAY_TRAITS), len(columns)), np.nan)
    counts = [["" for _ in columns] for _ in DISPLAY_TRAITS]
    stars = [["" for _ in columns] for _ in DISPLAY_TRAITS]
    for i, trait in enumerate(DISPLAY_TRAITS):
        for j, (kind, tail) in enumerate(columns):
            row = results[(results.analysis_type == kind) & np.isclose(results.cutoff, tail)
                          & (results.Trait == trait)].iloc[0]
            direction = 1 if row.observed_query_loci_overlapping >= row.null_mean_query_loci_overlapping else -1
            matrix[i, j] = direction * -np.log10(max(
                float(row.empirical_p_two_sided), 1 / (PERMUTATIONS + 1)))
            counts[i][j] = str(int(row.observed_query_loci_overlapping))
            stars[i][j] = star(float(row.empirical_p_two_sided))
    fig = plt.figure(figsize=(18 / 2.54, 7 / 2.54))
    ax = fig.add_axes([.10, .20, .885, .60])
    image = ax.imshow(matrix, cmap="RdBu_r",
                      norm=TwoSlopeNorm(vmin=-5, vcenter=0, vmax=5), aspect="auto")
    ax.set_xticks(range(len(columns)), [TAIL_LABEL[t] for _, t in columns], fontsize=5.5)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(DISPLAY_TRAITS)),
                  [TRAIT_LABEL.get(x, x) for x in DISPLAY_TRAITS], fontsize=5.5)
    ax.tick_params(length=0, pad=4)
    for i in range(len(DISPLAY_TRAITS)):
        for j in range(len(columns)):
            colour = "white" if abs(matrix[i, j]) >= 2 else "#202020"
            ax.text(j, i, counts[i][j], ha="center", va="center", fontsize=5.7,
                    color=colour)
            if stars[i][j]:
                ax.text(j + .43, i - .29, stars[i][j], ha="right", va="top",
                        fontsize=7.0, fontweight="bold", color=colour)
    axes_width_mm, axes_height_mm = 180 * .885, 70 * .60
    x_gap = len(columns) / axes_width_mm
    y_gap = len(DISPLAY_TRAITS) / axes_height_mm
    col_bounds = [-.5, *[x - .5 for x in range(3, len(columns), 3)], len(columns) - .5]
    row_bounds = [-.5, 4.5, 11.5, 12.5, len(DISPLAY_TRAITS) - .5]
    for boundary in col_bounds[1:-1]:
        ax.axvspan(boundary - x_gap / 2, boundary + x_gap / 2,
                   color="white", linewidth=0, zorder=4)
    for boundary in row_bounds[1:-1]:
        ax.axhspan(boundary - y_gap / 2, boundary + y_gap / 2,
                   color="white", linewidth=0, zorder=4)
    for ci, (left, right) in enumerate(zip(col_bounds[:-1], col_bounds[1:])):
        x0, x1 = left + (x_gap / 2 if ci else 0), right - (x_gap / 2 if ci < 6 else 0)
        for ri, (top, bottom) in enumerate(zip(row_bounds[:-1], row_bounds[1:])):
            y0, y1 = top + (y_gap / 2 if ri else 0), bottom - (y_gap / 2 if ri < 3 else 0)
            ax.add_patch(mpl.patches.Rectangle((x0, y0), x1 - x0, y1 - y0,
                         fill=False, edgecolor="#242424", linewidth=.55, zorder=6))
    for centre, kind in zip(range(1, len(columns), 3), ANALYSES):
        ax.text(centre, 1.08, LABEL[kind], transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=5.8, fontweight="bold")
    cax = fig.add_axes([(1 - 2 / 18) / 2, .145, 2 / 18, .025])
    cbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    cbar.set_ticks([-5, 0, 5]); cbar.set_label("signed -log10(P)", fontsize=4.8)
    cbar.ax.tick_params(labelsize=4.8, length=2)
    fig.text(.10, .012, "Blue = depletion; red = enrichment. Cell = observed callable QTLs. Stars use raw empirical P: * <= 0.05, ** <= 0.01, *** <= 0.001.", fontsize=4.3)
    stem = OUT / "Dole7_callable_QTL_D_pi_fDM_overlap_heatmap"
    for suffix, kwargs in ((".png", {"dpi": 600}), (".pdf", {}), (".svg", {})):
        fig.savefig(stem.with_suffix(suffix), facecolor="white", **kwargs)
    plt.close(fig)


def write_report(results, callability, comparisons):
    all_qtl = results[results.Trait == "ALL"].sort_values(["analysis_type", "cutoff"])
    lines = [
        "# Dole7-specific, callability-aware QTL overlap analysis", "",
        "The primary analysis compared the 305 French-population QTLs with Dole7-specific "
        "low Tajima's D and low nucleotide diversity (pi), and with high-fDM windows from "
        "contrasts in which Dole7 is the mapped recipient. Missing windows were treated as "
        "non-testable rather than as absence of signal.", "",
        "For D and pi, the callable universe was the intersection of Dole7 and Wild Serbian "
        "callable windows because both values are required to classify a signal as farm-specific. "
        "For fDM, the primary universe was the intersection of valid, lifted windows across "
        f"the Dole7-recipient contrasts: {', '.join(comparisons)}. Multipartite masks were "
        "intersections of their component masks.", "",
        f"Each null distribution used {PERMUTATIONS:,} placements per QTL. Random QTLs retained "
        "chromosome and width, were centred on SNPs in the complete GWAS scan, and were matched "
        "to the observed QTL's callable fraction. QTLs with no callable base were excluded.", "",
        "Benjamini-Hochberg correction was applied across all displayed primary Dole7 tests; "
        "panel-specific q values are also provided in the result table.", "",
        "## All-QTL results", "",
    ]
    for row in all_qtl.itertuples(index=False):
        lines.append(
            f"- {LABEL[row.analysis_type]}, {row.cutoff_label}: "
            f"{int(row.observed_query_loci_overlapping)} of {int(row.n_query_loci_callable)} "
            f"callable QTLs overlapped (null mean {row.null_mean_query_loci_overlapping:.2f}; "
            f"{row.direction}; P={row.empirical_p_two_sided:.5g}; "
            f"global BH q={row.fdr_q_all_primary:.5g})."
        )
    lines += ["", "## Callability", ""]
    for row in callability[callability.Trait == "ALL"].itertuples(index=False):
        lines.append(
            f"- {LABEL[row.analysis_type]}: {int(row.n_callable_QTLs)} of "
            f"{int(row.n_total_QTLs)} QTLs callable; median callable fraction "
            f"{row.median_callable_fraction:.3f}; mask {row.callable_mask_bp / 1e6:.2f} Mb."
        )
    lines += ["", "Complete outputs include QTL-level callability, null-matching diagnostics, "
              "overlap pairs, empirical P values, q values, and QC checks.", ""]
    (OUT / "ANALYSIS_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source.TAILS, source.TAIL_LABEL = TAILS, TAIL_LABEL
    core.TAILS, core.TAIL_LABEL = TAILS, TAIL_LABEL
    qtl, marker_windows = source.load_qtls_and_marker_windows(True)
    lengths, markers = core.chromosome_lengths(core.GFF), gwas_marker_universe()
    d_mask, d_targets, d_thresholds = metric_signal_and_mask("D")
    pi_mask, pi_targets, pi_thresholds = metric_signal_and_mask("pi")
    (fdm_mask, fdm_targets, fdm_thresholds, fdm_comparison_regions,
     fdm_callable_by_comparison) = fdm_signal_and_mask()
    comparisons = sorted(fdm_callable_by_comparison)
    masks = {
        "D": d_mask, "pi": pi_mask, "fDM": fdm_mask,
        "D_pi": intersect_many(d_mask, pi_mask),
        "D_fDM": intersect_many(d_mask, fdm_mask),
        "pi_fDM": intersect_many(pi_mask, fdm_mask),
        "D_pi_fDM": intersect_many(d_mask, pi_mask, fdm_mask),
    }
    targets = {kind: {} for kind in ANALYSES}
    for tail in TAILS:
        targets["D"][tail] = intersect_many(d_targets[tail], masks["D"])
        targets["pi"][tail] = intersect_many(pi_targets[tail], masks["pi"])
        targets["fDM"][tail] = intersect_many(fdm_targets[tail], masks["fDM"])
        targets["D_pi"][tail] = intersect_many(d_targets[tail], pi_targets[tail], masks["D_pi"])
        targets["D_fDM"][tail] = intersect_many(d_targets[tail], fdm_targets[tail], masks["D_fDM"])
        targets["pi_fDM"][tail] = intersect_many(pi_targets[tail], fdm_targets[tail], masks["pi_fDM"])
        targets["D_pi_fDM"][tail] = intersect_many(
            d_targets[tail], pi_targets[tail], fdm_targets[tail], masks["D_pi_fDM"])

    mask_parts = []
    for kind, mask in masks.items():
        z = mask.copy(); z.insert(0, "analysis_type", kind); mask_parts.append(z)
    pd.concat(mask_parts, ignore_index=True).to_csv(
        OUT / "Dole7_analysis_callable_masks.tsv", sep="\t", index=False)
    for comparison, mask in fdm_callable_by_comparison.items():
        z = mask.copy(); z.insert(0, "comparison", comparison)
        z.to_csv(OUT / f"fDM_callable_mask_{source.safe(comparison)}.tsv", sep="\t", index=False)
    pd.concat([d_thresholds, pi_thresholds, fdm_thresholds], ignore_index=True, sort=False).to_csv(
        OUT / "Dole7_selection_thresholds.tsv", sep="\t", index=False)
    fdm_comparison_regions.to_csv(
        OUT / "Dole7_fDM_selected_regions_by_comparison.tsv", sep="\t", index=False)

    result_parts, qtl_parts, count_rows, pair_parts, qc_rows = [], [], [], [], []
    for kind in ANALYSES:
        catalog, qtl_callability = build_catalog(qtl, masks[kind], markers, lengths, kind)
        qtl_parts.append(qtl_callability)
        for tail in TAILS:
            target = targets[kind][tail]
            result_parts.append(permutation_test(catalog, target, kind, tail))
            count_rows.append({"analysis_type": kind, "cutoff": tail,
                               "cutoff_label": TAIL_LABEL[tail],
                               "n_regions": len(merge(target)),
                               "selected_bp": interval_bp(target),
                               "callable_mask_bp": interval_bp(masks[kind])})
            pairs = core.interval_pairs(qtl, target)
            if len(pairs):
                pairs.insert(0, "analysis_type", kind); pairs.insert(1, "cutoff", tail)
                pairs.insert(2, "cutoff_label", TAIL_LABEL[tail]); pair_parts.append(pairs)
            outside_bp = interval_bp(source.subtract_interval_sets(target, masks[kind]))
            qc_rows.append({"check": "target within callable mask", "analysis_type": kind,
                            "cutoff": tail, "Trait": "ALL", "value": outside_bp,
                            "expected": 0, "pass": outside_bp == 0})

    results = pd.concat(result_parts, ignore_index=True)
    results["direction"] = np.where(
        results.observed_query_loci_overlapping >= results.null_mean_query_loci_overlapping,
        "enrichment", "depletion")
    results["log2_observed_expected"] = np.log2(
        (results.observed_query_loci_overlapping + .5) /
        (results.null_mean_query_loci_overlapping + .5))
    results["fdr_q_all_primary"] = core.bh(results.empirical_p_two_sided)
    results["panel"] = np.where(results.analysis_type.isin(PAIRWISE), "pairwise", "multipartite")
    results["fdr_q_panel"] = results.groupby("panel", group_keys=False)[
        "empirical_p_two_sided"].transform(core.bh)
    results.to_csv(OUT / "Dole7_callable_traitwise_permutation_tests.tsv", sep="\t", index=False)

    qtl_callability = pd.concat(qtl_parts, ignore_index=True)
    qtl_callability.to_csv(OUT / "Dole7_QTL_callability_and_null_matching.tsv", sep="\t", index=False)
    summary_rows = []
    for kind, group in qtl_callability.groupby("analysis_type", sort=False):
        for trait in DISPLAY_TRAITS:
            z = group if trait == "ALL" else group[group.Trait == trait]
            eligible = z[z.eligible]
            summary_rows.append({
                "analysis_type": kind, "Trait": trait, "n_total_QTLs": len(z),
                "n_callable_QTLs": len(eligible),
                "median_callable_fraction": eligible.callable_fraction.median(),
                "mean_callable_fraction": eligible.callable_fraction.mean(),
                "min_candidate_pool_size": int(eligible.candidate_pool_size.min()) if len(eligible) else 0,
                "median_candidate_pool_size": eligible.candidate_pool_size.median() if len(eligible) else 0,
                "callable_mask_bp": interval_bp(masks[kind]),
            })
    callability_summary = pd.DataFrame(summary_rows)
    callability_summary.to_csv(OUT / "Dole7_callability_summary.tsv", sep="\t", index=False)
    pd.DataFrame(count_rows).to_csv(OUT / "Dole7_selected_region_counts.tsv", sep="\t", index=False)
    all_pairs = pd.concat(pair_parts, ignore_index=True) if pair_parts else pd.DataFrame()
    all_pairs.to_csv(OUT / "Dole7_QTL_overlap_pairs.tsv", sep="\t", index=False)
    qtl.to_csv(OUT / "QTL_intervals_used.tsv", sep="\t", index=False)
    marker_windows.to_csv(OUT / "significant_SNP_50kb_windows.tsv", sep="\t", index=False)

    for row in results.itertuples(index=False):
        z = all_pairs[(all_pairs.analysis_type == row.analysis_type) &
                      np.isclose(all_pairs.cutoff, row.cutoff)] if len(all_pairs) else pd.DataFrame()
        if row.Trait != "ALL" and len(z):
            z = z[z.left_Trait == row.Trait]
        eligible = qtl_callability[(qtl_callability.analysis_type == row.analysis_type) &
                                   qtl_callability.eligible]
        if row.Trait != "ALL": eligible = eligible[eligible.Trait == row.Trait]
        actual = len(set(z.left_QTL_ID) & set(eligible.QTL_ID)) if len(z) else 0
        qc_rows.append({"check": "observed count reconciliation", "analysis_type": row.analysis_type,
                        "cutoff": row.cutoff, "Trait": row.Trait, "value": actual,
                        "expected": int(row.observed_query_loci_overlapping),
                        "pass": actual == int(row.observed_query_loci_overlapping)})
    qc = pd.DataFrame(qc_rows); qc.to_csv(OUT / "run_QC.tsv", sep="\t", index=False)
    if not qc["pass"].all():
        raise ValueError(f"QC failed in {(~qc['pass']).sum()} rows")
    if os.environ.get("SKIP_PLOTS", "0") != "1":
        make_heatmap(results)
    write_report(results, callability_summary, comparisons)
    print("Dole7 recipient comparisons:", "; ".join(comparisons))
    print("QC failures:", int((~qc["pass"]).sum()))
    print(results[results.Trait == "ALL"][["analysis_type", "cutoff_label",
          "n_query_loci_callable", "observed_query_loci_overlapping",
          "null_mean_query_loci_overlapping", "direction",
          "empirical_p_two_sided", "fdr_q_all_primary"]].to_string(index=False))


if __name__ == "__main__":
    main()
