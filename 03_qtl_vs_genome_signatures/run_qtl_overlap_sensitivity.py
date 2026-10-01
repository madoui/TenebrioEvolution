#!/usr/bin/env python3
"""Trait-wise physical-overlap sensitivity analysis for all merged QTLs.

Selection is defined from population-specific lower-tail Tajima's D windows.
Introgression is defined from comparison-specific upper-tail positive fDM windows,
retaining the established D < 0 direction requirement. Directly adjacent or
overlapping windows are merged before any overlap analysis. Coordinates are
treated as 1-based closed intervals.
"""
from __future__ import annotations

import os
import re
import sys
import zlib
from collections import defaultdict
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


ROOT = HERE.parents[0]
DATA = ROOT / "shared_data"
OUT = HERE / "results_core"
OUT.mkdir(parents=True, exist_ok=True)
QTL_FILE = DATA / "qtl" / "qtl_intervals.csv"
D_FILE = DATA / "selection" / "tajima_D_windows.tsv.gz"
INTRO_META = DATA / "Introgression" / "results_100kb" / "100kb_valid_windows_liftover_metadata.tsv"
INTRO_COORDS = DATA / "Introgression" / "introgression_100kb_candidate_windows.tsv"
GFF = DATA / "genome" / "chromosome_lengths.tsv"
TAILS = [0.05, 0.01, 0.001]
TAIL_LABEL = {0.05: "5%", 0.01: "1%", 0.001: "0.1%"}
PERMUTATIONS = int(os.environ.get("N_PERMUTATIONS", "10000"))
SEED = 20260907
RECIPIENT_MAP = {"Mystik Canada": "Mystik Canada", "Ynsect7 France": "Dole7"}
TRAIT_ORDER = ["EHR", "TELC", "DELR2", "DELR4", "DELR5", "IMM21", "IMM35",
               "IMM49", "IMM56", "IMM63", "AW", "PW", "DTEP categorical"]
TRAIT_LABEL = {"DTEP categorical": "DTEP"}


def chrom(value: object) -> str:
    x = str(value).strip()
    if x.lower().startswith("chr"):
        x = x[3:]
    return "10" if x == "X" else x


def safe(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")


def chromosome_lengths(path: Path) -> dict[str, int]:
    if path.suffix == ".tsv":
        table = pd.read_csv(path, sep="\t")
        return {chrom(c): int(n) for c, n in zip(table.chromosome, table.length)}
    ans: defaultdict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) == 9:
                ans[chrom(f[0])] = max(ans[chrom(f[0])], int(f[4]))
    return dict(ans)


def merge_simple(intervals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if intervals.empty:
        return pd.DataFrame(columns=["chromosome", "start", "end"])
    for c, group in intervals.groupby("chromosome", sort=False):
        current = None
        for r in group.sort_values(["start", "end"]).itertuples(index=False):
            s, e = int(r.start), int(r.end)
            if current is not None and s <= current[1] + 1:
                current[1] = max(current[1], e)
            else:
                if current is not None:
                    rows.append({"chromosome": str(c), "start": current[0], "end": current[1]})
                current = [s, e]
        if current is not None:
            rows.append({"chromosome": str(c), "start": current[0], "end": current[1]})
    return pd.DataFrame(rows)


def make_selection_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(D_FILE, sep="\t", na_values=["NA", "NaN"])
    raw["chromosome"] = raw["chr"].map(chrom)
    raw["start"] = pd.to_numeric(raw.pos, errors="coerce")
    raw["end"] = pd.to_numeric(raw.pos2, errors="coerce")
    raw["D"] = pd.to_numeric(raw.D, errors="coerce")
    raw = raw.dropna(subset=["start", "end", "D", "pop"]).copy()
    raw[["start", "end"]] = raw[["start", "end"]].astype(int)
    threshold_rows, region_rows = [], []
    for tail in TAILS:
        for pop, group in raw.groupby("pop", sort=True):
            threshold = float(group.D.quantile(tail))
            selected = group[group.D <= threshold].copy()
            threshold_rows.append({
                "selection_tail": tail, "selection_tail_label": TAIL_LABEL[tail],
                "population": pop, "D_threshold": threshold,
                "n_callable_windows": len(group), "n_outlier_windows": len(selected),
            })
            for c, cg in selected.groupby("chromosome", sort=False):
                block = None
                for r in cg.sort_values(["start", "end"]).itertuples(index=False):
                    if block is not None and int(r.start) <= block[1] + 1:
                        block[1] = max(block[1], int(r.end))
                        block[2] += 1
                        block[3] = min(block[3], float(r.D))
                    else:
                        if block is not None:
                            region_rows.append((tail, pop, str(c), *block))
                        block = [int(r.start), int(r.end), 1, float(r.D)]
                if block is not None:
                    region_rows.append((tail, pop, str(c), *block))
    regions = pd.DataFrame(region_rows, columns=[
        "selection_tail", "population", "chromosome", "start", "end", "n_10kb_windows", "min_D"])
    regions["selection_tail_label"] = regions.selection_tail.map(TAIL_LABEL)
    regions["selection_region_id"] = [
        f"D{TAIL_LABEL[t].replace('%','pct')}_{safe(p)}_{c}_R{i + 1:05d}"
        for i, (t, p, c) in enumerate(zip(regions.selection_tail, regions.population, regions.chromosome))]
    return pd.DataFrame(threshold_rows), regions


def make_introgression_regions() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(INTRO_META, sep="\t")
    meta = meta[(meta.valid_window == True) & (meta.nsnp >= 100)].copy()  # noqa: E712
    coords = pd.read_csv(INTRO_COORDS, sep="\t").drop_duplicates("window_id").copy()
    coords["chromosome"] = coords.chromosome.map(chrom)
    coords[["start", "end"]] = coords[["start", "end"]].astype(int)
    threshold_rows, region_rows = [], []
    for tail in TAILS:
        for comparison, group in meta.groupby("comparison", sort=True):
            positive = group.loc[group.FdM > 0, "FdM"].dropna()
            threshold = float(positive.quantile(1 - tail))
            threshold_rows.append({
                "fdm_tail": tail, "fdm_tail_label": TAIL_LABEL[tail],
                "comparison": comparison, "fdm_threshold": threshold,
                "n_valid_windows": len(group), "n_positive_fdm_windows": len(positive),
            })
            ids = set(group.loc[(group.D < 0) & (group.FdM > 0) & (group.FdM >= threshold), "window_id"])
            selected = coords[(coords.comparison == comparison) & coords.window_id.isin(ids)].copy()
            p2 = str(group.P2.iloc[0])
            recipient = RECIPIENT_MAP.get(p2, p2)
            for c, cg in selected.groupby("chromosome", sort=False):
                block = None
                block_ids = []
                for r in cg.sort_values(["start", "end"]).itertuples(index=False):
                    if block is not None and int(r.start) <= block[1] + 1:
                        block[1] = max(block[1], int(r.end))
                        block[2] += 1
                        block[3] = max(block[3], float(r.FdM))
                        block_ids.append(str(r.window_id))
                    else:
                        if block is not None:
                            region_rows.append((tail, comparison, p2, recipient, str(c), *block, ",".join(block_ids)))
                        block = [int(r.start), int(r.end), 1, float(r.FdM)]
                        block_ids = [str(r.window_id)]
                if block is not None:
                    region_rows.append((tail, comparison, p2, recipient, str(c), *block, ",".join(block_ids)))
    regions = pd.DataFrame(region_rows, columns=[
        "fdm_tail", "comparison", "P2", "recipient_population", "chromosome", "start", "end",
        "n_100kb_windows", "max_fdm", "window_ids"])
    regions["fdm_tail_label"] = regions.fdm_tail.map(TAIL_LABEL)
    regions["introgression_region_id"] = [
        f"FDM{TAIL_LABEL[t].replace('%','pct')}_{safe(comp)}_{c}_R{i + 1:04d}"
        for i, (t, comp, c) in enumerate(zip(regions.fdm_tail, regions.comparison, regions.chromosome))]
    return pd.DataFrame(threshold_rows), regions


def interval_pairs(left: pd.DataFrame, right: pd.DataFrame, matched_recipient: bool = False) -> pd.DataFrame:
    rows = []
    if left.empty or right.empty:
        return pd.DataFrame()
    left_cols, right_cols = list(left.columns), list(right.columns)
    for values in left.itertuples(index=False, name=None):
        a = dict(zip(left_cols, values))
        z = right[right.chromosome == a["chromosome"]]
        if matched_recipient:
            z = z[z.population == a["recipient_population"]]
        z = z[(z.start <= a["end"]) & (z.end >= a["start"])]
        for bvalues in z.itertuples(index=False, name=None):
            b = dict(zip(right_cols, bvalues))
            row = {f"left_{col}": a[col] for col in left_cols}
            row.update({f"right_{col}": b[col] for col in right_cols})
            row["overlap_start"] = max(int(a["start"]), int(b["start"]))
            row["overlap_end"] = min(int(a["end"]), int(b["end"]))
            row["overlap_bp"] = row["overlap_end"] - row["overlap_start"] + 1
            rows.append(row)
    return pd.DataFrame(rows)


def target_index(target: pd.DataFrame, key_cols: list[str] | None = None):
    key_cols = key_cols or []
    ans = {}
    group_cols = [*key_cols, "chromosome"]
    grouper = group_cols[0] if len(group_cols) == 1 else group_cols
    for key, group in target.groupby(grouper, sort=False):
        if not isinstance(key, tuple):
            key = (key,)
        merged = merge_simple(group)
        z = merged.sort_values("start")
        ans[tuple(str(x) for x in key)] = (z.start.to_numpy(np.int64), z.end.to_numpy(np.int64))
    return ans


def observed_hits(query: pd.DataFrame, index, key_cols: list[str] | None = None) -> int:
    key_cols = key_cols or []
    count = 0
    for r in query.itertuples(index=False):
        key = tuple(str(getattr(r, x)) for x in key_cols) + (str(r.chromosome),)
        ts, te = index.get(key, (np.array([], dtype=np.int64), np.array([], dtype=np.int64)))
        j = np.searchsorted(ts, int(r.end), side="right") - 1
        count += int(j >= 0 and te[j] >= int(r.start))
    return count


def permutation_test(query: pd.DataFrame, target: pd.DataFrame, lengths: dict[str, int],
                     key: str, key_cols: list[str] | None = None) -> dict[str, float]:
    key_cols = key_cols or []
    idx = target_index(target, key_cols)
    observed = observed_hits(query, idx, key_cols)
    rng = np.random.default_rng(SEED + zlib.crc32(key.encode("utf-8")))
    null = np.zeros(PERMUTATIONS, dtype=np.int32)
    for r in query.itertuples(index=False):
        c = str(r.chromosome)
        width = int(r.end) - int(r.start) + 1
        max_start = int(lengths[c]) - width + 1
        if max_start < 1:
            raise ValueError(f"Interval wider than chromosome: {c}:{r.start}-{r.end}")
        starts = rng.integers(1, max_start + 1, size=PERMUTATIONS)
        ends = starts + width - 1
        target_key = tuple(str(getattr(r, x)) for x in key_cols) + (c,)
        ts, te = idx.get(target_key, (np.array([], dtype=np.int64), np.array([], dtype=np.int64)))
        if len(ts):
            j = np.searchsorted(ts, ends, side="right") - 1
            valid = j >= 0
            jj = np.maximum(j, 0)
            null += valid & (te[jj] >= starts)
    p_upper = (1 + int((null >= observed).sum())) / (PERMUTATIONS + 1)
    p_lower = (1 + int((null <= observed).sum())) / (PERMUTATIONS + 1)
    mean = float(null.mean())
    return {
        "n_query_loci": len(query), "observed_query_loci_overlapping": observed,
        "percent_query_loci_overlapping": 100 * observed / len(query) if len(query) else np.nan,
        "null_mean_query_loci_overlapping": mean,
        "null_sd_query_loci_overlapping": float(null.std(ddof=1)),
        "fold_enrichment": observed / mean if mean > 0 else np.nan,
        "empirical_p_enrichment": p_upper, "empirical_p_depletion": p_lower,
        "empirical_p_two_sided": min(1.0, 2 * min(p_upper, p_lower)),
        "n_permutations": PERMUTATIONS,
    }


def bh(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(q, 1)
    return pd.Series(out, index=values.index)


def stars(q):
    if q <= 0.001: return "***"
    if q <= 0.01: return "**"
    if q <= 0.05: return "*"
    return ""


def plot_heatmaps(pairwise: pd.DataFrame, tripartite: pd.DataFrame):
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "pdf.fonttype": 42})
    traits = ["ALL", *TRAIT_ORDER]
    fig = plt.figure(figsize=(19 / 2.54, 13 / 2.54))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 3], height_ratios=[4, 1],
                          left=.09, right=.985, bottom=.10, top=.94, wspace=.18, hspace=.34)

    def trait_panel(ax, frame, title, col_key):
        mat = np.full((len(traits), 3), np.nan); txt = [["" for _ in TAILS] for _ in traits]
        for i, trait in enumerate(traits):
            for j, tail in enumerate(TAILS):
                z = frame[(frame.Trait == trait) & np.isclose(frame[col_key], tail)]
                if len(z):
                    r = z.iloc[0]; mat[i, j] = np.log2((r.observed_query_loci_overlapping + .5) /
                                                        (r.null_mean_query_loci_overlapping + .5))
                    txt[i][j] = f"{int(r.observed_query_loci_overlapping)}{stars(r.fdr_q_two_sided)}"
        ax.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2), aspect="auto")
        ax.set_title(title, fontsize=8, fontweight="bold", pad=5)
        ax.set_xticks(range(3), [TAIL_LABEL[x] for x in TAILS], fontsize=7)
        ax.set_xlabel("Outlier tail", fontsize=7)
        ax.set_yticks(range(len(traits)), [TRAIT_LABEL.get(x, x) for x in traits], fontsize=6.5)
        for i in range(len(traits)):
            for j in range(3):
                ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=6)
        ax.tick_params(length=0)
        for spine in ax.spines.values(): spine.set_visible(False)

    qsel = pairwise[pairwise.analysis == "QTL-selection"]
    qintro = pairwise[pairwise.analysis == "QTL-introgression"]
    trait_panel(fig.add_subplot(gs[0, 0]), qsel, "QTL–Tajima’s D", "selection_tail")
    ax2 = fig.add_subplot(gs[0, 1]); trait_panel(ax2, qintro, "QTL–fdm", "fdm_tail")

    ax3 = fig.add_subplot(gs[0, 2])
    combos = [(d, f) for d in TAILS for f in TAILS]
    mat = np.full((len(traits), len(combos)), np.nan); txt = [["" for _ in combos] for _ in traits]
    for i, trait in enumerate(traits):
        for j, (dt, ft) in enumerate(combos):
            z = tripartite[(tripartite.Trait == trait) & np.isclose(tripartite.selection_tail, dt) &
                            np.isclose(tripartite.fdm_tail, ft)]
            if len(z):
                r = z.iloc[0]; mat[i, j] = np.log2((r.observed_query_loci_overlapping + .5) /
                                                    (r.null_mean_query_loci_overlapping + .5))
                txt[i][j] = f"{int(r.observed_query_loci_overlapping)}{stars(r.fdr_q_two_sided)}"
    image = ax3.imshow(mat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2), aspect="auto")
    ax3.set_title("QTL–selection–introgression", fontsize=8, fontweight="bold", pad=5)
    ax3.set_xticks(range(len(combos)), [f"D{TAIL_LABEL[d]}\nf{TAIL_LABEL[f]}" for d, f in combos], fontsize=6)
    ax3.set_yticks(range(len(traits)), [])
    for i in range(len(traits)):
        for j in range(len(combos)):
            ax3.text(j, i, txt[i][j], ha="center", va="center", fontsize=5.5)
    ax3.tick_params(length=0)
    for spine in ax3.spines.values(): spine.set_visible(False)

    si = pairwise[pairwise.analysis == "selection-introgression"]
    ax4 = fig.add_subplot(gs[1, 0])
    smat = np.full((3, 3), np.nan); stxt = [["" for _ in TAILS] for _ in TAILS]
    for i, dt in enumerate(TAILS):
        for j, ft in enumerate(TAILS):
            z = si[np.isclose(si.selection_tail, dt) & np.isclose(si.fdm_tail, ft)]
            if len(z):
                r = z.iloc[0]; smat[i, j] = np.log2((r.observed_query_loci_overlapping + .5) /
                                                     (r.null_mean_query_loci_overlapping + .5))
                stxt[i][j] = f"{int(r.observed_query_loci_overlapping)}{stars(r.fdr_q_two_sided)}"
    ax4.imshow(smat, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2), aspect="auto")
    ax4.set_title("Recipient-matched selection–introgression", fontsize=7, fontweight="bold")
    ax4.set_xticks(range(3), [f"f{TAIL_LABEL[x]}" for x in TAILS], fontsize=6)
    ax4.set_yticks(range(3), [f"D{TAIL_LABEL[x]}" for x in TAILS], fontsize=6)
    for i in range(3):
        for j in range(3): ax4.text(j, i, stxt[i][j], ha="center", va="center", fontsize=6)
    ax4.tick_params(length=0)
    for spine in ax4.spines.values(): spine.set_visible(False)
    ax_note = fig.add_subplot(gs[1, 1:]); ax_note.axis("off")
    ax_note.text(0, .88, "Cell value = observed overlapping loci; stars = two-sided FDR (* q≤0.05, ** q≤0.01, *** q≤0.001).\n"
                 "Colour = log2 enrichment after a 0.5-count continuity correction. Three-way overlap requires a shared genomic base.",
                 ha="left", va="top", fontsize=7)
    cbar = fig.colorbar(image, ax=ax_note, orientation="horizontal", fraction=.20, pad=.28, aspect=35)
    cbar.set_label("log2 enrichment", fontsize=7); cbar.ax.tick_params(labelsize=6, length=2)
    for suffix, kwargs in [(".png", {"dpi": 600}), (".pdf", {}), (".svg", {})]:
        fig.savefig((OUT / "overlap_sensitivity_matrix").with_suffix(suffix), facecolor="white", **kwargs)
    plt.close(fig)


def write_report(pairwise, tripartite, sel_thresholds, intro_thresholds, selection, introgression):
    lines = [
        "# All merged-QTL overlap sensitivity analysis", "",
        f"This analysis used all 305 ±50-kb merged QTLs and {PERMUTATIONS:,} chromosome- and width-preserving permutations per test.", "",
        "Selection outliers were the population-specific lower 5%, 1%, and 0.1% tails of callable 10-kb Tajima's D windows. Directly adjacent outlier windows were merged within population and chromosome. Nucleotide diversity was not used.", "",
        "Introgression outliers were comparison-specific upper 5%, 1%, and 0.1% tails of positive fDM among valid 100-kb windows, with nsnp ≥ 100 and D < 0. Overlapping or adjacent windows were merged within comparison and chromosome. Selection–introgression and tripartite analyses were recipient-matched.", "",
        "QTL–selection tests use the genomic union across population-specific selection regions; interval-pair output retains the contributing population. A tripartite hit requires at least one genomic base shared simultaneously by a QTL, a recipient-matched Tajima's-D region, and an introgression region.", "",
        "## Numbers of merged scan regions", "",
    ]
    for tail in TAILS:
        lines.append(f"- D {TAIL_LABEL[tail]}: {len(selection[selection.selection_tail == tail]):,} population-specific merged regions")
    for tail in TAILS:
        lines.append(f"- fDM {TAIL_LABEL[tail]}: {len(introgression[introgression.fdm_tail == tail]):,} comparison-specific merged regions")
    lines += ["", "## Global pairwise and tripartite counts", ""]
    for tail in TAILS:
        r = pairwise[(pairwise.analysis == "QTL-selection") & (pairwise.Trait == "ALL") &
                     np.isclose(pairwise.selection_tail, tail)].iloc[0]
        lines.append(f"- QTL–D {TAIL_LABEL[tail]}: {int(r.observed_query_loci_overlapping)} of 305 QTLs overlap (null mean {r.null_mean_query_loci_overlapping:.2f}; enrichment P={r.empirical_p_enrichment:.4g})")
    for tail in TAILS:
        r = pairwise[(pairwise.analysis == "QTL-introgression") & (pairwise.Trait == "ALL") &
                     np.isclose(pairwise.fdm_tail, tail)].iloc[0]
        lines.append(f"- QTL–fdm {TAIL_LABEL[tail]}: {int(r.observed_query_loci_overlapping)} of 305 QTLs overlap (null mean {r.null_mean_query_loci_overlapping:.2f}; enrichment P={r.empirical_p_enrichment:.4g})")
    tri_total = int(tripartite[tripartite.Trait == "ALL"].observed_query_loci_overlapping.sum())
    lines.append(f"- Tripartite: {tri_total} QTL hits summed across the nine threshold combinations; every individual combination has zero exact three-way overlaps." if tri_total == 0 else f"- Tripartite: {tri_total} QTL hits summed across threshold combinations.")
    lines += ["", "## FDR-significant departures from the null", ""]
    sig_pair = pairwise[pairwise.fdr_q_two_sided <= .05]
    sig_tri = tripartite[tripartite.fdr_q_two_sided <= .05]
    if sig_pair.empty and sig_tri.empty:
        lines.append("No enrichment test remained significant at FDR 5%.")
    else:
        for r in pd.concat([sig_pair, sig_tri], ignore_index=True).itertuples():
            trait = getattr(r, "Trait", "ALL")
            dlab = getattr(r, "selection_tail_label", "")
            flab = getattr(r, "fdm_tail_label", "")
            threshold_text = ", ".join(x for x in [f"D {dlab}" if pd.notna(dlab) and dlab else "",
                                                        f"fdm {flab}" if pd.notna(flab) and flab else ""] if x)
            direction = "enrichment" if r.observed_query_loci_overlapping >= r.null_mean_query_loci_overlapping else "depletion"
            lines.append(f"- {r.analysis}, {trait}, {threshold_text}: {direction}; observed {r.observed_query_loci_overlapping}, expected {r.null_mean_query_loci_overlapping:.2f}, fold {r.fold_enrichment:.2f}, two-sided q={r.fdr_q_two_sided:.4g}")
    lines += ["", "Detailed observed counts, null expectations, empirical P values, and FDR values are in `pairwise_overlap_tests.tsv` and `tripartite_overlap_tests.tsv`.", ""]
    (OUT / "ANALYSIS_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    qtl = pd.read_csv(QTL_FILE).rename(columns={"CHR": "chromosome"})
    qtl["chromosome"] = qtl.chromosome.map(chrom)
    qtl[["start", "end"]] = qtl[["start", "end"]].astype(int)
    if qtl.QTL_ID.nunique() != 305:
        raise ValueError(f"Expected 305 all-merged QTLs, found {qtl.QTL_ID.nunique()}")

    sel_thresholds, selection = make_selection_regions()
    intro_thresholds, introgression = make_introgression_regions()
    qc_rows = []
    selection_lengths_ok = bool(((selection.end - selection.start + 1) ==
                                 selection.n_10kb_windows * 10_000).all())
    qc_rows.append({"check": "merged selection length equals n_windows × 10 kb", "failures": 0 if selection_lengths_ok else 1})
    for name, frame, tail_col, groups in [
        ("selection", selection, "selection_tail", ["population", "chromosome"]),
        ("introgression", introgression, "fdm_tail", ["comparison", "chromosome"]),
    ]:
        for strict, relaxed in [(0.01, 0.05), (0.001, 0.01)]:
            a, b = frame[np.isclose(frame[tail_col], strict)], frame[np.isclose(frame[tail_col], relaxed)]
            failures = 0
            for r in a.itertuples(index=False):
                mask = pd.Series(True, index=b.index)
                for col in groups:
                    mask &= b[col].astype(str) == str(getattr(r, col))
                mask &= (b.start <= int(r.start)) & (b.end >= int(r.end))
                failures += int(not mask.any())
            qc_rows.append({"check": f"{name} {TAIL_LABEL[strict]} regions nested within {TAIL_LABEL[relaxed]}", "failures": failures})
    pd.DataFrame(qc_rows).to_csv(OUT / "run_QC.tsv", sep="\t", index=False)
    if any(row["failures"] for row in qc_rows):
        raise ValueError(f"Sensitivity-set QC failed: {qc_rows}")
    sel_thresholds.to_csv(OUT / "tajima_D_population_thresholds.tsv", sep="\t", index=False)
    selection.to_csv(OUT / "tajima_D_merged_outlier_regions.tsv", sep="\t", index=False)
    intro_thresholds.to_csv(OUT / "fdm_comparison_thresholds.tsv", sep="\t", index=False)
    introgression.to_csv(OUT / "fdm_merged_outlier_regions.tsv", sep="\t", index=False)

    lengths = chromosome_lengths(GFF)
    for frame in [qtl, selection, introgression]:
        for c, end in frame.groupby("chromosome").end.max().items():
            lengths[str(c)] = max(lengths.get(str(c), 0), int(end))

    pairwise_tests, qsel_all, qintro_all, si_all, tri_all, trip_tests = [], [], [], [], [], []
    trait_groups = [("ALL", qtl), *[(t, qtl[qtl.Trait == t]) for t in TRAIT_ORDER]]

    for dt in TAILS:
        sel = selection[selection.selection_tail == dt]
        qsel = interval_pairs(qtl, sel)
        if len(qsel):
            qsel.insert(0, "selection_tail", dt); qsel.insert(1, "selection_tail_label", TAIL_LABEL[dt])
            qsel_all.append(qsel)
        target = merge_simple(sel)
        for trait, group in trait_groups:
            test = permutation_test(group, target, lengths, f"qsel|{dt}|{trait}")
            pairwise_tests.append({"analysis": "QTL-selection", "Trait": trait,
                                   "selection_tail": dt, "selection_tail_label": TAIL_LABEL[dt],
                                   "fdm_tail": np.nan, "fdm_tail_label": "", **test})

    for ft in TAILS:
        intro = introgression[introgression.fdm_tail == ft]
        qintro = interval_pairs(qtl, intro)
        if len(qintro):
            qintro.insert(0, "fdm_tail", ft); qintro.insert(1, "fdm_tail_label", TAIL_LABEL[ft])
            qintro_all.append(qintro)
        target = merge_simple(intro)
        for trait, group in trait_groups:
            test = permutation_test(group, target, lengths, f"qintro|{ft}|{trait}")
            pairwise_tests.append({"analysis": "QTL-introgression", "Trait": trait,
                                   "selection_tail": np.nan, "selection_tail_label": "",
                                   "fdm_tail": ft, "fdm_tail_label": TAIL_LABEL[ft], **test})

    for dt in TAILS:
        sel = selection[selection.selection_tail == dt]
        for ft in TAILS:
            intro = introgression[introgression.fdm_tail == ft]
            si = interval_pairs(intro, sel, matched_recipient=True)
            if len(si):
                si.insert(0, "selection_tail", dt); si.insert(1, "selection_tail_label", TAIL_LABEL[dt])
                si.insert(2, "fdm_tail", ft); si.insert(3, "fdm_tail_label", TAIL_LABEL[ft]); si_all.append(si)
                common = si.rename(columns={"overlap_start": "start", "overlap_end": "end",
                                             "left_chromosome": "chromosome"}).copy()
                common["common_region_id"] = [f"COMMON_D{TAIL_LABEL[dt]}_F{TAIL_LABEL[ft]}_{i+1:05d}" for i in range(len(common))]
            else:
                common = pd.DataFrame(columns=["chromosome", "start", "end"])
            matched_target = sel.rename(columns={"population": "recipient_population"}).copy()
            test = permutation_test(intro, matched_target, lengths, f"si|{dt}|{ft}",
                                    key_cols=["recipient_population"])
            pairwise_tests.append({"analysis": "selection-introgression", "Trait": "ALL",
                                   "selection_tail": dt, "selection_tail_label": TAIL_LABEL[dt],
                                   "fdm_tail": ft, "fdm_tail_label": TAIL_LABEL[ft], **test})

            common_target = merge_simple(common)
            tri = interval_pairs(qtl, common) if len(common) else pd.DataFrame()
            if len(tri):
                tri.insert(0, "selection_tail", dt); tri.insert(1, "selection_tail_label", TAIL_LABEL[dt])
                tri.insert(2, "fdm_tail", ft); tri.insert(3, "fdm_tail_label", TAIL_LABEL[ft]); tri_all.append(tri)
            for trait, group in trait_groups:
                ttest = permutation_test(group, common_target, lengths, f"tri|{dt}|{ft}|{trait}")
                trip_tests.append({"analysis": "tripartite", "Trait": trait,
                                   "selection_tail": dt, "selection_tail_label": TAIL_LABEL[dt],
                                   "fdm_tail": ft, "fdm_tail_label": TAIL_LABEL[ft], **ttest})

    pairwise = pd.DataFrame(pairwise_tests)
    pairwise["fdr_q_enrichment"] = pairwise.groupby("analysis", group_keys=False)["empirical_p_enrichment"].apply(bh)
    pairwise["fdr_q_depletion"] = pairwise.groupby("analysis", group_keys=False)["empirical_p_depletion"].apply(bh)
    pairwise["fdr_q_two_sided"] = pairwise.groupby("analysis", group_keys=False)["empirical_p_two_sided"].apply(bh)
    tripartite = pd.DataFrame(trip_tests)
    tripartite["fdr_q_enrichment"] = bh(tripartite.empirical_p_enrichment)
    tripartite["fdr_q_depletion"] = bh(tripartite.empirical_p_depletion)
    tripartite["fdr_q_two_sided"] = bh(tripartite.empirical_p_two_sided)
    pairwise.to_csv(OUT / "pairwise_overlap_tests.tsv", sep="\t", index=False)
    tripartite.to_csv(OUT / "tripartite_overlap_tests.tsv", sep="\t", index=False)
    pd.concat(qsel_all, ignore_index=True).to_csv(OUT / "QTL_selection_overlap_pairs.tsv", sep="\t", index=False)
    pd.concat(qintro_all, ignore_index=True).to_csv(OUT / "QTL_introgression_overlap_pairs.tsv", sep="\t", index=False)
    pd.concat(si_all, ignore_index=True).to_csv(OUT / "selection_introgression_recipient_matched_pairs.tsv", sep="\t", index=False)
    if tri_all:
        pd.concat(tri_all, ignore_index=True).to_csv(OUT / "tripartite_overlap_pairs.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(OUT / "tripartite_overlap_pairs.tsv", sep="\t", index=False)
    plot_heatmaps(pairwise, tripartite)
    write_report(pairwise, tripartite, sel_thresholds, intro_thresholds, selection, introgression)


if __name__ == "__main__":
    main()
