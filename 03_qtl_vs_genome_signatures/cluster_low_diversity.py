#!/usr/bin/env python3
"""Population-resolved low-diversity enrichment in STRING functional clusters.

This analysis deliberately uses the observed Tajima's D and nucleotide-diversity
values.  For every gene and population, the value retained is the lowest value
from a genomic window that physically overlaps the transcript.  STRING clusters
are then tested by covariate-stratified label permutation at PPI score 0.70.
"""
from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[0]
DATA = ROOT / "shared_data"
PPI = ROOT


def normalise_chrom(value: object) -> str:
    x = str(value).strip()
    if x.lower().startswith("chr"):
        x = x[3:]
    if x == "10":
        x = "X"
    if x.startswith("Tmol_scf_"):
        return x
    return f"chr{x}"


def parse_attributes(text: str) -> dict[str, str]:
    ans: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            ans[key] = value
    return ans


def read_genes(gff_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = []
    features = []
    with gff_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs = fields
            if feature not in {"mRNA", "CDS", "UTR"}:
                continue
            a = parse_attributes(attrs)
            gene_id = a.get("ID") if feature == "mRNA" else a.get("Parent")
            if not gene_id:
                continue
            row = {
                "gene_id": gene_id,
                "chrom": normalise_chrom(chrom),
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "feature": feature,
            }
            features.append(row)
            if feature == "mRNA":
                genes.append({k: row[k] for k in ("gene_id", "chrom", "start", "end", "strand")})
    g = pd.DataFrame(genes).drop_duplicates("gene_id")
    g["gene_length"] = g["end"] - g["start"] + 1
    return g, pd.DataFrame(features)


def read_metric(path: Path, metric: str) -> pd.DataFrame:
    # The source files contain an unnamed row-number field.  pandas correctly
    # treats it as the index because data rows have one more field than the header.
    dat = pd.read_csv(path, sep="\t", na_values=["NA", "NaN"])
    dat = dat.reset_index(drop=True)
    dat["chrom"] = dat["chr"].map(normalise_chrom)
    dat[metric] = pd.to_numeric(dat[metric], errors="coerce")
    dat["pos"] = pd.to_numeric(dat["pos"], errors="coerce").astype("Int64")
    dat["pos2"] = pd.to_numeric(dat["pos2"], errors="coerce").astype("Int64")
    dat = dat.dropna(subset=["pos", "pos2"]).copy()
    dat["pos"] = dat["pos"].astype(np.int64)
    dat["pos2"] = dat["pos2"].astype(np.int64)
    if metric == "D":
        dat["window_start"] = dat["pos"]
        dat["window_end"] = dat["pos2"]
    else:
        # pi is stored at successive 10-kb window midpoints (5,000, 15,000, ...).
        # Reconstruct their 1-based physical spans without altering observed pi.
        dat["window_start"] = np.maximum(1, dat["pos"] - 4_999)
        dat["window_end"] = dat["pos"] + 5_000
    return dat[["chrom", "window_start", "window_end", metric, "pop"]]


def window_thresholds(dat: pd.DataFrame, metric: str, tail: float) -> pd.DataFrame:
    rows = []
    for pop, group in dat.groupby("pop", sort=True):
        values = group[metric].dropna()
        rows.append({
            "metric": metric,
            "population": pop,
            "tail_probability": tail,
            "lower_tail_threshold": values.quantile(tail),
            "n_callable_windows": len(values),
        })
    return pd.DataFrame(rows)


def gene_population_minima(
    genes: pd.DataFrame, dat: pd.DataFrame, metric: str, thresholds: pd.DataFrame
) -> pd.DataFrame:
    threshold_map = thresholds.set_index("population")["lower_tail_threshold"].to_dict()
    rows: list[dict[str, object]] = []
    genes_by_chrom = {chrom: x for chrom, x in genes.groupby("chrom", sort=False)}
    for (pop, chrom), windows in dat.groupby(["pop", "chrom"], sort=False):
        if chrom not in genes_by_chrom:
            continue
        w = windows.sort_values("window_start")
        starts = w["window_start"].to_numpy(np.int64)
        ends = w["window_end"].to_numpy(np.int64)
        values = w[metric].to_numpy(float)
        for gene in genes_by_chrom[chrom].itertuples(index=False):
            hi = int(np.searchsorted(starts, gene.end, side="right"))
            lo = int(np.searchsorted(ends, gene.start, side="left"))
            if hi <= lo:
                continue
            local_values = values[lo:hi]
            finite = np.flatnonzero(np.isfinite(local_values))
            if not len(finite):
                continue
            relative = finite[int(np.argmin(local_values[finite]))]
            idx = lo + int(relative)
            minimum = float(values[idx])
            rows.append({
                "gene_id": gene.gene_id,
                "chrom": chrom,
                "gene_start": gene.start,
                "gene_end": gene.end,
                "population": pop,
                "metric": metric,
                "minimum_value": minimum,
                "n_overlapping_callable_windows": int(len(finite)),
                "minimum_window_start": int(starts[idx]),
                "minimum_window_end": int(ends[idx]),
                "lower_tail_threshold": float(threshold_map[pop]),
                "low_tail": bool(minimum <= threshold_map[pop]),
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["within_population_gene_percentile"] = out.groupby("population")["minimum_value"].rank(
            method="average", pct=True
        )
        out["continuous_lower_tail_score"] = -np.log10(
            out["within_population_gene_percentile"].clip(lower=1e-12)
        )
    return out


def semijoin(values) -> str:
    vals = sorted({str(x) for x in values if pd.notna(x) and str(x)})
    return "; ".join(vals)


def annotate_gene_summary(
    genes: pd.DataFrame,
    dmin: pd.DataFrame,
    pmin: pd.DataFrame,
    annotations: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def one_metric(x: pd.DataFrame, prefix: str) -> pd.DataFrame:
        idx = x.groupby("gene_id")["minimum_value"].idxmin()
        global_min = x.loc[idx, ["gene_id", "minimum_value", "population", "minimum_window_start", "minimum_window_end"]]
        global_min = global_min.rename(columns={
            "minimum_value": f"minimum_{prefix}",
            "population": f"minimum_{prefix}_population",
            "minimum_window_start": f"minimum_{prefix}_window_start",
            "minimum_window_end": f"minimum_{prefix}_window_end",
        })
        agg = x.groupby("gene_id", as_index=False).agg(
            **{
                f"n_callable_{prefix}_populations": ("population", "nunique"),
                f"n_low_{prefix}_populations": ("low_tail", "sum"),
                f"total_callable_{prefix}_windows": ("n_overlapping_callable_windows", "sum"),
            }
        )
        lowp = (x[x.low_tail].groupby("gene_id")["population"].agg(semijoin)
                .rename(f"low_{prefix}_populations").reset_index())
        top3 = (x.groupby("gene_id")["continuous_lower_tail_score"]
                .apply(lambda z: float(np.mean(np.sort(z.to_numpy(float))[-3:])))
                .rename(f"continuous_low_{prefix}_score_top3").reset_index())
        return (agg.merge(global_min, on="gene_id", how="left")
                .merge(lowp, on="gene_id", how="left").merge(top3, on="gene_id", how="left"))

    ds = one_metric(dmin, "D")
    ps = one_metric(pmin, "pi")
    common = dmin.merge(
        pmin, on=["gene_id", "chrom", "gene_start", "gene_end", "population"],
        suffixes=("_D", "_pi"), how="inner"
    )
    common["joint_low_D_pi"] = common["low_tail_D"] & common["low_tail_pi"]
    common["continuous_joint_tail_score"] = np.minimum(
        common["continuous_lower_tail_score_D"], common["continuous_lower_tail_score_pi"]
    )
    joint = common.groupby("gene_id", as_index=False).agg(
        n_joint_low_D_pi_populations=("joint_low_D_pi", "sum"),
        n_joint_callable_populations=("population", "nunique"),
    )
    jointp = (common[common.joint_low_D_pi].groupby("gene_id")["population"].agg(semijoin)
              .rename("joint_low_D_pi_populations").reset_index())
    joint = joint.merge(jointp, on="gene_id", how="left")
    joint_score = (common.groupby("gene_id")["continuous_joint_tail_score"]
                   .apply(lambda z: float(np.mean(np.sort(z.to_numpy(float))[-3:])))
                   .rename("continuous_joint_low_D_pi_score_top3").reset_index())
    joint = joint.merge(joint_score, on="gene_id", how="left")
    out = genes.merge(ds, on="gene_id", how="left").merge(ps, on="gene_id", how="left").merge(joint, on="gene_id", how="left")
    out = out.merge(annotations[["gene_id", "short_description"]], on="gene_id", how="left")
    for col in ["low_D_populations", "low_pi_populations", "joint_low_D_pi_populations", "short_description"]:
        out[col] = out[col].fillna("")
    for col in ["n_callable_D_populations", "n_low_D_populations", "total_callable_D_windows",
                "n_callable_pi_populations", "n_low_pi_populations", "total_callable_pi_windows",
                "n_joint_low_D_pi_populations", "n_joint_callable_populations"]:
        out[col] = out[col].fillna(0).astype(int)
    out["low_D_any"] = out["n_low_D_populations"] > 0
    out["low_pi_any"] = out["n_low_pi_populations"] > 0
    out["joint_low_D_pi_any"] = out["n_joint_low_D_pi_populations"] > 0
    return out, common


def overlap_introgression(genes: pd.DataFrame, path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    regions = pd.read_csv(path, sep="\t")
    regions = regions[regions["threshold_level"].eq("stringent")].copy()
    regions["chrom"] = regions["chromosome"].map(normalise_chrom)
    rows = []
    for chrom, reg in regions.groupby("chrom"):
        gs = genes[genes.chrom.eq(chrom)]
        if gs.empty:
            continue
        for r in reg.itertuples(index=False):
            hit = gs[(gs.start <= int(r.end)) & (gs.end >= int(r.start))]
            for g in hit.itertuples(index=False):
                rows.append({
                    "gene_id": g.gene_id,
                    "population": r.P2_population,
                    "comparison": r.comparison,
                    "introgression_region_id": r.introgression_region_id,
                    "region_chrom": chrom,
                    "region_start": int(r.start),
                    "region_end": int(r.end),
                    "gene_start": int(g.start),
                    "gene_end": int(g.end),
                })
    overlaps = pd.DataFrame(rows)
    if overlaps.empty:
        summary = pd.DataFrame(columns=["gene_id", "introgression_populations", "introgression_comparisons", "n_introgression_regions"])
    else:
        summary = overlaps.groupby("gene_id", as_index=False).agg(
            introgression_populations=("population", semijoin),
            introgression_comparisons=("comparison", semijoin),
            n_introgression_regions=("introgression_region_id", "nunique"),
        )
    return overlaps, summary


def qbins(values: np.ndarray, n: int) -> np.ndarray:
    if len(values) == 0:
        return np.array([], dtype=int)
    ranks = pd.Series(values).rank(method="average", pct=True).to_numpy()
    return np.minimum((ranks * n).astype(int), n - 1)


def make_strata(nodes: pd.DataFrame) -> np.ndarray:
    chrom_group = nodes["chrom"].where(nodes["chrom"].isin([f"chr{x}" for x in range(1, 10)] + ["chrX"]), "scaffold")
    length_bin = qbins(np.log1p(nodes["gene_length"].to_numpy(float)), 3)
    callable_bin = qbins(np.log1p((nodes["total_callable_D_windows"] + nodes["total_callable_pi_windows"]).to_numpy(float)), 3)
    degree_bin = qbins(nodes["degree"].to_numpy(float), 4)
    burden_bin = qbins(nodes["eligible_cluster_memberships"].to_numpy(float), 4)
    fine = np.array([f"{c}|{a}|{b}|{d}|{e}" for c, a, b, d, e in zip(chrom_group, length_bin, callable_bin, degree_bin, burden_bin)], dtype=object)
    counts = Counter(fine)
    medium_degree = qbins(nodes["degree"].to_numpy(float), 2)
    medium_burden = qbins(nodes["eligible_cluster_memberships"].to_numpy(float), 2)
    medium = np.array([f"{c}|{a}|{b}|{d}|{e}" for c, a, b, d, e in zip(chrom_group, length_bin, callable_bin, medium_degree, medium_burden)], dtype=object)
    out = fine.copy()
    for i, value in enumerate(out):
        if counts[value] < 8:
            out[i] = "M|" + medium[i]
    counts = Counter(out)
    coarse_len = qbins(np.log1p(nodes["gene_length"].to_numpy(float)), 2)
    coarse_call = qbins(np.log1p((nodes["total_callable_D_windows"] + nodes["total_callable_pi_windows"]).to_numpy(float)), 2)
    for i, value in enumerate(out):
        if counts[value] < 8:
            out[i] = f"C|{chrom_group.iloc[i]}|{coarse_len[i]}|{coarse_call[i]}"
    counts = Counter(out)
    for i, value in enumerate(out):
        if counts[value] < 8:
            out[i] = f"A|{coarse_len[i]}|{coarse_call[i]}"
    return pd.Categorical(out).codes.astype(np.int32)


def bh(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def cluster_counts(cluster_ix: np.ndarray, codes_at_membership: np.ndarray, n_clusters: int) -> tuple[np.ndarray, ...]:
    d = np.bincount(cluster_ix, weights=((codes_at_membership & 1) > 0), minlength=n_clusters).astype(int)
    p = np.bincount(cluster_ix, weights=((codes_at_membership & 2) > 0), minlength=n_clusters).astype(int)
    j = np.bincount(cluster_ix, weights=((codes_at_membership & 4) > 0), minlength=n_clusters).astype(int)
    both = np.bincount(cluster_ix, weights=((codes_at_membership & 3) == 3), minlength=n_clusters).astype(int)
    return d, p, j, both


def cluster_continuous_scores(
    cluster_ix: np.ndarray, score_at_membership: np.ndarray,
    n_clusters: int, denominators: np.ndarray,
) -> tuple[np.ndarray, ...]:
    values = []
    for column in range(score_at_membership.shape[1]):
        total = np.bincount(cluster_ix, weights=score_at_membership[:, column], minlength=n_clusters)
        values.append(total / denominators)
    return tuple(values)


def run_cluster_test(
    gene_summary: pd.DataFrame,
    score_min: float,
    min_cluster_size: int,
    max_cluster_size: int,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    info = pd.read_csv(PPI / "STRG0A51GAV.clusters.info.v12.0.txt.gz", sep="\t")
    membership = pd.read_csv(PPI / "STRG0A51GAV.clusters.proteins.v12.0.txt.gz", sep="\t")
    info = info[(info.cluster_size >= min_cluster_size) & (info.cluster_size <= max_cluster_size)].copy()
    membership = membership[membership.cluster_id.isin(info.cluster_id)].copy()
    membership["gene_id"] = membership["protein_id"].str.replace(r"^[^.]+\.", "", regex=True)

    evaluable = gene_summary[(gene_summary.n_callable_D_populations > 0) & (gene_summary.n_callable_pi_populations > 0)].copy()
    eval_ids = set(evaluable.gene_id)
    membership = membership[membership.gene_id.isin(eval_ids)].copy()
    tested_sizes = membership.groupby("cluster_id").gene_id.nunique()
    keep = tested_sizes[tested_sizes >= min_cluster_size].index
    info = info[info.cluster_id.isin(keep)].copy()
    membership = membership[membership.cluster_id.isin(keep)].drop_duplicates(["cluster_id", "protein_id"]).copy()

    clusters = info.cluster_id.to_numpy()
    cluster_index = {x: i for i, x in enumerate(clusters)}
    nodes = np.sort(membership.protein_id.unique())
    node_index = {x: i for i, x in enumerate(nodes)}
    membership["cluster_ix"] = membership.cluster_id.map(cluster_index).astype(np.int32)
    membership["node_ix"] = membership.protein_id.map(node_index).astype(np.int32)
    cluster_ix = membership.cluster_ix.to_numpy(np.int32)
    membership_node_ix = membership.node_ix.to_numpy(np.int32)

    node_table = pd.DataFrame({"node_id": nodes})
    node_table["gene_id"] = node_table.node_id.str.replace(r"^[^.]+\.", "", regex=True)
    node_table = node_table.merge(evaluable, on="gene_id", how="left")

    network = pd.read_csv(HERE / "project_input" / "network.tsv", sep="\t")
    score = pd.to_numeric(network.combined_score, errors="coerce")
    if score.max() > 1:
        score = score / 1000.0
    network = network[score >= score_min].copy()
    degree = np.zeros(len(nodes), dtype=np.int32)
    for col in ("protein1", "protein2"):
        idx = network[col].map(node_index).dropna().astype(int)
        degree += np.bincount(idx, minlength=len(nodes))
    node_table["degree"] = degree
    node_table["eligible_cluster_memberships"] = np.bincount(membership_node_ix, minlength=len(nodes))
    node_table["stratum"] = make_strata(node_table)

    codes = np.zeros(len(nodes), dtype=np.int8)
    codes[node_table.low_D_any.to_numpy(bool)] |= 1
    codes[node_table.low_pi_any.to_numpy(bool)] |= 2
    codes[node_table.joint_low_D_pi_any.to_numpy(bool)] |= 4
    observed = cluster_counts(cluster_ix, codes[membership_node_ix], len(clusters))
    exceed = [np.zeros(len(clusters), dtype=np.int32) for _ in observed]
    score_columns = ["continuous_low_D_score_top3", "continuous_low_pi_score_top3",
                     "continuous_joint_low_D_pi_score_top3"]
    node_scores = node_table[score_columns].fillna(0).to_numpy(float)
    denominators = np.bincount(cluster_ix, minlength=len(clusters)).astype(float)
    observed_continuous = cluster_continuous_scores(
        cluster_ix, node_scores[membership_node_ix], len(clusters), denominators
    )
    exceed_continuous = [np.zeros(len(clusters), dtype=np.int32) for _ in observed_continuous]
    strata_members = [np.flatnonzero(node_table.stratum.to_numpy() == x) for x in np.unique(node_table.stratum)]
    rng = np.random.default_rng(seed)
    perm_index = np.arange(len(nodes), dtype=np.int32)
    for _ in range(permutations):
        for ix in strata_members:
            perm_index[ix] = rng.permutation(ix)
        values = cluster_counts(cluster_ix, codes[perm_index][membership_node_ix], len(clusters))
        for j in range(len(values)):
            exceed[j] += values[j] >= observed[j]
        continuous_values = cluster_continuous_scores(
            cluster_ix, node_scores[perm_index][membership_node_ix], len(clusters), denominators
        )
        for j in range(len(continuous_values)):
            exceed_continuous[j] += continuous_values[j] >= observed_continuous[j]

    result = info.set_index("cluster_id").loc[clusters].reset_index()
    result["tested_genes"] = result.cluster_id.map(membership.groupby("cluster_id").gene_id.nunique()).astype(int)
    intro_ids = set(gene_summary.loc[gene_summary.n_introgression_regions > 0, "gene_id"])
    selection_path = HERE / "project_input" / "selection_genes.tsv"
    selection_ids = set(pd.read_csv(selection_path, sep="\t").gene_id) if selection_path.exists() else set()
    result["introgressed_genes"] = result.cluster_id.map(membership[membership.gene_id.isin(intro_ids)].groupby("cluster_id").gene_id.nunique()).fillna(0).astype(int)
    result["selection_genes"] = result.cluster_id.map(membership[membership.gene_id.isin(selection_ids)].groupby("cluster_id").gene_id.nunique()).fillna(0).astype(int)
    names = ["low_D_genes", "low_pi_genes", "joint_same_population_genes", "low_D_and_low_pi_genes"]
    for name, obs, exc in zip(names, observed, exceed):
        result[name] = obs
        result[f"p_{name}"] = (exc + 1) / (permutations + 1)
        result[f"FDR_{name}"] = bh(result[f"p_{name}"].to_numpy())
        result[f"neglog10p_{name}"] = -np.log10(result[f"p_{name}"].clip(lower=1 / (permutations + 1)))
    continuous_names = ["continuous_low_D", "continuous_low_pi", "continuous_joint_low_D_pi"]
    for name, obs, exc in zip(continuous_names, observed_continuous, exceed_continuous):
        result[f"mean_{name}_score"] = obs
        result[f"p_{name}"] = (exc + 1) / (permutations + 1)
        result[f"FDR_{name}"] = bh(result[f"p_{name}"].to_numpy())
    result["dual_nominal"] = (result.p_low_D_genes <= 0.05) & (result.p_low_pi_genes <= 0.05)
    result["dual_FDR10"] = (result.FDR_low_D_genes <= 0.10) & (result.FDR_low_pi_genes <= 0.10)
    result["joint_FDR10"] = result.FDR_joint_same_population_genes <= 0.10
    return result, membership, node_table, network


def nearest_parent_edges(membership: pd.DataFrame, cluster_table: pd.DataFrame, min_containment: float = 0.50) -> pd.DataFrame:
    cluster_sets = {k: set(v) for k, v in membership.groupby("cluster_id").gene_id}
    sizes = {k: len(v) for k, v in cluster_sets.items()}
    gene_clusters: dict[str, list[str]] = defaultdict(list)
    for cid, genes in cluster_sets.items():
        for gene in genes:
            gene_clusters[gene].append(cid)
    edges = []
    for cid, genes in cluster_sets.items():
        counts: Counter[str] = Counter()
        size = sizes[cid]
        for gene in genes:
            for other in gene_clusters[gene]:
                if other != cid and sizes[other] >= size:
                    counts[other] += 1
        if not counts:
            continue
        candidates = []
        for other, intersection in counts.items():
            containment = intersection / size
            union = size + sizes[other] - intersection
            jaccard = intersection / union
            if containment >= min_containment:
                candidates.append((containment, jaccard, -sizes[other], other, intersection))
        if candidates:
            containment, jaccard, _, parent, intersection = max(candidates)
            edges.append({
                "cluster1": cid, "cluster2": parent, "intersection_genes": intersection,
                "containment_smaller": containment, "jaccard": jaccard,
            })
    return pd.DataFrame(edges)


def collapse_redundant(result: pd.DataFrame, membership: pd.DataFrame, p_cutoff: float = 0.05) -> pd.DataFrame:
    formal = result.dual_FDR10 | result.joint_FDR10
    candidates = result[formal].copy()
    if candidates.empty:
        candidates = result[result.dual_nominal].copy()
        if candidates.empty:
            candidates = result.nsmallest(20, ["p_joint_same_population_genes", "p_low_D_genes", "p_low_pi_genes"]).copy()
        candidates["exploratory_only"] = True
    else:
        candidates["exploratory_only"] = False
    candidates["combined_evidence"] = (
        candidates.neglog10p_low_D_genes + candidates.neglog10p_low_pi_genes
        + candidates.neglog10p_joint_same_population_genes
    )
    # Prefer the smallest significant member of a nested STRING hierarchy as the
    # biological representative; retain broader parents in the output table.
    candidates = candidates.sort_values(["joint_FDR10", "dual_FDR10", "tested_genes", "combined_evidence"], ascending=[False, False, True, False])
    sets = {k: set(v) for k, v in membership[membership.cluster_id.isin(candidates.cluster_id)].groupby("cluster_id").gene_id}
    representatives: list[str] = []
    groups: dict[str, str] = {}
    for row in candidates.itertuples(index=False):
        cid = row.cluster_id
        best_rep = None
        best_j = 0.0
        for rep in representatives:
            inter = len(sets[cid] & sets[rep])
            union = len(sets[cid] | sets[rep])
            jac = inter / union if union else 0.0
            containment = inter / min(len(sets[cid]), len(sets[rep])) if min(len(sets[cid]), len(sets[rep])) else 0.0
            if max(jac, containment) >= 0.70 and jac > best_j:
                best_rep, best_j = rep, jac
        if best_rep is None:
            representatives.append(cid)
            groups[cid] = cid
        else:
            groups[cid] = best_rep
    candidates["redundancy_representative"] = candidates.cluster_id.map(groups)
    candidates["is_representative"] = candidates.cluster_id.eq(candidates.redundancy_representative)
    return candidates


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tail", type=float, default=0.01)
    ap.add_argument("--score-min", type=float, default=0.70)
    ap.add_argument("--min-cluster-size", type=int, default=5)
    ap.add_argument("--max-cluster-size", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260904)
    ap.add_argument("--outdir", default=str(HERE / "low_diversity_analysis"))
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    genes, features = read_genes(DATA / "gene" / "tmol_gene.gff")
    annotations = pd.read_csv(HERE / "project_input" / "gene_short_annotations.tsv", sep="\t", keep_default_na=False)
    d = read_metric(DATA / "Selection" / "allD_long_chrv3.txt", "D")
    pi = read_metric(DATA / "Selection" / "allpi_long_chrv3.txt", "pi")
    thresholds = pd.concat([window_thresholds(d, "D", args.tail), window_thresholds(pi, "pi", args.tail)], ignore_index=True)
    thresholds.to_csv(out / "population_lower_tail_thresholds.tsv", sep="\t", index=False)

    dmin = gene_population_minima(genes, d, "D", thresholds[thresholds.metric.eq("D")])
    pmin = gene_population_minima(genes, pi, "pi", thresholds[thresholds.metric.eq("pi")])
    minima = pd.concat([dmin, pmin], ignore_index=True)
    minima.to_csv(out / "gene_population_observed_minima.tsv", sep="\t", index=False)
    gene_summary, common = annotate_gene_summary(genes, dmin, pmin, annotations)

    intro_overlap, intro_summary = overlap_introgression(genes, DATA / "Introgression" / "introgression_100kb_candidate_regions.tsv")
    intro_overlap.to_csv(out / "gene_introgression_overlaps.tsv", sep="\t", index=False)
    gene_summary = gene_summary.merge(intro_summary, on="gene_id", how="left")
    gene_summary["introgression_populations"] = gene_summary["introgression_populations"].fillna("")
    gene_summary["introgression_comparisons"] = gene_summary["introgression_comparisons"].fillna("")
    gene_summary["n_introgression_regions"] = gene_summary["n_introgression_regions"].fillna(0).astype(int)
    gene_summary.to_csv(out / "gene_low_diversity_summary.tsv", sep="\t", index=False)

    cluster_result, membership, node_table, network = run_cluster_test(
        gene_summary, args.score_min, args.min_cluster_size, args.max_cluster_size,
        args.permutations, args.seed,
    )
    cluster_result.to_csv(out / "cluster_low_diversity_enrichment.tsv", sep="\t", index=False)
    membership.to_csv(out / "tested_cluster_memberships.tsv", sep="\t", index=False)
    node_table.to_csv(out / "node_permutation_strata.tsv", sep="\t", index=False)
    features.to_csv(out / "gene_model_features.tsv", sep="\t", index=False)
    edges = nearest_parent_edges(membership, cluster_result)
    edges.to_csv(out / "cluster_overlap_graph_edges.tsv", sep="\t", index=False)
    nonredundant = collapse_redundant(cluster_result, membership)
    nonredundant.to_csv(out / "nonredundant_low_diversity_pathways.tsv", sep="\t", index=False)

    common_joint = common[common.joint_low_D_pi].copy()
    intro_joint = intro_overlap.merge(common_joint[["gene_id", "population", "minimum_value_D", "minimum_value_pi",
                                                     "minimum_window_start_D", "minimum_window_end_D",
                                                     "minimum_window_start_pi", "minimum_window_end_pi"]],
                                      on=["gene_id", "population"], how="inner")
    intro_joint = intro_joint.merge(annotations[["gene_id", "short_description"]], on="gene_id", how="left")
    intro_joint.to_csv(out / "same_population_introgression_lowD_lowpi_genes.tsv", sep="\t", index=False)

    print(f"Genes: {len(genes):,}; D minima: {len(dmin):,}; pi minima: {len(pmin):,}")
    print(f"Eligible clusters: {len(cluster_result):,}; permutation strata: {node_table.stratum.nunique():,}")
    print("FDR <= 0.10:", {
        "low_D": int((cluster_result.FDR_low_D_genes <= .10).sum()),
        "low_pi": int((cluster_result.FDR_low_pi_genes <= .10).sum()),
        "joint_same_population": int((cluster_result.FDR_joint_same_population_genes <= .10).sum()),
        "dual": int(cluster_result.dual_FDR10.sum()),
    })
    print("Nominal dual clusters:", int(cluster_result.dual_nominal.sum()))
    print("Stringent-introgression genes jointly low in recipient population:", intro_joint.gene_id.nunique())
    show = cluster_result.sort_values(["p_joint_same_population_genes", "p_low_D_genes", "p_low_pi_genes"]).head(15)
    print(show[["cluster_id", "cluster_size", "best_described_by", "low_D_genes", "low_pi_genes",
                "joint_same_population_genes", "introgressed_genes", "p_low_D_genes", "p_low_pi_genes",
                "p_joint_same_population_genes", "FDR_joint_same_population_genes"]].to_string(index=False))


if __name__ == "__main__":
    main()
