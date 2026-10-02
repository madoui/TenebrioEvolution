#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np, pandas as pd

def overlaps(query,target):
    rows=[]; ids=set(); by={str(c):g for c,g in target.groupby("chromosome",sort=False)}
    for q in query.itertuples(index=False):
        g=by.get(str(q.chromosome)); hits=g[(g.start<=q.end)&(g.end>=q.start)] if g is not None else pd.DataFrame()
        if len(hits): ids.add(q.introgression_region_id)
        for t in hits.itertuples(index=False):
            s,e=max(q.start,t.start),min(q.end,t.end); rows.append({"introgression_region_id":q.introgression_region_id,"selection_region_id":t.selection_region_id,"chromosome":q.chromosome,"overlap_start":s,"overlap_end":e,"overlap_bp":e-s+1})
    return len(ids),pd.DataFrame(rows)

def main():
    here=Path(__file__).resolve().parent; root=here.parent; ap=argparse.ArgumentParser(); ap.add_argument("--permutations",type=int,default=10000); ap.add_argument("--seed",type=int,default=20260907); ap.add_argument("--outdir",type=Path,default=here/"results"); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    sf=root/"01_selective_sweep_signatures/results/overlapaware_selection_loci_5pct.tsv"
    if not sf.exists(): raise SystemExit("Run subsection 01 first")
    selection=pd.read_csv(sf,sep="\t",dtype={"chromosome":str}); selection=selection[selection.selection_class.eq("farmed_only")]
    intro=pd.read_csv(root/"shared_data/introgression/fdm_merged_regions.tsv",sep="\t",dtype={"chromosome":str}); intro=intro[np.isclose(intro.fdm_tail,.05)]
    lengths=pd.read_csv(root/"shared_data/genome/chromosome_lengths.tsv",sep="\t",dtype={"chromosome":str}); lm=dict(zip(lengths.chromosome.str.replace("chr","",regex=False),lengths.length.astype(int))); lm["10"]=lm.get("X",lm.get("10"))
    observed,pairs=overlaps(intro,selection); pairs.to_csv(a.outdir/"introgression_selection_overlap_pairs.tsv",sep="\t",index=False)
    rng=np.random.default_rng(a.seed); null=np.zeros(a.permutations,dtype=np.int32); target={str(c):g[["start","end"]].to_numpy(int) for c,g in selection.groupby("chromosome")}
    for q in intro.itertuples(index=False):
        chrom=str(q.chromosome); width=int(q.end-q.start+1); starts=rng.integers(1,lm[chrom]-width+2,size=a.permutations); ends=starts+width-1; hit=np.zeros(a.permutations,bool)
        for s,e in target.get(chrom,np.empty((0,2),int)): hit|=(starts<=e)&(ends>=s)
        null+=hit
    expected=float(null.mean()); p=(1+int((null>=observed).sum()))/(a.permutations+1)
    result=pd.DataFrame([{"n_introgression_regions":len(intro),"n_farm_only_selection_loci":len(selection),"observed_introgression_regions_overlapping":observed,"observed_selection_loci_overlapping":pairs.selection_region_id.nunique(),"null_mean_introgression_regions_overlapping":expected,"fold_enrichment":observed/expected if expected else np.nan,"empirical_p_enrichment":p,"n_permutations":a.permutations,"seed":a.seed}]); result.to_csv(a.outdir/"overlap_permutation_test.tsv",sep="\t",index=False); print(result.to_string(index=False))
if __name__=="__main__": main()
