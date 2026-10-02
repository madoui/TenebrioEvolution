#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
FARMS=["Dole7","Starfood Holland","Kingsect Belgium","Pronutrix+ Belgium","Papek Czechia","Mystik Canada"]; WILD="Wild Serbian"

def merge_population_runs(windows):
    individual=[]
    for r in windows.itertuples(index=False):
        supporters=[x for x in str(r.farmed_populations).split(";") if x]
        if str(r.wild_support).lower()=="true": supporters.append(WILD)
        individual += [(p,str(r.chromosome),int(r.start),int(r.end)) for p in supporters]
    runs=[]; frame=pd.DataFrame(individual,columns=["population","chromosome","start","end"])
    for (pop,chrom),g in frame.groupby(["population","chromosome"]):
        block=None
        for r in g.sort_values(["start","end"]).itertuples(index=False):
            if block and r.start<=block[3]+1: block[3]=max(block[3],r.end)
            else:
                if block: runs.append(block)
                block=[pop,chrom,r.start,r.end]
        if block: runs.append(block)
    rows=[]
    for chrom,g in pd.DataFrame(runs,columns=["population","chromosome","start","end"]).groupby("chromosome",sort=False):
        block=None
        for r in g.sort_values(["start","end"]).itertuples(index=False):
            if block and r.start<=block["end"]: block["end"]=max(block["end"],r.end); block["support"].add(r.population)
            else:
                if block: rows.append(block)
                block={"chromosome":chrom,"start":r.start,"end":r.end,"support":{r.population}}
        if block: rows.append(block)
    out=pd.DataFrame(rows); out["farmed_populations"]=out.support.map(lambda x:";".join(sorted(set(FARMS)&x))); out["wild_support"]=out.support.map(lambda x:WILD in x)
    out["selection_class"]=["farmed_and_wild" if f and w else "farmed_only" if f else "wild_only" for f,w in zip(out.farmed_populations.astype(bool),out.wild_support)]
    out=out.drop(columns="support").sort_values(["selection_class","chromosome","start"]).reset_index(drop=True); out.insert(0,"selection_region_id",[f"D5_{i+1:05d}" for i in range(len(out))]); return out

def main():
    here=Path(__file__).resolve().parent; root=here.parent; ap=argparse.ArgumentParser(); ap.add_argument("--outdir",type=Path,default=here/"results"); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    loci=pd.read_csv(root/"shared_data/selection/tajima_D_merged_regions.tsv",sep="\t",dtype={"chromosome":str}); loci=loci[loci.selection_tail.round(3).eq(.05)].copy()
    audit=pd.read_csv(root/"shared_data/selection/farm_region_serbia_spatial50_audit.tsv",sep="\t"); reclassified=set(audit.loc[audit.exclude_from_farm_specific,"selection_region_id"]); hit=loci.selection_region_id.isin(reclassified); loci.loc[hit,"selection_class"]="farmed_and_wild"; loci.loc[hit,"selection_class_label"]="Farmed + wild"; loci.loc[hit,"wild_support"]=True
    loci.to_csv(a.outdir/"selection_loci_5pct.tsv",sep="\t",index=False)
    windows=pd.read_csv(root/"shared_data/selection/tajima_D_outlier_windows.tsv",sep="\t",dtype={"chromosome":str}); windows=windows[windows.selection_tail.round(3).eq(.05)].copy(); windows["farmed_populations"]=windows.farmed_populations.fillna("")
    overlapaware=merge_population_runs(windows); overlapaware.to_csv(a.outdir/"overlapaware_selection_loci_5pct.tsv",sep="\t",index=False)
    summary=loci.groupby("selection_class",sort=False).size().rename("n_loci").reset_index(); summary["percent"]=100*summary.n_loci/len(loci); summary.to_csv(a.outdir/"selection_locus_counts.tsv",sep="\t",index=False)
    pd.read_csv(root/"shared_data/selection/population_tail_enrichment.tsv",sep="\t").to_csv(a.outdir/"coding_polymorphism_enrichment.tsv",sep="\t",index=False)
    print(summary.to_string(index=False)); print(f"total_loci={len(loci)} adaptive_introgression_farm_only_loci={(overlapaware.selection_class=='farmed_only').sum()}")
if __name__=="__main__": main()
