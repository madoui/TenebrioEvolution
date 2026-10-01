# Production-trait QTLs are not yet targeted by selection

`run_analysis.py` reproduces the manuscript's Dole7-specific, callability-aware
QTL analysis. It tests the 305 merged ±50-kb QTLs against Dole7-specific low
Tajima's D, low nucleotide diversity, and Dole7-recipient FdM signals at the
matched 10%, 5%, and 1% cutoffs. Null QTLs retain chromosome and width, are
centred on observed GWAS markers, and are matched for callable fraction.

Full manuscript run:

```bash
python run_analysis.py
```

Fast verification run:

```bash
N_PERMUTATIONS=100 SKIP_PLOTS=1 python run_analysis.py
```

The full analysis uses 100,000 placements per test and a fixed seed. The
checked full-run results reported in the manuscript are retained under
`results/`. The two helper modules in this folder are retained because the
primary script imports their interval and input-parsing functions.

`SKIP_PLOTS=1` skips only the PDF/PNG/SVG heatmap and is useful on headless
systems; all tables, tests, reports, and QC checks are still produced.
