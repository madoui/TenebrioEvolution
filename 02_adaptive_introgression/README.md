# Adaptive introgression in farmed populations

This analysis intersects merged upper-5% positive-FdM regions with the
overlap-aware farm-only lower-5% Tajima's-D loci produced by subsection 01. It then performs
10,000 chromosome- and interval-length-preserving random placements of the FdM
regions. The fixed seed makes the empirical result reproducible.

```bash
python ../01_selective_sweep_signatures/run_analysis.py
python run_analysis.py
```

Use `--permutations 100` for a quick smoke test. The manuscript run uses the
default 10,000 permutations.

The compact predigested coding layer is in `../shared_data/introgression/`:
direction-consistent SNV records (gzip-compressed), the 70 resolved
nonsynonymous variants/38 genes table, exact-overlap genes, direct GO results,
and the final manuscript statistic table. No SnpEff-annotated VCF is included.
