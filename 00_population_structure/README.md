# Chromosome-scale genome assembly and population structure

`run_population_structure.R` is the public, eight-population version of the
original `poolfstat8pop.Rmd` analysis. It reads the single unannotated Pool-seq
VCF in `../shared_data/8pop_poolseq.vcf.gz`; no private 25-population file is
needed.

Run from this directory:

```bash
Rscript run_population_structure.R
```

The script applies the manuscript MAF and coverage filters, then writes PCA
coordinates/plots, pairwise FST, expected heterozygosity, and farm-versus-wild
Wald contrasts. Each pool represents 40 diploid individuals (`poolsizes=80`).

Required R packages: `poolfstat` 3.1.0, `FactoMineR`, `factoextra`, `pheatmap`,
`RColorBrewer`, and `ggplot2`.
