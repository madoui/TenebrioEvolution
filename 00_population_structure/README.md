# Chromosome-scale genome assembly and population structure

`run_population_structure.R` is the eight-population version of the original
`poolfstat8pop.Rmd`. It reads `../shared_data/8pop_poolseq.vcf.gz` directly.
The VCF is not stored in this GitHub repository. Download
`8pop_poolseq.vcf.gz` from
[Zenodo (DOI: 10.5281/zenodo.23099934)](https://doi.org/10.5281/zenodo.23099934)
and place it in `../shared_data/` before running the script.

Population IDs and sample types are defined in the main `../README.md`.

```bash
Rscript run_population_structure.R
```

Required packages: poolfstat 3.1.0, FactoMineR, factoextra, pheatmap,
RColorBrewer, and ggplot2. Each pool contains 40 diploid individuals, hence
`poolsizes=80`.
