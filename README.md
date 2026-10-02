# Code and data availability

Analysis code and predigested data for the study of early domestication in
the yellow mealworm, *Tenebrio molitor*, organized by Results subsection.

## Directory guide

| Directory | Purpose |
|---|---|
| [00_population_structure](00_population_structure/) | PCA, pairwise FST, heterozygosity, and farm-versus-Serbia comparisons from the VCF. |
| [01_selective_sweep_signatures](01_selective_sweep_signatures/) | Merge low-Tajima's-D windows, classify farm and Serbia signals, and summarize coding-polymorphism enrichment. |
| [02_adaptive_introgression](02_adaptive_introgression/) | Test enrichment of overlaps between FdM regions and farm-only selection signals using permutations. |
| [03_qtl_vs_genome_signatures](03_qtl_vs_genome_signatures/) | Test Dole7 QTL overlaps with selection and introgression signals, accounting for callable regions and alternative thresholds. |
| [shared_data](shared_data/) | Input windows, intervals, QTLs, markers, gene summaries, and chromosome lengths. |

## Population IDs

Each pool contains 40 diploid beetles (80 chromosome copies).

| VCF ID | Population label | Country | Sample type |
|---|---|---|---|
| AAI | Dole7 | France | Industrial farm |
| AAK | Starfood Holland | Netherlands | Industrial farm |
| AAL | Kingsect Belgium | Belgium | Industrial farm (Belgium 1) |
| AAQ | Wild Serbian | Serbia | Synanthropic non-farmed reference |
| AAU | Pronutrix+ Belgium | Belgium | Industrial farm (Belgium 2) |
| AAW | Papek Czechia | Czechia | Industrial farm |
| AAX | Mystik Canada | Canada | Industrial farm |
| ABE | USA | United States | Laboratory line |

## Variant data

Download `8pop_poolseq.vcf.gz` from
[Zenodo: 10.5281/zenodo.23099934](https://doi.org/10.5281/zenodo.23099934)
and place it in `shared_data/` before running population structure. This
unannotated VCF contains the eight pools in the order above, with pooled
reference and alternate read counts (`RD:AD`). The other analyses use the
included predigested tables.

## Run the analyses

Run from the repository root (locally, `code_data_availability/`):

```bash
Rscript 00_population_structure/run_population_structure.R
python 01_selective_sweep_signatures/run_analysis.py
python 02_adaptive_introgression/run_analysis.py
python 03_qtl_vs_genome_signatures/run_analysis.py
python verify_bundle.py
```

Dependencies: Python with NumPy, pandas, and Matplotlib; R packages are listed
in [subsection 00](00_population_structure/README.md). Analyses use fixed
seeds and 10,000 introgression / 100,000 QTL permutations. Subdirectory READMEs
provide details; outputs are written to their respective `results/` folders.

`verify_bundle.py` checks the included result tables. Add `--full-vcf` to
validate the downloaded VCF completely. File checksums are in `SHA256SUMS.txt`.

## Manuscript

This code is part of the study by Rocha Ferreira et al., *Mosaic selection,
gene flow and the genetic basis of production traits at the onset of grain
beetle domestication* (in preparation).
