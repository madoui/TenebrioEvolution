# Code and predigested data for the manuscript

This directory is the compact public reproducibility package for the four
Results subsections. It contains one unannotated eight-population Pool-seq VCF
and only the predigested interval/statistical tables required by the analyses.
The private 25-population VCF, SnpEff VCFs, genome FASTA, large GWAS archive,
and redundant figure formats are intentionally excluded.

## Directory structure

```text
00_population_structure/             Chromosome-scale genome assembly and population structure
01_selective_sweep_signatures/        Selective-sweep signatures in farmed populations
02_adaptive_introgression/            Adaptive introgression in farmed populations
03_qtl_vs_genome_signatures/          Production-trait QTLs are not yet targeted by selection
shared_data/                          Inputs shared across Results subsections
```

Each subsection has its own README, runnable code, and `results/` directory.
Run subsection 01 before subsection 02 because the latter consumes its
overlap-aware low-D5 locus table.

## Eight-population VCF

`shared_data/8pop_poolseq.vcf.gz` contains, in this order:

| ID | Manuscript population | Type |
|---|---|---|
| AAI | Dole7 | farm |
| AAK | Starfood Holland | farm |
| AAL | Kingsect Belgium | farm |
| AAQ | Wild Serbian | wild |
| AAU | Pronutrix+ Belgium | farm |
| AAW | Papek Czechia | farm |
| AAX | Mystik Canada | farm |
| ABE | USA | laboratory |

Only chromosome, position, alleles, and pooled `RD:AD` fields are retained.
There are no SnpEff `ANN`/`EFF` fields. All pools contain DNA from 40 diploid
individuals, so the population-structure script uses pool size 80.

The VCF is suitable for Zenodo. Because it exceeds GitHub's ordinary 100-MB
file limit, either store it only in the linked Zenodo record and document its
DOI in the GitHub release, or track it with Git LFS. Do not split or duplicate
the VCF across subsection folders.

## Software

- Python 3.10 or later
- NumPy, pandas, and Matplotlib (see `requirements.txt`)
- R 4.3.1 and the packages listed in `00_population_structure/README.md`

All randomization scripts use fixed seeds. The adaptive-introgression default
is 10,000 permutations; the QTL default is 100,000 placements per test.
Reduced-permutation commands in the subsection READMEs are smoke tests only.

## Reproduction order

```bash
python 01_selective_sweep_signatures/run_analysis.py
python 02_adaptive_introgression/run_analysis.py
python 03_qtl_vs_genome_signatures/run_analysis.py
Rscript 00_population_structure/run_population_structure.R
python verify_bundle.py --full-vcf
```

The expected selection counts are 672 farm-only, 774 farm-and-wild, and 213
wild-only loci. The adaptive-introgression analysis expects 68 of 163 FdM5
regions to overlap 96 overlap-aware farm-only low-D5 loci (2.06-fold,
empirical P=0.0001). The checked 100,000-placement QTL results are retained in
`03_qtl_vs_genome_signatures/results/`.

## Public-release notes

- Add the final manuscript citation, GitHub URL, Zenodo DOI, authors, and an
  explicit software license before publication.
- Keep `SHA256SUMS.txt` with the Zenodo upload and regenerate it after any file
  changes.
- Cite both GitHub (versioned code) and Zenodo (immutable data/release archive)
  in the manuscript Data and Code Availability statement.
