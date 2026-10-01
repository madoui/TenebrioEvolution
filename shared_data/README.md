# Shared predigested data

- `8pop_poolseq.vcf.gz`: unannotated eight-population Pool-seq VCF with `RD:AD`.
- `selection/`: callable Tajima's-D and nucleotide-diversity windows, compact
  classified interval tables, and final gene-level coding-polymorphism metrics.
- `introgression/`: lifted 100-kb D/FdM windows and merged FdM outlier regions.
- `qtl/`: 305 merged QTL intervals, significant marker annotations, and the
  deduplicated GWAS marker-position universe used for null placement.
- `genome/chromosome_lengths.tsv`: chromosome/scaffold bounds used for interval
  randomization; the 287-MB reference FASTA is not duplicated in this deposit.

`make_8pop_vcf.py` documents the deterministic extraction of the public VCF
from the private headered 25-pool call set. The source file itself is not part
of the public package.
