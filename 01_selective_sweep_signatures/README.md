# Selective-sweep signatures in farmed populations

This subsection rebuilds the lower-5% Tajima's-D locus classes used in the
manuscript from compact predigested interval tables. Adjacent windows are
merged within their initial farm/wild class. A nominal farm-only locus is then
reclassified as farm-and-wild when a component farm window overlaps at least
50% of a callable Serbian low-D window. USA is retained in the VCF but excluded
from the farm-versus-wild locus classes.

```bash
python run_analysis.py
```

Expected primary counts are 1,659 loci: 672 farm-only, 774 farm-and-wild, and
213 wild-only. The copied coding-polymorphism enrichment table is the final
population-stratified result used for the manuscript subsection.
