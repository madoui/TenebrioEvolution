# Dole7-specific, callability-aware QTL overlap analysis

The primary analysis compared the 305 French-population QTLs with Dole7-specific low Tajima's D and low nucleotide diversity (pi), and with high-fDM windows from contrasts in which Dole7 is the mapped recipient. Missing windows were treated as non-testable rather than as absence of signal.

For D and pi, the callable universe was the intersection of Dole7 and Wild Serbian callable windows because both values are required to classify a signal as farm-specific. For fDM, the primary universe was the intersection of valid, lifted windows across the Dole7-recipient contrasts: France_Netherlands_vs_Belgium2, France_Netherlands_vs_Czechia. Multipartite masks were intersections of their component masks.

Each null distribution used 100,000 placements per QTL. Random QTLs retained chromosome and width, were centred on SNPs in the complete GWAS scan, and were matched to the observed QTL's callable fraction. QTLs with no callable base were excluded.

Benjamini-Hochberg correction was applied across all displayed primary Dole7 tests; panel-specific q values are also provided in the result table.

## All-QTL results

- D, 1%: 5 of 282 callable QTLs overlapped (null mean 8.49; depletion; P=0.2637; global BH q=1).
- D, 5%: 18 of 282 callable QTLs overlapped (null mean 28.89; depletion; P=0.01586; global BH q=1).
- D, 10%: 35 of 282 callable QTLs overlapped (null mean 48.88; depletion; P=0.01906; global BH q=1).
- D-fDM, 1%: 0 of 282 callable QTLs overlapped (null mean 0.26; depletion; P=1; global BH q=1).
- D-fDM, 5%: 0 of 282 callable QTLs overlapped (null mean 3.06; depletion; P=0.090719; global BH q=1).
- D-fDM, 10%: 2 of 282 callable QTLs overlapped (null mean 8.94; depletion; P=0.01116; global BH q=1).
- D-pi, 1%: 2 of 282 callable QTLs overlapped (null mean 3.67; depletion; P=0.54273; global BH q=1).
- D-pi, 5%: 9 of 282 callable QTLs overlapped (null mean 12.96; depletion; P=0.25642; global BH q=1).
- D-pi, 10%: 20 of 282 callable QTLs overlapped (null mean 18.21; enrichment; P=0.71997; global BH q=1).
- D-pi-fDM, 1%: 0 of 282 callable QTLs overlapped (null mean 0.00; enrichment; P=1; global BH q=1).
- D-pi-fDM, 5%: 0 of 282 callable QTLs overlapped (null mean 1.24; depletion; P=0.57073; global BH q=1).
- D-pi-fDM, 10%: 1 of 282 callable QTLs overlapped (null mean 3.56; depletion; P=0.2524; global BH q=1).
- fDM, 1%: 2 of 305 callable QTLs overlapped (null mean 3.87; depletion; P=0.51035; global BH q=1).
- fDM, 5%: 10 of 305 callable QTLs overlapped (null mean 15.20; depletion; P=0.20106; global BH q=1).
- fDM, 10%: 23 of 305 callable QTLs overlapped (null mean 33.30; depletion; P=0.058039; global BH q=1).
- pi, 1%: 20 of 282 callable QTLs overlapped (null mean 15.55; enrichment; P=0.28248; global BH q=1).
- pi, 5%: 54 of 282 callable QTLs overlapped (null mean 54.84; depletion; P=0.96359; global BH q=1).
- pi, 10%: 99 of 282 callable QTLs overlapped (null mean 79.26; enrichment; P=0.0071799; global BH q=1).
- pi-fDM, 1%: 0 of 282 callable QTLs overlapped (null mean 0.24; depletion; P=1; global BH q=1).
- pi-fDM, 5%: 2 of 282 callable QTLs overlapped (null mean 3.57; depletion; P=0.60889; global BH q=1).
- pi-fDM, 10%: 8 of 282 callable QTLs overlapped (null mean 11.51; depletion; P=0.35978; global BH q=1).

## Callability

- D: 282 of 305 QTLs callable; median callable fraction 1.000; mask 225.14 Mb.
- pi: 282 of 305 QTLs callable; median callable fraction 1.000; mask 221.94 Mb.
- fDM: 305 of 305 QTLs callable; median callable fraction 1.000; mask 261.18 Mb.
- D-pi: 282 of 305 QTLs callable; median callable fraction 1.000; mask 221.24 Mb.
- D-fDM: 282 of 305 QTLs callable; median callable fraction 1.000; mask 219.51 Mb.
- pi-fDM: 282 of 305 QTLs callable; median callable fraction 1.000; mask 217.02 Mb.
- D-pi-fDM: 282 of 305 QTLs callable; median callable fraction 1.000; mask 216.38 Mb.

Complete outputs include QTL-level callability, null-matching diagnostics, overlap pairs, empirical P values, q values, and QC checks.
