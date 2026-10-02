# Production-trait QTLs are not yet targeted by selection

`python run_analysis.py` reproduces the Dole7-specific, callability-aware QTL
analysis with 100,000 placements per test. For a smoke test use:

```bash
N_PERMUTATIONS=100 SKIP_PLOTS=1 python run_analysis.py
```

The checked full-run tables are retained under `results/`.

The other Python files are helper modules; use `run_analysis.py` as the entry
point for this package.
