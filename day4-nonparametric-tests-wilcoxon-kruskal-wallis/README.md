# Day 4 — Non-Parametric Tests: Wilcoxon, Kruskal-Wallis & Effect Sizes

**Series:** Applied Statistics for Microbiome & Genomics Data  
**Blog post:** [Read Day 4](https://jojyjohn28.github.io/blog/day4-nonparametric-tests-wilcoxon-kruskal-wallis/)

---

PERMANOVA told us communities differ by season. This day identifies _which specific MAGs_ are driving that difference — and whether the differences are biologically meaningful.

## Files

```
day4--nonparametric/
├── day4_nonparametric.R              # Full annotated R script
├── day4_nonparametric_figure.tiff   # Volcano plot + boxplots
│
└── day4_wilcoxon_results.csv     # All 61 MAGs × 3 comparisons

```

Data lives in `../day1/data/` — script points there.

## What the script does

1. Wilcoxon rank-sum tests — Season, Size Fraction, Salinity regime (<15 vs ≥15 psu)
2. FDR correction (Benjamini-Hochberg) across all 61 MAGs per comparison
3. Cliff's delta effect size for every test
4. Kruskal-Wallis across three salinity zones with post-hoc pairwise tests
5. Volcano plot — effect size vs significance for all MAGs
6. Boxplots of top 6 most significant MAGs

## Why these tests?

With 61 MAGs × 3 comparisons = 183 tests, uncorrected p-values give ~9 false positives by chance. FDR correction (BH method) keeps the false discovery rate at 5% across the full set.

## Quick start

```r
setwd("path/to/microbiome_stats/day4--nonparametric")
source("day4_nonparametric.R")
```

**Requires:** `vegan`, `ggplot2`, `ggrepel`, `effsize`, `dplyr`, `tidyr`, `patchwork`

---

_← [Day 3](../day3--permanova/) | [Series README](../README.md) | Day 5 →_
