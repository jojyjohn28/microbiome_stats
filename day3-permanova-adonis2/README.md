# Day 3 — PERMANOVA with adonis2

**Series:** Applied Statistics for Microbiome & Genomics Data  
**Blog post:** [Read Day 3](https://your-blog-url/day3-permanova-adonis2-microbiome)

---

You saw the groups separate in the RDA triplot. This day tests whether that separation is statistically real.

## Files

```
day3--permanova/
├── day3_PERMANOVA.R          # Full annotated R script
├── day3_permanova_figure.tiff   # PCoA + betadisper figure
└── results.md                # Real output with 1-2 line interpretations
```

Data lives in `../day1/data/` — script points there.

## What the script does

1. Calculates Bray-Curtis dissimilarity matrix
2. Single-factor PERMANOVA — Season, Bay, Size Fraction
3. Multi-factor additive model with marginal SS (`by = "margin"`)
4. Interaction models — Season × Bay, Season × SF
5. Betadisper — tests homogeneity of dispersion assumption
6. Pairwise PERMANOVA across Season × Bay combinations
7. Two-panel figure — PCoA + betadisper boxplot

## Key result

Season was the only significant driver (R² = 0.325, p = 0.001). Bay and size fraction were not significant. Betadisper confirmed homogeneous dispersion across all groups — the PERMANOVA result is trustworthy.

## Quick start

```r
setwd("path/to/microbiome_stats/day3--permanova")
source("day3_PERMANOVA.R")
```

**Requires:** `vegan`, `ggplot2`, `patchwork`, `pairwiseAdonis`

---

_← [Day 2](../day2--rda-dbrda/) | [Series README](../README.md) | Day 4 →_
