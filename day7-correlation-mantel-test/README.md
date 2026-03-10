# Day 7 — Correlation Analysis & Mantel Test

Part of the series **Applied Statistics for Microbiome & Genomics Data**
Dataset: 61 MAGs × 27 metagenomes, Chesapeake Bay + Delaware Bay

## Contents

```
day7-correlation/
├── day7_correlation_mantel.R     # full analysis script
├── Env_variable_corr.png         # Panel A: env–env Spearman heatmap
├── MAG_corr.png                  # Panel B: top 20 MAGs × env correlation heatmap
├── Mantel_test_bar_plot.png      # Panel C: Mantel r per env variable
├── mantel_distribution.png       # Panel D: permutation null distribution
├── Mantel_corr.png               # community vs environment distance scatter
└── Correlation_results.png       # combined results figure
```

## What this does

Two complementary analyses linking MAG abundance to environmental gradients:

**Correlation** — Spearman r between each of 61 MAGs and 9 environmental variables (Salinity, Temperature, Bacterial Production, Cell Count, Chl-a, Nitrate, Ammonium, Phosphate, Silicate). BH FDR correction applied across all 549 tests.

**Mantel test** — tests whether Bray-Curtis community dissimilarity scales with multivariate environmental dissimilarity across all sample pairs. Run for the full env matrix, per individual variable, and as a partial Mantel (Salinity | Temperature).

## Dependencies

```r
install.packages(c("vegan", "Hmisc", "corrplot", "ggcorrplot",
                   "ggplot2", "patchwork", "reshape2", "dplyr"))
```

## Blog post

[Day 7 — Correlation Analysis & Mantel Test](https://jojyjohn28.github.io/blog/day7-correlation-mantel-test-microbiome/)
