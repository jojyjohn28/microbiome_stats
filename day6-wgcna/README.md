# Day 6 — WGCNA: MAG Co-Abundance Networks

Part of the series **Applied Statistics for Microbiome & Genomics Data**
Dataset: 61 MAGs × 27 metagenomes, Chesapeake Bay + Delaware Bay

## Contents

```
day6-wgcna/
├── day6_wgcna.R                          # full analysis script
├── figures/
│   ├── day6_wgcna_figure.tiff            # 3-panel figure (dendrogram · TOM · module-trait heatmap)
│   └── scale_vs_connectivity.png         # soft threshold diagnostic plots
└── results/
    ├── day6_MAGInfo.csv                  # per-MAG module color, GS, and MM scores
    └── day6_networkConstruction.RData    # TaxaTree, dissTOM, moduleColors, mergedMEs
```

## What this does

Identifies MAG co-abundance modules using WGCNA (Weighted Gene Co-expression Network Analysis) and correlates each module eigengene with environmental variables.

Key parameters: Hellinger transformation · signed Pearson · β = 6 · TOM clustering · `minClusterSize = 5` · merge cutHeight = 0.25

## Dependencies

```r
install.packages(c("vegan", "WGCNA", "ggplot2", "ggdendro",
                   "patchwork", "reshape2", "dplyr"))
```

## Blog post

[Day 6 — WGCNA for Microbiome Co-Abundance Networks](https://jojyjohn28.github.io/blog/day6-wgcna-mag-coabundance-networks/)
