# Day 2 — Constrained Ordination: RDA & dbRDA

**Series:** Applied Statistics for Microbiome & Genomics Data  
**Blog post:** [Day 2 — Constrained Ordination: RDA & dbRDA](https://your-blog-url-here)  
**Series repo:** [github.com/jojyjohn28/microbiome_stats](https://github.com/jojyjohn28/microbiome_stats)

---

## What This Day Covers

On Day 1 we established why standard statistics fail for microbiome data. Today we run the first real analysis: **constrained ordination** — connecting community composition patterns directly to measured environmental drivers.

By the end of this day you will know how to:

- Transform MAG abundance data for RDA (Hellinger transformation)
- Run forward selection to find significant environmental predictors
- Build and test an RDA model with permutation-based significance tests
- Extract and plot a full triplot: sample points, MAG scores, environmental arrows
- Run dbRDA with Bray-Curtis as an alternative approach
- Read and interpret every element of the ordination figure

---

## Files in This Directory

```
day2/
├── RDA_analysis.R                # Full annotated R script
├── RDA_triplot.tiff              # Final publication figure
│
└── README.md                 # This file
```

Source of image : https://academic.oup.com/ismecommun/article/6/1/ycag021/8454622, see supplementary information

---

---

## How to Run

```r
# 1. Set working directory to the day2 folder
setwd("path/to/microbiome_stats/day2")

# 2. Install required packages if needed
install.packages(c("vegan", "ggplot2", "ggrepel"))

# 3. Run the script
source("day2_RDA.R")
```

---

## Reading the Triplot

Read detailed instruction in [Day 2 Blog](https://jojyjohn28.github.io/blog/day2-rda-dbrda-microbiome-ordination/)

---

_Part of the 7-day series: Applied Statistics for Microbiome & Genomics Data_  
_← [Day 1](../day1/README.md) | [Day 3](../day3/README.md) →_
