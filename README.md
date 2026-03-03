# Applied Statistics for Microbiome & Genomics Data

### A 7-Day Practical Series

> **Real data · Real R code · Real biological interpretation**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Series](https://img.shields.io/badge/Series-7%20Days-teal)](https://github.com/jojyjohn28/microbiome_stats)

This repository contains all data, R scripts, and figures for a 7-day blog series on applied statistics for microbiome and genomics data. Every method is demonstrated using real Metagenome-Assembled Genomes (MAGs) from estuarine environments — the same dataset runs through all 7 days.

📖 **Blog series:** [your-blog-url-here]  
🐦 **Author:** Jojy John

---

## The Dataset

All analyses use two files that run consistently across every day of the series:

### `mag_abundance.csv`

- **61 MAGs** recovered from estuarine metagenomes
- **27 metagenomes** from Chesapeake Bay and Delaware Bay
- Sampling design: 2 Bays × 2 Seasons × 2 Size Fractions
- MAG orders include: _Pelagibacterales_, _Flavobacteriales_, _Nanopelagicales_, _Pseudomonadales_, _Burkholderiales_, _Rhodobacterales_, _Actinomycetales_, _Cytophagales_

### `metadata.csv`

| Column                 | Type        | Description                         |
| ---------------------- | ----------- | ----------------------------------- |
| `Season`               | Categorical | Spring / Summer                     |
| `Bay`                  | Categorical | Chesapeake / Delaware               |
| `SF`                   | Categorical | Particle-attached / Free_living     |
| `Salinity`             | Continuous  | Practical salinity units            |
| `Temperature`          | Continuous  | °C                                  |
| `Depth`                | Continuous  | metres                              |
| `PAR`                  | Continuous  | Photosynthetically active radiation |
| `Attenuation`          | Continuous  | Light attenuation coefficient       |
| `Bacterial_Production` | Continuous  | µg C L⁻¹ h⁻¹                        |
| `nCells`               | Continuous  | Bacterial cell counts               |
| `ChlA`                 | Continuous  | Chlorophyll _a_ (µg L⁻¹)            |
| `Nitrate`              | Continuous  | µmol L⁻¹                            |
| `Ammonium`             | Continuous  | µmol L⁻¹                            |
| `Phosphate`            | Continuous  | µmol L⁻¹                            |
| `Silicate`             | Continuous  | µmol L⁻¹                            |

---

## Repository Structure

```
microbiome_stats/
│── day1/
├── data/                          # Shared data files (used across all days)
│   ├── mag_abundance.csv
│   └── metadata.csv                        # Why normal statistics fail
│   └── README.md
│
├── day2/                          # RDA & dbRDA
│   ├── RDA_analysis.R
│   └── README.md
│
├── day3/                          # PERMANOVA
│   ├── day3_PERMANOVA.R
│   └── README.md
│
├── day4/                          # Non-parametric tests
│   ├── day4_nonparametric.R
│   └── README.md
│
├── day5/                          # Multiple regression
│   ├── day5_regression.R
│   └── README.md
│
├── day6/                          # WGCNA co-abundance networks
│   ├── day6_WGCNA.R
│   └── README.md
│
├── day7/                          # Full workflow synthesis
│   ├── day7_full_workflow.R
│   └── README.md
│
└── README.md                      # This file
```

---

## The 7-Day Series

### 📅 Monday — Day 1

**Why Normal Statistics Fail for Microbiome Data**  
📖 [Blog](https://jojyjohn28.github.io/blog/day1-why-normal-stats-fail-microbiome/) · 📁 [day1/](https://github.com/jojyjohn28/microbiome_stats/tree/main/day1/data)

The series foundation — not a methods post, but an intuition-building post. Covers what makes microbiome data fundamentally different (compositional, sparse, multivariate, non-normal) and why t-tests, ANOVA, and Pearson correlation fail on this data. Ends with a first look at the MAG abundance and metadata files used throughout the series.

**Concepts:** compositionality, sparsity, multivariate structure, non-normality

---

### 📅 Day 2

**Constrained Ordination: RDA & dbRDA**  
📖 [Blog](https://jojyjohn28.github.io/blog/day2-rda-dbrda-microbiome-ordination/)
Runs the first real analysis: RDA with Hellinger-transformed MAG abundance data and a full set of environmental predictors. Covers forward variable selection, permutation-based significance testing, triplot construction, and biological interpretation. Includes dbRDA with Bray-Curtis as an alternative approach. Detailed section on how to read every element of an ordination triplot.

**R functions:** `rda()`, `decostand()`, `ordiR2step()`, `anova.cca()`, `scores()`, `dbrda()`, `vegdist()`  
**Key packages:** `vegan`, `ggplot2`, `ggrepel`

---

### 📅 Day 3

**PERMANOVA with adonis2 — Testing Real Community Differences**

The natural follow-up to ordination: we've _seen_ communities separate visually — now we formally test whether those differences are statistically real. Covers single-factor and multi-factor PERMANOVA models (Season, Bay, SF and their interactions), the critical homogeneity of dispersion assumption (`betadisper`), and correct interpretation of R² vs p-value.

**R functions:** `adonis2()`, `betadisper()`, `permutest()`, `pairwiseAdonis()`  
**Key packages:** `vegan`, `pairwiseAdonis`

---

### 📅 Day 4

**Non-Parametric Tests — Wilcoxon, Kruskal-Wallis & Effect Sizes**

Stepping back from community-level to taxon-level comparisons. Covers when to use Wilcoxon rank-sum vs Kruskal-Wallis, multiple testing correction (FDR/BH), Cliff's delta as effect size, and critically — when NOT to use a t-test on microbiome abundance data. Applied to real comparisons: Chesapeake vs Delaware, Spring vs Summer, Particle-attached vs Free-living.

**R functions:** `wilcox.test()`, `kruskal.test()`, `p.adjust()`, `cliff.delta()`  
**Key packages:** `effsize`, `ggplot2`

---

### 📅 Day 5

**Multiple Regression — Linking MAG Abundance to Environment**

Bridges community-level ecology to individual MAG/gene abundance. Covers linear models for MAG abundance as a function of environmental variables, interaction terms (Salinity × Season), model diagnostics (residual plots, QQ plots, Cook's distance), and biological interpretation of regression coefficients. Uses the same MAG × metadata dataset.

**R functions:** `lm()`, `summary()`, `anova()`, `plot()` (diagnostics), `AIC()`  
**Key packages:** `base R`, `ggplot2`, `broom`

---

### 📅 Day 6

**WGCNA — From Individual MAGs to Co-Abundance Network Modules**

The advanced day — stops asking "which MAG?" and starts asking "which module of co-occurring MAGs?". Covers WGCNA adapted for MAG co-abundance (not gene expression), soft threshold selection, module detection, relating module eigengenes to environmental traits, and visualizing the trait-module heatmap. This application of WGCNA to MAG co-abundance is underrepresented in the literature.

**R functions:** `pickSoftThreshold()`, `blockwiseModules()`, `moduleEigengenes()`, `cor()`  
**Key packages:** `WGCNA`

---

### 📅 Day 7

**From Raw Data to Paper Figures — A Complete Workflow**

The synthesis post. Walks through the complete analysis pipeline from raw MAG abundance + metadata through ordination → PERMANOVA → differential abundance → regression → network — showing how each method answers a different biological question and how they fit together into a results section. Includes the Mantel test (distance-environment correlation) and a reusable analysis checklist.

**R functions:** All methods from Days 1–6 + `mantel()`, `vegdist()`  
**Key packages:** All packages from the full series

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/jojyjohn28/microbiome_stats.git
cd microbiome_stats
```

### 2. Install required R packages

```r
install.packages(c(
  "vegan",       # ordination, PERMANOVA, diversity
  "ggplot2",     # plotting
  "ggrepel",     # non-overlapping text labels
  "effsize",     # Cliff's delta effect size
  "broom",       # tidy model outputs
  "WGCNA"        # co-abundance network analysis
))

# pairwiseAdonis from GitHub
install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
```

### 3. Run any day's script

```r
# Set working directory to the relevant day folder
setwd("path/to/microbiome_stats/day2")

# Source the script
source("RDA_analysis.R")
```

---

## What Every Script Includes

Each day's R script follows the same structure so you always know where to look:

```
1. Libraries
2. Load data
3. Data preparation (with explanation of WHY each step is needed)
4. The analysis (fully annotated)
5. Significance testing
6. Visualization
7. Output / saving figures
```

---

## Citation

If you use this code or data in your research, please cite:

```
Jojy John (2026). Applied Statistics for Microbiome & Genomics Data —
A 7-Day Practical Series. GitHub: github.com/jojyjohn28/microbiome_stats
```

---

Ongoing
Last updated March 3
