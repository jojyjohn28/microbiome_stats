# Day 5 — Multiple Regression for MAG Abundance: Linking Genomics to Environment

**Series:** Applied Statistics for Microbiome & Genomics Data  
**Blog post:** [Read Day 5](https://jojyjohn28.github.io/blog/day5-multiple-regression-mag-abundance/)

---

Days 2–4 focused on **discrete group comparisons**. Day 5 moves to **continuous environmental gradients**, asking how MAG abundance changes with **salinity**, **temperature**, and **season** using multiple regression.

This is the bridge from “which groups differ?” to **“how much does abundance change per unit environment?”**

## Files

```text
day5-regression/
├── day5_regression.R         # Full annotated R script
├── day5_all_MAGs.csv         # Regression summary for all 61 MAGs
├── summary_model.txt         # Main model summaries and notes
├── LR_scatter_plot.png       # Example scatter plots with linear fits
├── lr_1.svg                  # Vector version of regression figure
├── model_stat.png            # Model statistics / diagnostics figure
└── Sality-season.png         # Interaction-style plot (rename to Salinity-season.png if desired)
```

#### What the script does

1. Loads metadata and MAG abundance table

2. Aligns samples and applies log1p() transformation to MAG abundances

3. Fits simple and additive linear models for focal MAGs

4. Tests predictors: Salinity, Temperature, and Season

5. Checks multicollinearity using VIF

Evaluates model assumptions using:

● residuals vs fitted

● QQ plot

● scale-location plot

● residual histogram

● Shapiro-Wilk test

● Breusch-Pagan test

6.Tests an interaction model for selected MAGs

7. Loops regression across all 61 MAGs

8. Applies FDR correction to salinity, temperature, and season p-values

9. Exports regression results and figures

#### Quick start

```R
setwd("path/to/microbiome_stats/day5-regression")
source("day5_regression.R")
```

Read more at : **[Blog day5](https://jojyjohn28.github.io/blog/day5-multiple-regression-mag-abundance/)**
