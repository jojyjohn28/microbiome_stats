# Day 1 — Why Normal Statistics Fail for Microbiome Data

**Series:** Applied Statistics for Microbiome & Genomics Data  
**Blog post:** [Read Day 1](https://your-blog-url/day1-why-normal-stats-fail-microbiome)

---

No R code today — this is the conceptual foundation for the entire series.

Covers why standard tests (t-test, ANOVA, Pearson) break down for microbiome data: compositionality, sparsity, multivariate structure, and non-normality. Ends with a first look at the dataset used across all 7 days.

## Data

```
day1/data/
├── mag_abundance.csv   # 61 MAGs × 27 metagenomes (raw coverage)
└── metadata.csv        # 27 samples — Season, Bay, SF + 12 env variables
```

These same two files are used in every subsequent day.

---

_← [Series README](../README.md) | [Day 2 →](../day2-rda-dbrda/)_
