# Day 3 — PERMANOVA Results Summary

**Dataset:** 61 MAGs × 27 metagenomes | Chesapeake Bay & Delaware Bay | Spring & Summer  
**Distance metric:** Bray-Curtis dissimilarity  
**Test:** PERMANOVA (`adonis2`) with 999 permutations | `by = "margin"` (marginal SS)

---

## Bray-Curtis Distance Distribution

| Statistic | Value |
| --------- | ----- |
| Min       | 0.035 |
| Median    | 0.881 |
| Mean      | 0.733 |
| Max       | 1.000 |

Communities are highly dissimilar overall — median pairwise Bray-Curtis of 0.88 indicates that most sample pairs share very few MAGs, which is typical for estuarine metagenomes spanning different seasons and locations.

---

## Single-Factor PERMANOVA

| Factor             | R²    | F     | p-value      | Interpretation            |
| ------------------ | ----- | ----- | ------------ | ------------------------- |
| Season             | 0.321 | 11.81 | 0.001 \*\*\* | Strong significant effect |
| Bay                | 0.044 | 1.15  | 0.316        | Not significant           |
| Size Fraction (SF) | 0.019 | 0.49  | 0.771        | Not significant           |

**Season is the dominant driver**, explaining 32% of total community variance on its own. Bay and size fraction show no significant independent effect.

---

## Multi-Factor Additive Model (Marginal SS)

`adonis2(bc_dist ~ Season + Bay + SF, by = "margin")`

| Factor   | R²    | F     | p-value      | Interpretation                                      |
| -------- | ----- | ----- | ------------ | --------------------------------------------------- |
| Season   | 0.325 | 12.22 | 0.001 \*\*\* | Remains dominant after accounting for other factors |
| Bay      | 0.049 | 1.84  | 0.121        | Marginal, non-significant                           |
| SF       | 0.019 | 0.70  | 0.550        | Not significant                                     |
| Residual | 0.612 | —     | —            | Unexplained variance                                |

Season remains the only significant factor in the full additive model. Bay trends toward significance (p = 0.121) but does not cross the threshold, suggesting estuarine identity plays a minor role compared to seasonal dynamics. Approximately 61% of community variance remains unexplained — typical for field metagenomics data.

---

## Interaction Models

| Interaction       | R²     | F     | p-value | Interpretation  |
| ----------------- | ------ | ----- | ------- | --------------- |
| Season × Bay      | 0.032  | 1.22  | 0.303   | Not significant |
| Season × SF       | 0.035  | 1.29  | 0.257   | Not significant |
| Season × Bay × SF | 0.0001 | 0.005 | 1.000   | Not significant |

No significant interactions detected. The effect of Season on community composition does not depend on which Bay was sampled, nor on size fraction. The three-way interaction is effectively zero (R² = 0.0001), consistent with expectations given only 2–6 samples per cell.

---

## Homogeneity of Dispersion (betadisper)

| Factor        | F      | p-value | Interpretation |
| ------------- | ------ | ------- | -------------- |
| Season        | 1.42   | 0.277   | Homogeneous ✅ |
| Bay           | 0.023  | 0.882   | Homogeneous ✅ |
| Size Fraction | 0.0004 | 0.983   | Homogeneous ✅ |

All three factors show homogeneous dispersion — community spread is similar across groups. **This confirms that the significant Season PERMANOVA result reflects genuine differences in community composition (centroid shifts), not differences in variability.** The PERMANOVA results are fully interpretable.

---

## Pairwise PERMANOVA (Season × Bay combinations)

| Comparison                             | R²    | F    | p-value    | Significant? |
| -------------------------------------- | ----- | ---- | ---------- | ------------ |
| Spring Chesapeake vs Summer Chesapeake | 0.421 | 6.54 | 0.003 \*\* | Yes          |
| Spring Chesapeake vs Spring Delaware   | 0.239 | 1.89 | 0.129      | No           |
| Spring Chesapeake vs Summer Delaware   | 0.275 | 5.31 | 0.007 \*\* | Yes          |
| Summer Chesapeake vs Spring Delaware   | 0.497 | 8.91 | 0.004 \*\* | Yes          |
| Summer Chesapeake vs Summer Delaware   | 0.095 | 1.78 | 0.157      | No           |
| Spring Delaware vs Summer Delaware     | 0.345 | 7.38 | 0.004 \*\* | Yes          |

**Key pattern:** Seasonal differences are significant within both bays (Spring vs Summer Chesapeake: p = 0.003; Spring vs Summer Delaware: p = 0.004). Within-season comparisons between bays are non-significant (Spring Chesapeake vs Spring Delaware: p = 0.129; Summer Chesapeake vs Summer Delaware: p = 0.157), confirming that **Season structures community composition more strongly than geographic location.**

---

## Summary Statement

> Season is the primary driver of MAG community composition across Chesapeake and Delaware Bay metagenomes, explaining 32% of total community variance (R² = 0.321, p = 0.001). Bay identity and size fraction did not significantly explain community differences. All betadisper tests confirmed homogeneous dispersion, validating the PERMANOVA results. Pairwise comparisons show that spring and summer communities differ significantly within each estuary, while same-season communities between estuaries are not significantly different.

---

## Figure

Read more at [Day2 Blog](https://jojyjohn28.github.io/blog/day3-permanova-adonis2-microbiome/)

_Part of the 7-day series: Applied Statistics for Microbiome & Genomics Data_  
_← [Day 2 RDA](../day2--rda-dbrda/) | [Day 4](../day4/) →_
