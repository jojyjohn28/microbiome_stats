## ============================================================
## Day 3 — PERMANOVA with adonis2
## Applied Statistics for Microbiome & Genomics Data
## Author: Jojy John
## Series: github.com/jojyjohn28/microbiome_stats
## ============================================================

## ── Libraries ────────────────────────────────────────────────
library(vegan)
library(ggplot2)
library(ggrepel)
#install.packages("devtools")
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

## ── 1. Load data ─────────────────────────────────────────────
metadata <- read.csv("~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/metadata.csv",
                     header = TRUE, row.names = 1, sep = ",")

mag_raw  <- read.csv("~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/mag_abundance.csv",
                     header = TRUE, row.names = 1, sep = ",")

# Transpose: samples as rows, MAGs as columns
mag_t <- t(mag_raw)   # 27 × 61

# Reorder metadata rows to match the column order of mag_abundance
# (same 27 samples, different sort order between files)
metadata <- metadata[match(rownames(mag_t), rownames(metadata)), ]

# Now verify alignment is correct
stopifnot(all(rownames(mag_t) == rownames(metadata)))
# Should pass silently — if still errors, run the diagnostic below


# See exactly which names differ
mismatch <- which(rownames(mag_t) != rownames(metadata))
data.frame(
  mag_t    = rownames(mag_t)[mismatch],
  metadata = rownames(metadata)[mismatch]
)

# Verify alignment — rows must match between files
stopifnot(all(rownames(mag_t) == rownames(metadata)))

## ── 2. Extract group variables ───────────────────────────────
Season   <- metadata$Season    # Spring (n=8), Summer (n=19)
Bay      <- metadata$Bay       # Chesapeake (n=11), Delaware (n=16)
Fraction <- metadata$SF        # Particle-attached (n=14), Free_living (n=13)

# Design summary
table(Season, Bay)
table(Season, Fraction)

## ── 3. Bray-Curtis dissimilarity ─────────────────────────────
# Standard choice for community composition:
# - ignores shared absences
# - bounded 0 (identical) to 1 (no shared MAGs)
# - robust to sparsity and compositionality
bc_dist <- vegdist(mag_t, method = "bray")

# Distribution of pairwise distances
summary(as.vector(bc_dist))

## ── 4. Single-factor PERMANOVA ───────────────────────────────
# Always set seed for reproducibility
set.seed(123)
perm_season <- adonis2(bc_dist ~ Season,
                       data = metadata,
                       permutations = 999)
print(perm_season)
# Key: R² (effect size) and Pr(>F) (permutation p-value)

set.seed(123)
perm_bay <- adonis2(bc_dist ~ Bay,
                    data = metadata,
                    permutations = 999)
print(perm_bay)

set.seed(123)
perm_sf <- adonis2(bc_dist ~ SF,
                   data = metadata,
                   permutations = 999)
print(perm_sf)

## ── 5. Multi-factor additive model ───────────────────────────
# by = "margin" = Type III SS — each term tested after all others
# Essential for unbalanced designs (Spring n=8 vs Summer n=19)
set.seed(123)
perm_multi <- adonis2(bc_dist ~ Season + Bay + SF,
                      data = metadata,
                      permutations = 999,
                      by = "margin")
print(perm_multi)

## ── 6. Interaction models ─────────────────────────────────────
# Season × Bay: does seasonal effect differ between estuaries?
set.seed(123)
perm_season_bay <- adonis2(bc_dist ~ Season * Bay,
                           data = metadata,
                           permutations = 999,
                           by = "margin")
print(perm_season_bay)

# Season × SF: does fraction matter differently per season?
set.seed(123)
perm_season_sf <- adonis2(bc_dist ~ Season * SF,
                          data = metadata,
                          permutations = 999,
                          by = "margin")
print(perm_season_sf)

# Full three-way (interpret cautiously — only 2-6 samples per cell)
set.seed(123)
perm_full <- adonis2(bc_dist ~ Season * Bay * SF,
                     data = metadata,
                     permutations = 999,
                     by = "margin")
print(perm_full)

## ── 7. Betadisper — homogeneity of dispersion ────────────────
# CRITICAL: tests whether PERMANOVA result reflects centroid
# differences (real composition change) vs spread differences
# (one group is just more variable)

# ── Season ──
disp_season <- betadisper(bc_dist, metadata$Season)
set.seed(123)
permutest(disp_season, permutations = 999)
# p > 0.05 = homogeneous dispersion → PERMANOVA result valid
# p < 0.05 = heterogeneous → interpret with caution

# ── Bay ──
disp_bay <- betadisper(bc_dist, metadata$Bay)
set.seed(123)
permutest(disp_bay, permutations = 999)

# ── Size Fraction ──
disp_sf <- betadisper(bc_dist, metadata$SF)
set.seed(123)
permutest(disp_sf, permutations = 999)

## ── 8. Pairwise PERMANOVA ─────────────────────────────────────
# Needed when a factor has >2 levels, or to compare specific
# group combinations (e.g. Spring_Chesapeake vs Summer_Delaware)
metadata$SeasonBay <- paste(metadata$Season, metadata$Bay, sep = "_")

set.seed(123)
pairwise_result <- pairwise.adonis2(bc_dist ~ SeasonBay,
                                    data = metadata,
                                    permutations = 999)
print(pairwise_result)

## ── 9. Figure — Panel A: PCoA ────────────────────────────────
pcoa_result <- cmdscale(bc_dist, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  PC1    = pcoa_result$points[, 1],
  PC2    = pcoa_result$points[, 2],
  Season = metadata$Season,
  Bay    = metadata$Bay,
  SF     = metadata$SF
)

# % variance explained
eig_pct <- round(100 * pcoa_result$eig /
                   sum(pcoa_result$eig[pcoa_result$eig > 0]), 1)

# Grab R² values for annotation
r2_season <- round(perm_multi$R2[1], 3)
r2_bay    <- round(perm_multi$R2[2], 3)
p_season  <- perm_multi$`Pr(>F)`[1]
p_bay     <- perm_multi$`Pr(>F)`[2]

annotation_text <- paste0(
  "Season: R²=", r2_season, ", p=", p_season, "\n",
  "Bay: R²=", r2_bay, ", p=", p_bay
)

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                               colour = Season, shape = Bay)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Season), type = "t",
               linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = c("Spring" = "#22c55e",
                                  "Summer" = "#ef4444")) +
  scale_shape_manual(values = c("Chesapeake" = 16,
                                 "Delaware"   = 17)) +
  annotate("text", x = Inf, y = Inf,
           label = annotation_text,
           hjust = 1.05, vjust = 1.5,
           size = 3.5, colour = "grey30") +
  labs(
    x = paste0("PCoA1 (", eig_pct[1], "%)"),
    y = paste0("PCoA2 (", eig_pct[2], "%)"),
    title = "Community composition — Bray-Curtis PCoA",
    colour = "Season", shape = "Bay"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

print(p_pcoa)

## ── 10. Figure — Panel B: betadisper boxplot ─────────────────
disp_df <- data.frame(
  Distance = disp_season$distances,
  Season   = metadata$Season
)

p_disp <- ggplot(disp_df, aes(x = Season, y = Distance, fill = Season)) +
  geom_boxplot(width = 0.45, outlier.shape = 21,
               outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("Spring" = "#22c55e",
                                "Summer" = "#ef4444")) +
  labs(
    x     = NULL,
    y     = "Distance to centroid",
    title = "Beta-dispersion by Season",
    subtitle = "Homogeneity of multivariate dispersion"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

print(p_disp)

## ── 11. Save figures ─────────────────────────────────────────
# Combine panels using patchwork if installed
# install.packages("patchwork")
library(patchwork)

combined_fig <- p_pcoa + p_disp +
  plot_annotation(tag_levels = "A")

print(combined_fig)

ggsave("~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day3_permanova_figure.tiff",
       combined_fig,
       width = 12, height = 5.5,
       dpi = 300,
       compression = "lzw")

## ── End of Day 3 ─────────────────────────────────────────────
