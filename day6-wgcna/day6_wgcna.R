## ============================================================
## Day 6 — WGCNA: MAG Co-Abundance Networks
##          Soft thresholding · TOM · Module detection ·
##          Eigengenes · Module–trait correlations
## Applied Statistics for Microbiome & Genomics Data
## Author: Jojy John
## Series: github.com/jojyjohn28/microbiome_stats
## ============================================================

library(vegan)
library(WGCNA)
library(ggplot2)
library(ggdendro)
library(patchwork)
library(reshape2)
library(dplyr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

## ── 1. Load data ──────────────────────────────────────────────
metadata <- read.csv(
  "~day1/data/metadata.csv",
  header = TRUE, row.names = 1, sep = ","
)
mag_raw <- read.csv(
  "~day1/data/mag_abundance.csv",
  header = TRUE, row.names = 1, sep = ","
)

mag_t    <- t(mag_raw)                                       # samples × MAGs
metadata <- metadata[match(rownames(mag_t), rownames(metadata)), ]
stopifnot(all(rownames(mag_t) == rownames(metadata)))

cat("Samples:", nrow(mag_t), "| MAGs:", ncol(mag_t), "\n")

## ── 2. Hellinger transformation ───────────────────────────────
# Square root of relative abundance.
# Reduces influence of dominant taxa, handles compositionality and
# zero-inflation without arbitrary pseudocounts.
hell <- decostand(mag_t, method = "hellinger")

## ── 3. Check data quality ─────────────────────────────────────
# Flags samples or genes with too many missing values or zero variance.
# gsg$allOK = TRUE means all passed — safe to proceed.
gsg <- goodSamplesGenes(hell, verbose = 3)
if (!gsg$allOK) {
  hell <- hell[gsg$goodSamples, gsg$goodGenes]
  cat("Removed", sum(!gsg$goodSamples), "samples,",
      sum(!gsg$goodGenes), "MAGs\n")
}

## ── 4. Soft threshold selection ───────────────────────────────
# Goal: find the lowest power β where the network fits a scale-free
# topology (R² ≥ 0.80) while mean connectivity remains > 0.
# For n < 30 samples, R² may never reach 0.80 cleanly — accept the
# best local peak in the 6–12 range. Never force β = 30 to hit 0.80;
# you will collapse the network.
powers <- c(1:20)
sft    <- pickSoftThreshold(
  hell,
  powerVector = powers,
  networkType = "signed",
  verbose     = 5
)

# Inspect the output table
print(sft$fitIndices[, c("Power", "SFT.R.sq", "slope", "mean.k.")])

# ── Soft threshold diagnostic plots ──────────────────────────
sft_df <- sft$fitIndices

p_sft_r2 <- ggplot(sft_df, aes(x = Power,
                               y = -sign(slope) * SFT.R.sq)) +
  geom_point(size = 3, colour = "#2563eb") +
  geom_line(colour = "#9ca3af") +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             colour = "#dc2626", linewidth = 0.8) +
  geom_text(aes(label = Power), vjust = -0.8, size = 3,
            colour = "#374151") +
  labs(x = "Soft threshold (β)",
       y = "Scale-free topology fit R²",
       title = "Scale independence") +
  theme_classic(base_size = 11)

p_sft_k <- ggplot(sft_df, aes(x = Power, y = mean.k.)) +
  geom_point(size = 3, colour = "#2563eb") +
  geom_line(colour = "#9ca3af") +
  geom_text(aes(label = Power), vjust = -0.8, size = 3,
            colour = "#374151") +
  labs(x = "Soft threshold (β)",
       y = "Mean connectivity",
       title = "Mean connectivity") +
  theme_classic(base_size = 11)

p_sft_r2 | p_sft_k
# Choose β at the first R² ≥ 0.80 (or the local inflection point)
# AND where mean connectivity is still visibly above zero.

softPower <- 6    # ← update based on your plot

## ── 5. Build weighted adjacency ───────────────────────────────
# Signed network: adjacency = (0.5 + 0.5 × Pearson r)^β
# Signed preserves direction — positively and negatively co-abundant
# MAGs are treated differently. Use type = "unsigned" only if you
# have a specific reason to ignore direction.
adjacency <- adjacency(hell, power = softPower, type = "signed")

## ── 6. Topological Overlap Matrix (TOM) ───────────────────────
# TOM_ij = (Σ_u a_iu·a_uj) / (min(k_i, k_j) + 1 − a_ij)
# Two MAGs score highly not just for direct co-occurrence but for
# sharing many co-occurrence partners — much more robust to noise
# than pairwise correlation in small n datasets.
# This is the slow step. Save and reload rather than recompute.
TOM     <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

## ── 7. Hierarchical clustering on TOM dissimilarity ──────────
TaxaTree <- hclust(as.dist(dissTOM), method = "average")

## ── 8. Dynamic tree cut ───────────────────────────────────────
# deepSplit = 2: moderately sensitive — balances splitting fine
#   structure vs. merging noisy small clusters.
# minClusterSize = 5: do not create modules smaller than 5 MAGs.
#   With n = 27 samples, modules of 2–3 MAGs are not robust.
# pamRespectsDendro = FALSE: allows PAM stage to assign peripheral
#   MAGs to the nearest module regardless of dendrogram position.
minModuleSize <- 5

dynamicMods   <- cutreeDynamic(
  dendro            = TaxaTree,
  distM             = dissTOM,
  deepSplit         = 2,
  pamRespectsDendro = FALSE,
  minClusterSize    = minModuleSize
)
dynamicColors <- labels2colors(dynamicMods)
cat("\nInitial module sizes:\n")
print(table(dynamicColors))

## ── 9. Merge similar modules ──────────────────────────────────
# cutHeight = 0.25 → merge any two modules whose eigengenes
# correlate > 0.75 (1 − 0.25 = 0.75). This collapses noise-split
# modules into the correct smaller number of biological guilds.
# Never skip this step.
merge        <- mergeCloseModules(hell, dynamicColors,
                                  cutHeight = 0.25, verbose = 3)
moduleColors <- merge$colors
mergedMEs    <- merge$newMEs

cat("\nFinal module sizes after merging:\n")
print(table(moduleColors))
# grey = unassigned; these MAGs do not co-vary clearly with any module.
# Investigate separately — do not simply discard them.

## ── 10. Save network objects ──────────────────────────────────
dir.create(
  "day6-wgcna",
  recursive = TRUE, showWarnings = FALSE
)
save(
  TaxaTree, dissTOM, moduleColors, mergedMEs,
  file = "day6_networkConstruction.RData"
)

## ── 11. Recalculate module eigengenes ─────────────────────────
# Eigengene = 1st principal component of module abundance matrix.
# One number per module per sample — the representative abundance
# profile of the entire guild. orderMEs sorts by inter-eigengene
# correlation so the heatmap layout is visually clean.
MEs0 <- moduleEigengenes(hell, moduleColors)$eigengenes
MEs  <- orderMEs(MEs0)

## ── 12. Select and align numerical traits ─────────────────────
trait_cols <- c("Salinity", "Temperature", "ChlA",
                "Bacterial_Production", "Nitrate", "Phosphate")
datTraits  <- metadata[, intersect(trait_cols, names(metadata)),
                       drop = FALSE]
datTraits  <- data.frame(lapply(datTraits, as.numeric))
rownames(datTraits) <- rownames(metadata)

## ── 13. Module–trait Pearson correlation + p-values ───────────
nSamples       <- nrow(hell)
modTraitCor    <- cor(MEs, datTraits, use = "p")
modTraitPvalue <- corPvalueStudent(modTraitCor, nSamples)

cat("\nModule–trait correlations:\n")
print(round(modTraitCor, 3))
cat("\nModule–trait p-values:\n")
print(round(modTraitPvalue, 4))

## ── 14. Module membership (MM) and trait significance (GS) ────
# MM: Pearson r between each MAG and each module eigengene.
#   High MM (> 0.8) = hub MAG, central to the module.
# GS: Pearson r between each MAG and a trait of interest.
#   High GS (> 0.5) = the MAG is strongly linked to that trait.
# Hub MAGs = high MM AND high GS — best candidates for follow-up.

sal      <- as.data.frame(datTraits$Salinity)
names(sal) <- "Salinity"
modNames <- substring(names(MEs), 3)

ModuleMembership <- as.data.frame(cor(hell, MEs, use = "p"))
MMPvalue         <- as.data.frame(
  corPvalueStudent(as.matrix(ModuleMembership),
                   nSamples))
names(ModuleMembership) <- paste0("MM", modNames)
names(MMPvalue)         <- paste0("p.MM", modNames)

TraitSignificance <- as.data.frame(cor(hell, sal, use = "p"))
TSPvalue          <- as.data.frame(
  corPvalueStudent(as.matrix(TraitSignificance),
                   nSamples))
names(TraitSignificance) <- "GS.Salinity"
names(TSPvalue)          <- "p.GS.Salinity"

## ── 15. Export per-MAG results table ──────────────────────────
# hell is a matrix → use colnames(), not names(), for MAG names.
# TraitSignificance and TSPvalue have MAG names as rownames —
# colnames(hell) aligns them correctly.
MAGInfo <- data.frame(
  MAG         = colnames(hell),
  moduleColor = moduleColors,
  TraitSignificance,
  TSPvalue,
  row.names   = NULL
)
for (mod in 1:ncol(ModuleMembership)) {
  old     <- names(MAGInfo)
  MAGInfo <- data.frame(MAGInfo,
                        ModuleMembership[, mod],
                        MMPvalue[, mod])
  names(MAGInfo) <- c(old,
                      paste0("MM.", modNames[mod]),
                      paste0("p.MM.", modNames[mod]))
}
MAGOrder <- order(MAGInfo$moduleColor,
                  -abs(MAGInfo$GS.Salinity))
MAGInfo  <- MAGInfo[MAGOrder, ]

dir.create(
  "/day6-wgcna/results",
  recursive = TRUE, showWarnings = FALSE
)
write.csv(
  MAGInfo,
  "day6-wgcna/results/day6_MAGInfo.csv",
  row.names = FALSE
)
cat("\nExported MAGInfo.csv — columns: moduleColor, GS.Salinity,",
    "MM per module\n")

## ── 16. Panel A: MAG dendrogram + module color bar ────────────
# Leaf order from TaxaTree so color bar aligns with dendrogram leaves.
leaf_order <- TaxaTree$order

cbar_df <- data.frame(
  x     = seq_along(leaf_order),
  color = moduleColors[leaf_order]
)

p_dend <- ggdendrogram(TaxaTree, rotate = FALSE, labels = FALSE) +
  labs(
    title = "A   MAG Dendrogram (TOM dissimilarity)",
    x     = NULL,
    y     = "Dissimilarity (1 \u2212 TOM)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    plot.title    = element_text(face = "bold", size = 11),
    plot.margin   = margin(5, 5, 0, 5)
  )
p_dend
p_cbar <- ggplot(cbar_df, aes(x = x, y = 1, fill = color)) +
  geom_tile(height = 1) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 5, 5))

p_cbar

# Stack dendrogram above color bar; color bar takes ~10% of height
p_A <- p_dend / p_cbar +
  plot_layout(heights = c(10, 1))

p_A

## ── 17. Panel B: Topological Overlap Matrix (TOM) ─────────────
# Reorder by dendrogram leaf order so module block structure
# appears along the diagonal. Raise to the 7th power to increase
# contrast between low (off-module) and high (within-module) TOM.
tom_vals  <- 1 - dissTOM
tom_plot  <- tom_vals[leaf_order, leaf_order]
tom_plot7 <- tom_plot^7
diag(tom_plot7) <- NA

tom_df       <- melt(tom_plot7,
                     varnames   = c("MAG_x", "MAG_y"),
                     value.name = "TOM7")
tom_df$MAG_x <- as.numeric(tom_df$MAG_x)
tom_df$MAG_y <- as.numeric(tom_df$MAG_y)

p_B <- ggplot(tom_df, aes(x = MAG_x, y = MAG_y, fill = TOM7)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradientn(
    colours  = c("#f8fafc", "#fde68a", "#f59e0b",
                 "#dc2626", "#7f1d1d"),
    na.value = "white",
    name     = "TOM\u2077",
    guide    = guide_colorbar(barwidth = 0.8, barheight = 6)
  ) +
  labs(
    title = "B   Topological Overlap Matrix (TOM)",
    x     = "MAGs",
    y     = "MAGs"
  ) +
  coord_fixed() +
  theme_classic(base_size = 11) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    legend.position  = "right"
  )

p_B

## ── 18. Panel C: Module–trait heatmap ─────────────────────────
# Each cell = Pearson r (eigengene vs trait) with significance stars.
# Diverging red–white–blue: red = positive, blue = negative.
cor_df  <- melt(modTraitCor,
                varnames   = c("Module", "Trait"),
                value.name = "r")
pval_df <- melt(modTraitPvalue,
                varnames   = c("Module", "Trait"),
                value.name = "p")
heatmap_df <- merge(cor_df, pval_df)

heatmap_df$star <- with(heatmap_df,
                        ifelse(p < 0.001, "***",
                               ifelse(p < 0.01,  "**",
                                      ifelse(p < 0.05,  "*", ""))))

heatmap_df$label <- with(heatmap_df,
                         paste0(sprintf("%.2f", r),
                                ifelse(star != "", paste0("\n", star),
                                       paste0("\np=", sprintf("%.2f", p)))))

# Remove "ME" prefix from module names for a cleaner y-axis
heatmap_df$Module <- gsub("^ME", "", as.character(heatmap_df$Module))

# Human-readable trait labels
trait_display <- c(
  Salinity             = "Salinity\n(psu)",
  Temperature          = "Temperature\n(\u00b0C)",
  ChlA                 = "Chl-a\n(\u03bcg/L)",
  Bacterial_Production = "Bacterial\nProduction",
  Nitrate              = "Nitrate\n(\u03bcM)",
  Phosphate            = "Phosphate\n(\u03bcM)"
)
heatmap_df$Trait_label <- dplyr::recode(
  as.character(heatmap_df$Trait),
  !!!trait_display
)

# Order modules by salinity correlation for consistent layout
module_order <- heatmap_df %>%
  filter(Trait == "Salinity") %>%
  arrange(r) %>%
  pull(Module)
heatmap_df$Module <- factor(heatmap_df$Module, levels = module_order)

p_C <- ggplot(heatmap_df,
              aes(x = Trait_label, y = Module, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label),
            size       = 3.2,
            lineheight = 0.9,
            colour     = ifelse(abs(heatmap_df$r) > 0.65,
                                "white", "#111827")) +
  scale_fill_gradientn(
    colours = c("#2166ac", "#92c5de", "#f7f7f7",
                "#f4a582", "#b2182b"),
    limits  = c(-1, 1),
    name    = "Pearson r",
    guide   = guide_colorbar(barwidth = 0.8, barheight = 6)
  ) +
  labs(
    title    = "C   Module\u2013Trait Relationships",
    subtitle = "Pearson r  \u00b7  * p<0.05  ** p<0.01  *** p<0.001",
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 9, lineheight = 0.9),
    axis.text.y   = element_text(size = 10, face = "bold"),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey50"),
    legend.position = "right"
  )

p_C
## ── 19. Assemble and save figure ──────────────────────────────
combined_fig <- (p_A | p_B | p_C) +
  plot_annotation(
    title    = "Day 6 \u2014 WGCNA: MAG Co-Abundance Modules & Environmental Drivers",
    subtitle = paste0("61 MAGs  \u00b7  27 metagenomes  \u00b7  Signed Pearson",
                      "  \u00b7  \u03b2=6  \u00b7  TOM clustering",
                      "  \u00b7  Hellinger transformation"),
    tag_levels = "A",
    theme      = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey50")
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 13))

combined_fig


dir.create(
  "/day6-wgcna/figures",
  recursive = TRUE, showWarnings = FALSE
)
ggsave(
  "figures/day6_wgcna_figure.tiff",
  combined_fig,
  width       = 18,
  height      = 7,
  dpi         = 300,
  compression = "lzw"
)

cat("\nDone.\n")
cat("Network:  day6-wgcna/day6_networkConstruction.RData\n")
cat("Results:  day6-wgcna/results/day6_MAGInfo.csv\n")
cat("Figure:   day6-wgcna/figures/day6_wgcna_figure.tiff\n")
