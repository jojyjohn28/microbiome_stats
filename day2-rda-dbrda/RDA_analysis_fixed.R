# =========================================================
# RDA of MAG abundance vs environmental variables
# Fixed version: handles NA values and builds final model safely
# =========================================================

library(vegan)
library(ggplot2)
library(ggrepel)

# =========================
# Step 1 — Load data
# =========================

# Load metadata — rows = samples, columns = environmental variables
metadata <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/metadata.csv",
  header = TRUE,
  row.names = 1,
  sep = ",",
  check.names = FALSE
)

# Load MAG abundance — rows = MAGs, columns = samples
mag_raw <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/mag_abundance.csv",
  header = TRUE,
  row.names = 1,
  sep = ",",
  check.names = FALSE
)

# =========================
# Step 2 — Match sample names
# =========================

cat("Checking sample name matching...\n")
cat("Samples in abundance table but not metadata:\n")
print(setdiff(colnames(mag_raw), rownames(metadata)))

cat("Samples in metadata but not abundance table:\n")
print(setdiff(rownames(metadata), colnames(mag_raw)))

# Reorder metadata to match abundance table
metadata <- metadata[colnames(mag_raw), , drop = FALSE]

cat("Do sample names now match exactly?\n")
print(identical(colnames(mag_raw), rownames(metadata)))

# =========================
# Step 3 — Prepare community matrix
# =========================

# Transpose abundance table: samples as rows, MAGs as columns
mag_t <- t(mag_raw)

# Hellinger transform
spe.hel <- decostand(mag_t, method = "hellinger")

# =========================
# Step 4 — Prepare environmental matrix
# =========================

# Keep only continuous numeric variables for initial RDA
env <- metadata[, c("Salinity", "Temperature", "Depth", "PAR",
                    "Attenuation", "Bacterial_Production", "nCells",
                    "ChlA", "Nitrate", "Ammonium", "Phosphate", "Silicate"),
                drop = FALSE]

# Force numeric in case any columns were imported as character/factor
env <- as.data.frame(env)
env[] <- lapply(env, function(x) as.numeric(as.character(x)))

cat("\nMissing values per environmental variable:\n")
print(colSums(is.na(env)))

cat("\nSamples with incomplete environmental data:\n")
print(rownames(env)[!complete.cases(env)])

# Keep only complete samples
keep <- complete.cases(env)

env2      <- env[keep, , drop = FALSE]
metadata2 <- metadata[keep, , drop = FALSE]
spe.hel2  <- spe.hel[keep, , drop = FALSE]

# Standardize environmental variables
env.z <- as.data.frame(decostand(env2, method = "stand"))

# Optional: add categorical variables separately if needed later
metadata2$Season <- factor(metadata2$Season)
metadata2$Bay    <- factor(metadata2$Bay)
metadata2$SF     <- factor(metadata2$SF)

# =========================
# Step 5 — Full RDA model
# =========================

spe.rda.full <- rda(spe.hel2 ~ ., data = env.z)

cat("\nFull model summary:\n")
print(summary(spe.rda.full))

# =========================
# Step 6 — Forward selection
# =========================

set.seed(123)

fwd.sel <- ordiR2step(
  rda(spe.hel2 ~ 1, data = env.z),   # null model
  scope = formula(spe.rda.full),     # full model
  direction = "forward",
  R2scope = TRUE,
  permutations = 1000,
  trace = FALSE
)

cat("\nForward-selected model:\n")
print(fwd.sel)
cat("\nSelected formula:\n")
print(formula(fwd.sel))

# =========================
# Step 7 — Final RDA model
# =========================

# Use the selected model directly
spe.rda.signif <- fwd.sel

cat("\nAdjusted R2:\n")
print(RsquareAdj(spe.rda.signif))

cat("\nOverall model significance:\n")
print(anova.cca(spe.rda.signif, permutations = 1000))

cat("\nSignificance of each term:\n")
print(anova.cca(spe.rda.signif, permutations = 1000, by = "term"))

cat("\nSignificance of each axis:\n")
print(anova.cca(spe.rda.signif, permutations = 1000, by = "axis"))

# =========================
# Step 8 — Extract scores
# =========================

# % variance explained by first 2 constrained axes
perc <- round(100 * summary(spe.rda.signif)$cont$importance[2, 1:2], 2)

# Site scores
sc_sites <- scores(spe.rda.signif, display = "sites", choices = c(1, 2), scaling = 1)
sc_sites <- as.data.frame(sc_sites)
sc_sites$Sample  <- rownames(sc_sites)
sc_sites$Season  <- metadata2$Season
sc_sites$Bay     <- metadata2$Bay
sc_sites$SF      <- metadata2$SF

# Species (MAG) scores
sc_sp <- scores(spe.rda.signif, display = "species", choices = c(1, 2), scaling = 1)
sc_sp <- as.data.frame(sc_sp)
sc_sp$MAG <- rownames(sc_sp)

# Environmental arrows
sc_bp <- scores(spe.rda.signif, display = "bp", choices = c(1, 2), scaling = 1)
sc_bp <- as.data.frame(sc_bp)
sc_bp$Variable <- rownames(sc_bp)

# Label only strong MAG contributors
high_score_mags <- sc_sp$MAG[apply(abs(sc_sp[, 1:2]), 1, max) > 0.35]

# =========================
# Step 9 — Base R triplot
# =========================

plot(spe.rda.signif,
     scaling = 1,
     type = "none",
     frame = FALSE,
     xlab = paste0("RDA1 (", perc[1], "%)"),
     ylab = paste0("RDA2 (", perc[2], "%)"),
     main = "RDA of MAG community vs environment")

# Samples colored by Season
points(sc_sites[sc_sites$Season == "Summer", 1:2],
       pch = 21, col = "black", bg = "red2", cex = 1.2)

points(sc_sites[sc_sites$Season == "Spring", 1:2],
       pch = 21, col = "black", bg = "green3", cex = 1.2)

# MAGs
points(sc_sp[, 1:2], pch = 22, col = "black", bg = "#f2bd33", cex = 1.1)

# Label selected MAGs
sel_sp <- sc_sp[sc_sp$MAG %in% high_score_mags, ]
text(sel_sp$RDA1 + 0.03, sel_sp$RDA2 + 0.03,
     labels = sel_sp$MAG, col = "grey30", cex = 0.6, font = 3)

# Environmental arrows
arrows(0, 0, sc_bp$RDA1, sc_bp$RDA2, col = "blue", lwd = 1.2, length = 0.08)
text(sc_bp$RDA1 * 1.08, sc_bp$RDA2 * 1.08,
     labels = sc_bp$Variable, col = "blue", cex = 0.9, font = 2)

legend("topright",
       legend = c("Summer", "Spring", "MAGs"),
       pt.bg = c("red2", "green3", "#f2bd33"),
       pch = c(21, 21, 22),
       pt.cex = 1.2,
       bty = "n")

# =========================
# Step 10 — ggplot version
# =========================

p <- ggplot() +
  geom_point(data = sc_sites,
             aes(x = RDA1, y = RDA2, fill = Season, shape = SF),
             size = 3.5, color = "black", stroke = 0.4) +
  geom_point(data = sc_sp,
             aes(x = RDA1, y = RDA2),
             shape = 22, size = 2.4, fill = "#f2bd33", color = "black", stroke = 0.3) +
  geom_segment(data = sc_bp,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.18, "cm")),
               linewidth = 0.7, color = "blue") +
  geom_text_repel(data = subset(sc_sp, MAG %in% high_score_mags),
                  aes(x = RDA1, y = RDA2, label = MAG),
                  size = 3, color = "grey20", fontface = "italic", max.overlaps = 50) +
  geom_text_repel(data = sc_bp,
                  aes(x = RDA1, y = RDA2, label = Variable),
                  size = 3.5, color = "blue", fontface = "bold") +
  labs(
    x = paste0("RDA1 (", perc[1], "%)"),
    y = paste0("RDA2 (", perc[2], "%)"),
    title = "RDA of MAG community and environmental variables"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(p)

# =========================
# Step 11 — Save outputs
# =========================

write.csv(sc_sites, "RDA_site_scores.csv", row.names = FALSE)
write.csv(sc_sp, "RDA_MAG_scores.csv", row.names = FALSE)
write.csv(sc_bp, "RDA_env_arrows.csv", row.names = FALSE)

ggsave("RDA_MAG_environment_triplot.pdf", p, width = 9, height = 7)
ggsave("RDA_MAG_environment_triplot.png", p, width = 9, height = 7, dpi = 300)

# =========================
# Optional Step 12 — If you want to test SF too
# =========================
# SF is categorical, so do not include it inside env.z directly.
# Instead combine numeric env + SF like this:

env.plus <- cbind(env.z, SF = metadata2$SF)
env.plus$SF <- factor(env.plus$SF)

# Example model including SF explicitly
# spe.rda.sf <- rda(spe.hel2 ~ SF + Salinity + Temperature + ChlA, data = env.plus)
# summary(spe.rda.sf)
# anova.cca(spe.rda.sf, permutations = 1000)

