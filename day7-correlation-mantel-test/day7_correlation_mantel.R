## ============================================================
## Day 7 — Correlation Analysis & Mantel Test
## Linking MAG abundance and environmental variables
## Applied Statistics for Microbiome & Genomics Data
## Author: Jojy John
## github.com/jojyjohn28/microbiome_stats
## ============================================================

library(vegan)
library(Hmisc)
library(corrplot)
library(ggcorrplot)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)

## ── 1. Load data ──────────────────────────────────────────────
metadata  <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/metadata.csv",
  header = TRUE, row.names = 1, sep = ","
)
mag_raw <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/mag_abundance.csv",
  header = TRUE, row.names = 1, sep = ","
)
mag_t <- t(mag_raw)   # samples × MAGs

# Align
metadata <- metadata[match(rownames(mag_t), rownames(metadata)), ]
stopifnot(all(rownames(mag_t) == rownames(metadata)))
cat("Samples:", nrow(mag_t), "| MAGs:", ncol(mag_t), "\n")

## ── 2. Select numeric environmental variables ─────────────────
env_cols <- c("Salinity", "Temperature", "Bacterial_Production",
              "nCells", "ChlA", "Nitrate",
              "Ammonium", "Phosphate", "Silicate")
env_df <- metadata[, env_cols]

# Replace negative values (instrument artefacts) with NA
env_df[env_df < 0] <- NA

# Remove columns with >30% missing
keep <- colMeans(is.na(env_df)) < 0.3
env_df <- env_df[, keep]
cat("Environmental variables retained:", ncol(env_df), "\n")
cat(names(env_df), "\n")

## ── 3. PART A: Environmental variable correlation ─────────────
# Spearman — robust to non-normality and outliers (better for
# field environmental data than Pearson)

corr_env     <- rcorr(as.matrix(env_df), type = "spearman")
corr_mat_env <- corr_env$r
pval_mat_env <- corr_env$P
diag(pval_mat_env) <- 1   # set diagonal p to 1 (self-correlation)

cat("\nEnvironment–environment Spearman correlations:\n")
print(round(corr_mat_env, 2))

## ── 4. Panel A: Environmental correlation heatmap ─────────────
dir.create(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day7-correlation/figures",
  recursive = TRUE, showWarnings = FALSE
)

# clean variable labels
colnames(corr_mat_env) <- rownames(corr_mat_env) <- c(
  "Salinity", "Temperature", "Bact. Prod.",
  "Cell Count", "Chl-a", "Nitrate",
  "Ammonium", "Phosphate", "Silicate"
)[seq_len(ncol(corr_mat_env))]

colnames(pval_mat_env) <- rownames(pval_mat_env) <- colnames(corr_mat_env)

p_A <- ggcorrplot(
  corr_mat_env,
  hc.order    = TRUE,
  type        = "lower",
  lab         = TRUE,
  lab_size    = 3,
  p.mat       = pval_mat_env,
  sig.level   = 0.05,
  insig       = "blank",
  colors      = c("#2166ac", "#f7f7f7", "#b2182b"),
  outline.col = "white",
  tl.cex      = 9,
  tl.col      = "grey20"
) +
  labs(
    title    = "A   Environmental Variable Correlations",
    subtitle = "Spearman r · blank = p > 0.05"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey50"),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 7)
  )
p_A
## ── 5. PART B: MAG abundance vs environment correlation ───────
# Log-transform abundances
mag_log <- log1p(mag_t)

# Spearman r: each MAG vs each env variable
# rcorr requires a single matrix — bind them together
combined <- cbind(mag_log, env_df)
corr_all     <- rcorr(as.matrix(combined), type = "spearman")
corr_all_r   <- corr_all$r
corr_all_p   <- corr_all$P

n_mags <- ncol(mag_log)
n_env  <- ncol(env_df)

# Extract the MAG × env block only
mag_env_r <- corr_all_r[1:n_mags, (n_mags + 1):(n_mags + n_env)]
mag_env_p <- corr_all_p[1:n_mags, (n_mags + 1):(n_mags + n_env)]

# BH correction across all MAG–env pairs
mag_env_q <- matrix(
  p.adjust(as.vector(mag_env_p), method = "BH"),
  nrow = n_mags, ncol = n_env,
  dimnames = dimnames(mag_env_p)
)

# Summary: how many significant pairs per env variable
sig_per_env <- colSums(mag_env_q < 0.05, na.rm = TRUE)
cat("\nSignificant MAG–env pairs (q < 0.05) per variable:\n")
print(sig_per_env)

# Top 20 MAGs by max |r| across env variables (for plotting)
max_r <- apply(abs(mag_env_r), 1, max, na.rm = TRUE)
top20 <- names(sort(max_r, decreasing = TRUE))[1:20]

mag_env_r_top <- mag_env_r[top20, ]
mag_env_q_top <- mag_env_q[top20, ]

# Clean MAG labels (remove long prefixes)
rownames(mag_env_r_top) <- gsub("_MAG_", " MAG ", rownames(mag_env_r_top))

# Melt for ggplot
r_df <- melt(mag_env_r_top, varnames = c("MAG", "Env"), value.name = "r")
q_df <- melt(mag_env_q_top, varnames = c("MAG", "Env"), value.name = "q")
r_df$MAG <- gsub("_MAG_", " MAG ", as.character(r_df$MAG))
q_df$MAG <- gsub("_MAG_", " MAG ", as.character(q_df$MAG))
hm_df <- merge(r_df, q_df)

hm_df$star <- with(hm_df,
  ifelse(q < 0.001, "***",
  ifelse(q < 0.01,  "**",
  ifelse(q < 0.05,  "*", ""))))

p_B <- ggplot(hm_df, aes(x = Env, y = MAG, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = star), size = 2.8, colour = "white",
            vjust = 0.75) +
  scale_fill_gradientn(
    colours = c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#b2182b"),
    limits  = c(-1, 1),
    name    = "Spearman r"
  ) +
  labs(
    title    = "B   Top 20 MAGs \u2013 Environment Correlations",
    subtitle = "Spearman r · BH-corrected  * q<0.05  ** q<0.01  *** q<0.001",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y   = element_text(size   = 8),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey50"),
    legend.position = "right"
  )
p_B
## ── 6. PART C: Mantel test ────────────────────────────────────
# Question: does community dissimilarity correlate with
# environmental dissimilarity across samples?
#
# Community distance: Bray-Curtis on Hellinger-transformed abundances
# Environmental distance: Euclidean on scaled numeric env variables

# Make sure sample order matches community matrix
stopifnot(all(rownames(env_df) == rownames(mag_t)))

# Keep only numeric environmental variables
env_num <- env_df[, sapply(env_df, is.numeric), drop = FALSE]

# Remove columns with any NA or zero variance for global Mantel
env_num_clean <- env_num[, sapply(env_num, function(x) all(!is.na(x)) && sd(x) > 0), drop = FALSE]

# Community distance
hell    <- decostand(mag_t, method = "hellinger")
bc_dist <- vegdist(hell, method = "bray")

# Environmental distance
env_scaled <- scale(env_num_clean)
env_dist   <- dist(env_scaled, method = "euclidean")

# ── Full Mantel ───────────────────────────────────────────────
set.seed(42)
mantel_full <- mantel(bc_dist, env_dist, method = "pearson", permutations = 9999)

cat("\nMantel test (community ~ full environment):\n")
print(mantel_full)

# ── Mantel per individual env variable ────────────────────────
bc_mat <- as.matrix(bc_dist)

mantel_results <- do.call(rbind, lapply(names(env_num), function(var) {
  
  x <- env_num[[var]]
  
  if (any(is.na(x))) {
    return(data.frame(Variable = var, r = NA, p = NA, Note = "contains NA"))
  }
  
  if (sd(x) == 0) {
    return(data.frame(Variable = var, r = NA, p = NA, Note = "zero variance"))
  }
  
  d <- dist(scale(x), method = "euclidean")
  m <- mantel(bc_dist, d, method = "pearson", permutations = 9999)
  
  data.frame(
    Variable = var,
    r = unname(m$statistic),
    p = m$signif,
    Note = "ok"
  )
}))

# add significance labels for plotting
mantel_results$sig <- ifelse(
  is.na(mantel_results$p), "",
  ifelse(mantel_results$p < 0.001, "***",
         ifelse(mantel_results$p < 0.01, "**",
                ifelse(mantel_results$p < 0.05, "*", "ns")))
)

print(mantel_results)

# ── Partial Mantel: Salinity after controlling for Temperature ─
if (!"Salinity" %in% names(env_num)) stop("Salinity column not found")
if (!"Temperature" %in% names(env_num)) stop("Temperature column not found")

if (any(is.na(env_num$Salinity)) || any(is.na(env_num$Temperature))) {
  stop("Salinity or Temperature contains NA")
}
if (sd(env_num$Salinity) == 0 || sd(env_num$Temperature) == 0) {
  stop("Salinity or Temperature has zero variance")
}

sal_dist  <- dist(scale(env_num[["Salinity"]]), method = "euclidean")
temp_dist <- dist(scale(env_num[["Temperature"]]), method = "euclidean")

set.seed(42)
partial_mantel <- mantel.partial(
  bc_dist, sal_dist, temp_dist,
  method = "pearson",
  permutations = 9999
)

cat("\nPartial Mantel (community ~ Salinity | Temperature):\n")
print(partial_mantel)

## ── 7. Panel C: Mantel r bar chart ────────────────────────────
plot_df <- subset(mantel_results, !is.na(r))

plot_df$Variable <- factor(
  plot_df$Variable,
  levels = plot_df$Variable[order(plot_df$r)]
)

plot_df$fill_col <- ifelse(plot_df$p < 0.05, "#2166ac", "#aac4dc")

xmin <- min(0, min(plot_df$r, na.rm = TRUE) * 1.2)
xmax <- max(0, max(plot_df$r, na.rm = TRUE) * 1.4)

p_C <- ggplot(plot_df, aes(x = r, y = Variable, fill = fill_col)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(sig, "  r=", round(r, 2))),
    hjust = ifelse(plot_df$r >= 0, -0.1, 1.1),
    size = 3,
    colour = "grey30"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60") +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(xmin, xmax),
    expand = c(0.02, 0.02)
  ) +
  labs(
    title = "C   Mantel Test: Community vs Environment",
    subtitle = "Pearson r · 9999 permutations · blue = p < 0.05",
    x = "Mantel r",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey50"),
    axis.text.y = element_text(size = 9)
  )

p_C

## ── 8. Mantel permutation distribution (Panel D) ──────────────
# Visualise the null distribution from the full Mantel test
perm_r <- mantel_full$perm
obs_r  <- mantel_full$statistic

perm_df <- data.frame(r = perm_r)

p_D <- ggplot(perm_df, aes(x = r)) +
  geom_histogram(bins = 50, fill = "#aac4dc",
                 colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_r, colour = "#b2182b",
             linewidth = 1.2, linetype = "dashed") +
  annotate("text",
           x = obs_r + 0.01,
           y = Inf, vjust = 1.6, hjust = 0,
           label = paste0("Observed r = ", round(obs_r, 3),
                          "\np = ", mantel_full$signif),
           size = 3.2, colour = "#b2182b") +
  labs(
    title    = "D   Mantel Permutation Distribution",
    subtitle = "Full environment · 9999 permutations",
    x        = "Permuted Mantel r",
    y        = "Count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5, colour = "grey50")
  )
p_D
## ── 9. Assemble and save ──────────────────────────────────────
combined_fig <- (p_A | p_B) / (p_C | p_D) +
  plot_annotation(
    title    = "Day 7 \u2014 Correlation Analysis & Mantel Test",
    subtitle = paste0(
      "61 MAGs  \u00b7  27 metagenomes  \u00b7  Spearman r  \u00b7  ",
      "BH FDR  \u00b7  Chesapeake + Delaware Bay"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey50")
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day7-correlation/figures/day7_correlation_mantel_figure.tiff",
  combined_fig,
  width = 16, height = 12, dpi = 300, compression = "lzw"
)

## ── 10. Save results tables ───────────────────────────────────
dir.create(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day7-correlation/results",
  recursive = TRUE, showWarnings = FALSE
)

# Full MAG–env correlation table
mag_env_results <- data.frame(
  MAG       = rep(rownames(mag_env_r), ncol(mag_env_r)),
  Env_var   = rep(colnames(mag_env_r), each = nrow(mag_env_r)),
  r         = as.vector(mag_env_r),
  p         = as.vector(mag_env_p),
  q_BH      = as.vector(mag_env_q)
)
mag_env_results <- mag_env_results[order(mag_env_results$q_BH), ]

write.csv(
  mag_env_results,
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day7-correlation/results/day7_MAG_env_correlations.csv",
  row.names = FALSE
)
write.csv(
  mantel_results,
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day7-correlation/results/day7_mantel_results.csv",
  row.names = FALSE
)

cat("\nDone.\n")
cat("Figure:  day7-correlation/figures/day7_correlation_mantel_figure.tiff\n")
cat("Results: day7-correlation/results/day7_MAG_env_correlations.csv\n")
cat("         day7-correlation/results/day7_mantel_results.csv\n")

###########################################################
#Mantel relationships between a biological matrix and environmental variables/groups
#laoding required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(linkET)

#Preparing data

# keep only numeric environmental variables
env_num <- env_df[, sapply(env_df, is.numeric), drop = FALSE]

# make sure row order matches community matrix
env_num <- env_num[rownames(mag_t), , drop = FALSE]

# optional: remove NA-containing columns
env_num <- env_num[, colSums(is.na(env_num)) == 0, drop = FALSE]

# optional: remove zero-variance columns
env_num <- env_num[, sapply(env_num, function(x) sd(x) > 0), drop = FALSE]

#Correlation matrix
cor_env <- correlate(env_num, method = "spearman")

#mantel test
hell <- decostand(mag_t, method = "hellinger")

mantel_df <- do.call(rbind, lapply(names(env_num), function(v) {
  d_env <- dist(scale(env_num[[v]]))
  m <- mantel(vegdist(hell, method = "bray"), d_env, permutations = 9999)
  
  data.frame(
    from = "Community",
    to = v,
    r = unname(m$statistic),
    p = m$signif
  )
}))

mantel_df <- mantel_df %>%
  mutate(
    p_cat = case_when(
      p < 0.05 ~ "< 0.05",
      p < 0.10 ~ "0.05 - 0.1",
      TRUE ~ ">= 0.1"
    ),
    r_cat = case_when(
      r < 0.2 ~ "< 0.2",
      r < 0.4 ~ "0.2 - 0.4",
      TRUE ~ ">= 0.4"
    )
  )
#plot
p <- qcorrplot(cor_env,
               type = "lower",
               diag = FALSE) +
  geom_square(aes(fill = r), size = 5) +
  scale_fill_gradient2(
    low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0
  ) +
  geom_couple(
    aes(colour = p_cat, size = r_cat),
    data = mantel_df,
    curvature = nice_curvature()
  ) +
  scale_size_manual(
    values = c("< 0.2" = 0.6, "0.2 - 0.4" = 1.0, ">= 0.4" = 1.4)
  ) +
  scale_colour_manual(
    values = c("< 0.05" = "limegreen", "0.05 - 0.1" = "orange", ">= 0.1" = "grey70")
  ) +
  guides(
    size = guide_legend(title = "Mantel's r"),
    colour = guide_legend(title = "Mantel's p"),
    fill = guide_colorbar(title = "Spearman's r")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )

p
################################################333
End