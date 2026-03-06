## ============================================================
## Day 5 — Multiple Regression for MAG Abundance
##          Linking Genomics to Environment
## Applied Statistics for Microbiome & Genomics Data
## Author: Jojy John
## Series: github.com/jojyjohn28/microbiome_stats
## ============================================================
##
## Models covered:
##   1. Simple regression:   MAG ~ Salinity
##   2. Additive model:      MAG ~ Salinity + Temperature + Season
##   3. Interaction model:   MAG ~ Salinity * Temperature
##   4. Loop across 61 MAGs with FDR correction
##
## Diagnostics: residuals vs fitted, Q-Q, scale-location,
##              residual histogram, Shapiro-Wilk, Breusch-Pagan
##
## Outputs:
##   results/day5_regression_all_MAGs.csv
##   figures/day5_regression_figure.tiff
## ============================================================

library(vegan)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(broom)    # tidy(), glance() — clean model summaries
library(car)      # vif() — variance inflation factors
library(lmtest)   # bptest() — Breusch-Pagan heteroscedasticity test
library(dplyr)
library(tidyr)


## ── 1. Load and align data ───────────────────────────────────

metadata <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/metadata.csv",
  header = TRUE, row.names = 1, sep = ","
)
mag_raw <- read.csv(
  "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/mag_abundance.csv",
  header = TRUE, row.names = 1, sep = ","
)

mag_t    <- t(mag_raw)                                       # samples × MAGs
metadata <- metadata[match(rownames(mag_t), rownames(metadata)), ]
stopifnot(all(rownames(mag_t) == rownames(metadata)))

cat("Samples:", nrow(metadata), "\n")
cat("MAGs:   ", ncol(mag_t),    "\n")


## ── 2. Log-transform & build predictors ──────────────────────

# CRITICAL: log1p before any regression — raw counts are right-skewed
mag_log <- log1p(mag_t)

# Encode Season as binary (Spring = 0, Summer = 1)
metadata$Season_num <- ifelse(metadata$Season == "Summer", 1, 0)

# Explicit interaction term (used in Model 3)
metadata$Sal_x_Temp <- metadata$Salinity * metadata$Temperature

# Salinity groups for reference
metadata$Sal_group <- ifelse(metadata$Salinity < 15, "Low (<15 psu)", "High (≥15 psu)")

cat("\nSalinity range:", range(metadata$Salinity), "\n")
cat("Temperature range:", range(metadata$Temperature), "\n")
cat("Season breakdown:\n"); print(table(metadata$Season))


## ── 3. Build combined modelling dataframe ────────────────────

# Using Salinity, Temperature, ChlA — all 27 samples, no NAs
# (Bacterial_Production has 2 NAs; Secchi/PAR/Attenuation have 10 NAs)
focal_mags <- c(
  "Pelagibacterales_MAG_18",   # SAR11 — marine specialist, salinity+
  "Nanopelagicales_MAG_18",    # acI   — freshwater specialist, salinity-
  "Flavobacteriales_MAG_6",    # particle-degrader, season-linked
  "Rhodobacterales_MAG_15"     # Roseobacter, moderate salinity response
)

reg_df <- cbind(
  metadata[, c("Season", "Season_num", "Bay", "SF",
               "Salinity", "Temperature", "ChlA", "nCells",
               "Sal_x_Temp")],
  mag_log[, focal_mags]
)
reg_df <- as.data.frame(reg_df)

# Confirm no NAs (n should remain 27)
reg_df <- reg_df[complete.cases(reg_df[, c("Salinity","Temperature","Season_num",
                                             focal_mags)]), ]
cat("Complete cases for modelling: n =", nrow(reg_df), "\n")


## ── 4. Step 1 — Always scatter-plot before modelling ─────────

p_A <- ggplot(reg_df, aes(Salinity, Pelagibacterales_MAG_18, colour = Season)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, colour = "#1d4ed8", linewidth = 1.4,
              fill = "#bfdbfe", alpha = 0.15) +
  scale_colour_manual(values = c("Spring" = "#16a34a", "Summer" = "#dc2626")) +
  annotate("text", x = 1, y = 11.5,
           label = "r² = 0.654\nβ = +0.327 per psu\np < 0.001",
           hjust = 0, size = 3.5, family = "mono", colour = "#1d4ed8") +
  labs(x = "Salinity (psu)", y = "log(abundance + 1)",
       title = "Pelagibacterales MAG 18 ~ Salinity",
       subtitle = "Simple OLS — n = 27") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

p_B <- ggplot(reg_df, aes(Salinity, Nanopelagicales_MAG_18, colour = Season)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, colour = "#7c3aed", linewidth = 1.4,
              fill = "#ddd6fe", alpha = 0.12) +
  scale_colour_manual(values = c("Spring" = "#16a34a", "Summer" = "#dc2626")) +
  annotate("text", x = 16, y = 9.5,
           label = "r² = 0.769\nβ = −0.325 per psu\np < 0.001",
           hjust = 0, size = 3.5, family = "mono", colour = "#7c3aed") +
  labs(x = "Salinity (psu)", y = "log(abundance + 1)",
       title = "Nanopelagicales MAG 18 ~ Salinity",
       subtitle = "Opposite salinity response") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))


## ── 5. Additive multiple regression models ───────────────────

# ── Model 1: Pelagibacterales ~ Sal + Temp + Season ──────────
m1 <- lm(Pelagibacterales_MAG_18 ~ Salinity + Temperature + Season_num,
          data = reg_df)
cat("\n============================================================\n")
cat("Model 1: Pelagibacterales_MAG_18 ~ Sal + Temp + Season\n")
print(summary(m1))
print(tidy(m1, conf.int = TRUE))   # clean table with 95% CI
print(glance(m1))                  # R², AIC, etc.

# Check multicollinearity
cat("\nVIF for Model 1:\n"); print(vif(m1))
# Values < 5 are acceptable

# Biological interpretation of Salinity coefficient
beta_sal <- coef(m1)["Salinity"]
cat(sprintf(
  "\nBiological interpretation — Salinity:\n",
  "  β = %.4f\n",
  "  Per 1 psu increase: ×%.2f (%.0f%% change in abundance)\n",
  "  Across full range (0.18–30.52 psu = 30.3 units): %.1f log-units\n"
))
cat(sprintf("  β = %.4f\n", beta_sal))
cat(sprintf("  Per 1 psu: ×%.2f = %.0f%% change\n",
            exp(beta_sal), (exp(beta_sal) - 1) * 100))
cat(sprintf("  Full salinity range: %.1f log-units predicted change\n",
            beta_sal * 30.3))

# ── Model 2: Nanopelagicales ~ Sal + Temp + Season ───────────
m2 <- lm(Nanopelagicales_MAG_18 ~ Salinity + Temperature + Season_num,
          data = reg_df)
cat("\n============================================================\n")
cat("Model 2: Nanopelagicales_MAG_18 ~ Sal + Temp + Season\n")
print(summary(m2))
cat("\nVIF:\n"); print(vif(m2))

# ── Model 3: Flavobacteriales ~ Sal + Temp + Season ──────────
m3 <- lm(Flavobacteriales_MAG_6 ~ Salinity + Temperature + Season_num,
          data = reg_df)
cat("\n============================================================\n")
cat("Model 3: Flavobacteriales_MAG_6 ~ Sal + Temp + Season\n")
cat("Note: High R² (0.879) but no significant individual predictors\n")
cat("      — collinear predictors with jointly strong signal\n")
print(summary(m3))

# ── Model 4: Rhodobacterales ~ Sal + Temp + Season ───────────
m4 <- lm(Rhodobacterales_MAG_15 ~ Salinity + Temperature + Season_num,
          data = reg_df)
cat("\n============================================================\n")
cat("Model 4: Rhodobacterales_MAG_15 ~ Sal + Temp + Season\n")
print(summary(m4))


## ── 6. Interaction model ─────────────────────────────────────

# Tests whether the salinity effect depends on temperature
# (Salinity × Temperature interaction)
m5 <- lm(Pelagibacterales_MAG_18 ~ Salinity * Temperature,
          data = reg_df)
cat("\n============================================================\n")
cat("Model 5: Pelagibacterales_MAG_18 ~ Salinity * Temperature\n")
print(summary(m5))
# Sal:Temp  β = +0.0148  p = 0.015 * → significant interaction

# Compare additive vs interaction with AIC
cat("\nModel comparison (lower AIC = better):\n")
cat("  Additive  AIC:", AIC(m1), "\n")
cat("  Interaction AIC:", AIC(m5), "\n")
cat("  ΔAIC:", AIC(m1) - AIC(m5), "\n")
# ΔAIC > 2 favours the interaction model


## ── 7. Season-split slopes (manual interaction visualisation) ─

cat("\n── Season-split regression slopes — Pelagibacterales ──\n")
for (seas in c("Spring", "Summer")) {
  d   <- reg_df[reg_df$Season == seas, ]
  fit <- lm(Pelagibacterales_MAG_18 ~ Salinity, data = d)
  cat(sprintf("  %s (n=%d):  β_sal = %+.4f  r² = %.3f  p = %.4f\n",
              seas, nrow(d),
              coef(fit)["Salinity"],
              summary(fit)$r.squared,
              summary(fit)$coefficients["Salinity", "Pr(>|t|)"]))
}
# Spring: β ≈ +0.020  r² = 0.036  p = 0.655 ns  — flat
# Summer: β = +0.309  r² = 0.784  p < 0.001    — steep


## ── 8. Model diagnostics for Model 1 ────────────────────────

# Run all four base-R diagnostic plots together
dev.new()
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
plot(m1, which = 1:4,
     col    = "#1d4ed8",
     pch    = 19,
     cex    = 0.8,
     sub.caption = "Model 1: Pelagibacterales ~ Sal + Temp + Season")
par(mfrow = c(1, 1))

# Formal statistical tests for assumptions
cat("\n── Assumption tests — Model 1 ──\n")
sw <- shapiro.test(residuals(m1))
cat(sprintf("  Shapiro-Wilk (normality):  W = %.3f  p = %.4f  %s\n",
            sw$statistic, sw$p.value,
            ifelse(sw$p.value > 0.05, "✓ cannot reject normality",
                   "✗ normality violated")))

bp <- bptest(m1)
cat(sprintf("  Breusch-Pagan (homoscedast.): BP = %.2f  p = %.4f  %s\n",
            bp$statistic, bp$p.value,
            ifelse(bp$p.value > 0.05, "✓ no heteroscedasticity",
                   "✗ heteroscedasticity detected")))

# Check outlier influence
cd <- cooks.distance(m1)
cat(sprintf("  Max Cook's distance: %.3f  (threshold ≈ 1.0)\n", max(cd)))
cat(sprintf("  Influential points (Cook's > 0.5): %d\n", sum(cd > 0.5)))


## ── 9. ggplot2 diagnostic panels ─────────────────────────────

diag_df <- data.frame(
  fitted    = fitted(m1),
  resid     = residuals(m1),
  std_resid = rstandard(m1),
  sqrt_abs  = sqrt(abs(rstandard(m1))),
  cooks     = cooks.distance(m1),
  sample_id = rownames(reg_df)
)

# Panel D — Residuals vs Fitted
p_D <- ggplot(diag_df, aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = "#dc2626",
             linetype = "dashed", linewidth = 1.2) +
  geom_point(colour = "#1d4ed8", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE,
              colour = "#d97706", linewidth = 1.5) +
  labs(x = "Fitted values", y = "Residuals",
       title = "D  Residuals vs Fitted",
       subtitle = "Want: random scatter around 0") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Panel E — Normal Q-Q
p_E <- ggplot(diag_df, aes(sample = std_resid)) +
  stat_qq(colour = "#1d4ed8", size = 2.5, alpha = 0.8) +
  stat_qq_line(colour = "#dc2626",
               linetype = "dashed", linewidth = 1.2) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "E  Normal Q-Q Plot",
       subtitle = "Want: points near the line") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Panel F — Scale-Location (homoscedasticity)
p_F <- ggplot(diag_df, aes(fitted, sqrt_abs)) +
  geom_point(colour = "#1d4ed8", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE,
              colour = "#d97706", linewidth = 1.5) +
  labs(x = "Fitted values", y = "√|Standardised residuals|",
       title = "F  Scale-Location",
       subtitle = "Want: flat horizontal band") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Panel G — Residual histogram with normal curve overlay
p_G <- ggplot(diag_df, aes(resid)) +
  geom_histogram(bins = 11, fill = "#1d4ed8",
                 colour = "white", alpha = 0.65) +
  stat_function(
    fun  = function(x) dnorm(x, mean(diag_df$resid), sd(diag_df$resid)) *
                       nrow(diag_df) * diff(range(diag_df$resid)) / 11,
    colour    = "#dc2626",
    linetype  = "dashed",
    linewidth = 1.8
  ) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "#374151", linewidth = 0.8) +
  labs(x = "Residuals", y = "Count",
       title = "G  Residual Distribution",
       subtitle = "Want: bell-shaped, centred at 0") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))


## ── 10. Interaction visualisation ───────────────────────────

# Separate regression line per season
p_H <- ggplot(reg_df, aes(Salinity, Pelagibacterales_MAG_18,
                            colour = Season)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 2.2) +
  scale_colour_manual(values = c("Spring" = "#16a34a",
                                  "Summer" = "#dc2626")) +
  annotate("text", x = 0.5, y = 1.2,
           label = paste0(
             "Spring  β = +0.020  (r² = 0.036, p = 0.655 ns)\n",
             "Summer β = +0.309  (r² = 0.784, p < 0.001)\n",
             "Δβ = −0.289  →  interaction term needed"
           ),
           hjust = 0, size = 3.3, family = "mono",
           colour = "#111827") +
  labs(
    x        = "Salinity (psu)",
    y        = "log(Pelagibacterales MAG 18 + 1)",
    title    = "H  Salinity × Season Interaction",
    subtitle = "Different slopes per season → additive model is insufficient"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "#6b7280"))


## ── 11. Coefficient forest plot ─────────────────────────────

# Extract tidy coefficients for both models
tidy_m1 <- tidy(m1, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(model = "Pelagib. MAG 18",
         term  = recode(term,
                        "Salinity"    = "Salinity",
                        "Temperature" = "Temperature",
                        "Season_num"  = "Season (Summer=1)"))

tidy_m2 <- tidy(m2, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(model = "Nanopel. MAG 18",
         term  = recode(term,
                        "Salinity"    = "Salinity",
                        "Temperature" = "Temperature",
                        "Season_num"  = "Season (Summer=1)"))

coef_df <- bind_rows(tidy_m1, tidy_m2)

p_C <- ggplot(coef_df,
              aes(estimate, term,
                  colour    = model,
                  shape     = p.value < 0.05,
                  xmin      = conf.low,
                  xmax      = conf.high)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "#9ca3af", linewidth = 0.9) +
  geom_errorbarh(height = 0.15, linewidth = 0.9,
                 position = position_dodge(width = 0.45)) +
  geom_point(size = 3.5,
             position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Pelagib. MAG 18" = "#1d4ed8",
                                  "Nanopel. MAG 18" = "#7c3aed")) +
  scale_shape_manual(values     = c("TRUE" = 18, "FALSE" = 16),
                     labels     = c("TRUE" = "p < 0.05", "FALSE" = "ns"),
                     name       = NULL) +
  labs(x      = "Coefficient value",
       y      = NULL,
       title  = "C  Multiple Regression Coefficients with 95% CI",
       colour = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position   = "right",
        plot.title        = element_text(face = "bold"),
        axis.text.y       = element_text(face = "bold", size = 10))


## ── 12. Model summary table ──────────────────────────────────

# Build summary table as a data frame for display / export
model_summary <- data.frame(
  MAG           = c("Pelagib. MAG 18", "Pelagib. MAG 18",
                    "Nanopel. MAG 18", "Pelagib. MAG 18",
                    "Flavob. MAG 6",   "Rhodb. MAG 15"),
  Predictors    = c("Salinity (simple)",
                    "Sal + Temp + Season",
                    "Sal + Temp + Season",
                    "Sal × Temp (interaction)",
                    "Sal + Temp + Season",
                    "Sal + Temp + Season"),
  R2            = c(0.654, 0.817, 0.815, 0.857, 0.879, 0.655),
  adj_R2        = c(0.640, 0.793, 0.790, 0.838, 0.864, 0.610),
  Key_predictor = c("Salinity",  "Salinity",  "Salinity",
                    "Sal × Temp", "Season",   "Salinity"),
  p_value       = c("< 0.001 ***", "< 0.001 ***", "< 0.001 ***",
                    "= 0.015 *",   "= 0.306 ns",  "< 0.001 ***")
)
print(model_summary)


## ── 13. Loop regression across all 61 MAGs ──────────────────

mag_names <- colnames(mag_log)

results_reg <- do.call(rbind, lapply(mag_names, function(mag) {
  df_m <- cbind(metadata, abundance = mag_log[, mag])
  df_m <- df_m[complete.cases(df_m[, c("Salinity","Temperature",
                                         "Season_num","abundance")]), ]
  if (nrow(df_m) < 10) return(NULL)

  m  <- lm(abundance ~ Salinity + Temperature + Season_num, data = df_m)
  sm <- summary(m)
  cf <- coef(sm)

  data.frame(
    MAG        = mag,
    n          = nrow(df_m),
    R2         = round(sm$r.squared,      4),
    adj_R2     = round(sm$adj.r.squared,  4),
    beta_sal   = round(cf["Salinity",    "Estimate"], 4),
    se_sal     = round(cf["Salinity",    "Std. Error"], 4),
    p_sal      = cf["Salinity",          "Pr(>|t|)"],
    beta_temp  = round(cf["Temperature", "Estimate"], 4),
    p_temp     = cf["Temperature",       "Pr(>|t|)"],
    beta_seas  = round(cf["Season_num",  "Estimate"], 4),
    p_seas     = cf["Season_num",        "Pr(>|t|)"]
  )
}))

# FDR correction (BH) across all 61 MAGs per predictor
results_reg$q_sal  <- p.adjust(results_reg$p_sal,  method = "BH")
results_reg$q_temp <- p.adjust(results_reg$p_temp, method = "BH")
results_reg$q_seas <- p.adjust(results_reg$p_seas, method = "BH")

# Summary
cat("\n── All 61 MAGs — additive model results ──\n")
cat("Significant Salinity effect   (q < 0.05):", sum(results_reg$q_sal  < 0.05), "\n")
cat("Significant Temperature effect (q < 0.05):", sum(results_reg$q_temp < 0.05), "\n")
cat("Significant Season effect      (q < 0.05):", sum(results_reg$q_seas < 0.05), "\n")

# Top MAGs by salinity effect
cat("\nTop 10 MAGs by |β_Salinity| (q < 0.05):\n")
top_sal <- results_reg %>%
  filter(q_sal < 0.05) %>%
  arrange(desc(abs(beta_sal)))
print(head(top_sal[, c("MAG","R2","beta_sal","q_sal")], 10))

# Save full results
dir.create("/day5-regression/results",
           showWarnings = FALSE, recursive = TRUE)
write.csv(
  results_reg,
  "day5_regression_all_MAGs.csv",
  row.names = FALSE
)
cat("\nResults saved.\n")


## ── 14. Volcano-style plot: β_Sal vs −log10(q) ──────────────

results_reg$neg_log10_q <- -log10(results_reg$q_sal + 1e-10)
results_reg$direction   <- ifelse(results_reg$beta_sal > 0,
                                   "Salinity+", "Salinity-")
results_reg$sig         <- results_reg$q_sal < 0.05

p_volcano <- ggplot(results_reg,
                    aes(beta_sal, neg_log10_q,
                        colour = interaction(sig, direction),
                        size   = sig)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_text_repel(
    data   = results_reg[results_reg$sig, ],
    aes(label = MAG),
    size   = 2.5, max.overlaps = 20,
    colour = "grey20"
  ) +
  scale_colour_manual(
    values = c("TRUE.Salinity+"  = "#1d4ed8",
               "TRUE.Salinity-"  = "#7c3aed",
               "FALSE.Salinity+" = "grey75",
               "FALSE.Salinity-" = "grey75"),
    name   = NULL,
    labels = c("TRUE.Salinity+"  = "Salinity+ (q<0.05)",
               "TRUE.Salinity-"  = "Salinity- (q<0.05)",
               "FALSE.Salinity+" = "ns",
               "FALSE.Salinity-" = "ns")
  ) +
  scale_size_manual(values = c("TRUE" = 2.8, "FALSE" = 1.5),
                    guide = "none") +
  labs(
    x        = "β_Salinity  (log-abundance change per 1 psu)",
    y        = "-log10(q-value, BH corrected)",
    title    = "All 61 MAGs — Salinity Effect Size vs Significance",
    subtitle = "Model: log(abundance+1) ~ Salinity + Temperature + Season"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))


## ── 15. Assemble and save figure ────────────────────────────

# Row 1: scatter panels + coefficient forest
row1 <- (p_A | p_B | p_C) +
  plot_layout(widths = c(1, 1, 1.4))

# Row 2: four diagnostic panels
row2 <- (p_D | p_E | p_F | p_G)

# Row 3: interaction + volcano
row3 <- (p_H | p_volcano) +
  plot_layout(widths = c(1.2, 1))

combined <- row1 / row2 / row3 +
  plot_annotation(
    title    = "Day 5 — Multiple Regression: Linking MAG Abundance to Environment",
    subtitle = paste0("61 MAGs  ·  27 metagenomes  ·  OLS linear regression  ·  ",
                      "log(abundance + 1)  ·  Chesapeake Bay + Delaware Bay"),
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size =  9, colour = "#6b7280")
    )
  )

dir.create(
  "day5-regression/figures",
  showWarnings = FALSE, recursive = TRUE
)

ggsave(
  "day5_regression_figure.tiff",
  combined,
  width = 18, height = 16, dpi = 300, compression = "lzw"
)
cat("Figure saved.\n")


## ── 16. Optional: Robust regression for comparison ──────────
## If diagnostics show heavy-tailed residuals or influential outliers,
## compare OLS with M-estimation via MASS::rlm()

# library(MASS)
# m1_robust <- rlm(Pelagibacterales_MAG_18 ~ Salinity + Temperature + Season_num,
#                  data = reg_df, method = "M")
# summary(m1_robust)
# Compare coefficients:
# cbind(OLS = coef(m1), Robust = coef(m1_robust))


## ── 17. Optional: Generalised linear model (Gamma/log-link) ─
## For MAGs with heavy zero-inflation, GLM may fit better than OLS + log1p

# library(glmmTMB)  # or use base glm() for simple cases
# m_glm <- glm(Pelagibacterales_MAG_18_raw ~ Salinity + Temperature + Season_num,
#              data   = reg_df,
#              family = Gamma(link = "log"),
#              subset = Pelagibacterales_MAG_18_raw > 0)
# summary(m_glm)
# # Coefficients now on log-count scale — same interpretation as OLS on log1p


## ── Session info ─────────────────────────────────────────────
cat("\n── Session info ──\n")
sessionInfo()

## ── End of Day 5 script ──────────────────────────────────────
