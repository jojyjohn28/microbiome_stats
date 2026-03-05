## ============================================================
## Day 4 — Non-Parametric Tests: Wilcoxon, Kruskal-Wallis
##          & Cliff's delta effect size
## Applied Statistics for Microbiome & Genomics Data
## Author: Jojy John
## Series: github.com/jojyjohn28/microbiome_stats
## ============================================================

library(vegan)  
library(ggplot2)
library(ggrepel)
library(effsize)
library(dplyr) 
library(tidyr)
library(patchwork)

## ── 1. Load and align data ───────────────────────────────────
metadata <- read.csv("~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/metadata.csv",  header=TRUE, row.names=1, sep=",")
mag_raw  <- read.csv("~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day1/data/mag_abundance.csv", header=TRUE, row.names=1, sep=",")
mag_t    <- t(mag_raw)
metadata <- metadata[match(rownames(mag_t), rownames(metadata)), ]
stopifnot(all(rownames(mag_t) == rownames(metadata)))

## ── 2. Group variables ───────────────────────────────────────
Season   <- metadata$Season
Fraction <- metadata$SF
Salinity_group <- ifelse(metadata$Salinity < 15, "Low_salinity", "High_salinity")
mag_names <- colnames(mag_t)

cat("Total MAGs:", length(mag_names),
    "| Total tests:", length(mag_names)*3,
    "| Expected false positives (uncorrected):", length(mag_names)*3*0.05, "\n")

## ── 3. Single MAG walkthrough ────────────────────────────────
mag_abund    <- mag_t[, "Pelagibacterales_MAG_18"]
spring_abund <- mag_abund[Season == "Spring"]
summer_abund <- mag_abund[Season == "Summer"]

boxplot(list(Spring = spring_abund, Summer = summer_abund),
        main="Pelagibacterales_MAG_18 by Season", ylab="Abundance",
        col=c("#22c55e","#ef4444"))

wt <- wilcox.test(spring_abund, summer_abund, exact=FALSE)
cd <- cliff.delta(summer_abund, spring_abund)
cat("p =", wt$p.value, "| Cliff's d =", round(cd$estimate,3), "(", cd$magnitude, ")\n")

## ── 4. All MAGs: Season ──────────────────────────────────────
results_season <- do.call(rbind, lapply(mag_names, function(mag) {
  wt <- wilcox.test(mag_t[Season=="Spring", mag],
                    mag_t[Season=="Summer", mag], exact=FALSE)
  cd <- cliff.delta(mag_t[Season=="Summer", mag],
                    mag_t[Season=="Spring", mag])
  data.frame(MAG=mag, p_value=wt$p.value,
             cliffs_d=cd$estimate, magnitude=as.character(cd$magnitude))
}))
results_season$p_adj <- p.adjust(results_season$p_value, method="BH")
results_season <- results_season[order(results_season$p_adj), ]
sig_season <- results_season[results_season$p_adj < 0.05, ]
cat("\nSeason — significant MAGs:", nrow(sig_season), "\n"); print(sig_season)

## ── 5. All MAGs: Size fraction ───────────────────────────────
results_sf <- do.call(rbind, lapply(mag_names, function(mag) {
  wt <- wilcox.test(mag_t[Fraction=="Particle-attached", mag],
                    mag_t[Fraction=="Free_living",        mag], exact=FALSE)
  cd <- cliff.delta(mag_t[Fraction=="Particle-attached", mag],
                    mag_t[Fraction=="Free_living",        mag])
  data.frame(MAG=mag, p_value=wt$p.value,
             cliffs_d=cd$estimate, magnitude=as.character(cd$magnitude))
}))
results_sf$p_adj <- p.adjust(results_sf$p_value, method="BH")
results_sf <- results_sf[order(results_sf$p_adj), ]
sig_sf <- results_sf[results_sf$p_adj < 0.05, ]
cat("\nSize Fraction — significant MAGs:", nrow(sig_sf), "\n"); print(sig_sf)

## ── 6. All MAGs: Salinity ────────────────────────────────────
results_sal <- do.call(rbind, lapply(mag_names, function(mag) {
  wt <- wilcox.test(mag_t[Salinity_group=="Low_salinity",  mag],
                    mag_t[Salinity_group=="High_salinity", mag], exact=FALSE)
  cd <- cliff.delta(mag_t[Salinity_group=="High_salinity", mag],
                    mag_t[Salinity_group=="Low_salinity",  mag])
  data.frame(MAG=mag, p_value=wt$p.value,
             cliffs_d=cd$estimate, magnitude=as.character(cd$magnitude))
}))
results_sal$p_adj <- p.adjust(results_sal$p_value, method="BH")
results_sal <- results_sal[order(results_sal$p_adj), ]
sig_sal <- results_sal[results_sal$p_adj < 0.05, ]
cat("\nSalinity — significant MAGs:", nrow(sig_sal), "\n"); print(sig_sal)

## ── 7. Kruskal-Wallis: 3 salinity zones ─────────────────────
Salinity_3 <- cut(metadata$Salinity, breaks=c(-Inf,5,18,Inf),
                  labels=c("Oligohaline","Mesohaline","Polyhaline"))

results_kw <- do.call(rbind, lapply(mag_names, function(mag) {
  kw <- kruskal.test(mag_t[, mag] ~ Salinity_3)
  data.frame(MAG=mag, p_value=kw$p.value, chi_sq=as.numeric(kw$statistic))
}))
results_kw$p_adj <- p.adjust(results_kw$p_value, method="BH")
sig_kw <- results_kw[results_kw$p_adj < 0.05, ]

for (mag in sig_kw$MAG) {
  cat("\nPost-hoc:", mag, "\n")
  print(pairwise.wilcox.test(mag_t[, mag], Salinity_3,
                              p.adjust.method="BH", exact=FALSE)$p.value)
}

## ── 8. Save results table ────────────────────────────────────
all_results <- data.frame(
  MAG         = mag_names,
  season_padj = results_season$p_adj[match(mag_names, results_season$MAG)],
  season_d    = results_season$cliffs_d[match(mag_names, results_season$MAG)],
  sf_padj     = results_sf$p_adj[match(mag_names, results_sf$MAG)],
  sf_d        = results_sf$cliffs_d[match(mag_names, results_sf$MAG)],
  sal_padj    = results_sal$p_adj[match(mag_names, results_sal$MAG)],
  sal_d       = results_sal$cliffs_d[match(mag_names, results_sal$MAG)]
)
all_results$any_sig <- with(all_results,
  season_padj < 0.05 | sf_padj < 0.05 | sal_padj < 0.05)

dir.create("results", showWarnings=FALSE)
write.csv(all_results, "~/Jojy_Research_Sync/website_assets/projects/microbiome_stats/day4-nonparametric-tests-wilcoxon-kruskal-wallis/day4_wilcoxon_results.csv", row.names=FALSE)

## ── 9. Figure A: Volcano plot ────────────────────────────────
vdf <- results_season
vdf$neg_log10_padj <- -log10(vdf$p_adj + 1e-10)
vdf$sig   <- vdf$p_adj < 0.05
vdf$label <- ifelse(vdf$sig & abs(vdf$cliffs_d) >= 0.474,
                     as.character(vdf$MAG), NA)

p_volcano <- ggplot(vdf, aes(cliffs_d, neg_log10_padj, colour=sig, size=sig)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey50") +
  geom_vline(xintercept=c(-0.474, 0.474), linetype="dashed", colour="grey50") +
  geom_text_repel(aes(label=label), size=3, max.overlaps=15, colour="grey20") +
  scale_colour_manual(values=c("FALSE"="grey75","TRUE"="#ef4444"),
                       labels=c("Not significant","q < 0.05")) +
  scale_size_manual(values=c("FALSE"=1.5,"TRUE"=2.5), guide="none") +
  annotate("text", x=-0.85, y=0.5, label="Spring-enriched",
           colour="#22c55e", fontface="italic", size=3.5) +
  annotate("text", x= 0.85, y=0.5, label="Summer-enriched",
           colour="#ef4444", fontface="italic", size=3.5) +
  labs(x="Cliff's delta  (negative = Spring, positive = Summer)",
       y="-log10(adjusted p-value)",
       title="Differential abundance — Spring vs Summer", colour=NULL) +
  theme_classic(base_size=13) +
  theme(legend.position="top", plot.title=element_text(face="bold"))

print(p_volcano)

## ── 10. Figure B: Top MAG boxplots ───────────────────────────
top_mags  <- head(sig_season$MAG, 6)
plot_long <- pivot_longer(
  cbind(as.data.frame(mag_t[, top_mags, drop=FALSE]), Season=Season),
  -Season, names_to="MAG", values_to="Abundance")

p_box <- ggplot(plot_long, aes(Season, Abundance, fill=Season)) +
  geom_boxplot(width=0.5, outlier.shape=21, outlier.size=1.5, alpha=0.85) +
  geom_jitter(width=0.1, size=1.5, alpha=0.5, colour="grey30") +
  scale_fill_manual(values=c("Spring"="#22c55e","Summer"="#ef4444")) +
  facet_wrap(~MAG, scales="free_y", ncol=3) +
  labs(x=NULL, y="Abundance",
       title="Top differentially abundant MAGs — Spring vs Summer") +
  theme_classic(base_size=11) +
  theme(legend.position="none",
        strip.text=element_text(face="italic", size=8),
        plot.title=element_text(face="bold"))
p_box
## ── 11. Save figure ──────────────────────────────────────────
dir.create("figures", showWarnings=FALSE)
ggsave("figures/day4_nonparametric_figure.tiff",
       p_volcano / p_box + plot_annotation(tag_levels="A"),
       width=12, height=10, dpi=300, compression="lzw")

cat("\nDone. Results in results/ | Figure in figures/\n")

#################################################################
# =========================
# Heatmap with q-value side bars
# =========================
library(dplyr)
library(tidyr)
library(ggplot2)

# choose MAGs to show (limit to keep it readable)
show_mags <- all_results %>%
  filter(any_sig) %>%
  arrange(pmin(season_padj, sf_padj, sal_padj, na.rm = TRUE)) %>%
  slice_head(n = 30) %>%              # show top 30; adjust to 20/40 if needed
  pull(MAG)

# build abundance long table
hm <- as.data.frame(mag_t[, show_mags, drop = FALSE])
hm$Sample <- rownames(hm)
hm$Season <- Season
hm <- hm %>%
  pivot_longer(cols = all_of(show_mags), names_to = "MAG", values_to = "Abundance")

# scale per MAG (z-score) for display
hm <- hm %>%
  group_by(MAG) %>%
  mutate(z = as.numeric(scale(Abundance))) %>%
  ungroup()

# q-value annotation tiles (convert q to -log10(q) for visibility)
ann <- all_results %>%
  filter(MAG %in% show_mags) %>%
  mutate(
    season_sig = -log10(season_padj + 1e-12),
    sf_sig     = -log10(sf_padj + 1e-12),
    sal_sig    = -log10(sal_padj + 1e-12)
  ) %>%
  select(MAG, season_sig, sf_sig, sal_sig) %>%
  pivot_longer(-MAG, names_to = "Test", values_to = "neglog10_q") %>%
  mutate(Test = recode(Test,
                       season_sig = "Season q",
                       sf_sig     = "SF q",
                       sal_sig    = "Salinity q"))

# reorder MAGs by strongest signal
mag_order <- ann %>%
  group_by(MAG) %>%
  summarise(maxsig = max(neglog10_q, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(maxsig)) %>%
  pull(MAG)

hm$MAG  <- factor(hm$MAG, levels = mag_order)
ann$MAG <- factor(ann$MAG, levels = mag_order)

# Panel B: abundance heatmap
p_heat <- ggplot(hm, aes(x = Sample, y = MAG, fill = z)) +
  geom_tile() +
  labs(title = "Top significant MAGs (scaled abundance)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title = element_text(face = "bold"))
p_heat
# Panel C: q-value bars (left-side idea)
p_q <- ggplot(ann, aes(x = Test, y = MAG, fill = neglog10_q)) +
  geom_tile() +
  labs(title = "Significance (-log10 q)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# Combine Heatmap + q tiles side-by-side
p_heatmap_combo <- p_q + p_heat + plot_layout(widths = c(1.2, 4))
p_heatmap_combo

combined <- ( p_heat / p_box) +
  plot_annotation(tag_levels = "A")

combined 

#significant only table
sig_table <- all_results %>%
  filter(any_sig) %>%
  mutate(
    season = sprintf("q=%.3g, d=%.2f", season_padj, season_d),
    sf     = sprintf("q=%.3g, d=%.2f", sf_padj, sf_d),
    sal    = sprintf("q=%.3g, d=%.2f", sal_padj, sal_d)
  ) %>%
  select(MAG, season, sf, sal) %>%
  arrange(MAG)

sig_table

write.csv(sig_table, "results/day4_significant_summary_table.csv", row.names = FALSE)


