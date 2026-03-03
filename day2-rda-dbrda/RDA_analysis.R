library(vegan)
library(ggplot2)
library(ggrepel)

#Step1 Load data

# Load metadata — rows = samples, columns = environmental variables
metadata <- read.csv("data/metadata.csv", header = TRUE, row.names = 1, sep = ",")

# Load MAG abundance — rows = MAGs, columns = samples
mag_raw <- read.csv("data/mag_abundance.csv", header = TRUE, row.names = 1, sep = ",")
```

# Step 2 — Prepare the community matrix

# The abundance table has MAGs as rows and samples as columns
# RDA needs samples as rows — so we transpose
mag_t <- t(mag_raw)   # now 27 rows (samples) × 61 columns (MAGs)

# Hellinger transformation:
# Divides each value by the row sum (converts to relative abundance)
# then takes the square root — reduces dominance of highly abundant MAGs
spe.hel <- decostand(mag_t, method = "hellinger")

# Step 3 — Prepare the environmental matrix
# Select continuous environmental variables only
# metadata now has: Season, SF, Bay (categorical) — exclude these
# keep only numeric predictors for the RDA
env <- metadata[, c("Salinity", "Temperature", "Depth", "PAR",
                    "Attenuation", "Bacterial_Production", "nCells",
                    "ChlA", "Nitrate", "Ammonium", "Phosphate", "Silicate")]

# Standardize to mean = 0, sd = 1 so all variables are comparable
# (Salinity in ppt and nCells in millions are on completely different scales)
env.z <- decostand(env, method = "stand")
```

### Step 4 — Forward selection of significant variables

Rather than throwing all 13 variables into the model, we use forward selection to find only those that significantly and independently explain community variance:
  
  ```r
# Full model (upper bound)
spe.rda.full <- rda(spe.hel ~ ., data = env.z)

# Forward selection — adds variables one by one, keeps only significant ones
# stops when adding more variables would exceed the full model's R²
fwd.sel <- ordiR2step(
  rda(spe.hel ~ 1, data = env.z),   # null model (lower bound)
  scope = formula(spe.rda.full),     # upper bound
  direction = "forward",
  R2scope = TRUE,
  pstep = 1000,
  trace = FALSE
)

# See which variables were selected
fwd.sel$call

# Step 5 — Build the final model
# Final RDA model with forward-selected variables
spe.rda.signif <- rda(spe.hel ~ Temperature + Salinity + Phosphate +
                        ChlA + Silicate + Nitrate + SF + nCells,
                      data = env.z)

# Adjusted R² — how much community variance is explained?
Rsquareadj(spe.rda.signif)

# Test overall model significance (permutation test)
anova.cca(spe.rda.signif, step = 1000)

# Test significance of each individual term
anova.cca(spe.rda.signif, step = 1000, by = "term")

# Test significance of each axis
anova.cca(spe.rda.signif, step = 1000, by = "axis")

# Step 6 — Extract scores and build the triplot
# Extract % variance explained by RDA1 and RDA2
perc <- round(100 * (summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

# Season, Bay and size fraction — now clean columns in metadata
Season   <- metadata$Season   # "Spring" or "Summer"
Bay      <- metadata$Bay      # "Chesapeake" or "Delaware"
Fraction <- metadata$SF       # "Particle-attached" or "Free_living"

# Site scores — split by season for colour-coding
sc_si_summer <- scores(spe.rda.signif, display = "sites")[Season == "Summer", ]
sc_si_spring <- scores(spe.rda.signif, display = "sites")[Season == "Spring", ]

# Species (MAG) scores
sc_sp <- scores(spe.rda.signif, display = "species", choices = c(1, 2), scaling = 1)

# Biplot arrows — environmental variables
sc_bp <- scores(spe.rda.signif, display = "bp", choices = c(1, 2), scaling = 1)

# Only label MAGs with strong contributions (|score| > 0.35 on either axis)
high_score_mags <- rownames(sc_sp)[rowSums(abs(sc_sp) > 0.35) > 0]

# Build the triplot
plot(spe.rda.signif,
     scaling = 1, type = "none", frame = FALSE,
     xlim = c(-1.0, 1.5), ylim = c(-0.5, 0.5),
     xlab = paste0("RDA1 (", perc[1], "%)"),
     ylab = paste0("RDA2 (", perc[2], "%)"))

# Summer samples — red
points(sc_si_summer, pch = 21, col = "black", bg = "red2", cex = 1.2)

# Spring samples — green
points(sc_si_spring, pch = 21, col = "black", bg = "green3", cex = 1.2)

# MAG positions — yellow squares
points(sc_sp, pch = 22, col = "black", bg = "#f2bd33", cex = 1.2)

# Labels for most influential MAGs only
text(sc_sp[high_score_mags, ] + c(0.20, 0.15),
     labels = high_score_mags, col = "grey4", font = 3, cex = 0.6)

# Environmental arrows
arrows(0, 0, sc_bp[, 1], sc_bp[, 2], col = "blue", lwd = 1.2)
text(x = sc_bp[, 1] - 0.1, y = sc_bp[, 2] - 0.03,
     labels = rownames(sc_bp), col = "blue", cex = 1, font = 2)