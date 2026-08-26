# ── Load & clean ──────────────────────────────────────────────────────────────
df <- read.csv("pop_wings.csv", stringsAsFactors = FALSE, na.strings = c("", "NA", "na"))
df$Species <- trimws(df$Species)

keep <- c("Lycaeides anna", "Lycaeides anna lotis", "Lycaeides idas")
df <- df[df$Species %in% keep, ]

# ── Assign populations ────────────────────────────────────────────────────────
# L. anna split by site (YG = Yuga Gap, LS = Leek Springs).
# L. anna lotis and L. idas kept as single groups.
# Individuals with no site or ID are dropped.
assign_pop <- function(species, site, id) {
  site <- ifelse(is.na(site), "", site)
  id   <- ifelse(is.na(id),   "", id)
  ifelse(
    species == "Lycaeides anna lotis", "L. anna lotis",
    ifelse(
      species == "Lycaeides idas", "L. idas",
      ifelse(
        grepl("YG", site) | grepl("^YG", id), "L. anna (YG)",
        ifelse(
          grepl("LS", site) | grepl("^LS", id), "L. anna (LS)",
          NA
        )
      )
    )
  )
}

df$Population <- assign_pop(df$Species, df$Site, df$ID..)
df <- df[!is.na(df$Population), ]

# ── Feature columns ───────────────────────────────────────────────────────────
spot_cols <- c("Spot.1", "Spot.2", "Spot.3", "Spot.5", "Spot.6", "Spot.7",
               "Spot.8", "Spot.9", "Spot.10", "Spot.11", "Spot.12", "Spot.13")
feat_cols <- c("Total.wing.area", spot_cols)

for (col in feat_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

# ── Per-population mean imputation ────────────────────────────────────────────
for (col in feat_cols) {
  for (pop in unique(df$Population)) {
    idx   <- df$Population == pop
    pmean <- mean(df[[col]][idx], na.rm = TRUE)
    df[[col]][idx & is.na(df[[col]])] <- pmean
  }
}

# ── PCA ───────────────────────────────────────────────────────────────────────
X      <- as.matrix(df[, feat_cols])
labels <- df$Population

pc  <- prcomp(X, center = TRUE, scale. = TRUE)
o   <- summary(pc)
pct <- round(o$importance[2, 1:3] * 100, 1)

# Export scores
scores            <- as.data.frame(pc$x[, 1:3])
scores$Population <- labels
scores$ID         <- df$ID..
write.csv(scores, "pca_scores.csv", row.names = FALSE)

# ── Colours & symbols ─────────────────────────────────────────────────────────
pop_levels <- c("L. anna (YG)", "L. anna (LS)", "L. anna lotis", "L. idas")
pop_cols   <- c("L. anna (YG)"  = "#185FA5",
                "L. anna (LS)"  = "#73B8F2",
                "L. anna lotis" = "#0F6E56",
                "L. idas"       = "#993C1D")
pop_pch    <- c("L. anna (YG)"  = 19,
                "L. anna (LS)"  = 19,
                "L. anna lotis" = 17,
                "L. idas"       = 19)

pt_col <- pop_cols[labels]
pt_pch <- pop_pch[labels]

# ── Plot ──────────────────────────────────────────────────────────────────────
cm <- 1.4; ca <- 1.1; cl <- 1.2

pdf("pca_wings.pdf", width = 5, height = 5)
par(mar = c(4.5, 5, 2.5, 1.5))

pad  <- 1.5
xlim <- range(pc$x[, 1]) + c(-pad, pad)
ylim <- range(pc$x[, 2]) + c(-pad, pad)

plot(pc$x[, 1], pc$x[, 2],
     pch     = pt_pch,
     col     = pt_col,
     xlim    = xlim,
     ylim    = ylim,
     xlab    = paste0("PC1 (", pct[1], "%)"),
     ylab    = paste0("PC2 (", pct[2], "%)"),
     cex     = 1.2,
     cex.lab = cl,
     cex.axis= ca,
     las     = 1)

title(main = "PCA wing measurements", cex.main = cm)

abline(h = 0, v = 0, col = "gray80", lwd = 0.8)

# ── Confidence ellipses (1 SD) ────────────────────────────────────────────────
draw_ellipse <- function(x, y, n_sd = 1, col, lty = 1, n_pts = 200) {
  if (length(x) < 3) return(invisible(NULL))
  S        <- cov(cbind(x, y))
  mu       <- c(mean(x), mean(y))
  ev       <- eigen(S)
  theta    <- seq(0, 2 * pi, length.out = n_pts)
  circ     <- rbind(cos(theta), sin(theta))
  ell      <- t(ev$vectors %*% diag(n_sd * sqrt(ev$values)) %*% circ)
  ell[, 1] <- ell[, 1] + mu[1]
  ell[, 2] <- ell[, 2] + mu[2]
  rgb_vals <- col2rgb(col) / 255
  fill_col <- rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha = 0.12)
  polygon(ell, border = col, col = fill_col, lty = lty, lwd = 1.5)
}

for (pop in pop_levels) {
  idx <- labels == pop
  lty <- ifelse(pop == "L. anna (LS)", 2, 1)
  draw_ellipse(pc$x[idx, 1], pc$x[idx, 2],
               n_sd = 1, col = pop_cols[pop], lty = lty)
}

# redraw points on top of ellipses
points(pc$x[, 1], pc$x[, 2], pch = pt_pch, col = pt_col, cex = 1.2)

legend("topright",
       legend = pop_levels,
       pch    = pop_pch[pop_levels],
       col    = pop_cols[pop_levels],
       ncol   = 1,
       cex    = 0.75,
       bty    = "n")

dev.off()
cat("Saved to pca_wings.pdf\n")
