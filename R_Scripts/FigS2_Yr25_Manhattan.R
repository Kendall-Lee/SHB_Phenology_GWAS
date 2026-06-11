#!/usr/bin/env Rscript
# FigS2_Yr25_Manhattan.R
# Supplementary Figure S2 — Year 2025 stress GWAS genome-wide Manhattan
# Four traits (DTFlower, DTFruit, Flow2Fruit, FruitWeight) from yr.25 GWASpoly LOCO scan.
# The Plant Genome specs: double-column 180 mm (7.087 in), 300 DPI, min 8 pt at print.
#
# Data: LINEAR/ALL_TRAITS/GWASpoly_scans.RData  →  data.loco.scan (scan_all)
# Output: Linear_MS/Figures/FigS2_Yr25_Manhattan.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(GWASpoly)
  library(scales)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS/Figures")
THRESH  <- 5.52   # Meff LOCO threshold (same as main analysis)

# Trait colours matching main manuscript
SIG_COLS <- c(
  DTFlower    = "#8B0000",
  DTFruit     = "#1565C0",
  Flow2Fruit  = "#1A5276",
  FruitWeight = "#4A235A"
)
CHR_COLS <- c("#5D6D7E", "#AEB6BF")

pub_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text    = element_text(color = "black"),
    axis.title   = element_text(color = "black", face = "bold"),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
    plot.margin  = margin(4, 6, 4, 4)
  )

# ── Load yr.25 scan ────────────────────────────────────────────────────────────
cat("Loading scan_all (yr.25 traits)...\n")
load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"))
scan_all <- data.loco.scan
cat("Loaded. Traits available:", paste(names(scan_all@scores), collapse=", "), "\n")

# ── Cumulative position map ────────────────────────────────────────────────────
build_cumpos <- function(scan) {
  map <- as.data.table(scan@map)
  names(map)[1:3] <- c("Marker", "Chrom", "Position")
  map[, Position := as.numeric(Position)]
  map <- map[!grepl("[.][0-9]+$", Marker)]
  map <- unique(map, by = "Marker")
  chr_order <- paste0("Chr.", sprintf("%02d", 1:12))
  chr_order <- chr_order[chr_order %in% unique(map$Chrom)]
  map[, Chrom := factor(Chrom, levels = chr_order)]
  setorder(map, Chrom, Position)
  chr_max <- map[, .(ChrLen = max(Position)), by = Chrom]
  chr_max <- chr_max[order(Chrom)]
  chr_max[, CumStart := c(0, cumsum(as.numeric(ChrLen))[-.N])]
  chr_max[, CumStart := CumStart + seq_len(.N) * 5e6]
  map <- merge(map, chr_max[, .(Chrom, CumStart)], by = "Chrom")
  map[, CumPos := Position + CumStart]
  chr_mids <- map[, .(Mid = mean(CumPos)), by = Chrom]
  list(map = map, chr_mids = chr_mids)
}

cp <- build_cumpos(scan_all)

# ── Extract scores — additive model only for yr.25 (avoids rare-dosage inflation) ──
extract_scores <- function(scan, trait, cumpos) {
  s <- scan@scores[[trait]]
  if (is.null(s)) { cat("  Trait not found:", trait, "\n"); return(NULL) }
  markers <- rownames(s)
  markers <- markers[!grepl("[.][0-9]+$", markers)]
  s_sub <- s[markers, , drop = FALSE]
  # Use additive model only; if absent fall back to max across models
  best <- if ("additive" %in% colnames(s_sub)) {
    as.numeric(s_sub[, "additive"])
  } else {
    apply(s_sub, 1, max, na.rm = TRUE)
  }
  dt <- data.table(Marker = markers, Score = best)
  merge(dt, cumpos$map[, .(Marker, Chrom, CumPos)], by = "Marker")
}

# ── Single-trait Manhattan panel ──────────────────────────────────────────────
make_panel <- function(dt, cumpos, trait_label, sig_col,
                       n_sig, show_xlab = FALSE, show_ylab = TRUE) {
  chr_mids <- cumpos$chr_mids
  dt[, ChrIdx := as.integer(Chrom)]
  dt[, ChrCol := ifelse(ChrIdx %% 2 == 1, CHR_COLS[1], CHR_COLS[2])]
  dt[, Sig    := Score >= THRESH]

  y_max <- max(max(dt$Score, na.rm = TRUE) * 1.10, THRESH * 1.4)
  xlabs <- gsub("^Chr\\.0?", "", as.character(chr_mids$Chrom))

  ggplot(dt) +
    geom_point(data = dt[Sig == FALSE],
               aes(x = CumPos, y = Score, colour = ChrCol),
               size = 0.4, alpha = 0.45, stroke = 0) +
    geom_point(data = dt[Sig == TRUE],
               aes(x = CumPos, y = Score),
               colour = sig_col, size = 1.3, alpha = 0.85, stroke = 0) +
    geom_hline(yintercept = THRESH, linetype = "dashed",
               colour = "grey50", linewidth = 0.55) +
    scale_colour_identity() +
    scale_x_continuous(breaks = chr_mids$Mid, labels = xlabs,
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = pretty_breaks(n = 4),
                       expand = c(0, 0)) +
    labs(
      title = trait_label,
      subtitle = sprintf("yr.2025  |  n=536  |  %d QTL (−log₁₀p ≥ %.2f)", n_sig, THRESH),
      x = if (show_xlab) "Chromosome" else NULL,
      y = if (show_ylab) expression(bold(-log[10](italic(p)))) else NULL
    ) +
    pub_theme +
    theme(
      panel.border  = element_rect(color = sig_col, fill = NA, linewidth = 1.1),
      plot.title    = element_text(size = 11, face = "bold", colour = sig_col),
      plot.subtitle = element_text(size = 8,  colour = "grey45"),
      axis.text.x   = element_text(size = 8),
      axis.text.y   = element_text(size = 8),
      axis.title    = element_text(size = 9, face = "bold")
    )
}

# ── Build panels ──────────────────────────────────────────────────────────────
traits_yr25 <- c(
  DTFlower    = "DTFlower_yr.25",
  DTFruit     = "DTFruit_yr.25",
  Flow2Fruit  = "Flow2Fruit_yr.25",
  FruitWeight = "FruitWeight_yr.25"
)

panels <- list()
for (nm in names(traits_yr25)) {
  cat("Extracting", traits_yr25[nm], "...\n")
  dt <- extract_scores(scan_all, traits_yr25[nm], cp)
  if (is.null(dt) || nrow(dt) == 0) next
  n_sig  <- sum(dt$Score >= THRESH, na.rm = TRUE)
  is_bot <- nm %in% c("Flow2Fruit", "FruitWeight")
  panels[[nm]] <- make_panel(
    dt, cp,
    trait_label = c(DTFlower    = "Days to 50% Flowering",
                    DTFruit     = "Days to 50% Ripe Fruit",
                    Flow2Fruit  = "Flower-to-Fruit Interval",
                    FruitWeight = "25-Fruit Weight")[nm],
    sig_col   = SIG_COLS[nm],
    n_sig     = n_sig,
    show_xlab = is_bot,
    show_ylab = nm %in% c("DTFlower", "Flow2Fruit")
  )
}

# ── Title strip ───────────────────────────────────────────────────────────────
title_strip <- ggdraw() +
  draw_label(
    "Year 2025 stress GWAS — genome-wide Manhattan (GWASpoly LOCO, additive model)",
    fontface = "bold", size = 11, x = 0.02, hjust = 0
  ) +
  draw_label(
    paste0("yr.25 individual-year BLUEs (n=536); additive model displayed; Meff threshold = 5.52. ",
           "281 total significant QTL across all dosage models (see Methods)."),
    size = 8, colour = "grey40", x = 0.02, y = 0.25, hjust = 0
  )

# ── Assemble 2×2 grid ─────────────────────────────────────────────────────────
top_row <- if (!is.null(panels$DTFlower) && !is.null(panels$DTFruit))
  plot_grid(panels$DTFlower, panels$DTFruit, nrow = 1, labels = c("A","B"),
            label_size = 13, label_fontface = "bold")
bot_row <- if (!is.null(panels$Flow2Fruit) && !is.null(panels$FruitWeight))
  plot_grid(panels$Flow2Fruit, panels$FruitWeight, nrow = 1, labels = c("C","D"),
            label_size = 13, label_fontface = "bold")

grid <- plot_grid(top_row, bot_row, ncol = 1, rel_heights = c(1, 1))

fig_final <- plot_grid(
  title_strip, grid,
  ncol = 1, rel_heights = c(0.055, 1)
)

# ── Save — TPG double-column 180 mm, 300 DPI ──────────────────────────────────
W <- 7.087; H <- 7.5
ggsave(file.path(OUT_DIR, "FigS2_Yr25_Manhattan.pdf"), fig_final,
       width = W, height = H, units = "in", device = "pdf")
ggsave(file.path(OUT_DIR, "FigS2_Yr25_Manhattan.png"), fig_final,
       width = W, height = H, units = "in", dpi = 300)
cat(sprintf("Saved FigS2_Yr25_Manhattan (%.3f x %.1f in, 300 dpi)\n", W, H))
