#!/usr/bin/env Rscript
# Fig2_Manhattan.R
# Publication Figure 2 — Linear GWAS Manhattan plots
# Layout: 4 trait rows
#   DTFlower / DTFruit / Flow2Fruit: primary BLUE_exc.25 (large) + yr.23 / yr.24 / yr.25 (small)
#   FruitWeight: single full-width BLUE panel (not confounded by yr.25)
# Output: Fig2_Manhattan.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(GWASpoly)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION")
THRESH  <- 5.52

# ── Colours ───────────────────────────────────────────────────────────────────
CHR_COLS   <- c("#5D6D7E", "#AEB6BF")   # alternating chr colours
SIG_COLS <- c(
  DTFlower   = "#8B0000",
  DTFruit    = "#1565C0",
  Flow2Fruit = "#1A5276",
  FruitWeight= "#4A235A"
)
BORDER_COLS <- c(
  primary = "#8B0000",
  support = "grey60"
)

pub_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black", face = "bold"),
    plot.title       = element_text(size = 9.5, face = "bold"),
    plot.subtitle    = element_text(size = 7.5, color = "grey40"),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
    plot.margin      = margin(3, 4, 2, 4)
  )

# ── Load scans ────────────────────────────────────────────────────────────────
cat("Loading scan objects...\n")
load(file.path(BASE, "LINEAR/RedoNo25/GWASpoly_scans.reBLUE.RData"))
scan_exc25 <- data.loco.scan

load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"))
scan_all <- data.loco.scan
cat("Loaded\n")

# ── Build cumulative position map ─────────────────────────────────────────────
build_cumpos <- function(scan) {
  map <- as.data.table(scan@map)
  names(map)[1:3] <- c("Marker","Chrom","Position")
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
  chr_max[, CumStart := CumStart + seq_len(.N) * 5e6]  # gap between chrs

  map <- merge(map, chr_max[, .(Chrom, CumStart)], by = "Chrom")
  map[, CumPos := Position + CumStart]

  chr_mids <- map[, .(Mid = mean(CumPos)), by = Chrom]
  list(map = map, chr_mids = chr_mids, chr_max = chr_max)
}

cp_exc25 <- build_cumpos(scan_exc25)
cp_all   <- build_cumpos(scan_all)

# ── Extract best score per marker per trait ───────────────────────────────────
# additive_only = TRUE: use only the additive column (for yr.25 support panels)
extract_scores <- function(scan, trait, cumpos_map, additive_only = FALSE) {
  s <- scan@scores[[trait]]
  if (is.null(s)) return(NULL)
  markers <- rownames(s)
  markers <- markers[!grepl("[.][0-9]+$", markers)]
  s_sub <- s[markers, , drop = FALSE]
  if (additive_only && "additive" %in% colnames(s_sub)) {
    best <- as.numeric(s_sub[, "additive"])
  } else {
    best <- apply(s_sub, 1, max, na.rm = TRUE)
  }
  dt <- data.table(Marker = markers, Score = best)
  merge(dt, cumpos_map$map[, .(Marker, Chrom, CumPos)], by = "Marker")
}

# ── Manhattan plot builder ────────────────────────────────────────────────────
make_manhattan <- function(dt, cumpos_info, title_txt, subtitle_txt = "",
                           sig_col, is_primary = FALSE,
                           y_max = NULL, show_xlab = FALSE,
                           show_ylab = TRUE) {

  chr_mids <- cumpos_info$chr_mids
  chr_lvls <- levels(dt$Chrom)
  n_chr    <- length(chr_lvls)

  dt[, ChrIdx := as.integer(Chrom)]
  dt[, ChrCol := ifelse(ChrIdx %% 2 == 1, CHR_COLS[1], CHR_COLS[2])]
  dt[, Sig    := Score >= THRESH]

  if (is.null(y_max)) y_max <- max(dt$Score, na.rm = TRUE) * 1.08
  y_max <- max(y_max, THRESH * 1.3)

  # x-axis labels: show chr number only
  xlabs <- gsub("Chr.", "", chr_mids$Chrom)
  xlabs <- gsub("^0", "", xlabs)

  p <- ggplot(dt) +
    # non-significant points
    geom_point(data = dt[Sig == FALSE],
               aes(x = CumPos, y = Score, colour = ChrCol),
               size = 0.45, alpha = 0.55, stroke = 0) +
    # significant points
    geom_point(data = dt[Sig == TRUE],
               aes(x = CumPos, y = Score),
               colour = sig_col, size = if (is_primary) 1.6 else 1.1,
               alpha = 0.9, stroke = 0) +
    # threshold
    geom_hline(yintercept = THRESH, linetype = "dashed",
               colour = "grey50", linewidth = 0.55) +
    scale_colour_identity() +
    scale_x_continuous(
      breaks = chr_mids$Mid,
      labels = xlabs,
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = scales::pretty_breaks(n = if (is_primary) 4 else 3),
      expand = c(0, 0)
    ) +
    labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x        = if (show_xlab) "Chromosome" else NULL,
      y        = if (show_ylab) expression(bold(-log[10](italic(p)))) else NULL
    ) +
    pub_theme +
    theme(
      panel.border   = element_rect(
        color     = if (is_primary) sig_col else "grey55",
        fill      = NA,
        linewidth = if (is_primary) 1.1 else 0.5
      ),
      axis.text.x    = element_text(size = if (is_primary) 8 else 6.5),
      axis.text.y    = element_text(size = if (is_primary) 8 else 6.5),
      axis.title.x   = element_text(size = 8),
      axis.title.y   = element_text(size = if (is_primary) 9 else 7.5),
      axis.ticks.x   = element_line(linewidth = 0.3),
      plot.title     = element_text(
        size   = if (is_primary) 10 else 8,
        face   = "bold",
        colour = if (is_primary) sig_col else "grey20"
      ),
      plot.subtitle  = element_text(size = 7, colour = "grey45")
    )
  p
}

# ── Helper: build one trait row ───────────────────────────────────────────────
# include_yr25: set FALSE to drop yr.25 panel entirely
build_trait_row <- function(trait_base, primary_trait, primary_scan,
                            primary_cp, support_traits, support_scan,
                            support_cp, sig_col,
                            primary_title, show_xlab = FALSE,
                            include_yr25 = FALSE) {

  dt_primary <- extract_scores(primary_scan, primary_trait, primary_cp)
  if (is.null(dt_primary) || nrow(dt_primary) == 0) {
    cat("  WARNING: no data for", primary_trait, "\n"); return(NULL)
  }

  global_ymax <- max(dt_primary$Score, na.rm = TRUE) * 1.08

  p_main <- make_manhattan(dt_primary, primary_cp,
    title_txt    = primary_title,
    subtitle_txt = paste0("BLUE exc. yr.2025  |  n=626  |  threshold: -log10p=", THRESH),
    sig_col      = sig_col,
    is_primary   = TRUE,
    y_max        = global_ymax,
    show_xlab    = show_xlab,
    show_ylab    = TRUE)

  yr_labels <- c(yr.23 = "2023", yr.24 = "2024", yr.25 = "2025")

  # Drop yr.25 if not requested
  if (!include_yr25) support_traits <- support_traits[names(support_traits) != "yr.25"]

  support_plots <- lapply(names(support_traits), function(yr_key) {
    tr  <- support_traits[[yr_key]]
    use_additive <- yr_key == "yr.25"
    dt_s <- extract_scores(support_scan, tr, support_cp, additive_only = use_additive)
    if (is.null(dt_s) || nrow(dt_s) == 0) return(NULL)
    ymax_s <- max(dt_s$Score, na.rm = TRUE) * 1.08
    sub <- if (use_additive) "additive model only" else ""
    make_manhattan(dt_s, support_cp,
      title_txt    = paste0("yr.", yr_labels[yr_key]),
      subtitle_txt = sub,
      sig_col      = sig_col,
      is_primary   = FALSE,
      y_max        = ymax_s,
      show_xlab    = show_xlab,
      show_ylab    = FALSE)
  })
  support_plots <- Filter(Negate(is.null), support_plots)

  plot_grid(plotlist = c(list(p_main), support_plots),
            nrow = 1,
            rel_widths = c(2.2, rep(1, length(support_plots))),
            align = "h", axis = "tb")
}

# ── Build all rows ────────────────────────────────────────────────────────────
cat("Building DTFlower row...\n")
row_flower <- build_trait_row(
  trait_base     = "DTFlower",
  primary_trait  = "DTFlower_BLUE",
  primary_scan   = scan_exc25,
  primary_cp     = cp_exc25,
  support_traits = list(yr.23 = "DTFlower_yr.23",
                        yr.24 = "DTFlower_yr.24",
                        yr.25 = "DTFlower_yr.25"),
  support_scan   = scan_all,
  support_cp     = cp_all,
  sig_col        = SIG_COLS["DTFlower"],
  primary_title  = "Days to 50% Flowering",
  show_xlab      = FALSE,
  include_yr25   = TRUE
)

cat("Building DTFruit row...\n")
row_fruit <- build_trait_row(
  trait_base     = "DTFruit",
  primary_trait  = "DTFruit_BLUE",
  primary_scan   = scan_exc25,
  primary_cp     = cp_exc25,
  support_traits = list(yr.23 = "DTFruit_yr.23",
                        yr.24 = "DTFruit_yr.24",
                        yr.25 = "DTFruit_yr.25"),
  support_scan   = scan_all,
  support_cp     = cp_all,
  sig_col        = SIG_COLS["DTFruit"],
  primary_title  = "Days to 50% Ripe Fruit",
  show_xlab      = FALSE,
  include_yr25   = TRUE
)

cat("Building Flow2Fruit row...\n")
row_f2f <- build_trait_row(
  trait_base     = "Flow2Fruit",
  primary_trait  = "Flow2Fruit_BLUE",
  primary_scan   = scan_exc25,
  primary_cp     = cp_exc25,
  support_traits = list(yr.23 = "Flow2Fruit_yr.23",
                        yr.24 = "Flow2Fruit_yr.24",
                        yr.25 = "Flow2Fruit_yr.25"),
  support_scan   = scan_all,
  support_cp     = cp_all,
  sig_col        = SIG_COLS["Flow2Fruit"],
  primary_title  = "Flower-to-Fruit Interval",
  show_xlab      = TRUE,
  include_yr25   = TRUE
)

# ── Panel labels ──────────────────────────────────────────────────────────────
add_label <- function(row, lbl) {
  plot_grid(
    ggdraw() + draw_label(lbl, fontface = "bold", size = 16,
                          x = 0.005, hjust = 0, vjust = 0.5),
    row,
    ncol = 1, rel_heights = c(0.08, 1)
  )
}

fig2 <- plot_grid(
  add_label(row_flower, "A"),
  add_label(row_fruit,  "B"),
  add_label(row_f2f,    "C"),
  ncol = 1,
  rel_heights = c(1, 1, 1)
)

# ── Save ──────────────────────────────────────────────────────────────────────
cat("Saving...\n")
ggsave(file.path(OUT_DIR, "Fig2_Manhattan_withYr25.pdf"), fig2,
       width = 7.09, height = 6.5, units = "in", device = "pdf")
ggsave(file.path(OUT_DIR, "Fig2_Manhattan_withYr25.png"), fig2,
       width = 7.09, height = 6.5, units = "in", dpi = 300, device = "png")

cat("\nFig2 saved to", OUT_DIR, "\n")
