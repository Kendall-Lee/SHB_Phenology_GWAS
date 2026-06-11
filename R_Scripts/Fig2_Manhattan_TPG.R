#!/usr/bin/env Rscript
# Fig2_Manhattan_TPG.R
# The Plant Genome-ready version — minimum 8pt text at 178mm print width

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(GWASpoly)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS")
THRESH  <- 5.52

CHR_COLS <- c("#5D6D7E", "#AEB6BF")
SIG_COLS <- c(
  DTFlower    = "#8B0000",
  DTFruit     = "#1565C0",
  Flow2Fruit  = "#1A5276",
  FruitWeight = "#4A235A"
)

# ── Theme: base_size 13 so all derived text >= 8pt at 178mm ──────────────────
pub_theme <- theme_classic(base_size = 13) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black", face = "bold"),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
    plot.margin      = margin(4, 5, 6, 5)
  )

# ── Load scans ────────────────────────────────────────────────────────────────
cat("Loading scan objects...\n")
load(file.path(BASE, "LINEAR/RedoNo25/GWASpoly_scans.reBLUE.RData"))
scan_exc25 <- data.loco.scan
load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"))
scan_all <- data.loco.scan
cat("Loaded\n")

# ── Cumulative position map ───────────────────────────────────────────────────
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
  chr_max[, CumStart := CumStart + seq_len(.N) * 5e6]
  map <- merge(map, chr_max[, .(Chrom, CumStart)], by = "Chrom")
  map[, CumPos := Position + CumStart]
  chr_mids <- map[, .(Mid = mean(CumPos)), by = Chrom]
  list(map = map, chr_mids = chr_mids, chr_max = chr_max)
}

cp_exc25 <- build_cumpos(scan_exc25)
cp_all   <- build_cumpos(scan_all)

# ── Extract scores ────────────────────────────────────────────────────────────
extract_scores <- function(scan, trait, cumpos_map, additive_only = FALSE) {
  s <- scan@scores[[trait]]
  if (is.null(s)) return(NULL)
  markers <- rownames(s)
  markers <- markers[!grepl("[.][0-9]+$", markers)]
  s_sub <- s[markers, , drop = FALSE]
  best <- if (additive_only && "additive" %in% colnames(s_sub)) {
    as.numeric(s_sub[, "additive"])
  } else {
    apply(s_sub, 1, max, na.rm = TRUE)
  }
  dt <- data.table(Marker = markers, Score = best)
  merge(dt, cumpos_map$map[, .(Marker, Chrom, CumPos)], by = "Marker")
}

# ── Manhattan builder — enlarged text ─────────────────────────────────────────
make_manhattan <- function(dt, cumpos_info, title_txt, subtitle_txt = "",
                           sig_col, is_primary = FALSE,
                           y_max = NULL, show_xlab = FALSE,
                           show_ylab = TRUE) {

  chr_mids <- cumpos_info$chr_mids
  dt[, ChrIdx := as.integer(Chrom)]
  dt[, ChrCol := ifelse(ChrIdx %% 2 == 1, CHR_COLS[1], CHR_COLS[2])]
  dt[, Sig    := Score >= THRESH]

  if (is.null(y_max)) y_max <- max(dt$Score, na.rm = TRUE) * 1.08
  y_max <- max(y_max, THRESH * 1.3)

  xlabs <- gsub("^Chr\\.0?", "", as.character(chr_mids$Chrom))

  p <- ggplot(dt) +
    geom_point(data = dt[Sig == FALSE],
               aes(x = CumPos, y = Score, colour = ChrCol),
               size = 0.5, alpha = 0.55, stroke = 0) +
    geom_point(data = dt[Sig == TRUE],
               aes(x = CumPos, y = Score),
               colour = sig_col,
               size   = if (is_primary) 2.0 else 1.4,
               alpha  = 0.9, stroke = 0) +
    geom_hline(yintercept = THRESH, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    scale_colour_identity() +
    scale_x_continuous(breaks = chr_mids$Mid, labels = xlabs,
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = scales::pretty_breaks(n = if (is_primary) 4 else 3),
                       expand = c(0, 0)) +
    labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x        = if (show_xlab) "Chromosome" else NULL,
      y        = if (show_ylab) expression(bold(-log[10](italic(p)))) else NULL
    ) +
    pub_theme +
    theme(
      panel.border  = element_rect(
        color     = if (is_primary) sig_col else "grey55",
        fill      = NA,
        linewidth = if (is_primary) 1.2 else 0.5
      ),
      # ── text sizes: primary vs support ──
      axis.text.x  = element_text(size = if (is_primary) 10 else  8.5),
      axis.text.y  = element_text(size = if (is_primary) 10 else  8.5),
      axis.title.x = element_text(size = if (is_primary) 11 else  9.5, face = "bold"),
      axis.title.y = element_text(size = if (is_primary) 11 else  9.5, face = "bold"),
      plot.title   = element_text(
        size   = if (is_primary) 13 else 10,
        face   = "bold",
        colour = if (is_primary) sig_col else "grey20"
      ),
      plot.subtitle = element_text(size = if (is_primary) 9 else 8, colour = "grey45")
    )
  p
}

# ── Trait row builder ─────────────────────────────────────────────────────────
build_trait_row <- function(primary_trait, primary_scan, primary_cp,
                            support_traits, support_scan, support_cp,
                            sig_col, primary_title, show_xlab = FALSE,
                            include_yr25 = TRUE) {

  dt_primary <- extract_scores(primary_scan, primary_trait, primary_cp)
  if (is.null(dt_primary) || nrow(dt_primary) == 0) return(NULL)

  global_ymax <- max(dt_primary$Score, na.rm = TRUE) * 1.08

  p_main <- make_manhattan(dt_primary, primary_cp,
    title_txt    = primary_title,
    subtitle_txt = paste0("BLUE exc. yr.2025  |  n=626  |  −log₁₀p threshold = ", THRESH),
    sig_col      = sig_col, is_primary = TRUE,
    y_max        = global_ymax, show_xlab = show_xlab, show_ylab = TRUE)

  if (!include_yr25) support_traits <- support_traits[names(support_traits) != "yr.25"]

  support_plots <- lapply(names(support_traits), function(yr_key) {
    tr  <- support_traits[[yr_key]]
    use_add <- yr_key == "yr.25"
    dt_s <- extract_scores(support_scan, tr, support_cp, additive_only = use_add)
    if (is.null(dt_s) || nrow(dt_s) == 0) return(NULL)
    ymax_s <- max(dt_s$Score, na.rm = TRUE) * 1.08
    yr_lab <- c(yr.23 = "yr.2023", yr.24 = "yr.2024", yr.25 = "yr.2025")[yr_key]
    sub    <- if (use_add) "additive model only" else ""
    make_manhattan(dt_s, support_cp,
      title_txt = yr_lab, subtitle_txt = sub,
      sig_col = sig_col, is_primary = FALSE,
      y_max = ymax_s, show_xlab = show_xlab, show_ylab = FALSE)
  })
  support_plots <- Filter(Negate(is.null), support_plots)

  plot_grid(plotlist = c(list(p_main), support_plots),
            nrow = 1, rel_widths = c(2.2, rep(1, length(support_plots))),
            align = "h", axis = "tb")
}

# ── Build rows ────────────────────────────────────────────────────────────────
cat("Building rows...\n")
row_flower <- build_trait_row(
  primary_trait  = "DTFlower_BLUE",   primary_scan = scan_exc25, primary_cp = cp_exc25,
  support_traits = list(yr.23 = "DTFlower_yr.23", yr.24 = "DTFlower_yr.24",
                        yr.25 = "DTFlower_yr.25"),
  support_scan = scan_all, support_cp = cp_all,
  sig_col = SIG_COLS["DTFlower"], primary_title = "Days to 50% Flowering",
  show_xlab = FALSE, include_yr25 = TRUE)

row_fruit <- build_trait_row(
  primary_trait  = "DTFruit_BLUE",    primary_scan = scan_exc25, primary_cp = cp_exc25,
  support_traits = list(yr.23 = "DTFruit_yr.23", yr.24 = "DTFruit_yr.24",
                        yr.25 = "DTFruit_yr.25"),
  support_scan = scan_all, support_cp = cp_all,
  sig_col = SIG_COLS["DTFruit"], primary_title = "Days to 50% Ripe Fruit",
  show_xlab = FALSE, include_yr25 = TRUE)

row_f2f <- build_trait_row(
  primary_trait  = "Flow2Fruit_BLUE", primary_scan = scan_exc25, primary_cp = cp_exc25,
  support_traits = list(yr.23 = "Flow2Fruit_yr.23", yr.24 = "Flow2Fruit_yr.24",
                        yr.25 = "Flow2Fruit_yr.25"),
  support_scan = scan_all, support_cp = cp_all,
  sig_col = SIG_COLS["Flow2Fruit"], primary_title = "Flower-to-Fruit Interval",
  show_xlab = TRUE, include_yr25 = TRUE)

# ── Assemble panels with built-in labels (no separate strip to clip) ──────────
fig2 <- plot_grid(
  row_flower, row_fruit, row_f2f,
  ncol          = 1,
  rel_heights   = c(1, 1, 1),
  labels        = c("A", "B", "C"),
  label_size    = 18,
  label_fontface = "bold",
  hjust         = -0.15,
  vjust         =  1.1
)

# ── Legend panel ─────────────────────────────────────────────────────────────
legend_df <- data.frame(
  x     = c(1, 2, 3, 5, 6),
  y     = c(1, 1, 1, 1, 1),
  label = c("Days to 50% Flowering", "Days to 50% Ripe Fruit",
            "Flower-to-Fruit Interval", "Significant (above threshold)",
            "Non-significant"),
  col   = c(SIG_COLS["DTFlower"], SIG_COLS["DTFruit"], SIG_COLS["Flow2Fruit"],
            "black", CHR_COLS[1]),
  sz    = c(3, 3, 3, 3, 2),
  stringsAsFactors = FALSE
)

legend_panel <- ggplot() +
  # trait colour dots
  geom_point(data = legend_df[1:3,],
             aes(x = x, y = y), colour = legend_df$col[1:3], size = 4) +
  # sig / non-sig dots
  geom_point(aes(x = 5, y = 1), colour = "black",        size = 4) +
  geom_point(aes(x = 6, y = 1), colour = CHR_COLS[1],    size = 2.5, alpha = 0.6) +
  # labels below dots
  annotate("text", x = 1,   y = 0.55, label = "DTFlower",
           size = 3.5, colour = SIG_COLS["DTFlower"],   fontface = "bold", hjust = 0.5) +
  annotate("text", x = 2,   y = 0.55, label = "DTFruit",
           size = 3.5, colour = SIG_COLS["DTFruit"],    fontface = "bold", hjust = 0.5) +
  annotate("text", x = 3,   y = 0.55, label = "Flow2Fruit",
           size = 3.5, colour = SIG_COLS["Flow2Fruit"], fontface = "bold", hjust = 0.5) +
  annotate("text", x = 5,   y = 0.55, label = "Significant",
           size = 3.5, colour = "black", hjust = 0.5) +
  annotate("text", x = 6,   y = 0.55, label = "Non-significant",
           size = 3.5, colour = "grey40", hjust = 0.5) +
  annotate("text", x = 5.0, y = 0.25, label = "Threshold (−log₁₀p = 5.52)",
           size = 3.5, colour = "grey40", hjust = 0.5) +
  xlim(0.4, 6.8) + ylim(0.1, 1.4) +
  theme_void() +
  theme(plot.margin = margin(2, 5, 4, 5))

# ── Assemble with legend at bottom ───────────────────────────────────────────
fig2_final <- plot_grid(
  fig2,
  legend_panel,
  ncol = 1,
  rel_heights = c(1, 0.10)
)

# ── Save — 178mm wide (7.01 in), taller for legibility, 600 dpi ──────────────
cat("Saving...\n")
W <- 7.01   # 178 mm — The Plant Genome max 2-col width
H <- 13.0   # increased panel height so axis text is not clipped

ggsave(file.path(OUT_DIR, "Fig2_Manhattan_TPG.pdf"), fig2_final,
       width = W, height = H, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig2_Manhattan_TPG.png"), fig2_final,
       width = W, height = H, units = "in", dpi = 600, device = "png")

cat("Saved to", OUT_DIR, "\n")
