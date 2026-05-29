#!/usr/bin/env Rscript
# Fig3_Chr05_Hotspot.R
# Publication Figure 3 — Chr.05:47–48 Mb regional association + allele effects + QTL-seq
# Linear GWAS only (BLUE_exc.25 primary + yr.23/yr.24 support)
# Output: Fig3_Chr05_Hotspot.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(GWASpoly)
  library(grid)
  library(ggtext)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION")
THRESH  <- 5.52

# ── Colours ───────────────────────────────────────────────────────────────────
COL_PRIMARY  <- "#8B0000"   # dark red  — BLUE_exc.25
COL_YR23     <- "#1A6B3C"   # dark green — yr.23
COL_YR24     <- "#5B2C8D"   # dark purple — yr.24
COL_DTFRUIT  <- "#1565C0"   # dark blue — DTFruit (secondary trait)
COL_QTLSEQ   <- "#B45309"   # amber — QTL-seq
COL_THRESH   <- "#888888"
GENE_FILL    <- "#2C3E50"

pub_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text        = element_text(color = "black", size = 10),
    axis.title       = element_text(color = "black", size = 11, face = "bold"),
    plot.title       = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

cat("Loading scan objects...\n")

# ── 1. Load GWAS scans ────────────────────────────────────────────────────────
load(file.path(BASE, "LINEAR/RedoNo25/GWASpoly_scans.reBLUE.RData"))
scan_exc25 <- data.loco.scan

load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"))
scan_all <- data.loco.scan

cat("Scans loaded\n")

# ── 2. Extract Chr.05 46–49.5 Mb scores ──────────────────────────────────────
WIN_START <- 46e6; WIN_END <- 49.5e6

extract_chr05 <- function(scan, traits, label_suffix = "") {
  map <- as.data.table(scan@map)
  # Standardise column names regardless of exact count
  names(map)[1:3] <- c("Marker","Chrom","Position")
  map[, Position := as.numeric(Position)]
  chr5_markers <- map[Chrom == "Chr.05" & Position >= WIN_START & Position <= WIN_END, Marker]
  chr5_markers <- chr5_markers[!grepl("[.][0-9]+$", chr5_markers)]  # drop dup artifacts

  rbindlist(lapply(traits, function(tr) {
    s <- scan@scores[[tr]]
    if (is.null(s)) return(NULL)
    keep <- intersect(chr5_markers, rownames(s))
    if (length(keep) == 0) return(NULL)
    sub <- as.data.table(s[keep, , drop = FALSE], keep.rownames = "Marker")
    sub[, BestScore := apply(.SD, 1, max, na.rm = TRUE), .SDcols = colnames(s)]
    pos_dt <- map[Marker %in% keep, .(Marker, Position)]
    sub <- merge(sub, pos_dt, by = "Marker")
    sub[, .(Marker, Position, Score = BestScore, Trait = tr)]
  }))
}

# BLUE_exc.25 traits
exc25_scores <- extract_chr05(scan_exc25,
  c("DTFlower_BLUE", "DTFruit_BLUE", "Flow2Fruit_BLUE"))

# yr.23 and yr.24 individual scans for support
yr_scores <- extract_chr05(scan_all,
  c("DTFlower_yr.23", "DTFlower_yr.24",
    "DTFruit_yr.23",  "DTFruit_yr.24"))

gwas_dt <- rbind(exc25_scores, yr_scores)
gwas_dt[, Pos_Mb := Position / 1e6]

# Clean display labels and assign colours + shapes
gwas_dt[, TraitLabel := fcase(
  Trait == "DTFlower_BLUE",     "DTFlower BLUE exc.25",
  Trait == "DTFruit_BLUE",      "DTFruit BLUE exc.25",
  Trait == "Flow2Fruit_BLUE",   "Flow2Fruit BLUE exc.25",
  Trait == "DTFlower_yr.23",    "DTFlower yr.23",
  Trait == "DTFlower_yr.24",    "DTFlower yr.24",
  Trait == "DTFruit_yr.23",     "DTFruit yr.23",
  Trait == "DTFruit_yr.24",     "DTFruit yr.24",
  default = Trait
)]

# Order for legend
label_order <- c("DTFlower BLUE exc.25","DTFruit BLUE exc.25","Flow2Fruit BLUE exc.25",
                  "DTFlower yr.23","DTFlower yr.24","DTFruit yr.23","DTFruit yr.24")
gwas_dt[, TraitLabel := factor(TraitLabel, levels = label_order)]

trait_colours <- c(
  "DTFlower BLUE exc.25"  = COL_PRIMARY,
  "DTFruit BLUE exc.25"   = COL_DTFRUIT,
  "Flow2Fruit BLUE exc.25"= "#2E86C1",
  "DTFlower yr.23"        = COL_YR23,
  "DTFlower yr.24"        = COL_YR24,
  "DTFruit yr.23"         = "#5DADE2",
  "DTFruit yr.24"         = "#A569BD"
)
trait_alpha <- c(
  "DTFlower BLUE exc.25"  = 1,
  "DTFruit BLUE exc.25"   = 1,
  "Flow2Fruit BLUE exc.25"= 0.85,
  "DTFlower yr.23"        = 0.65,
  "DTFlower yr.24"        = 0.65,
  "DTFruit yr.23"         = 0.55,
  "DTFruit yr.24"         = 0.55
)
trait_size <- c(
  "DTFlower BLUE exc.25"  = 2.2,
  "DTFruit BLUE exc.25"   = 2.2,
  "Flow2Fruit BLUE exc.25"= 1.8,
  "DTFlower yr.23"        = 1.5,
  "DTFlower yr.24"        = 1.5,
  "DTFruit yr.23"         = 1.4,
  "DTFruit yr.24"         = 1.4
)

cat(sprintf("GWAS points extracted: %d\n", nrow(gwas_dt)))

# ── 3. Gene track (46–49.5 Mb) ────────────────────────────────────────────────
genes_bed <- fread("/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/All_Markers/genes.bed",
                   col.names = c("Chrom","Start","End","Gene"), header = FALSE)
genes_chr5 <- genes_bed[Chrom == "Chr.05" & Start >= WIN_START & End <= WIN_END]
genes_chr5[, Start_Mb := Start / 1e6]
genes_chr5[, End_Mb   := End   / 1e6]
genes_chr5[, Mid_Mb   := (Start_Mb + End_Mb) / 2]

# Candidate genes to label
candidates <- data.table(
  Gene  = c("g54686", "g54697"),
  Label = c("g54686\n(PR-10/MLP)", "g54697\n(PR-10/MLP)"),
  Pos_Mb = c(47.711, 47.899)
)

cat(sprintf("Genes in region: %d\n", nrow(genes_chr5)))

# ── 4. QTL-seq Chr.05 data ────────────────────────────────────────────────────
cat("Loading QTL-seq data...\n")
QTLSEQ_FILE <- paste0(
  "/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/",
  "LRLP_QTLVar/DIY_QTLSeq/PopStructureCorrection/",
  "results_Pangenome_PCAcorr/DTFlower/All_markers_stats_DTFlower.txt.gz")

qtlseq_raw <- fread(QTLSEQ_FILE,
                    select = c("Chrom","Position","R2","R2_smooth","PVE",
                               "padj","candidate_QTL"))
qtlseq_chr5 <- qtlseq_raw[Chrom == "Chr.05" &
                            Position >= WIN_START & Position <= WIN_END]
qtlseq_chr5[, Pos_Mb := Position / 1e6]
qtlseq_chr5[, Sig := padj < 0.05]

cat(sprintf("QTL-seq Chr.05 markers: %d\n", nrow(qtlseq_chr5)))

# ── 5. Allele effect data for Panel B ─────────────────────────────────────────
cat("Loading geno + phenotype data for allele effects...\n")

load(file.path(BASE, "LINEAR/RedoNo25/data_loco.tet.reBLUE.RData"))
geno_mat   <- data.loco@geno
sample_ids <- rownames(geno_mat)
rm(data.loco); gc()

TOP_MARKER <- "Chr.05_48252083"
dos_vec <- as.numeric(geno_mat[, TOP_MARKER])
dos_dt  <- data.table(Sample = sample_ids, Dosage = dos_vec)
dos_dt[, DosClass := factor(floor(pmin(Dosage + 0.5, 4)), levels = 0:4)]
rm(geno_mat); gc()

pheno_exc25 <- fread(file.path(BASE, "Phenotypes/BLUE_exc.25/DTFlower_2324only_BLUE.txt"))
setnames(pheno_exc25, 1:2, c("Sample","Pheno"))

pheno_all <- fread(file.path(BASE, "Phenotypes/SHB_BLUE_all_pheno.csv"))
setnames(pheno_all, "DNA ID", "Sample")

make_dose_df <- function(pheno_dt, pheno_col = "Pheno") {
  df <- merge(dos_dt, pheno_dt[, .(Sample, Pheno = get(pheno_col))], by = "Sample")
  df[!is.na(Dosage) & !is.na(Pheno)]
}

df_exc <- make_dose_df(pheno_exc25)
df_23  <- make_dose_df(pheno_all, "DTFlower_yr.23")
df_24  <- make_dose_df(pheno_all, "DTFlower_yr.24")

# Tukey HSD CLD
get_cld <- function(df) {
  tryCatch({
    df2 <- copy(df); df2[, DosChar := as.character(DosClass)]
    grp_ok <- df2[, .N, by = DosChar][N >= 2]
    if (nrow(grp_ok) < 2) return(NULL)
    df_sub <- df2[DosChar %in% grp_ok$DosChar]
    df_sub[, DosChar := factor(DosChar)]
    mod <- aov(Pheno ~ DosChar, data = df_sub)
    tukey <- TukeyHSD(mod, "DosChar")$DosChar
    pvals <- tukey[, "p adj"]
    groups <- levels(df_sub$DosChar); n <- length(groups)
    adj <- matrix(TRUE, n, n, dimnames = list(groups, groups))
    for (nm in names(pvals)) {
      pts <- strsplit(nm, "-")[[1]]
      if (length(pts) == 2 && all(pts %in% groups)) {
        adj[pts[1], pts[2]] <- pvals[nm] >= 0.05
        adj[pts[2], pts[1]] <- pvals[nm] >= 0.05
      }
    }
    letter_sets <- setNames(vector("list", n), groups); cur_letter <- 1L
    for (g in groups) {
      added <- FALSE
      for (let in unique(unlist(letter_sets))) {
        members <- groups[sapply(groups, function(x) let %in% letter_sets[[x]])]
        if (all(adj[g, members])) {
          letter_sets[[g]] <- c(letter_sets[[g]], let); added <- TRUE
        }
      }
      if (!added || length(letter_sets[[g]]) == 0) {
        letter_sets[[g]] <- c(letter_sets[[g]], letters[cur_letter])
        cur_letter <- cur_letter + 1L
      }
    }
    cld_letters <- sapply(groups, function(g) paste(sort(unique(letter_sets[[g]])), collapse = ""))
    grp_means <- df_sub[, .(emmean = mean(Pheno, na.rm = TRUE)), by = DosChar]
    w_tops <- df_sub[, .(w_top = boxplot.stats(Pheno)$stats[5],
                         DosChar = as.character(DosClass[1])), by = DosClass][
                       , .(DosChar = as.character(DosClass), w_top)]
    cld_dt <- data.table(DosClass = groups, cld_letter = cld_letters)
    cld_dt <- merge(cld_dt, w_tops, by.x = "DosClass", by.y = "DosChar", all.x = TRUE)
    y_rng <- diff(range(df$Pheno, na.rm = TRUE))
    cld_dt[, label_y := w_top + y_rng * 0.08]
    cld_dt
  }, error = function(e) NULL)
}

make_box_panel <- function(df, panel_label, fill_col, score = NA, pve = NA) {
  counts <- df[, .N, by = DosClass][order(DosClass)]
  y_min <- min(df$Pheno, na.rm = TRUE); y_max <- max(df$Pheno, na.rm = TRUE)
  y_rng <- y_max - y_min
  y_bot  <- y_min - y_rng * 0.18
  cld_dt <- get_cld(df)
  y_top <- if (!is.null(cld_dt)) max(cld_dt$label_y, na.rm = TRUE) + y_rng * 0.05 else y_max + y_rng * 0.1

  sub_txt <- paste(Filter(nchar, c(
    if (!is.na(score)) sprintf("-log10p = %.2f", score),
    if (!is.na(pve) && is.finite(pve)) sprintf("PVE = %.1f%%", pve)
  )), collapse = "  |  ")

  p <- ggplot(df, aes(x = DosClass, y = Pheno)) +
    geom_boxplot(fill = fill_col, colour = "#1a1a1a", alpha = 0.75,
                 outlier.size = 0.6, outlier.alpha = 0.35, width = 0.52) +
    geom_jitter(width = 0.16, alpha = 0.18, size = 0.55, colour = "#1a1a1a") +
    stat_summary(fun = mean, geom = "point", shape = 23, fill = "white",
                 colour = fill_col, size = 2.5, stroke = 1.2) +
    geom_text(data = counts,
              aes(x = DosClass, y = y_bot + y_rng * 0.02,
                  label = paste0("n=", N)),
              size = 2.6, colour = "grey40", inherit.aes = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(y_bot, y_top)) +
    labs(title = panel_label,
         subtitle = sub_txt,
         x = "Dosage (ALT copies)", y = "Days to 50% Flowering") +
    pub_theme +
    theme(
      plot.title    = element_text(size = 9.5, face = "bold", colour = fill_col),
      plot.subtitle = element_text(size = 7.5, colour = "grey40"),
      axis.title.y  = element_text(size = 9),
      panel.border  = element_rect(color = fill_col, fill = NA,
                                   linewidth = if (fill_col == COL_PRIMARY) 1.4 else 0.8)
    )
  if (!is.null(cld_dt) && nrow(cld_dt) > 0)
    p <- p + geom_text(data = cld_dt,
                       aes(x = DosClass, y = label_y, label = cld_letter),
                       size = 3.6, fontface = "bold", colour = fill_col,
                       inherit.aes = FALSE)
  p
}

# Scores for Chr.05_48252083
get_score <- function(scan, trait) {
  s <- scan@scores[[trait]]
  if (is.null(s) || !TOP_MARKER %in% rownames(s)) return(NA)
  max(as.numeric(s[TOP_MARKER, , drop = FALSE]), na.rm = TRUE)
}

score_exc <- get_score(scan_exc25, "DTFlower_BLUE")
score_23  <- get_score(scan_all,   "DTFlower_yr.23")
score_24  <- get_score(scan_all,   "DTFlower_yr.24")
pve_exc   <- 7.79  # from TopMarkers table

cat(sprintf("Scores — BLUE_exc25: %.2f  yr23: %.2f  yr24: %.2f\n",
            score_exc, score_23, score_24))

# ── 6. Build panels ───────────────────────────────────────────────────────────

cat("Building Panel A — regional association plot...\n")

# ---- Panel A1: GWAS scores ----
pA1 <- ggplot(gwas_dt[Score > 0],
              aes(x = Pos_Mb, y = Score,
                  colour = TraitLabel, size = TraitLabel, alpha = TraitLabel)) +
  geom_point(stroke = 0) +
  geom_hline(yintercept = THRESH, linetype = "dashed",
             colour = COL_THRESH, linewidth = 0.7) +
  annotate("text", x = 46.15, y = THRESH + 0.25, label = "Meff threshold",
           size = 2.8, colour = COL_THRESH, hjust = 0) +
  # Shade hotspot
  annotate("rect", xmin = 47.5, xmax = 48.5,
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "#8B0000") +
  scale_colour_manual(values = trait_colours, name = NULL) +
  scale_size_manual(values = trait_size, name = NULL) +
  scale_alpha_manual(values = trait_alpha, name = NULL) +
  scale_x_continuous(limits = c(46, 49.5), breaks = seq(46, 49.5, 0.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 14.5), breaks = seq(0, 14, 2)) +
  labs(x = NULL, y = expression(bold(-log[10](italic(p))))) +
  pub_theme +
  theme(
    legend.position   = c(0.17, 0.80),
    legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.4),
    legend.key.size   = unit(0.4, "cm"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.title.x      = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)),
         size = "none", alpha = "none")

# ---- Panel A2: gene track ----
pA2 <- ggplot(genes_chr5) +
  geom_rect(aes(xmin = Start_Mb, xmax = End_Mb, ymin = 0, ymax = 1),
            fill = "#BDC3C7", colour = NA, alpha = 0.7) +
  # Highlight candidates — wider bar + dot above for visibility
  geom_rect(data = genes_chr5[Gene %in% c("g54686","g54697")],
            aes(xmin = Start_Mb, xmax = End_Mb, ymin = 0, ymax = 1),
            fill = GENE_FILL, colour = GENE_FILL, linewidth = 0.4) +
  geom_point(data = genes_chr5[Gene %in% c("g54686","g54697")],
             aes(x = Mid_Mb, y = 1.0),
             shape = 25, fill = GENE_FILL, colour = GENE_FILL, size = 2.5,
             inherit.aes = FALSE) +
  annotate("rect", xmin = 47.5, xmax = 48.5,
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "#8B0000") +
  scale_x_continuous(limits = c(46, 49.5), breaks = seq(46, 49.5, 0.5),
                     expand = c(0, 0), name = "Chr.05 Position (Mb)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Chr.05 Position (Mb)", y = NULL) +
  pub_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x  = element_text(size = 9),
    plot.margin  = margin(0, 5.5, 2, 5.5)
  )

cat("Building Panel B — allele effects...\n")

pB1 <- make_box_panel(df_exc,
  panel_label = "BLUE exc. yr.25 (primary)",
  fill_col = COL_PRIMARY, score = score_exc, pve = pve_exc)

pB2 <- make_box_panel(df_23,
  panel_label = "yr.23 (support)",
  fill_col = COL_YR23, score = score_23)

pB3 <- make_box_panel(df_24,
  panel_label = "yr.24 (support)",
  fill_col = COL_YR24, score = score_24)

pB1 <- pB1 + labs(y = "Days to 50% Flowering")
pB2 <- pB2 + labs(y = NULL)
pB3 <- pB3 + labs(y = NULL)

cat("Building Panel C — QTL-seq...\n")

# ---- Panel C: QTL-seq R² ----
pC <- ggplot(qtlseq_chr5, aes(x = Pos_Mb)) +
  annotate("rect", xmin = 47.5, xmax = 48.5,
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "#8B0000") +
  geom_point(aes(y = R2 * 100), colour = "grey75", size = 0.5, alpha = 0.5) +
  geom_line(aes(y = R2_smooth * 100), colour = COL_QTLSEQ,
            linewidth = 1.1, alpha = 0.9) +
  geom_point(data = qtlseq_chr5[Sig == TRUE],
             aes(y = R2 * 100), colour = COL_QTLSEQ,
             size = 2.0, shape = 21, fill = "white", stroke = 1.2) +
  # Mark QTL-seq peaks at 47.7 and 47.9
  annotate("segment", x = 47.712, xend = 47.712,
           y = 29, yend = 27,
           arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
           colour = COL_QTLSEQ, linewidth = 0.7) +
  annotate("text", x = 47.3, y = 32,
           label = "R²=26.6%", size = 2.8, colour = COL_QTLSEQ, hjust = 0.5) +
  annotate("segment", x = 47.3, xend = 47.67,
           y = 31.2, yend = 29.2,
           colour = COL_QTLSEQ, linewidth = 0.4, linetype = "dotted") +
  annotate("segment", x = 47.900, xend = 47.900,
           y = 29, yend = 27,
           arrow = arrow(length = unit(0.11, "cm"), type = "closed"),
           colour = COL_QTLSEQ, linewidth = 0.7) +
  annotate("text", x = 48.3, y = 32,
           label = "R²=27.3%", size = 2.8, colour = COL_QTLSEQ, hjust = 0.5) +
  annotate("segment", x = 48.3, xend = 47.95,
           y = 31.2, yend = 29.2,
           colour = COL_QTLSEQ, linewidth = 0.4, linetype = "dotted") +
  scale_x_continuous(limits = c(46, 49.5), breaks = seq(46, 49.5, 0.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 35), breaks = seq(0, 35, 5)) +
  labs(x = "Chr.05 Position (Mb)", y = "QTL-seq R² (%)") +
  pub_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9)
  )

# ── 7. Assemble figure ────────────────────────────────────────────────────────
cat("Assembling figure...\n")

# Combine regional plot (GWAS scores + gene track)
pA_combined <- plot_grid(pA1, pA2,
  ncol = 1, rel_heights = c(5.5, 1),
  align = "v", axis = "lr")

# Allele effects row
pB_row <- plot_grid(pB1, pB2, pB3,
  nrow = 1, rel_widths = c(1.08, 0.96, 0.96),
  labels = c("", "", ""), align = "h", axis = "tb")

# QTL-seq
pC_full <- pC

# Panel labels
pA_labelled <- plot_grid(
  ggdraw() + draw_label("A", fontface = "bold", size = 14, x = 0.01, hjust = 0),
  pA_combined,
  ncol = 1, rel_heights = c(0.05, 1)
)
pBC_label <- ggdraw() + draw_label("B", fontface = "bold", size = 14, x = 0.01, hjust = 0)
pC_label  <- ggdraw() + draw_label("C", fontface = "bold", size = 14, x = 0.01, hjust = 0)

pB_labelled <- plot_grid(pBC_label, pB_row, ncol = 1, rel_heights = c(0.07, 1))
pC_labelled <- plot_grid(pC_label,  pC_full, ncol = 1, rel_heights = c(0.07, 1))

bottom_row <- plot_grid(pB_labelled, pC_labelled,
  nrow = 1, rel_widths = c(1.3, 0.85))

fig3 <- plot_grid(pA_labelled, bottom_row,
  ncol = 1, rel_heights = c(1, 1.1))

# ── 8. Save ───────────────────────────────────────────────────────────────────
cat("Saving...\n")

pdf_path <- file.path(OUT_DIR, "Fig3_Chr05_Hotspot.pdf")
png_path <- file.path(OUT_DIR, "Fig3_Chr05_Hotspot.png")

ggsave(pdf_path, fig3, width = 7.09, height = 6.0, units = "in", device = "pdf")
ggsave(png_path, fig3, width = 7.09, height = 6.0, units = "in", dpi = 300, device = "png")

cat(sprintf("\nFig3 saved:\n  %s\n  %s\n", pdf_path, png_path))
