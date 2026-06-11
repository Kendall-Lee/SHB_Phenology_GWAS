#!/usr/bin/env Rscript
# FigS_Allele_Effects_Secondary.R
# Allele dosage effect plots for all secondary association loci
# One row per locus: BLUE_exc.25 (primary) + yr.23 + yr.24 (support)
# Sorted: Table S4 locus first (Chr.06:22-23Mb), then Table S5 loci by chromosome

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(GWASpoly)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS")
THRESH  <- 5.52

# ── Locus list: (locus_label, marker, primary_trait, score_BLUE, PVE) ─────────
locus_list <- list(
  # Table S4 — independently replicated
  list(label="Chr.06:22–23 Mb†", marker="Chr.06_22292639",  trait="DTFlower", score=6.988, pve=4.47),
  # Table S5 — candidate loci
  list(label="Chr.02:3–4 Mb",    marker="Chr.02_3615748",   trait="DTFlower", score=8.207, pve=5.39),
  list(label="Chr.02:37–38 Mb",  marker="Chr.02_37452394",  trait="DTFlower", score=7.718, pve=4.95),
  list(label="Chr.02:47–48 Mb",  marker="Chr.02_47854803",  trait="DTFlower", score=6.100, pve=3.85),
  list(label="Chr.03:16–17 Mb",  marker="Chr.03_16371025",  trait="DTFruit",  score=7.721, pve=4.94),
  list(label="Chr.05:20–21 Mb",  marker="Chr.05_20189714",  trait="DTFruit",  score=7.537, pve=4.82),
  list(label="Chr.05:44–45 Mb",  marker="Chr.05_44735183",  trait="DTFruit",  score=7.359, pve=4.69),
  list(label="Chr.05:45–46 Mb",  marker="Chr.05_45296330",  trait="DTFruit",  score=8.370, pve=5.40),
  list(label="Chr.05:46–47 Mb",  marker="Chr.05_46175030",  trait="DTFlower", score=6.874, pve=4.37),
  list(label="Chr.06:29–30 Mb",  marker="Chr.06_29469225",  trait="DTFlower", score=6.533, pve=4.19),
  list(label="Chr.06:42–43 Mb",  marker="Chr.06_42369200",  trait="DTFlower", score=6.778, pve=4.32),
  list(label="Chr.08:18–19 Mb",  marker="Chr.08_18849701",  trait="DTFlower", score=6.591, pve=4.22),
  list(label="Chr.08:27–28 Mb",  marker="Chr.08_27127327",  trait="DTFlower", score=6.157, pve=3.89),
  list(label="Chr.08:28–29 Mb",  marker="Chr.08_28429653",  trait="DTFlower", score=6.186, pve=4.00),
  list(label="Chr.08:29–30 Mb",  marker="Chr.08_29470828",  trait="DTFlower", score=6.783, pve=4.36),
  list(label="Chr.09:1–2 Mb",    marker="Chr.09_1693527",   trait="DTFlower", score=6.205, pve=3.93),
  list(label="Chr.11:19–20 Mb",  marker="Chr.11_19498635",  trait="DTFlower", score=5.852, pve=3.74)
)

# ── Load scans and genotypes ───────────────────────────────────────────────────
cat("Loading GWAS scans...\n")
load(file.path(BASE, "LINEAR/RedoNo25/GWASpoly_scans.reBLUE.RData"))
scan_exc25 <- data.loco.scan
load(file.path(BASE, "LINEAR/ALL_TRAITS/GWASpoly_scans.RData"))
scan_all <- data.loco.scan

cat("Loading genotype matrix...\n")
load(file.path(BASE, "LINEAR/RedoNo25/data_loco.tet.reBLUE.RData"))
geno_mat   <- data.loco@geno
sample_ids <- rownames(geno_mat)
rm(data.loco); gc()

# ── Load phenotypes ────────────────────────────────────────────────────────────
cat("Loading phenotypes...\n")
pheno_fl_exc <- fread(file.path(BASE, "Phenotypes/BLUE_exc.25/DTFlower_2324only_BLUE.txt"))
pheno_fr_exc <- fread(file.path(BASE, "Phenotypes/BLUE_exc.25/DTFruit_2324only_BLUE.txt"))
setnames(pheno_fl_exc, 1:2, c("Sample","Pheno"))
setnames(pheno_fr_exc, 1:2, c("Sample","Pheno"))

pheno_all <- fread(file.path(BASE, "Phenotypes/SHB_BLUE_all_pheno.csv"))
setnames(pheno_all, "DNA ID", "Sample")

# ── Helpers ────────────────────────────────────────────────────────────────────
get_score <- function(scan, marker, trait) {
  s <- scan@scores[[trait]]
  if (is.null(s) || !marker %in% rownames(s)) return(NA_real_)
  round(max(as.numeric(s[marker, , drop=FALSE]), na.rm=TRUE), 2)
}

get_dosage <- function(marker) {
  if (!marker %in% colnames(geno_mat)) return(NULL)
  dos <- as.numeric(geno_mat[, marker])
  data.table(Sample = sample_ids, Dosage = dos,
             DosClass = factor(floor(pmin(dos + 0.5, 4)), levels=0:4))
}

get_cld <- function(df) {
  tryCatch({
    df2 <- copy(df); df2[, DosChar := as.character(DosClass)]
    grp_ok <- df2[, .N, by=DosChar][N >= 2]
    if (nrow(grp_ok) < 2) return(NULL)
    df_sub <- df2[DosChar %in% grp_ok$DosChar]
    df_sub[, DosChar := factor(DosChar)]
    mod    <- aov(Pheno ~ DosChar, data=df_sub)
    tukey  <- TukeyHSD(mod, "DosChar")$DosChar
    groups <- levels(df_sub$DosChar); n <- length(groups)
    adj    <- matrix(TRUE, n, n, dimnames=list(groups,groups))
    for (nm in names(tukey[,"p adj"])) {
      pts <- strsplit(nm,"-")[[1]]
      if (length(pts)==2 && all(pts %in% groups)) {
        adj[pts[1],pts[2]] <- tukey[nm,"p adj"] >= 0.05
        adj[pts[2],pts[1]] <- tukey[nm,"p adj"] >= 0.05
      }
    }
    letter_sets <- setNames(vector("list",n), groups); cur <- 1L
    for (g in groups) {
      added <- FALSE
      for (let in unique(unlist(letter_sets))) {
        members <- groups[sapply(groups, function(x) let %in% letter_sets[[x]])]
        if (all(adj[g, members])) { letter_sets[[g]] <- c(letter_sets[[g]], let); added <- TRUE }
      }
      if (!added || length(letter_sets[[g]])==0) {
        letter_sets[[g]] <- c(letter_sets[[g]], letters[cur]); cur <- cur+1L
      }
    }
    cld <- sapply(groups, function(g) paste(sort(unique(letter_sets[[g]])),collapse=""))
    w_tops <- df_sub[, .(w_top=boxplot.stats(Pheno)$stats[5],
                         DosChar=as.character(DosClass[1])), by=DosClass][,.(DosChar=as.character(DosClass),w_top)]
    cld_dt <- data.table(DosClass=groups, cld_letter=cld)
    cld_dt <- merge(cld_dt, w_tops, by.x="DosClass", by.y="DosChar", all.x=TRUE)
    y_rng  <- diff(range(df$Pheno, na.rm=TRUE))
    cld_dt[, label_y := w_top + y_rng*0.1]
    cld_dt
  }, error=function(e) NULL)
}

pub_theme <- theme_classic(base_size=11) +
  theme(axis.text      = element_text(color="black", size=9),
        axis.title     = element_text(size=10, face="bold"),
        plot.title     = element_text(size=10, face="bold"),
        plot.subtitle  = element_text(size=8,  color="grey40"),
        panel.border   = element_rect(color="black", fill=NA, linewidth=0.6))

make_box <- function(dos_dt, pheno_dt, pheno_col="Pheno",
                     fill_col, title_txt, sub_txt, is_primary=FALSE, y_label=NULL) {
  df <- merge(dos_dt, pheno_dt[, .(Sample, Pheno=get(pheno_col))], by="Sample")
  df <- df[!is.na(Dosage) & !is.na(Pheno)]
  if (nrow(df) < 10) return(NULL)
  counts <- df[, .N, by=DosClass][order(DosClass)]
  y_min  <- min(df$Pheno, na.rm=TRUE); y_max <- max(df$Pheno, na.rm=TRUE)
  y_rng  <- y_max - y_min
  y_bot  <- y_min - y_rng*0.22
  cld_dt <- get_cld(df)
  y_top  <- if (!is.null(cld_dt)) max(cld_dt$label_y, na.rm=TRUE)+y_rng*0.05 else y_max+y_rng*0.1

  p <- ggplot(df, aes(x=DosClass, y=Pheno)) +
    geom_boxplot(fill=fill_col, colour="#1a1a1a", alpha=0.75,
                 outlier.size=0.5, outlier.alpha=0.3, width=0.52) +
    geom_jitter(width=0.15, alpha=0.15, size=0.45, colour="#1a1a1a") +
    stat_summary(fun=mean, geom="point", shape=23, fill="white",
                 colour=fill_col, size=2.2, stroke=1.1) +
    geom_text(data=counts,
              aes(x=DosClass, y=y_bot+y_rng*0.08, label=paste0("n=",N)),
              size=1.9, colour="grey40", angle=90, hjust=0, vjust=0.5,
              inherit.aes=FALSE) +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(limits=c(y_bot, y_top)) +
    labs(title=title_txt, subtitle=sub_txt,
         x="ALT copies",
         y=if(!is.null(y_label)) y_label else NULL) +
    pub_theme +
    theme(plot.title   = element_text(size=9.5, face="bold", colour=fill_col),
          plot.subtitle= element_text(size=7.5, colour="grey40"),
          panel.border = element_rect(color=fill_col, fill=NA,
                                      linewidth=if(is_primary) 1.3 else 0.7))
  if (!is.null(cld_dt) && nrow(cld_dt)>0)
    p <- p + geom_text(data=cld_dt, aes(x=DosClass, y=label_y, label=cld_letter),
                       size=3.2, fontface="bold", colour=fill_col, inherit.aes=FALSE)
  p
}

# ── Colors ─────────────────────────────────────────────────────────────────────
COL_FL_PRIMARY <- "#8B0000"   # dark red
COL_FL_YR23    <- "#1A6B3C"   # dark green
COL_FL_YR24    <- "#5B2C8D"   # dark purple
COL_FR_PRIMARY <- "#1565C0"   # dark blue
COL_FR_YR23    <- "#0D7A3C"
COL_FR_YR24    <- "#4A148C"

# ── Build all panels ───────────────────────────────────────────────────────────
cat("Building allele effect panels...\n")
all_rows <- list()

for (loc in locus_list) {
  cat(sprintf("  %s (%s)...\n", loc$label, loc$marker))
  dos_dt <- get_dosage(loc$marker)
  if (is.null(dos_dt)) { cat("    Marker not found in geno matrix\n"); next }

  trait       <- loc$trait
  is_fl       <- trait == "DTFlower"
  pheno_exc   <- if (is_fl) pheno_fl_exc else pheno_fr_exc
  pheno_col23 <- if (is_fl) "DTFlower_yr.23" else "DTFruit_yr.23"
  pheno_col24 <- if (is_fl) "DTFlower_yr.24" else "DTFruit_yr.24"
  y_lbl       <- if (is_fl) "Days to 50% Flowering" else "Days to 50% Ripe Fruit"
  col_pri     <- if (is_fl) COL_FL_PRIMARY else COL_FR_PRIMARY
  col_23      <- if (is_fl) COL_FL_YR23    else COL_FR_YR23
  col_24      <- if (is_fl) COL_FL_YR24    else COL_FR_YR24

  sc_exc <- loc$score
  sc_23  <- get_score(scan_all, loc$marker, paste0(trait,"_yr.23"))
  sc_24  <- get_score(scan_all, loc$marker, paste0(trait,"_yr.24"))

  sub_exc <- sprintf("-log₁₀p = %.2f  |  PVE = %.1f%%", sc_exc, loc$pve)
  sub_23  <- if (!is.na(sc_23)) sprintf("-log₁₀p = %.2f", sc_23) else ""
  sub_24  <- if (!is.na(sc_24)) sprintf("-log₁₀p = %.2f", sc_24) else ""

  p1 <- make_box(dos_dt, pheno_exc,  "Pheno",    col_pri, "BLUE exc. yr.25", sub_exc, TRUE,  y_lbl)
  p2 <- make_box(dos_dt, pheno_all,  pheno_col23, col_23,  "yr.23",           sub_23,  FALSE, NULL)
  p3 <- make_box(dos_dt, pheno_all,  pheno_col24, col_24,  "yr.24",           sub_24,  FALSE, NULL)

  panels <- Filter(Negate(is.null), list(p1, p2, p3))
  if (length(panels) == 0) next

  row_plot <- plot_grid(plotlist=panels, nrow=1,
                        rel_widths=c(1.1, rep(1, length(panels)-1)),
                        align="h", axis="tb")

  # Add locus label on left
  locus_label <- ggdraw() +
    draw_label(loc$label, fontface="bold", size=9, angle=90, x=0.5, y=0.5)

  row_final <- plot_grid(locus_label, row_plot, nrow=1, rel_widths=c(0.08, 1))
  all_rows[[length(all_rows)+1]] <- row_final
}

cat(sprintf("Built %d locus rows\n", length(all_rows)))

# ── Assemble: 4 loci per page ─────────────────────────────────────────────────
n_per_page <- 4
pages      <- split(all_rows, ceiling(seq_along(all_rows)/n_per_page))

out_pdf <- file.path(OUT_DIR, "FigS_Allele_Effects_Secondary_Loci.pdf")
out_png <- file.path(OUT_DIR, "FigS_Allele_Effects_Secondary_Loci_p1.png")

cat("Saving PDF...\n")
pdf(out_pdf, width=10, height=11, onefile=TRUE)
for (pg in pages) {
  fig <- plot_grid(plotlist=pg, ncol=1,
                   rel_heights=rep(1, length(pg)))
  print(fig)
}
dev.off()

# Save page 1 as PNG preview
cat("Saving PNG preview (page 1)...\n")
fig_p1 <- plot_grid(plotlist=pages[[1]], ncol=1,
                    rel_heights=rep(1, length(pages[[1]])))
ggsave(out_png, fig_p1, width=10, height=11, dpi=300)

cat(sprintf("\nDone.\n  PDF: %s\n  PNG (p1): %s\n", out_pdf, out_png))
