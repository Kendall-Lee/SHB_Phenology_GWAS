# ============================================================================
# Linear_QTLSeq.R
# QTL-seq analysis using SNP + SV dosage files from Linear/
#
# WORKFLOW:
#   1. Compute per-marker ΔDS, R², P-value using BLUE phenotype bulks
#   2. Smooth ΔDS along each chromosome
#   3. Manhattan-style genome-wide plot of smoothed |ΔDS|
#   4. Extract QTL regions (peak ΔDS, R², P-value)
#   5. Validate each QTL using year-specific bulk comparisons
#   6. Per-QTL plots showing BLUE + year-by-year ΔDS
#   7. Write QTL summary table
#
# Usage:
#   Rscript Linear_QTLSeq.R FruitWT
#   Rscript Linear_QTLSeq.R DTFruit
#   Rscript Linear_QTLSeq.R DTFlower
#   Rscript Linear_QTLSeq.R Flow2Fruit
#   Rscript Linear_QTLSeq.R ALL          # runs all traits
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(zoo)
})

# ─── CONFIGURATION ───────────────────────────────────────────────────────────

BASE_DIR   <- "/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq"
LINEAR_DIR <- file.path(BASE_DIR, "Linear")

SV_FILE  <- file.path(LINEAR_DIR, "LRLP_suzihap1_SVs.mindep2.minminor.5.maxmiss0.8.dosage.file")
SNP_FILE <- file.path(LINEAR_DIR, "LRLP_suzihap1_SNPs.mindep2.minminor.5.maxmiss0.8.dosage.file")

BULK_SIZE      <- 10      # individuals per bulk (for year-specific bulk re-ranking)
WINDOW_SIZE    <- 201     # markers for rolling mean smoothing (odd number)
QTL_PERCENTILE <- 0.995   # top 0.5% of smoothed |ΔDS| defines QTL candidates
QTL_MERGE_BP   <- 5e6     # merge candidate peaks within 5 Mb
MIN_N_VALID    <- 15      # min samples with both geno+pheno for R²/P-value

# Trait definitions: phenotype file, pre-defined BLUE bulk lists, plot label
TRAIT_CONFIG <- list(
  FruitWT = list(
    pheno = "FruitWT_SHB_allPheno.txt",
    high  = "FruitWT.High.list",
    low   = "FruitWT.Low.list",
    label = "Fruit Weight (g)"
  ),
  DTFruit = list(
    pheno = "DTFruit_SHB_allPheno.txt",
    high  = "D2Fruit.high.list",
    low   = "D2Fruit.low.list",
    label = "Days to Fruit"
  ),
  DTFlower = list(
    pheno = "DTFlower_SHB_allPheno.txt",
    high  = "D2Flower.High.list",
    low   = "D2Flower.Low.list",
    label = "Days to Flower"
  ),
  Flow2Fruit = list(
    pheno = "Flow2Fruit_SHB_allPheno.txt",
    high  = "Flow2Fruit.High.list",
    low   = "Flow2Fruit.Low.list",
    label = "Flower to Fruit Interval"
  )
)

# ─── COMMAND LINE ─────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript Linear_QTLSeq.R <TRAIT|ALL> [--snp]\n")
  cat("  --snp  : use SNP dosage file (PVE >= 10%), outputs to results_SNP_<TRAIT>/\n")
  cat("Traits:", paste(names(TRAIT_CONFIG), collapse = ", "), "\n")
  stop("Provide a trait name", call. = FALSE)
}
RUN_TRAIT  <- args[1]
SNP_MODE   <- "--snp" %in% args
DOSAGE_FILE <- if (SNP_MODE) SNP_FILE else SV_FILE
OUT_PREFIX  <- if (SNP_MODE) "results_SNP_" else "results_"
MIN_R2_USE  <- if (SNP_MODE) 0.10 else 0.15   # PVE >= 10% for SNPs

if (RUN_TRAIT == "ALL") {
  TRAITS_TO_RUN <- names(TRAIT_CONFIG)
} else {
  if (!RUN_TRAIT %in% names(TRAIT_CONFIG)) {
    stop(sprintf("Unknown trait '%s'. Valid: %s",
                 RUN_TRAIT, paste(names(TRAIT_CONFIG), collapse = ", ")))
  }
  TRAITS_TO_RUN <- RUN_TRAIT
}
message(sprintf("Mode: %s | Min R2: %.2f | File: %s",
                if(SNP_MODE) "SNP" else "SV", MIN_R2_USE, basename(DOSAGE_FILE)))

setwd(BASE_DIR)

# ─── LOAD DOSAGE DATA (once, shared across traits) ────────────────────────────

message("\n=== Loading dosage data ===")
message("  Reading: ", basename(DOSAGE_FILE))
dos <- fread(DOSAGE_FILE, showProgress = FALSE, fill = Inf)
# Drop any extra columns beyond the expected header count (caused by malformed rows)
expected_ncol <- ncol(fread(DOSAGE_FILE, nrows = 0, showProgress = FALSE))
if (ncol(dos) > expected_ncol) dos <- dos[, seq_len(expected_ncol), with = FALSE]
# Normalize column names (SNP file uses Ref/Alt, SV file uses REF/ALT)
setnames(dos, old = intersect(c("Ref","Alt"), colnames(dos)),
              new = c("REF","ALT")[seq_along(intersect(c("Ref","Alt"), colnames(dos)))])
# Coerce all sample columns to numeric (malformed rows can corrupt column types)
meta_cols <- c("Marker","Chrom","Position","REF","ALT")
for (col in setdiff(colnames(dos), meta_cols)) {
  if (!is.numeric(dos[[col]])) set(dos, j = col, value = suppressWarnings(as.numeric(dos[[col]])))
}
invisible(gc())

dos <- dos[order(Chrom, Position)]
SAMPLE_COLS <- setdiff(colnames(dos), c("Marker", "Chrom", "Position", "REF", "ALT"))
CHROMOSOMES <- unique(dos$Chrom)

# Chromosome factor order (Chr.01 → Chr.12)
chr_levels <- CHROMOSOMES[order(as.integer(sub(".*\\.", "", CHROMOSOMES)))]
dos[, Chrom := factor(Chrom, levels = chr_levels)]

message(sprintf("  Total: %d markers, %d chromosomes, %d samples",
                nrow(dos), length(CHROMOSOMES), length(SAMPLE_COLS)))

# ─── HELPER FUNCTIONS ─────────────────────────────────────────────────────────

# Vectorized R² and P-value via point-biserial correlation
# Processes one chromosome's marker matrix at a time to control memory
# dos_sub: data.table subset for one chromosome
# pheno_vec: named numeric vector (name = sample ID, value = phenotype)
compute_r2_pval_chrom <- function(dos_sub, pheno_vec) {
  samps <- intersect(SAMPLE_COLS, names(pheno_vec))
  if (length(samps) < MIN_N_VALID) {
    return(data.table(R2 = rep(NA_real_, nrow(dos_sub)),
                      Pvalue = rep(NA_real_, nrow(dos_sub))))
  }

  M <- as.matrix(dos_sub[, ..samps])   # n_markers × n_samples
  y <- pheno_vec[samps]
  y2 <- y^2

  W  <- !is.na(M) + 0L            # 0/1 indicator matrix
  X  <- M; X[is.na(X)] <- 0       # NA → 0 for sum operations

  nv  <- rowSums(W)                # n valid per marker
  sx  <- rowSums(X)                # Σx per marker
  sy  <- as.vector(W  %*% y)       # Σy (only valid samples)
  sxy <- as.vector(X  %*% y)       # Σxy
  sx2 <- rowSums(X^2)              # Σx²
  sy2 <- as.vector(W  %*% y2)     # Σy² (only valid samples)

  # Pearson r (bivariate regression equivalent)
  cov_xy <- sxy - sx * sy / nv
  var_x  <- sx2 - sx^2 / nv
  var_y  <- sy2 - sy^2 / nv

  denom <- sqrt(var_x * var_y)
  r <- ifelse(denom > 0, cov_xy / denom, NA_real_)
  r <- pmin(pmax(r, -1), 1)

  r2   <- r^2
  tval <- r * sqrt(nv - 2) / sqrt(pmax(1 - r^2, 1e-12))
  pval <- 2 * pt(-abs(tval), df = nv - 2)

  # Mask low-coverage markers
  bad <- nv < MIN_N_VALID | var_x <= 0
  r2[bad]   <- NA_real_
  pval[bad] <- NA_real_

  data.table(R2 = r2, Pvalue = pval)
}

# ΔDS for a given set of high/low sample IDs
compute_delta <- function(dos_sub, high_ids, low_ids) {
  h <- intersect(high_ids, SAMPLE_COLS)
  l <- intersect(low_ids,  SAMPLE_COLS)
  rowMeans(as.matrix(dos_sub[, ..h]), na.rm = TRUE) -
  rowMeans(as.matrix(dos_sub[, ..l]), na.rm = TRUE)
}

# Rolling smooth per chromosome
smooth_chrom <- function(dt, col_in, col_out, w = WINDOW_SIZE) {
  dt[order(Chrom, Position),
     (col_out) := rollapply(get(col_in), width = w, FUN = mean,
                             align = "center", fill = NA, na.rm = TRUE),
     by = Chrom]
}

# Extract and merge QTL peaks
extract_qtl <- function(dt, delta_smooth_col, r2_col, pval_col) {
  thresh <- quantile(abs(dt[[delta_smooth_col]]), QTL_PERCENTILE, na.rm = TRUE)
  cands  <- dt[abs(get(delta_smooth_col)) >= thresh][order(Chrom, Position)]
  if (nrow(cands) == 0) return(data.table())

  res <- lapply(levels(cands$Chrom), function(chr) {
    sub <- cands[Chrom == chr]
    if (nrow(sub) == 0) return(NULL)
    sub[, grp := cumsum(c(0L, diff(Position) > QTL_MERGE_BP)) + 1L]
    sub[, {
      best <- which.max(abs(get(delta_smooth_col)))
      .(Chrom       = chr,
        Start       = min(Position),
        End         = max(Position),
        Width_Mb    = round((max(Position) - min(Position)) / 1e6, 3),
        PeakPos     = Position[best],
        PeakDelta   = round(get(delta_smooth_col)[best], 4),
        PeakR2      = round(max(get(r2_col),   na.rm = TRUE), 4),
        MinPval     = signif(min(get(pval_col), na.rm = TRUE), 3),
        NMarkers    = .N)
    }, by = grp][, grp := NULL]
  })
  rbindlist(res, fill = TRUE)
}

# ─── MANHATTAN PLOT (genome-wide, all chromosomes) ────────────────────────────

plot_manhattan <- function(dt, y_col, qtl_dt, trait_label, out_pdf) {
  # Compute cumulative x-axis offsets per chromosome
  chr_sizes <- dt[, .(MaxPos = max(Position, na.rm = TRUE)), by = Chrom]
  chr_sizes <- chr_sizes[order(as.integer(sub(".*\\.", "", as.character(Chrom))))]
  chr_sizes[, Offset := c(0, cumsum(as.numeric(MaxPos))[-.N])]

  dt2 <- merge(dt, chr_sizes[, .(Chrom, Offset)], by = "Chrom")
  dt2[, GenPos := Position + Offset]

  # Mid-point of each chromosome for x-axis labels
  chr_mids <- dt2[, .(Mid = mean(range(GenPos, na.rm = TRUE))), by = Chrom]

  # Alternating colors for chromosomes
  n_chr <- nlevels(dt2$Chrom)
  chr_cols <- setNames(rep(c("#2166AC", "#4DAF4A"), length.out = n_chr),
                       levels(dt2$Chrom))

  p <- ggplot(dt2[!is.na(get(y_col))],
              aes(x = GenPos, y = abs(get(y_col)), color = Chrom)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = chr_cols, guide = "none") +
    scale_x_continuous(
      breaks = chr_mids$Mid,
      labels = sub("Chr\\.", "", as.character(chr_mids$Chrom))
    ) +
    labs(
      title  = sprintf("%s - Smoothed |DeltaDS| (window = %d markers)", trait_label, WINDOW_SIZE),
      x      = "Chromosome",
      y      = expression("|" * Delta * "DS|  (smoothed)")
    ) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 9))

  # Add QTL threshold line
  thresh <- quantile(abs(dt[[y_col]]), QTL_PERCENTILE, na.rm = TRUE)
  p <- p + geom_hline(yintercept = thresh, linetype = "dashed",
                       color = "red", linewidth = 0.6) +
    annotate("text", x = max(dt2$GenPos, na.rm = TRUE) * 0.02,
             y = thresh * 1.05, label = sprintf("Top %.1f%%", (1 - QTL_PERCENTILE) * 100),
             color = "red", hjust = 0, size = 3)

  # Highlight QTL peak positions
  if (!is.null(qtl_dt) && nrow(qtl_dt) > 0) {
    qtl_pos <- merge(qtl_dt[, .(Chrom, PeakPos)],
                     chr_sizes[, .(Chrom, Offset)], by = "Chrom")
    qtl_pos[, GenPos := PeakPos + Offset]
    # Get y-value at peak from dt2
    qtl_y <- dt2[, .(GenPos, y_val = abs(get(y_col)))]
    qtl_pos2 <- qtl_pos[, {
      idx <- which.min(abs(dt2$GenPos - GenPos))
      .(GenPos = GenPos, y_val = dt2[[y_col]][idx])
    }, by = .(Chrom, PeakPos)]
    p <- p + geom_point(data = qtl_pos2,
                         aes(x = GenPos, y = abs(y_val)),
                         color = "red", size = 2, shape = 23,
                         fill = "red", inherit.aes = FALSE)
  }

  ggsave(out_pdf, plot = p, width = 14, height = 5)
  message("  Saved: ", out_pdf)
}

# ─── PER-QTL YEAR SUPPORT PLOTS ──────────────────────────────────────────────

# For each QTL region, plot ΔDS from BLUE bulks vs year-specific bulks
plot_qtl_year_support <- function(dos_region, qtl_row, pheno_dt, trait_label, out_pdf) {
  # BLUE bulk ΔDS is already in dos_region$DeltaDS_BLUE
  # Compute year-specific ΔDS by re-ranking individuals per year

  year_cols <- c("yr.23", "yr.24", "yr.25")
  year_labels <- c("Year 2023", "Year 2024", "Year 2025")

  plot_data <- list()

  # BLUE
  plot_data[["BLUE"]] <- data.table(
    Position = dos_region$Position,
    DeltaDS  = dos_region$DeltaDS_BLUE,
    Source   = "BLUE (pre-defined bulks)"
  )

  # Per-year: re-rank by year phenotype, take top/bottom BULK_SIZE
  for (i in seq_along(year_cols)) {
    yr <- year_cols[i]
    if (!yr %in% colnames(pheno_dt)) next
    yr_pheno <- pheno_dt[!is.na(get(yr)), .(sample = sample, val = get(yr))]
    yr_pheno <- yr_pheno[sample %in% SAMPLE_COLS]
    yr_pheno <- yr_pheno[order(-val)]
    h_yr <- head(yr_pheno$sample, BULK_SIZE)
    l_yr <- tail(yr_pheno$sample, BULK_SIZE)
    delta_yr <- compute_delta(dos_region, h_yr, l_yr)
    plot_data[[yr]] <- data.table(
      Position = dos_region$Position,
      DeltaDS  = delta_yr,
      Source   = year_labels[i]
    )
  }

  pdat <- rbindlist(plot_data)
  pdat[, Source := factor(Source,
                           levels = c("BLUE (pre-defined bulks)",
                                      "Year 2023", "Year 2024", "Year 2025"))]

  # Smoothed per-source
  pdat[order(Source, Position),
       DeltaDS_sm := rollapply(DeltaDS, width = 51, FUN = mean,
                                align = "center", fill = NA, na.rm = TRUE),
       by = Source]

  region_label <- sprintf("%s  %s: %s–%s Mb",
                          trait_label,
                          qtl_row$Chrom,
                          round(qtl_row$Start / 1e6, 2),
                          round(qtl_row$End / 1e6, 2))

  p <- ggplot(pdat[!is.na(DeltaDS_sm)],
              aes(x = Position / 1e6, y = DeltaDS_sm, color = Source)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = qtl_row$PeakPos / 1e6,
               linetype = "dotted", color = "red", linewidth = 0.7) +
    scale_color_manual(values = c("BLUE (pre-defined bulks)" = "black",
                                  "Year 2023" = "#E41A1C",
                                  "Year 2024" = "#377EB8",
                                  "Year 2025" = "#4DAF4A")) +
    facet_wrap(~Source, ncol = 2) +
    labs(title  = region_label,
         x      = "Position (Mb)",
         y      = expression(Delta * "DS  (smoothed)"),
         color  = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "grey90"))

  ggsave(out_pdf, plot = p, width = 10, height = 7)
  message("  Saved: ", out_pdf)
}

# ─── ALLELE EFFECT BOXPLOT ───────────────────────────────────────────────────
#
# For a single marker: group all individuals by dosage class (0–4),
# boxplot phenotype vs. class, pairwise Wilcoxon + BH correction,
# draw significance brackets for all significant pairs.

plot_allele_effect <- function(marker_id, dos_dt, pheno_vec, trait_label) {

  row <- dos_dt[Marker == marker_id]
  if (nrow(row) == 0) return(NULL)

  samps <- intersect(SAMPLE_COLS, names(pheno_vec))
  dosage_vals <- as.numeric(row[1, ..samps])
  names(dosage_vals) <- samps

  df <- data.table(
    sample    = samps,
    dosage    = dosage_vals,
    phenotype = pheno_vec[samps]
  )
  df <- df[!is.na(dosage) & !is.na(phenotype)]
  df[, dosage_class := factor(round(dosage))]

  # Keep only classes with ≥ 3 individuals
  keep <- df[, .N, by = dosage_class][N >= 3, dosage_class]
  df   <- df[dosage_class %in% keep]

  if (nrow(df) < 10 || length(unique(df$dosage_class)) < 2) return(NULL)

  # ── Pairwise Wilcoxon tests (BH-corrected) ──────────────────────────────
  classes <- sort(as.integer(as.character(unique(df$dosage_class))))
  pairs   <- combn(classes, 2, simplify = FALSE)

  raw_pvals <- sapply(pairs, function(pr) {
    g1 <- df[dosage_class == pr[1], phenotype]
    g2 <- df[dosage_class == pr[2], phenotype]
    if (length(g1) < 3 || length(g2) < 3) return(NA_real_)
    suppressWarnings(wilcox.test(g1, g2)$p.value)
  })
  adj_pvals <- p.adjust(raw_pvals, method = "BH")

  sig_label <- function(p) {
    if (is.na(p) || p >= 0.05) return(NA_character_)
    if (p < 0.001) "***" else if (p < 0.01) "**" else "*"
  }

  sig_pairs <- Filter(Negate(is.null), lapply(seq_along(pairs), function(i) {
    lbl <- sig_label(adj_pvals[i])
    if (is.na(lbl)) return(NULL)
    list(x1 = pairs[[i]][1], x2 = pairs[[i]][2], label = lbl)
  }))

  # ── Base plot ────────────────────────────────────────────────────────────
  n_per_class  <- df[, .N, by = dosage_class][order(as.integer(as.character(dosage_class)))]
  x_labels     <- paste0(n_per_class$dosage_class, "\n(n=", n_per_class$N, ")")

  pal <- colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#D6604D"))(5)
  fill_map <- setNames(pal, as.character(0:4))
  use_fills <- fill_map[as.character(sort(as.integer(as.character(keep))))]

  # Pearson r for subtitle annotation
  r_val  <- cor(as.integer(as.character(df$dosage_class)), df$phenotype, use = "complete.obs")
  r_sign <- ifelse(r_val >= 0, "+", "-")

  p <- ggplot(df, aes(x = dosage_class, y = phenotype, fill = dosage_class)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.55, linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 1.4, alpha = 0.55, color = "grey30") +
    scale_fill_manual(values = use_fills, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    labs(
      title    = marker_id,
      subtitle = sprintf("%s:%s  |  r = %s%.3f  |  %s",
                         as.character(row$Chrom[1]),
                         format(row$Position[1], big.mark = ","),
                         r_sign, abs(r_val), trait_label),
      x = "Dosage Class  (0 = hom-ref, 4 = hom-alt)",
      y = trait_label
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  # ── Significance brackets ────────────────────────────────────────────────
  if (length(sig_pairs) > 0) {
    y_max   <- max(df$phenotype, na.rm = TRUE)
    y_range <- diff(range(df$phenotype, na.rm = TRUE))
    step    <- y_range * 0.09

    # Stack brackets shortest-span first so they don't overlap
    spans       <- sapply(sig_pairs, function(s) s$x2 - s$x1)
    sig_ordered <- sig_pairs[order(spans)]

    # Map dosage class value → x-axis position index
    class_levels <- levels(df$dosage_class)
    pos_of <- function(cls) which(class_levels == as.character(cls))

    for (bi in seq_along(sig_ordered)) {
      s  <- sig_ordered[[bi]]
      x1 <- pos_of(s$x1);  x2 <- pos_of(s$x2)
      yb <- y_max + step * bi
      tick_h <- step * 0.2

      p <- p +
        annotate("segment", x = x1, xend = x2, y = yb, yend = yb,
                 color = "black", linewidth = 0.45) +
        annotate("segment", x = x1, xend = x1, y = yb - tick_h, yend = yb,
                 color = "black", linewidth = 0.45) +
        annotate("segment", x = x2, xend = x2, y = yb - tick_h, yend = yb,
                 color = "black", linewidth = 0.45) +
        annotate("text", x = (x1 + x2) / 2, y = yb + step * 0.1,
                 label = s$label, size = 3.5, hjust = 0.5)
    }
    p <- p + expand_limits(y = y_max + step * (length(sig_ordered) + 2))
  }

  p
}

# ─── MAIN ANALYSIS LOOP ───────────────────────────────────────────────────────

all_qtl_tables <- list()

for (TRAIT in TRAITS_TO_RUN) {
  cfg <- TRAIT_CONFIG[[TRAIT]]
  message("\n", strrep("=", 60))
  message("Trait: ", TRAIT, "  (", cfg$label, ")")
  message(strrep("=", 60))

  out_dir <- file.path(LINEAR_DIR, paste0(OUT_PREFIX, TRAIT))
  dir.create(out_dir, showWarnings = FALSE)

  # ── Load phenotype data ──────────────────────────────────────────────────
  pheno <- fread(file.path(BASE_DIR, cfg$pheno))
  setnames(pheno, "DNA ID", "sample")

  # Identify BLUE column
  blue_col <- intersect(c("BLUE", "BLUES", "Data_BLUE", "blue"), colnames(pheno))[1]
  if (is.na(blue_col)) stop("Cannot find BLUE column in ", cfg$pheno)
  message("  BLUE column: ", blue_col)

  blue_pheno <- setNames(pheno[[blue_col]], pheno$sample)
  blue_pheno <- blue_pheno[!is.na(blue_pheno)]

  # Year phenotype vectors
  yr_cols <- c("yr.23", "yr.24", "yr.25")
  yr_cols <- yr_cols[yr_cols %in% colnames(pheno)]

  # ── Load pre-defined bulk lists ───────────────────────────────────────────
  high_bulk <- scan(file.path(BASE_DIR, cfg$high), what = "character", quiet = TRUE)
  low_bulk  <- scan(file.path(BASE_DIR, cfg$low),  what = "character", quiet = TRUE)
  high_bulk <- intersect(high_bulk, SAMPLE_COLS)
  low_bulk  <- intersect(low_bulk,  SAMPLE_COLS)
  message(sprintf("  BLUE bulks: High=%d, Low=%d samples", length(high_bulk), length(low_bulk)))

  # ── Per-chromosome analysis ───────────────────────────────────────────────
  result_list <- vector("list", length(levels(dos$Chrom)))
  names(result_list) <- levels(dos$Chrom)

  for (chr in levels(dos$Chrom)) {
    message("  Processing ", chr, " ...", appendLF = FALSE)
    sub <- dos[Chrom == chr]

    # ΔDS (BLUE pre-defined bulks)
    delta_blue <- compute_delta(sub, high_bulk, low_bulk)

    # R² and P-value (regression of BLUE phenotype on dosage, all individuals)
    rp <- compute_r2_pval_chrom(sub, blue_pheno)

    result_list[[chr]] <- data.table(
      Marker   = sub$Marker,
      Chrom    = sub$Chrom,
      Position = sub$Position,
      DeltaDS  = delta_blue,
      R2       = rp$R2,
      Pvalue   = rp$Pvalue
    )
    message(sprintf(" %d markers, mean |DeltaDS|=%.3f", nrow(sub),
                    mean(abs(delta_blue), na.rm = TRUE)))
  }

  res <- rbindlist(result_list)
  res[, Chrom := factor(Chrom, levels = levels(dos$Chrom))]

  # ── Multiple testing correction ───────────────────────────────────────────
  res[, Pvalue_BH := p.adjust(Pvalue, method = "BH")]

  # ── Smooth ΔDS per chromosome ─────────────────────────────────────────────
  message("  Smoothing DeltaDS (window=", WINDOW_SIZE, " markers) ...")
  smooth_chrom(res, "DeltaDS", "DeltaDS_smooth")

  # ── Manhattan plot ────────────────────────────────────────────────────────
  message("  Generating Manhattan plot ...")
  qtl_prelim <- extract_qtl(res, "DeltaDS_smooth", "R2", "Pvalue")

  man_pdf <- file.path(out_dir, paste0(TRAIT, "_Manhattan_DeltaDS.pdf"))
  plot_manhattan(res, "DeltaDS_smooth", qtl_prelim, cfg$label, man_pdf)

  # Also: R² Manhattan plot
  smooth_chrom(res, "R2", "R2_smooth")
  res_r2 <- res[!is.na(R2_smooth)]
  r2_man_pdf <- file.path(out_dir, paste0(TRAIT, "_Manhattan_R2.pdf"))

  # Compute cumulative offsets for R² plot
  chr_off <- res[, .(MaxPos = max(Position, na.rm = TRUE)), by = Chrom]
  chr_off[, Offset := c(0, cumsum(as.numeric(MaxPos))[-.N])]
  res_plot <- merge(res[!is.na(R2_smooth)], chr_off[, .(Chrom, Offset)], by = "Chrom")
  res_plot[, GenPos := Position + Offset]
  chr_mids_r2 <- res_plot[, .(Mid = mean(range(GenPos, na.rm = TRUE))), by = Chrom]
  n_chr <- nlevels(res_plot$Chrom)
  r2_cols <- setNames(rep(c("#D95F02", "#7570B3"), length.out = n_chr), levels(res_plot$Chrom))

  p_r2 <- ggplot(res_plot, aes(x = GenPos, y = R2_smooth, color = Chrom)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = r2_cols, guide = "none") +
    scale_x_continuous(breaks = chr_mids_r2$Mid,
                       labels = sub("Chr\\.", "", as.character(chr_mids_r2$Chrom))) +
    labs(title = sprintf("%s - Smoothed R2 (window = %d markers)", cfg$label, WINDOW_SIZE),
         x = "Chromosome", y = expression("R²  (smoothed)")) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size = 9))
  ggsave(r2_man_pdf, plot = p_r2, width = 14, height = 5)
  message("  Saved: ", r2_man_pdf)

  # ── Extract QTL regions ───────────────────────────────────────────────────
  qtl <- extract_qtl(res, "DeltaDS_smooth", "R2", "Pvalue")
  if (nrow(qtl) == 0) {
    message("  No QTL regions found above threshold.")
    next
  }
  qtl[, Trait := TRAIT]
  message(sprintf("  Found %d QTL region(s)", nrow(qtl)))
  print(qtl[, .(Chrom, Start, End, Width_Mb, PeakDelta, PeakR2, MinPval)])

  # ── Year-by-year support ──────────────────────────────────────────────────
  message("  Generating year support plots for each QTL ...")
  for (qi in seq_len(nrow(qtl))) {
    q <- qtl[qi]
    # Extract markers in QTL window (± 2.5 Mb around peak)
    flank   <- 2.5e6
    region  <- res[Chrom == q$Chrom &
                   Position >= (q$PeakPos - flank) &
                   Position <= (q$PeakPos + flank)]
    region_dos <- dos[Chrom == q$Chrom &
                      Position >= (q$PeakPos - flank) &
                      Position <= (q$PeakPos + flank)]
    region_dos[, DeltaDS_BLUE := region$DeltaDS]

    out_pdf <- file.path(out_dir,
                         sprintf("%s_QTL%02d_%s_%sMb_YearSupport.pdf",
                                 TRAIT, qi, q$Chrom,
                                 round(q$PeakPos / 1e6, 1)))
    plot_qtl_year_support(region_dos, q, pheno, cfg$label, out_pdf)
  }

  # ── Allele effect plots: one per QTL region peak ─────────────────────────
  message("  Generating allele effect plots for QTL peaks ...")
  blue_pheno_named <- setNames(pheno[[blue_col]], pheno$sample)
  blue_pheno_named <- blue_pheno_named[!is.na(blue_pheno_named)]

  for (qi in seq_len(nrow(qtl))) {
    q       <- qtl[qi]
    peak_id <- res[Chrom == q$Chrom &
                   Position == q$PeakPos, Marker][1]
    if (is.na(peak_id)) next

    ae_plot <- plot_allele_effect(peak_id, dos, blue_pheno_named, cfg$label)
    if (is.null(ae_plot)) {
      message("    QTL", qi, " - insufficient dosage classes, skipping")
      next
    }
    ae_pdf <- file.path(out_dir,
                        sprintf("%s_QTL%02d_%s_%sMb_AlleleEffect.pdf",
                                TRAIT, qi, q$Chrom,
                                round(q$PeakPos / 1e6, 1)))
    ggsave(ae_pdf, plot = ae_plot, width = 6, height = 6)
    message("  Saved: ", ae_pdf)
  }

  # ── Top-25 allele effect multi-page PDF ───────────────────────────────────
  message("  Generating top-25 allele effect multi-page PDF ...")
  top25_ae_pdf <- file.path(out_dir, paste0(TRAIT, "_Top25_AlleleEffect.pdf"))
  pdf(top25_ae_pdf, width = 6.5, height = 6.5)

  # ── Top-25 QTL markers per trait (with REF/ALT) ──────────────────────────
  # Join stats back to dosage table to recover REF/ALT columns
  # Use match() to avoid cartesian join from duplicate Marker keys in dos
  ref_alt_lu <- dos[!duplicated(Marker), .(Marker, REF, ALT)]
  res_annot  <- res[, .(Marker, Chrom, Position, DeltaDS, DeltaDS_smooth,
                         R2, Pvalue, Pvalue_BH)]
  idx <- match(res_annot$Marker, ref_alt_lu$Marker)
  res_annot[, REF := ref_alt_lu$REF[idx]]
  res_annot[, ALT := ref_alt_lu$ALT[idx]]
  top25 <- res_annot[order(-abs(DeltaDS_smooth))][1:min(25, .N)]
  top25[, Rank := seq_len(.N)]
  top25 <- top25[, .(Rank, Marker, Chrom, Position, REF, ALT,
                     DeltaDS       = round(DeltaDS, 4),
                     DeltaDS_smooth = round(DeltaDS_smooth, 4),
                     R2            = round(R2, 4),
                     Pvalue        = signif(Pvalue, 3),
                     Pvalue_BH     = signif(Pvalue_BH, 3))]

  top25_txt  <- file.path(out_dir, paste0(TRAIT, "_Top25_QTL_markers.txt"))
  top25_xlsx <- file.path(out_dir, paste0(TRAIT, "_Top25_QTL_markers.xlsx"))
  fwrite(top25, top25_txt, sep = "\t")
  message("  Saved top-25 markers (txt): ", top25_txt)

  # Write Excel if writexl is available
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(as.data.frame(top25), top25_xlsx)
    message("  Saved top-25 markers (xlsx): ", top25_xlsx)
  } else {
    message("  (Install 'writexl' for Excel output: install.packages('writexl'))")
  }

  # Print allele effect plots for each top-25 marker into the open PDF device
  n_plotted <- 0
  for (mid in top25$Marker) {
    ae <- plot_allele_effect(mid, dos, blue_pheno_named, cfg$label)
    if (!is.null(ae)) { print(ae); n_plotted <- n_plotted + 1 }
  }
  dev.off()
  message(sprintf("  Saved top-25 allele effect PDF (%d plots): %s",
                  n_plotted, top25_ae_pdf))

  # ── Save full results table ───────────────────────────────────────────────
  full_out <- file.path(out_dir, paste0(TRAIT, "_all_markers_stats.csv.gz"))
  fwrite(res[, .(Marker, Chrom, Position, DeltaDS, DeltaDS_smooth,
                 R2, R2_smooth, Pvalue, Pvalue_BH)],
         full_out)
  message("  Saved full stats: ", full_out)

  qtl_out <- file.path(out_dir, paste0(TRAIT, "_QTL_regions.txt"))
  fwrite(qtl, qtl_out, sep = "\t")
  message("  Saved QTL table: ", qtl_out)

  all_qtl_tables[[TRAIT]] <- qtl
}

# ─── COMBINED QTL SUMMARY ─────────────────────────────────────────────────────

if (length(all_qtl_tables) > 0) {
  combined <- rbindlist(all_qtl_tables, fill = TRUE)
  combined_out <- file.path(LINEAR_DIR, "QTL_Summary_AllTraits.txt")
  fwrite(combined[order(Chrom, Start)], combined_out, sep = "\t")
  message("\n=== Combined QTL summary saved: ", combined_out, " ===")
  message(sprintf("Total QTL across all traits: %d", nrow(combined)))
  print(combined[, .(Trait, Chrom, Start, End, PeakDelta, PeakR2, MinPval)])
}

message("\nDone.\n")
