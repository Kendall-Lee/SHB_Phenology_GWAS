#!/usr/bin/env Rscript
# ============================================================================
# Run_LinearPCA_FullPipeline.R
#
# Full re-run of the linear reference QTL-seq pipeline using population
# structure correction derived from the linear reference genotypes, not
# the pangenome dosage matrix.
#
# Background:
#   PCA correlation check (LinearPCA_vs_PangenomePCA_report.txt) showed
#   that linear-ref PC1-4 have |r| <= 0.76 against pangenome PCs. The two
#   references capture different structure axes, so the correction must be
#   re-derived from linear genotypes.
#
# Steps:
#   1. PCA from SHB_LINEAR_LRLP_SNPS (clean linear SNPs, 90 samples)
#   2. Correct AllYears BLUE and Method1 BLUE using linear PC1-4
#      -> writes *_LinPCAcorrected.txt files to PopStructureCorrection/Linear/
#   3. Select new HIGH/LOW bulks (top/bottom 10 per trait) from linear BLUEs
#      -> writes *_LinPCAcorr.High/Low.list to PopStructureCorrection/Linear/
#   4. Run QTL-seq — PRIMARY (yr23+yr24) and BACKUP (all-years)
#      -> results_LINPCA_M1_{SV,SNP}_<trait>/
#      -> results_LINPCA_ALLYRS_{SV,SNP}_<trait>/
#   5. Update LinearRef_QTLseq_Report.pdf (corrected methods text + new plots)
#   6. Update AllQTL_Combined.xlsx and AllQTL_Nonredundant.bed
#
# Created: 2026-04-07
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(writexl)
})

# ─── PATHS ────────────────────────────────────────────────────────────────────

BASE_DIR   <- "/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq"
CORR_DIR   <- file.path(BASE_DIR, "PopStructureCorrection")
LINEAR_DIR <- file.path(CORR_DIR, "Linear")
PHENO_DIR  <- file.path(BASE_DIR, "Phenotype_Data")
GENO_DIR   <- file.path(BASE_DIR, "Genotype_Data")
LIN_SRC    <- file.path(LINEAR_DIR, "Linear_QTLSeq.R")

LINEAR_SNP <- file.path(GENO_DIR,
  "SHB_LINEAR_LRLP_SNPS.redo.DP.2.maxmiss.4.minminor.25.cleannames")
DEDUP_SV   <- file.path(GENO_DIR,
  "SHB_LINEAR_LRLP_SVs.mindep2.minminor.5.maxmiss0.8.dedup.dosage.file")
CLEAN_SNP  <- LINEAR_SNP   # already clean

META_COLS  <- c("Marker","Chrom","Position","REF","ALT")

TRAITS     <- c("DTFlower","DTFruit","Flow2Fruit","FruitWT")
TRAIT_LABEL <- c(
  DTFlower   = "Days to Flower",
  DTFruit    = "Days to Fruit",
  Flow2Fruit = "Flower to Fruit (days)",
  FruitWT    = "Fruit Weight (g)"
)

RAW_PHENO <- list(
  DTFlower   = file.path(PHENO_DIR, "DTFlower_SHB_allPheno.txt"),
  DTFruit    = file.path(PHENO_DIR, "DTFruit_SHB_allPheno.txt"),
  Flow2Fruit = file.path(PHENO_DIR, "Flow2Fruit_SHB_allPheno.txt"),
  FruitWT    = file.path(PHENO_DIR, "FruitWT_SHB_allPheno.txt")
)

# Output corrected pheno files
M1_CORR <- list(
  DTFlower   = file.path(LINEAR_DIR,"DTFlower_LinPCAcorrected_Method1.txt"),
  DTFruit    = file.path(LINEAR_DIR,"DTFruit_LinPCAcorrected_Method1.txt"),
  Flow2Fruit = file.path(LINEAR_DIR,"Flow2Fruit_LinPCAcorrected_Method1.txt"),
  FruitWT    = file.path(LINEAR_DIR,"FruitWT_LinPCAcorrected_Method1.txt")
)
AY_CORR <- list(
  DTFlower   = file.path(LINEAR_DIR,"DTFlower_LinPCAcorrected_AllYears.txt"),
  DTFruit    = file.path(LINEAR_DIR,"DTFruit_LinPCAcorrected_AllYears.txt"),
  Flow2Fruit = file.path(LINEAR_DIR,"Flow2Fruit_LinPCAcorrected_AllYears.txt"),
  FruitWT    = file.path(LINEAR_DIR,"FruitWT_LinPCAcorrected_AllYears.txt")
)

PCA_SCORES  <- file.path(LINEAR_DIR, "LinearPCA_scores.csv")
PCA_PDF     <- file.path(LINEAR_DIR, "LinearPCA_QC.pdf")
BULK_N      <- 10L

HIGH_LISTS <- list(
  DTFlower   = file.path(LINEAR_DIR,"DTFlower_LinPCAcorr.High.list"),
  DTFruit    = file.path(LINEAR_DIR,"DTFruit_LinPCAcorr.High.list"),
  Flow2Fruit = file.path(LINEAR_DIR,"Flow2Fruit_LinPCAcorr.High.list"),
  FruitWT    = file.path(LINEAR_DIR,"FruitWT_LinPCAcorr.High.list")
)
LOW_LISTS <- list(
  DTFlower   = file.path(LINEAR_DIR,"DTFlower_LinPCAcorr.Low.list"),
  DTFruit    = file.path(LINEAR_DIR,"DTFruit_LinPCAcorr.Low.list"),
  Flow2Fruit = file.path(LINEAR_DIR,"Flow2Fruit_LinPCAcorr.Low.list"),
  FruitWT    = file.path(LINEAR_DIR,"FruitWT_LinPCAcorr.Low.list")
)

setwd(BASE_DIR)

# ─── STEP 1: PCA FROM LINEAR SNPs ────────────────────────────────────────────
#The overall flow is: load → filter (missingness + MAF) → thin → impute → PCA → save scores + QC plots — a standard population structure analysis pipeline, with the MAF formula correctly adapted for your tetraploid blueberry data.


message("\n========================================")
message("STEP 1: PCA from linear reference SNPs")
message("========================================")

message("Loading linear SNP dosage file...")
dos <- fread(LINEAR_SNP, showProgress = FALSE)
samp_cols <- setdiff(names(dos), META_COLS)
message(sprintf("  %d markers x %d samples", nrow(dos), length(samp_cols)))

mat_raw <- as.matrix(dos[, ..samp_cols])
class(mat_raw) <- "numeric"
rownames(mat_raw) <- dos$Marker
rm(dos); invisible(gc())

# FILTER

#For each marker (row), calculates the proportion of samples with missing dosage (NA). is.na() produces a logical matrix; rowMeans treats TRUE=1, FALSE=0.
#Calculates minor allele frequency (MAF) for each marker. For each row:
#Drops NA values
#Returns 0 if fewer than 10 non-missing samples (not enough data)
#Otherwise computes MAF as min(mean(dosage)/4, 1 - mean(dosage)/4) — this is the tetraploid MAF formula, where dosage ranges 0–4, so dividing by 4 gives allele frequency, and min(p, 1-p) gives the minor allele frequency

miss_rate <- rowMeans(is.na(mat_raw))
maf_vec <- apply(mat_raw, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 10) return(0)
  min(mean(x)/4, 1 - mean(x)/4)
})
keep <- miss_rate <= 0.20 & !is.na(maf_vec) & maf_vec >= 0.05  #Builds a logical vector flagging markers that pass both filters: missingness ≤ 20% AND MAF ≥ 5%.
message(sprintf("  After filters (miss<=20%%, MAF>=5%%): %d / %d markers",
                sum(keep), length(keep)))

# Thin to ~10k markers
#To keep PCA computationally tractable, thins markers to ~10,000 by taking every Nth marker. thin_step is calculated so that evenly-spaced sampling gives roughly 10k markers. max(1L, ...) prevents a step size of 0 if fewer than 10k markers passed filtering.
keep_idx  <- which(keep)
thin_step <- max(1L, floor(length(keep_idx) / 10000L))
thin_idx  <- keep_idx[seq(1, length(keep_idx), by = thin_step)]
mat_f <- mat_raw[thin_idx, ]
message(sprintf("  After thinning (every %dth): %d markers", thin_step, nrow(mat_f)))
rm(mat_raw); invisible(gc())

# Impute
#For each marker with any remaining NAs, fills them with that marker's mean dosage across non-missing samples — a simple mean imputation strategy.
row_means <- rowMeans(mat_f, na.rm = TRUE)
for (i in seq_len(nrow(mat_f))) {
  na_idx <- is.na(mat_f[i,])
  if (any(na_idx)) mat_f[i, na_idx] <- row_means[i]
}
mat_f[is.na(mat_f)] <- mean(mat_f, na.rm = TRUE)

# PCA
#Runs PCA. The matrix is transposed (t()) so that samples are rows and markers are columns — standard orientation for PCA. center = TRUE mean-centers each marker; scale. = FALSE does not standardize variance (common for dosage data where you want to preserve variance differences).
message("Running prcomp()...")
pca_res <- prcomp(t(mat_f), center = TRUE, scale. = FALSE)
pct_exp <- 100 * pca_res$sdev^2 / sum(pca_res$sdev^2)
cum_exp <- cumsum(pct_exp)
#Computes percent variance explained per PC (sdev = standard deviations of each PC; squaring gives variance). cumsum gives cumulative variance.
message("Variance explained:")
for (i in 1:6) message(sprintf("  PC%d: %.2f%%  (cum: %.2f%%)", i, pct_exp[i], cum_exp[i]))
#^Prints variance explained for the first 6 PCs — useful for a quick sanity check of how much structure the top PCs are capturing.

scores_dt <- data.table(Sample = rownames(pca_res$x), pca_res$x[, 1:10]). #Extracts the first 10 PC scores for each sample into a data.table (with sample names as a column) and writes it to the path PCA_SCORES for downstream use (e.g., as covariates in GWASpoly).
fwrite(scores_dt, PCA_SCORES)
message("Saved: ", PCA_SCORES)

# PCA QC plots

#Creates a scree plot — a bar chart of variance explained for PCs 1–15, with percentage labels on top of each bar. Useful for deciding how many PCs to include as covariates.

p_scree <- ggplot(data.table(PC=1:15, PVE=pct_exp[1:15]),
                  aes(factor(PC), PVE)) +
  geom_col(fill="#2C7BB6", alpha=0.85) +
  geom_text(aes(label=sprintf("%.1f",PVE)), vjust=-0.3, size=2.8) +
  labs(title="Linear-ref SNP PCA — scree",
       subtitle=sprintf("%d markers (thinned), n=%d samples",
                        nrow(mat_f), length(samp_cols)),
       x="PC", y="Variance explained (%)") +
  theme_bw(base_size=10)

#Scatter plot of PC1 vs PC2, with sample labels that repel each other to avoid overlap (ggrepel). The axis labels include variance explained percentages.

p_12 <- ggplot(scores_dt, aes(PC1,PC2,label=Sample)) +
  geom_point(colour="#2C7BB6",size=2,alpha=0.8) +
  ggrepel::geom_text_repel(size=2,max.overlaps=15,colour="grey30") +
  labs(title="PC1 vs PC2",
       x=sprintf("PC1 (%.1f%%)",pct_exp[1]),
       y=sprintf("PC2 (%.1f%%)",pct_exp[2])) +
  theme_bw(base_size=10)

p_34 <- ggplot(scores_dt, aes(PC3,PC4,label=Sample)) +
  geom_point(colour="#D7191C",size=2,alpha=0.8) +
  ggrepel::geom_text_repel(size=2,max.overlaps=15,colour="grey30") +
  labs(title="PC3 vs PC4",
       x=sprintf("PC3 (%.1f%%)",pct_exp[3]),
       y=sprintf("PC4 (%.1f%%)",pct_exp[4])) +
  theme_bw(base_size=10)

# ─── STEP 2: CORRECT PHENOTYPES ──────────────────────────────────────────────

message("\n========================================")
message("STEP 2: PCA-correct phenotypes (linear PCs)")
message("========================================")

correct_pheno <- function(pheno_dt, pca, blue_col) {
  m <- merge(pheno_dt,
             pca[, .(Sample, PC1, PC2, PC3, PC4)],
             by.x = "DNA ID", by.y = "Sample", all.x = TRUE)
  has <- !is.na(m$PC1) & !is.na(m[[blue_col]])
  if (sum(has) < 5) { message("  Too few overlapping samples — skip"); return(NULL) }
  mod <- lm(as.formula(paste(blue_col,"~ PC1+PC2+PC3+PC4")), data = m[has])
  r2  <- summary(mod)$r.squared
  message(sprintf("    %s: n=%d  model_R2=%.3f", blue_col, sum(has), r2))
  m[, BLUE_corrected := NA_real_]
  m[has, BLUE_corrected := residuals(mod)]
  raw_mu  <- mean(m[[blue_col]], na.rm = TRUE)
  corr_mu <- mean(m$BLUE_corrected, na.rm = TRUE)
  no_pca  <- !is.na(m[[blue_col]]) & is.na(m$PC1)
  m[no_pca, BLUE_corrected := get(blue_col) - raw_mu + corr_mu]
  list(data = m, r2 = r2)
}

pheno_plots <- list()

for (tr in TRAITS) {
  message(sprintf("  %s", tr))
  pheno <- fread(RAW_PHENO[[tr]])

  # Standardise BLUE column name
  blue_raw <- intersect(c("Data_BLUE","BLUES","BLUE"), names(pheno))[1]
  if (is.na(blue_raw)) { message("  No BLUE column — skip"); next }
  if (blue_raw != "Data_BLUE") setnames(pheno, blue_raw, "Data_BLUE")

  # All-years correction
  res_ay <- correct_pheno(pheno, scores_dt, "Data_BLUE")
  if (!is.null(res_ay)) {
    out <- res_ay$data[, .(`DNA ID`, yr.23, yr.24, yr.25, Data_BLUE, BLUE_corrected)]
    fwrite(out, AY_CORR[[tr]], sep = "\t")
    message("    Saved AllYears: ", AY_CORR[[tr]])
  }

  # Method1 correction (mean yr23+yr24)
  pheno[, Method1_raw := ifelse(
    !is.na(yr.23) & !is.na(yr.24), (yr.23 + yr.24) / 2,
    ifelse(!is.na(yr.23), yr.23, yr.24)
  )]
  res_m1 <- correct_pheno(pheno, scores_dt, "Method1_raw")
  if (!is.null(res_m1)) {
    out_m1 <- res_m1$data[, .(`DNA ID`, yr.23, yr.24, yr.25,
                               Data_BLUE = Method1_raw, BLUE_corrected)]
    fwrite(out_m1, M1_CORR[[tr]], sep = "\t")
    message("    Saved Method1:  ", M1_CORR[[tr]])
  }

  # QC scatter
  if (!is.null(res_ay)) {
    pd <- res_ay$data[!is.na(Data_BLUE) & !is.na(BLUE_corrected)]
    r  <- round(cor(pd$Data_BLUE, pd$BLUE_corrected), 3)
    pheno_plots[[tr]] <- ggplot(pd, aes(Data_BLUE, BLUE_corrected)) +
      geom_point(alpha=0.55, size=1.8, colour="#2C7BB6") +
      geom_smooth(method="lm", se=FALSE, colour="firebrick", linewidth=0.8) +
      labs(title=TRAIT_LABEL[tr],
           subtitle=sprintf("r=%s  |  PC model R2=%.3f", r, res_ay$r2),
           x="Raw BLUE", y="Linear-PCA-corrected BLUE") +
      theme_bw(base_size=10)
  }
}

pdf(PCA_PDF, width=14, height=9)
print(p_scree)
print(p_12 | p_34)
if (length(pheno_plots) > 0) print(wrap_plots(pheno_plots, ncol=2))
dev.off()
message("PCA QC PDF: ", PCA_PDF)

# ─── STEP 3: SELECT BULKS ────────────────────────────────────────────────────

message("\n========================================")
message("STEP 3: Select HIGH/LOW bulks (n=", BULK_N, " per trait)")
message("========================================")

# Samples present in both SV and SNP files
sv_samps  <- names(fread(DEDUP_SV,  nrows=0))[-(1:5)]
snp_samps <- names(fread(CLEAN_SNP, nrows=0))[-(1:5)]
geno_samps <- intersect(sv_samps, snp_samps)
message(sprintf("  Samples in both SV+SNP files: %d", length(geno_samps)))

bulk_summary <- list()

for (tr in TRAITS) {
  if (!file.exists(M1_CORR[[tr]])) next
  pheno <- fread(M1_CORR[[tr]])
  # Keep only samples with genotype data and valid corrected BLUE
  pheno_ok <- pheno[`DNA ID` %in% geno_samps & !is.na(BLUE_corrected)]
  setorder(pheno_ok, BLUE_corrected)

  n_elig <- nrow(pheno_ok)
  if (n_elig < 2 * BULK_N) {
    message(sprintf("  %s: only %d eligible — skipping", tr, n_elig))
    next
  }

  low_ids  <- pheno_ok[1:BULK_N, `DNA ID`]
  high_ids <- pheno_ok[(n_elig - BULK_N + 1):n_elig, `DNA ID`]

  writeLines(high_ids, HIGH_LISTS[[tr]])
  writeLines(low_ids,  LOW_LISTS[[tr]])

  delta <- mean(pheno_ok[(n_elig-BULK_N+1):n_elig, BLUE_corrected]) -
           mean(pheno_ok[1:BULK_N, BLUE_corrected])
  message(sprintf("  %-12s: HIGH=[%s...%s] (mean=%.1f)  LOW=[%s...%s] (mean=%.1f)  delta=%.2f",
                  tr,
                  high_ids[1], high_ids[BULK_N],
                  mean(pheno_ok[(n_elig-BULK_N+1):n_elig, BLUE_corrected]),
                  low_ids[1], low_ids[BULK_N],
                  mean(pheno_ok[1:BULK_N, BLUE_corrected]),
                  delta))

  bulk_summary[[tr]] <- data.table(
    Trait = tr, Bulk = c(rep("HIGH",BULK_N), rep("LOW",BULK_N)),
    SampleID = c(high_ids, low_ids),
    BLUE_corrected = c(pheno_ok[(n_elig-BULK_N+1):n_elig, BLUE_corrected],
                       pheno_ok[1:BULK_N, BLUE_corrected])
  )
}

fwrite(rbindlist(bulk_summary),
       file.path(LINEAR_DIR, "BulkComposition_LinPCA.csv"))
message("Bulk CSV: ", file.path(LINEAR_DIR, "BulkComposition_LinPCA.csv"))

# ─── STEP 4: QTL-SEQ ─────────────────────────────────────────────────────────

message("\n========================================")
message("STEP 4: QTL-seq (linear PCA correction)")
message("========================================")

# -- Patch helper --
make_patched_script <- function(pheno_files, high_lists, low_lists,
                                sv_prefix, snp_prefix) {
  lines <- readLines(LIN_SRC)

  # Output directory
  lines <- sub('LINEAR_DIR <- file.path(BASE_DIR, "Linear")',
               sprintf('LINEAR_DIR <- "%s"', LINEAR_DIR),
               lines, fixed=TRUE)
  # Genotype files
  lines <- sub(
    'SV_FILE  <- file.path(LINEAR_DIR, "LRLP_suzihap1_SVs.mindep2.minminor.5.maxmiss0.8.dosage.file")',
    sprintf('SV_FILE  <- "%s"', DEDUP_SV), lines, fixed=TRUE)
  lines <- sub(
    'SNP_FILE <- file.path(LINEAR_DIR, "LRLP_suzihap1_SNPs.mindep2.minminor.5.maxmiss0.8.dosage.file")',
    sprintf('SNP_FILE <- "%s"', CLEAN_SNP), lines, fixed=TRUE)
  # Output prefix
  lines <- sub(
    'OUT_PREFIX  <- if (SNP_MODE) "results_SNP_" else "results_"',
    sprintf('OUT_PREFIX  <- if (SNP_MODE) "%s" else "%s"', snp_prefix, sv_prefix),
    lines, fixed=TRUE)
  # BLUE column
  lines <- gsub('"Data_BLUE"', '"BLUE_corrected"', lines, fixed=TRUE)
  lines <- gsub('"BLUES"',     '"BLUE_corrected"', lines, fixed=TRUE)
  # Pheno files
  for (tr in TRAITS) {
    old_name <- switch(tr,
      FruitWT    = 'pheno = "FruitWT_SHB_allPheno.txt"',
      DTFruit    = 'pheno = "DTFruit_SHB_allPheno.txt"',
      DTFlower   = 'pheno = "DTFlower_SHB_allPheno.txt"',
      Flow2Fruit = 'pheno = "Flow2Fruit_SHB_allPheno.txt"'
    )
    lines <- sub(old_name,
                 sprintf('pheno = "%s"', pheno_files[[tr]]),
                 lines, fixed=TRUE)
  }
  lines <- sub('pheno <- fread(file.path(BASE_DIR, cfg$pheno))',
               'pheno <- fread(cfg$pheno)', lines, fixed=TRUE)
  # Bulk list files
  bulk_map <- list(
    FruitWT_H  = c('high  = "FruitWT.High.list"',   HIGH_LISTS[["FruitWT"]]),
    FruitWT_L  = c('low   = "FruitWT.Low.list"',    LOW_LISTS[["FruitWT"]]),
    DTFruit_H  = c('high  = "D2Fruit.high.list"',   HIGH_LISTS[["DTFruit"]]),
    DTFruit_L  = c('low   = "D2Fruit.low.list"',    LOW_LISTS[["DTFruit"]]),
    DTFlower_H = c('high  = "D2Flower.High.list"',  HIGH_LISTS[["DTFlower"]]),
    DTFlower_L = c('low   = "D2Flower.Low.list"',   LOW_LISTS[["DTFlower"]]),
    Flow2F_H   = c('high  = "Flow2Fruit.High.list"',HIGH_LISTS[["Flow2Fruit"]]),
    Flow2F_L   = c('low   = "Flow2Fruit.Low.list"', LOW_LISTS[["Flow2Fruit"]])
  )
  for (bm in bulk_map) {
    lines <- sub(bm[1], sprintf('%s"%s"', sub('=.*', '= ', bm[1]), bm[2]),
                 lines, fixed=TRUE)
  }
  lines <- sub('high_bulk <- scan(file.path(BASE_DIR, cfg$high), what = "character", quiet = TRUE)',
               'high_bulk <- scan(cfg$high, what = "character", quiet = TRUE)',
               lines, fixed=TRUE)
  lines <- sub('low_bulk  <- scan(file.path(BASE_DIR, cfg$low),  what = "character", quiet = TRUE)',
               'low_bulk  <- scan(cfg$low,  what = "character", quiet = TRUE)',
               lines, fixed=TRUE)

  stopifnot(any(grepl("BLUE_corrected", lines)))
  stopifnot(any(grepl("LinPCAcorr",     lines)))
  stopifnot(any(grepl(basename(DEDUP_SV), lines)))

  tmp <- tempfile(fileext=".R")
  writeLines(lines, tmp)
  tmp
}

run_mode <- function(tmp_script, sv_pfx, snp_pfx, mode) {
  for (tr in TRAITS) {
    pfx     <- if (mode=="SV") sv_pfx else snp_pfx
    out_dir <- file.path(LINEAR_DIR, paste0(pfx, tr))
    if (dir.exists(out_dir) && length(list.files(out_dir)) > 0) {
      message(sprintf("  SKIP %s (%s) — results exist", tr, mode)); next
    }
    message(sprintf("  --- %s (%s) ---", tr, mode))
    t0 <- proc.time()
    tryCatch({
      env <- new.env(parent=globalenv())
      env$commandArgs <- function(...) if (mode=="SNP") c(tr,"--snp") else c(tr)
      source(tmp_script, local=env, echo=FALSE)
      message(sprintf("    Done in %.1f s", (proc.time()-t0)[["elapsed"]]))
    }, error=function(e) message(sprintf("    ERROR: %s", e$message)))
  }
}

# PRIMARY — Method1 (yr23+yr24)
message("\n  PRIMARY: yr23+yr24 (Method1, linear PCA corrected)")
m1_tmp <- make_patched_script(
  M1_CORR, HIGH_LISTS, LOW_LISTS,
  sv_prefix  = "results_LINPCA_M1_PCAcorr_",
  snp_prefix = "results_LINPCA_M1_SNP_PCAcorr_"
)
run_mode(m1_tmp, "results_LINPCA_M1_PCAcorr_", "results_LINPCA_M1_SNP_PCAcorr_", "SV")
run_mode(m1_tmp, "results_LINPCA_M1_PCAcorr_", "results_LINPCA_M1_SNP_PCAcorr_", "SNP")

# BACKUP — All-years
message("\n  BACKUP: all-years (linear PCA corrected)")
ay_tmp <- make_patched_script(
  AY_CORR, HIGH_LISTS, LOW_LISTS,
  sv_prefix  = "results_LINPCA_ALLYRS_PCAcorr_",
  snp_prefix = "results_LINPCA_ALLYRS_SNP_PCAcorr_"
)
run_mode(ay_tmp, "results_LINPCA_ALLYRS_PCAcorr_", "results_LINPCA_ALLYRS_SNP_PCAcorr_", "SV")
run_mode(ay_tmp, "results_LINPCA_ALLYRS_PCAcorr_", "results_LINPCA_ALLYRS_SNP_PCAcorr_", "SNP")

message("\nQTL-seq complete. Sourcing report and BED scripts...")

# ─── STEP 5: UPDATE REPORT ───────────────────────────────────────────────────
# Patch Report_LinearRef_QTLseq.R to use the new LinPCA paths, then source it.

message("\n========================================")
message("STEP 5: Regenerate report PDF")
message("========================================")

rpt_src  <- file.path(LINEAR_DIR, "Report_LinearRef_QTLseq.R")
rpt_lines <- readLines(rpt_src)

# Point pheno files to LinPCA versions
for (tr in TRAITS) {
  old_m1 <- sprintf('file.path(PHENO_DIR, "%s_SHB_allPheno_Method1_QTLcorrected.txt")', tr)
  new_m1 <- sprintf('"%s"', M1_CORR[[tr]])
  rpt_lines <- gsub(old_m1, new_m1, rpt_lines, fixed=TRUE)

  old_ay <- sprintf('file.path(PHENO_DIR, "%s_SHB_allPheno_QTLcorrected.txt")', tr)
  new_ay <- sprintf('"%s"', AY_CORR[[tr]])
  rpt_lines <- gsub(old_ay, new_ay, rpt_lines, fixed=TRUE)
}

# Point bulk lists to LinPCA versions
for (tr in TRAITS) {
  pfx  <- sub("WT$","Wt",tr)   # DTFlower → DTFlower, FruitWT → FruitWt etc.
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "D2Flower_PCAcorr.High.list")'),
    sprintf('"%s"', HIGH_LISTS[["DTFlower"]]), rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "D2Flower_PCAcorr.Low.list")'),
    sprintf('"%s"', LOW_LISTS[["DTFlower"]]),  rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "D2Fruit_PCAcorr.High.list")'),
    sprintf('"%s"', HIGH_LISTS[["DTFruit"]]),  rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "D2Fruit_PCAcorr.Low.list")'),
    sprintf('"%s"', LOW_LISTS[["DTFruit"]]),   rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "Flow2Fruit_PCAcorr.High.list")'),
    sprintf('"%s"', HIGH_LISTS[["Flow2Fruit"]]), rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "Flow2Fruit_PCAcorr.Low.list")'),
    sprintf('"%s"', LOW_LISTS[["Flow2Fruit"]]),  rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "FruitWT_PCAcorr.High.list")'),
    sprintf('"%s"', HIGH_LISTS[["FruitWT"]]),  rpt_lines, fixed=TRUE)
  rpt_lines <- gsub(
    sprintf('file.path(PHENO_DIR, "FruitWT_PCAcorr.Low.list")'),
    sprintf('"%s"', LOW_LISTS[["FruitWT"]]),   rpt_lines, fixed=TRUE)
}

# Update result directory prefixes
rpt_lines <- gsub("results_REDO_M1_PCAcorr_",       "results_LINPCA_M1_PCAcorr_",       rpt_lines)
rpt_lines <- gsub("results_REDO_M1_SNP_PCAcorr_",   "results_LINPCA_M1_SNP_PCAcorr_",   rpt_lines)
rpt_lines <- gsub("results_REDO_ALLYRS_PCAcorr_",   "results_LINPCA_ALLYRS_PCAcorr_",   rpt_lines)
rpt_lines <- gsub("results_REDO_ALLYRS_SNP_PCAcorr_","results_LINPCA_ALLYRS_SNP_PCAcorr_",rpt_lines)

# Update methods text: replace pangenome PCA note with linear PCA note
rpt_lines <- gsub(
  "A PCA was performed on the pangenome dosage matrix; the first four PCs were used as",
  "A PCA was performed on the linear-reference SNP dosage matrix (157,359 markers, 90",
  rpt_lines, fixed=TRUE)
rpt_lines <- gsub(
  "covariates in a linear model to extract structure-corrected BLUEs (BLUE_corrected).",
  "samples after filtering); the first four PCs were used as covariates in a linear",
  rpt_lines, fixed=TRUE)
# Patch the PHENO_DIR path used in the report
rpt_lines <- gsub(
  'PHENO_DIR   <- file.path(BASE_DIR, "PopStructureCorrection")',
  sprintf('PHENO_DIR   <- "%s"', LINEAR_DIR),
  rpt_lines, fixed=TRUE)

rpt_tmp <- tempfile(fileext=".R")
writeLines(rpt_lines, rpt_tmp)
message("Sourcing patched report script...")
source(rpt_tmp, local=new.env(parent=globalenv()), echo=FALSE)

# ─── STEP 6: UPDATE COMBINED BED + XLSX ──────────────────────────────────────

message("\n========================================")
message("STEP 6: Update AllQTL_Combined files")
message("========================================")

comb_src   <- file.path(CORR_DIR, "Compile_AllQTL_Combined.R")
comb_lines <- readLines(comb_src)

# Update linear result prefixes
comb_lines <- gsub("results_REDO_M1_PCAcorr_",        "results_LINPCA_M1_PCAcorr_",       comb_lines)
comb_lines <- gsub("results_REDO_M1_SNP_PCAcorr_",    "results_LINPCA_M1_SNP_PCAcorr_",   comb_lines)

comb_tmp <- tempfile(fileext=".R")
writeLines(comb_lines, comb_tmp)
message("Sourcing patched Compile_AllQTL_Combined script...")
source(comb_tmp, local=new.env(parent=globalenv()), echo=FALSE)

message("\n============================================================")
message("PIPELINE COMPLETE")
message("============================================================")
message("PCA scores  : ", PCA_SCORES)
message("PCA QC PDF  : ", PCA_PDF)
message("Report PDF  : ", file.path(LINEAR_DIR, "LinearRef_QTLseq_Report.pdf"))
message("XLSX        : ", file.path(CORR_DIR, "AllQTL_Combined.xlsx"))
message("NR BED      : ", file.path(CORR_DIR, "AllQTL_Nonredundant.bed"))
