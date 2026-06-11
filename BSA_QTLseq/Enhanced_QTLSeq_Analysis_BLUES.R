#!/usr/bin/env Rscript
# Enhanced QTL-Seq Analysis with R² and PVE Calculation
# For tetraploid blueberry samples mapped on pangenome
# UPDATED: Uses BLUE phenotypes from [trait]_SHB_allPheno.txt files

# ============================================================================
# COMMAND-LINE ARGUMENTS
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Define valid traits
valid_traits <- c("DTflower", "DTfruit", "Flow2Fruit", "FruitWt")

if(length(args) < 1) {
  cat("\nUsage: Rscript Enhanced_QTLSeq_Analysis_BLUES.R <TRAIT>\n\n")
  cat("Available traits:\n")
  cat("  - DTflower    (Days to flowering)\n")
  cat("  - DTfruit     (Days to fruiting)\n")
  cat("  - Flow2Fruit  (Flowering to fruit interval)\n")
  cat("  - FruitWt     (Fruit weight)\n\n")
  cat("Example:\n")
  cat("  Rscript Enhanced_QTLSeq_Analysis_BLUES.R Flow2Fruit\n\n")
  stop("Please provide a trait name as argument", call. = FALSE)
}

TRAIT <- args[1]

# Validate the trait
if(!(TRAIT %in% valid_traits)) {
  cat(sprintf("\nError: '%s' is not a valid trait!\n\n", TRAIT))
  cat("Available traits: ", paste(valid_traits, collapse=", "), "\n\n")
  stop("Please use one of the listed traits", call. = FALSE)
}

cat(sprintf("Analyzing trait: %s\n", TRAIT))

setwd("/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq")

# Load required libraries
library(data.table)
library(ggplot2)
library(zoo)      # for rolling mean smoothing
library(dplyr)
library(tidyr)

cat("============================================\n")
cat("Enhanced QTL-Seq Analysis with BLUEs\n")
cat("============================================\n\n")
cat(sprintf("Trait: %s\n\n", TRAIT))

# ============================================================================
# 1. READ GENOTYPE AND PHENOTYPE DATA
# ============================================================================

cat("1. Reading genotype dosage data...\n")
dos <- fread("SHB_PAN_LR.DP.1.maxmiss.4.minminor.3")
cat(sprintf("   - Loaded %d markers for %d samples\n", nrow(dos), ncol(dos)-1))

# Read phenotype data
cat("\n2. Reading BLUE phenotype data...\n")
pheno_file <- paste0(TRAIT, "_SHB_allPheno.txt")

if(!file.exists(pheno_file)) {
  stop(sprintf("ERROR: Phenotype file not found: %s\n", pheno_file))
}

pheno <- fread(pheno_file)
cat(sprintf("   - Loaded phenotype data from: %s\n", pheno_file))
cat(sprintf("   - Phenotype data for %d samples\n", nrow(pheno)))

# Check for BLUE column (handle different naming conventions)
blue_col <- NULL
possible_names <- c("BLUE", "BLUES", "Data_BLUE", "blue", "blues")

for(col_name in possible_names) {
  if(col_name %in% colnames(pheno)) {
    blue_col <- col_name
    break
  }
}

if(is.null(blue_col)) {
  stop(sprintf("ERROR: BLUE column not found! Available columns: %s\n",
               paste(colnames(pheno), collapse=", ")))
}

cat(sprintf("   - Using BLUE column: '%s'\n", blue_col))

# Use BLUE column directly (no averaging needed!)
pheno$mean_pheno <- pheno[[blue_col]]

# Create a clean phenotype vector
pheno_clean <- pheno[!is.na(mean_pheno) & !is.infinite(mean_pheno),
                      .(`DNA ID`, mean_pheno)]
setnames(pheno_clean, "DNA ID", "sample")

cat(sprintf("   - %d samples with valid BLUE phenotype data\n", nrow(pheno_clean)))
cat(sprintf("   - BLUE range: %.2f - %.2f\n",
            min(pheno_clean$mean_pheno), max(pheno_clean$mean_pheno)))

# ============================================================================
# 2. ENSURE HIGH AND LOW BULKS ARE VALID
# ============================================================================

cat("\n3. Loading bulk sample lists...\n")

# Find bulk list files (handle inconsistent naming)
# Try different patterns
bulk_patterns <- c(
  paste0(TRAIT, ".High.list"),
  paste0(TRAIT, ".high.list"),
  paste0("D2", TRAIT, ".High.list"),
  paste0("D2", TRAIT, ".high.list"),
  # Handle special case for FruitWT vs FruitWt
  paste0(gsub("Wt", "WT", TRAIT), ".High.list"),
  paste0(gsub("Wt", "WT", TRAIT), ".high.list"),
  # Handle special case for D2Fruit vs DTfruit
  paste0("D2", gsub("DT", "", TRAIT), ".high.list"),
  paste0("D2", gsub("DT", "", TRAIT), ".High.list")
)

high_file <- NULL
for(pattern in bulk_patterns) {
  if(file.exists(pattern)) {
    high_file <- pattern
    break
  }
}

if(is.null(high_file)) {
  stop(sprintf("ERROR: Could not find HIGH bulk list file for trait: %s\n", TRAIT))
}

# Same for LOW bulk
bulk_patterns_low <- gsub("High|high", "Low", bulk_patterns)
bulk_patterns_low <- gsub("High|high", "low", bulk_patterns_low)

low_file <- NULL
for(pattern in bulk_patterns_low) {
  if(file.exists(pattern)) {
    low_file <- pattern
    break
  }
}

if(is.null(low_file)) {
  stop(sprintf("ERROR: Could not find LOW bulk list file for trait: %s\n", TRAIT))
}

# Load bulk lists
high <- scan(high_file, what = "character", quiet = TRUE)
low  <- scan(low_file, what = "character", quiet = TRUE)

cat(sprintf("   - HIGH bulk file: %s\n", high_file))
cat(sprintf("   - LOW bulk file: %s\n", low_file))

# Filter to samples present in genotype data
high <- intersect(high, colnames(dos))
low  <- intersect(low, colnames(dos))

cat(sprintf("   - High bulk: %d samples\n", length(high)))
cat(sprintf("   - Low bulk: %d samples\n", length(low)))

if(length(high) < 5 || length(low) < 5) {
  stop("ERROR: Fewer than 5 samples in high or low bulk after filtering!")
}

# ============================================================================
# 3. COMPUTE BULK MEANS AND ΔDS PER MARKER
# ============================================================================

cat("\n4. Computing bulk statistics...\n")
dos$high_mean <- rowMeans(dos[, ..high], na.rm = TRUE)
dos$low_mean  <- rowMeans(dos[, ..low], na.rm = TRUE)
dos$deltaDS   <- dos$high_mean - dos$low_mean
dos$abs_deltaDS <- abs(dos$deltaDS)

cat(sprintf("   - Mean deltaDS: %.4f\n", mean(dos$deltaDS, na.rm=TRUE)))
cat(sprintf("   - Max |deltaDS|: %.4f\n", max(dos$abs_deltaDS, na.rm=TRUE)))

# ============================================================================
# 4. CALCULATE R² AND PVE FOR EACH MARKER
# ============================================================================

cat("\n5. Calculating R² and PVE for each marker...\n")

# Get common samples between genotype and phenotype data
common_samples <- intersect(colnames(dos), pheno_clean$sample)
cat(sprintf("   - Found %d samples with both genotype and phenotype\n", length(common_samples)))

# Extract genotype matrix for common samples
geno_cols <- which(colnames(dos) %in% common_samples)
marker_col <- which(colnames(dos) == "Marker")

# Match phenotypes to genotype column order
pheno_matched <- pheno_clean[match(colnames(dos)[geno_cols], pheno_clean$sample)]

# Calculate R² for each marker
cat(sprintf("   - Computing regression for %d markers...\n", nrow(dos)))

r2_results <- lapply(1:nrow(dos), function(i) {
  if(i %% 5000 == 0) cat(sprintf("     Progress: %d / %d\n", i, nrow(dos)))

  # Get genotype for this marker
  geno <- as.numeric(dos[i, ..geno_cols])
  pheno_vec <- pheno_matched$mean_pheno

  # Remove missing data
  valid <- !is.na(geno) & !is.na(pheno_vec)

  if(sum(valid) < 10) {
    return(list(r2 = NA, pvalue = NA, beta = NA, n_samples = sum(valid)))
  }

  geno_clean <- geno[valid]
  pheno_clean_vec <- pheno_vec[valid]

  # Check for variance
  if(sd(geno_clean) == 0) {
    return(list(r2 = NA, pvalue = NA, beta = NA, n_samples = sum(valid)))
  }

  # Linear regression: phenotype ~ genotype
  tryCatch({
    fit <- lm(pheno_clean_vec ~ geno_clean)
    r2 <- summary(fit)$r.squared
    pval <- summary(fit)$coefficients[2, 4]  # p-value for genotype effect
    beta <- coef(fit)[2]  # effect size

    list(r2 = r2, pvalue = pval, beta = beta, n_samples = sum(valid))
  }, error = function(e) {
    list(r2 = NA, pvalue = NA, beta = NA, n_samples = sum(valid))
  })
})

# Add results to dos table
dos$R2 <- sapply(r2_results, function(x) x$r2)
dos$pvalue <- sapply(r2_results, function(x) x$pvalue)
dos$beta <- sapply(r2_results, function(x) x$beta)
dos$n_samples <- sapply(r2_results, function(x) x$n_samples)

# Calculate PVE (same as R² in simple linear regression)
dos$PVE <- dos$R2 * 100  # Convert to percentage

cat(sprintf("   - Mean R² across all markers: %.4f\n", mean(dos$R2, na.rm=TRUE)))
cat(sprintf("   - Max R² found: %.4f (%.2f%% PVE)\n",
            max(dos$R2, na.rm=TRUE), max(dos$PVE, na.rm=TRUE)))

# ============================================================================
# 5. APPLY MULTIPLE TESTING CORRECTION
# ============================================================================

cat("\n6. Applying multiple testing correction...\n")
dos$pvalue_adj <- p.adjust(dos$pvalue, method = "BH")
dos$significant_005 <- !is.na(dos$pvalue_adj) & dos$pvalue_adj < 0.05
dos$significant_001 <- !is.na(dos$pvalue_adj) & dos$pvalue_adj < 0.01

cat(sprintf("   - Significant at FDR < 0.05: %d markers\n", sum(dos$significant_005, na.rm=TRUE)))
cat(sprintf("   - Significant at FDR < 0.01: %d markers\n", sum(dos$significant_001, na.rm=TRUE)))

# ============================================================================
# 6. SMOOTHING ΔDS AND R² USING SLIDING WINDOW
# ============================================================================

cat("\n7. Smoothing statistics with sliding window...\n")
window_size <- 101  # odd number for symmetric window

dos$deltaDS_smooth <- rollapply(dos$deltaDS,
                                width = window_size,
                                FUN = mean,
                                align = "center",
                                fill = NA,
                                na.rm = TRUE)

dos$R2_smooth <- rollapply(dos$R2,
                          width = window_size,
                          FUN = mean,
                          align = "center",
                          fill = NA,
                          na.rm = TRUE)

# ============================================================================
# 7. PARSE CHROMOSOME AND POSITION
# ============================================================================

cat("\n8. Parsing genomic coordinates...\n")
if(!("Chrom" %in% colnames(dos))) {
  dos$Chrom <- sub("_.*","", dos$Marker)
}
if(!("Position" %in% colnames(dos))) {
  dos$Position <- as.integer(sub(".*_","", dos$Marker))
}

# Order chromosomes properly
dos$Chrom <- factor(dos$Chrom, levels = unique(dos$Chrom[order(as.numeric(gsub("\\D", "", dos$Chrom)))]))

# ============================================================================
# 8. IDENTIFY QTL REGIONS BASED ON MULTIPLE CRITERIA
# ============================================================================

cat("\n9. Identifying QTL regions...\n")

# Define candidate QTL using combined criteria:
# 1. Top 1% of smoothed |ΔDS|
# 2. R² > 0.1 (explains >10% variance)
# 3. FDR-adjusted p-value < 0.05

threshold_deltaDS <- quantile(abs(dos$deltaDS_smooth), 0.99, na.rm = TRUE)
threshold_R2 <- 0.1

dos$candidate_deltaDS <- abs(dos$deltaDS_smooth) >= threshold_deltaDS
dos$candidate_R2 <- !is.na(dos$R2) & dos$R2 >= threshold_R2
dos$candidate_pval <- !is.na(dos$significant_005) & dos$significant_005

# Combined QTL candidates (any of the criteria)
dos$candidate_QTL <- dos$candidate_deltaDS | dos$candidate_R2 | dos$candidate_pval

cat(sprintf("   - Candidate markers by ΔDS: %d\n", sum(dos$candidate_deltaDS, na.rm=TRUE)))
cat(sprintf("   - Candidate markers by R² > 0.1: %d\n", sum(dos$candidate_R2, na.rm=TRUE)))
cat(sprintf("   - Candidate markers by p-value: %d\n", sum(dos$candidate_pval, na.rm=TRUE)))
cat(sprintf("   - Total unique candidate markers: %d\n", sum(dos$candidate_QTL, na.rm=TRUE)))

# ============================================================================
# 9. DEFINE QTL REGIONS
# ============================================================================

setDT(dos)
setorder(dos, Chrom, Position)

dos[, qtl_region := rleid(candidate_QTL)]

qtl_regions <- dos[candidate_QTL == TRUE, .(
  start = min(Position),
  end = max(Position),
  n_markers = .N,
  peak_deltaDS = deltaDS_smooth[which.max(abs(deltaDS_smooth))],
  mean_deltaDS = mean(deltaDS_smooth, na.rm=TRUE),
  max_R2 = max(R2, na.rm=TRUE),
  mean_R2 = mean(R2, na.rm=TRUE),
  max_PVE = max(PVE, na.rm=TRUE),
  mean_PVE = mean(PVE, na.rm=TRUE),
  min_pvalue = min(pvalue_adj, na.rm=TRUE)
), by = .(Chrom, qtl_region)]

# Calculate QTL region size
qtl_regions$size_bp <- qtl_regions$end - qtl_regions$start + 1
qtl_regions$size_mb <- qtl_regions$size_bp / 1e6

# Sort by R² and ΔDS
qtl_regions <- qtl_regions[order(-max_R2, -abs(peak_deltaDS))]

cat("\n10. QTL regions identified:\n")
if(nrow(qtl_regions) > 0) {
  print(qtl_regions[, .(Chrom, start, end, size_mb, n_markers,
                        peak_deltaDS, max_R2, max_PVE, min_pvalue)])
} else {
  cat("   No QTL regions identified with current criteria\n")
}

# Save QTL regions to file
output_dir <- paste0("results_", TRAIT)
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

fwrite(qtl_regions, file.path(output_dir, paste0("QTL_regions_summary_", TRAIT, ".txt")), sep="\t")
cat(sprintf("\n    Saved to: %s/QTL_regions_summary_%s.txt\n", output_dir, TRAIT))

# ============================================================================
# 10. VISUALIZATION: R² MANHATTAN PLOT
# ============================================================================

cat("\n11. Creating visualization plots...\n")

# Calculate genome position for Manhattan plot
dos <- dos[order(Chrom, Position)]
dos[, chrom_num := as.numeric(factor(Chrom, levels = unique(Chrom)))]
dos[, genome_pos := Position + cumsum(c(0, as.numeric(diff(chrom_num) != 0))) * 1e8]

# Chromosome boundaries for x-axis
chrom_bounds <- dos[, .(
  start = min(genome_pos),
  end = max(genome_pos),
  mid = (min(genome_pos) + max(genome_pos)) / 2
), by = Chrom]

# Plot 1: R² Manhattan Plot
pdf(file.path(output_dir, paste0("R2_manhattan_plot_", TRAIT, ".pdf")), width = 14, height = 6)
p1 <- ggplot(dos, aes(x = genome_pos, y = R2)) +
  geom_point(aes(color = Chrom), size = 0.8, alpha = 0.6) +
  geom_point(data = dos[candidate_QTL == TRUE],
             aes(x = genome_pos, y = R2),
             color = "red", size = 1.5) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray50") +
  scale_x_continuous(breaks = chrom_bounds$mid,
                     labels = chrom_bounds$Chrom) +
  labs(x = "Chromosome",
       y = expression(R^2),
       title = paste0("QTL Analysis: ", TRAIT, " - Variance Explained by Each Marker")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "gray90")
  )
print(p1)
dev.off()
cat(sprintf("    - Saved: %s/R2_manhattan_plot_%s.pdf\n", output_dir, TRAIT))

# Plot 2: Smoothed R² by chromosome
pdf(file.path(output_dir, paste0("R2_smoothed_by_chromosome_", TRAIT, ".pdf")), width = 14, height = 8)
p2 <- ggplot(dos, aes(x = Position, y = R2_smooth, color = Chrom)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "blue") +
  facet_wrap(~Chrom, scales = "free_x", ncol = 4) +
  labs(x = "Position (bp)",
       y = expression(paste("Smoothed ", R^2)),
       title = paste0(TRAIT, " - Smoothed R² along genome (window = ", window_size, " markers)")) +
  theme_classic() +
  theme(legend.position = "none")
print(p2)
dev.off()
cat(sprintf("    - Saved: %s/R2_smoothed_by_chromosome_%s.pdf\n", output_dir, TRAIT))

# Plot 3: ΔDS vs R² scatter plot
pdf(file.path(output_dir, paste0("DeltaDS_vs_R2_scatter_", TRAIT, ".pdf")), width = 8, height = 6)

# Try hex plot if hexbin is available, otherwise use point density
if(requireNamespace("hexbin", quietly = TRUE)) {
  p3 <- ggplot(dos[!is.na(R2)], aes(x = abs(deltaDS_smooth), y = R2)) +
    geom_hex(bins = 50) +
    geom_point(data = dos[candidate_QTL == TRUE & !is.na(R2)],
               aes(x = abs(deltaDS_smooth), y = R2),
               color = "red", size = 2, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma") +
    labs(x = expression(paste("|", Delta, "DS| (smoothed)")),
         y = expression(R^2),
         title = paste0(TRAIT, " - Relationship between DeltaDS and Variance Explained"),
         fill = "Count") +
    theme_classic() +
    theme(legend.position = "right")
} else {
  p3 <- ggplot(dos[!is.na(R2)], aes(x = abs(deltaDS_smooth), y = R2)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_point(data = dos[candidate_QTL == TRUE & !is.na(R2)],
               aes(x = abs(deltaDS_smooth), y = R2),
               color = "red", size = 2, alpha = 0.7) +
    labs(x = expression(paste("|", Delta, "DS| (smoothed)")),
         y = expression(R^2),
         title = paste0(TRAIT, " - Relationship between DeltaDS and Variance Explained")) +
    theme_classic()
}
print(p3)
dev.off()
cat(sprintf("    - Saved: %s/DeltaDS_vs_R2_scatter_%s.pdf\n", output_dir, TRAIT))

# Plot 4: -log10(p-value) Manhattan plot
pdf(file.path(output_dir, paste0("Pvalue_manhattan_plot_", TRAIT, ".pdf")), width = 14, height = 6)
dos$neglog10p <- -log10(dos$pvalue_adj)
dos$neglog10p[is.infinite(dos$neglog10p)] <- max(dos$neglog10p[!is.infinite(dos$neglog10p)], na.rm=TRUE)

p4 <- ggplot(dos[!is.na(neglog10p)], aes(x = genome_pos, y = neglog10p)) +
  geom_point(aes(color = Chrom), size = 0.8, alpha = 0.6) +
  geom_point(data = dos[candidate_QTL == TRUE & !is.na(neglog10p)],
             aes(x = genome_pos, y = neglog10p),
             color = "red", size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = chrom_bounds$mid,
                     labels = chrom_bounds$Chrom) +
  labs(x = "Chromosome",
       y = expression(-log[10](p[adj])),
       title = paste0(TRAIT, " - QTL Analysis: Statistical Significance")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "gray90")
  )
print(p4)
dev.off()
cat(sprintf("    - Saved: %s/Pvalue_manhattan_plot_%s.pdf\n", output_dir, TRAIT))

# Plot 5: Combined plot (ΔDS and R² together)
pdf(file.path(output_dir, paste0("Combined_DeltaDS_R2_plot_", TRAIT, ".pdf")), width = 14, height = 10)

# Normalize both to 0-1 scale for comparison
dos$deltaDS_norm <- abs(dos$deltaDS_smooth) / max(abs(dos$deltaDS_smooth), na.rm=TRUE)
dos$R2_norm <- dos$R2_smooth / max(dos$R2_smooth, na.rm=TRUE)

# Create long format manually using rbind instead of pivot_longer
dos_delta <- dos[!is.na(deltaDS_norm), .(genome_pos, Chrom, value = deltaDS_norm,
                                          metric = "deltaDS", candidate_QTL)]
dos_r2 <- dos[!is.na(R2_norm), .(genome_pos, Chrom, value = R2_norm,
                                 metric = "R2", candidate_QTL)]
dos_long <- rbind(dos_delta, dos_r2)

p5 <- ggplot(dos_long, aes(x = genome_pos, y = value, color = metric)) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  geom_point(data = dos_long[candidate_QTL == TRUE],
             aes(x = genome_pos, y = value),
             color = "red", size = 1) +
  scale_color_manual(values = c("deltaDS" = "blue", "R2" = "darkgreen"),
                     labels = c("|DeltaDS|", "R2")) +
  scale_x_continuous(breaks = chrom_bounds$mid,
                     labels = chrom_bounds$Chrom) +
  labs(x = "Chromosome",
       y = "Normalized Value (0-1)",
       title = paste0(TRAIT, " - Combined DeltaDS and R2 QTL Analysis"),
       color = "Metric") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )
print(p5)
dev.off()
cat(sprintf("    - Saved: %s/Combined_DeltaDS_R2_plot_%s.pdf\n", output_dir, TRAIT))

# ============================================================================
# 11. SUMMARY STATISTICS AND TOP QTL
# ============================================================================

cat("\n============================================\n")
cat("SUMMARY STATISTICS\n")
cat("============================================\n\n")

cat(sprintf("Top 10 QTL regions by R² for %s:\n", TRAIT))
if(nrow(qtl_regions) > 0) {
  print(qtl_regions[1:min(10, nrow(qtl_regions)),
                    .(Chrom, start, end, size_mb, n_markers,
                      max_R2, max_PVE, peak_deltaDS, min_pvalue)])
}

cat(sprintf("\n\nTop 20 individual markers by R² for %s:\n", TRAIT))
top_markers <- dos[order(-R2)][1:20, .(Marker, Chrom, Position, R2, PVE,
                                        deltaDS_smooth, pvalue_adj, beta)]
print(top_markers)
fwrite(top_markers, file.path(output_dir, paste0("Top_markers_by_R2_", TRAIT, ".txt")), sep="\t")

cat("\n\nOverall statistics:\n")
cat(sprintf("  Trait: %s\n", TRAIT))
cat(sprintf("  Total markers analyzed: %d\n", nrow(dos)))
cat(sprintf("  Mean R²: %.4f\n", mean(dos$R2, na.rm=TRUE)))
cat(sprintf("  Median R²: %.4f\n", median(dos$R2, na.rm=TRUE)))
cat(sprintf("  Max R²: %.4f (%.2f%% PVE)\n", max(dos$R2, na.rm=TRUE), max(dos$PVE, na.rm=TRUE)))
cat(sprintf("  Markers with R² > 0.05: %d (%.2f%%)\n",
            sum(dos$R2 > 0.05, na.rm=TRUE),
            100 * sum(dos$R2 > 0.05, na.rm=TRUE) / sum(!is.na(dos$R2))))
cat(sprintf("  Markers with R² > 0.10: %d (%.2f%%)\n",
            sum(dos$R2 > 0.10, na.rm=TRUE),
            100 * sum(dos$R2 > 0.10, na.rm=TRUE) / sum(!is.na(dos$R2))))
cat(sprintf("  Significant markers (FDR < 0.05): %d\n", sum(dos$significant_005, na.rm=TRUE)))
cat(sprintf("  QTL regions identified: %d\n", nrow(qtl_regions)))

# Save complete results
fwrite(dos, file.path(output_dir, paste0("Complete_QTL_results_with_R2_PVE_", TRAIT, ".txt")), sep="\t")

cat("\n============================================\n")
cat("Analysis complete!\n")
cat("============================================\n")
cat(sprintf("\nOutput directory: %s/\n", output_dir))
cat("\nOutput files:\n")
cat(sprintf("  - Complete_QTL_results_with_R2_PVE_%s.txt\n", TRAIT))
cat(sprintf("  - QTL_regions_summary_%s.txt\n", TRAIT))
cat(sprintf("  - Top_markers_by_R2_%s.txt\n", TRAIT))
cat(sprintf("  - R2_manhattan_plot_%s.pdf\n", TRAIT))
cat(sprintf("  - R2_smoothed_by_chromosome_%s.pdf\n", TRAIT))
cat(sprintf("  - DeltaDS_vs_R2_scatter_%s.pdf\n", TRAIT))
cat(sprintf("  - Pvalue_manhattan_plot_%s.pdf\n", TRAIT))
cat(sprintf("  - Combined_DeltaDS_R2_plot_%s.pdf\n", TRAIT))
cat("\n")
