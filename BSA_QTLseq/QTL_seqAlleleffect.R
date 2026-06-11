# =========================
# Load libraries
# =========================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# =========================
# Load data
# =========================
setwd("/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq")

geno <- fread("SHB_PAN_LR.DP.1.maxmiss.4.minminor.3")
pheno <- fread("DTFruit_SHB_allPheno.txt")
top_markers <- fread("results_MultiPheno/DTFruit_SHB/Top50_markers_DTFruit_SHB.txt")
setnames(pheno, old = "DNA ID", new = "DNA_ID")
str(pheno)
# =========================
# Select marker
# =========================
marker_of_interest <- top_markers$Marker[2]

geno_sub <- geno[Marker == marker_of_interest]
# Extract REF / ALT from genotype
ref_allele <- geno_sub$REF[1]
alt_allele <- geno_sub$ALT[1]

# Extract stats from top_markers
marker_stats <- top_markers[Marker == marker_of_interest]

deltaDS <- marker_stats$deltaDS_smooth
R2_val  <- marker_stats$R2
# Keep only DNA_IDs that exist in phenotype
dna_cols <- intersect(pheno$DNA_ID, colnames(geno_sub))

# =========================
# Reshape genotype → long
# =========================
geno_long <- melt(
  geno_sub,
  id.vars = c("Marker","Chrom","Position","REF","ALT"),
  measure.vars = dna_cols,
  variable.name = "DNA_ID",
  value.name = "Allele"
)

# =========================
# Merge with phenotype
# =========================
df <- merge(geno_long, pheno, by = "DNA_ID", all.x = TRUE)

# =========================
# Bin allele dosage + clean
# =========================
df <- df %>%
  mutate(
    Allele_bin = cut(
      Allele,
      breaks = c(-Inf, 0.5, 1.5, 2.5, 3.5, 4.5, Inf),
      labels = c(0, 1, 2, 3, 4, 5)
    )
  ) %>%
  filter(!is.na(Allele_bin), !is.na(Data_BLUE))

# =========================
# Compute N per bin
# =========================
N_df <- df %>%
  group_by(Allele_bin) %>%
  summarise(
    N = n(),
    y_pos = max(Data_BLUE) + 1,
    .groups = "drop"
  )

# =========================
# Pairwise comparisons (only existing bins)
# =========================
existing_bins <- sort(unique(df$Allele_bin))
comparisons <- combn(existing_bins, 2, simplify = FALSE)

# =========================
# Plot
# =========================
annot_text <- paste0(
  "REF: ", ref_allele,
  " | ALT: ", alt_allele,
  "\nΔDS: ", round(deltaDS, 3),
  " | R²: ", round(R2_val, 3)
)

p <- ggplot(df, aes(x = Allele_bin, y = Data_BLUE, fill = Allele_bin)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_text(
    data = N_df,
    aes(x = Allele_bin, y = y_pos, label = paste0("N=", N)),
    inherit.aes = FALSE,
    size = 3
  ) +
  stat_compare_means(comparisons = comparisons, label = "p.signif") +
  
  # 🔥 ADD THIS BLOCK
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = annot_text,
    hjust = -0.1,
    vjust = 1.5,
    size = 4
  ) +
  
  labs(
    x = "Binned Allele Dosage",
    y = "Data_BLUE",
    title = marker_of_interest
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

p



#####################################
####################################

plot_marker <- function(marker_of_interest, geno, pheno, top_markers) {
  
  # Subset genotype
  geno_sub <- geno[Marker == marker_of_interest]
  if (nrow(geno_sub) == 0) return(NULL)
  
  # Match DNA IDs
  dna_cols <- intersect(pheno$DNA_ID, colnames(geno_sub))
  if (length(dna_cols) == 0) return(NULL)
  
  # Melt to long format
  geno_long <- melt(
    geno_sub,
    id.vars = c("Marker","Chrom","Position","REF","ALT"),
    measure.vars = dna_cols,
    variable.name = "DNA_ID",
    value.name = "Allele"
  )
  
  # Merge phenotype
  df <- merge(geno_long, pheno, by = "DNA_ID", all.x = TRUE)
  
  # Bin + clean
  df <- df %>%
    mutate(
      Allele_bin = cut(
        Allele,
        breaks = c(-Inf, 0.5, 1.5, 2.5, 3.5, 4.5, Inf),
        labels = c(0, 1, 2, 3, 4, 5)
      )
    ) %>%
    filter(!is.na(Allele_bin), !is.na(Data_BLUE))
  
  if (nrow(df) == 0) return(NULL)
  
  # N per bin
  N_df <- df %>%
    group_by(Allele_bin) %>%
    summarise(
      N = n(),
      y_pos = max(Data_BLUE) + 1,
      .groups = "drop"
    )
  
  # Comparisons
  bins <- sort(unique(df$Allele_bin))
  if (length(bins) < 2) return(NULL)
  comparisons <- combn(bins, 2, simplify = FALSE)
  
  # Marker stats
  ref_allele <- geno_sub$REF[1]
  alt_allele <- geno_sub$ALT[1]
  
  marker_stats <- top_markers[Marker == marker_of_interest]
  deltaDS <- marker_stats$deltaDS_smooth
  R2_val  <- marker_stats$R2
  
  annot_text <- sprintf(
    "REF: %s   ALT: %s\nΔDS = %.3f   R² = %.3f",
    ref_allele, alt_allele, deltaDS, R2_val
  )
  
  # Plot
  p <- ggplot(df, aes(x = Allele_bin, y = Data_BLUE, fill = Allele_bin)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_text(
      data = N_df,
      aes(x = Allele_bin, y = y_pos, label = paste0("N=", N)),
      inherit.aes = FALSE,
      size = 3
    ) +
    stat_compare_means(comparisons = comparisons, label = "p.signif") +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annot_text,
      hjust = -0.1,
      vjust = 1.5,
      size = 4
    ) +
    labs(
      x = "Binned Allele Dosage",
      y = "Data_BLUE",
      title = marker_of_interest
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}


dir.create("Marker_Plots", showWarnings = FALSE)



for (m in top_markers$Marker) {
  
  p <- tryCatch(
    plot_marker(m, geno, pheno, top_markers),
    error = function(e) NULL
  )
  
  if (!is.null(p)) {
    ggsave(
      filename = paste0("Marker_Plots/", m, ".png"),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
}
