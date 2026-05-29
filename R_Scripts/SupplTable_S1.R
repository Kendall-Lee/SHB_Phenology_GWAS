## SupplTable_S1.R
## Build Supplementary Table S1: Full 88-marker GWAS table
## BLUE_exc.25 (yr.23+yr.24 only), Linear reference, GWASpoly LOCO
## Meff significance threshold: -log10(p) >= 5.52 (alpha = 0.05)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

MEFF <- 5.52

# ── Load raw marker table ─────────────────────────────────────────────────────
d <- fread("/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/All_Markers/LINEAR/BLUE_exc.25/TopMarkers_LINEAR_BLUE_exc25.csv")

# ── Derived columns ───────────────────────────────────────────────────────────

# 1. Locus window (1-Mb bins, label as "Chr.05: 47-48 Mb")
d[, Locus := paste0(Chrom, ": ", floor(Position_Mb), "-", floor(Position_Mb)+1, " Mb")]

# 2. Clean primary trait name
d[, Primary_Trait := sub("_BLUE.*", "", Trait_exc25)]

# 3. Dosage distribution string
d[, Dosage_Dist := paste(N_D0, N_D1, N_D2, N_D3, N_D4, sep = ":")]

# 4. Per-year scores for the primary trait (for within-study replication)
d[, yr23_primary := fcase(
  Primary_Trait == "DTFlower",    DTFlower_yr.23,
  Primary_Trait == "DTFruit",     DTFruit_yr.23,
  Primary_Trait == "Flow2Fruit",  Flow2Fruit_yr.23,
  Primary_Trait == "FruitWeight", FruitWeight_yr.23
)]
d[, yr24_primary := fcase(
  Primary_Trait == "DTFlower",    DTFlower_yr.24,
  Primary_Trait == "DTFruit",     DTFruit_yr.24,
  Primary_Trait == "Flow2Fruit",  Flow2Fruit_yr.24,
  Primary_Trait == "FruitWeight", FruitWeight_yr.24
)]

# 5. Within-study replication flag
d[, Replicated := fcase(
  yr23_primary >= MEFF & yr24_primary >= MEFF, "yr.23 + yr.24",
  yr23_primary >= MEFF & yr24_primary <  MEFF, "yr.23 only",
  yr23_primary <  MEFF & yr24_primary >= MEFF, "yr.24 only",
  default = "Neither"
)]

# ── Sort: by chromosome, then by Score descending ────────────────────────────
setorder(d, Chrom, -Score_exc25)

# ── Select and rename columns for journal output ─────────────────────────────
out <- d[, .(
  Locus,
  Marker,
  Chr            = Chrom,
  `Position (Mb)` = round(Position_Mb, 3),
  `Position (bp)` = Position,
  `Variant type`  = VarType,
  `Primary trait` = Primary_Trait,
  `Dosage model`  = Model_exc25,
  `-log10(p) [BLUE exc.25]` = round(Score_exc25, 3),
  `PVE (%) [BLUE exc.25]`   = round(PVE_exc25, 2),
  `Alt allele freq`         = round(AF, 3),
  `N genotyped`             = N_called,
  `Dosage dist (D0:D1:D2:D3:D4)` = Dosage_Dist,
  `DTFlower score [BLUE exc.25]`  = round(DTFlower_BLUE, 3),
  `DTFruit score [BLUE exc.25]`   = round(DTFruit_BLUE, 3),
  `Flow2Fruit score [BLUE exc.25]`= round(Flow2Fruit_BLUES, 3),
  `FruitWeight score [BLUE exc.25]` = round(FruitWeight_BLUE, 3),
  `yr.23 score (primary trait)` = round(yr23_primary, 3),
  `yr.24 score (primary trait)` = round(yr24_primary, 3),
  `Within-study replication`    = Replicated
)]

# ── Write plain CSV ───────────────────────────────────────────────────────────
out_dir <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Figures"
csv_path <- file.path(out_dir, "SupplTable_S1_TopMarkers_LINEAR_BLUEexc25.csv")
fwrite(out, csv_path)
cat("CSV saved:", csv_path, "\n")

# ── Write formatted Excel ─────────────────────────────────────────────────────
wb <- createWorkbook()

# ---- Data sheet ----
addWorksheet(wb, "Table S1")

# Header style
hdr_style <- createStyle(
  fontSize = 10, fontColour = "white", fgFill = "#0D3B6E",
  halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE,
  border = "Bottom", borderColour = "white"
)
# Body styles
body_style <- createStyle(fontSize = 9, halign = "left", valign = "center")
num_style  <- createStyle(fontSize = 9, halign = "right", numFmt = "0.000")
pve_style  <- createStyle(fontSize = 9, halign = "right", numFmt = "0.00")
int_style  <- createStyle(fontSize = 9, halign = "right", numFmt = "0")

# Replication highlight styles
rep_both <- createStyle(fontSize = 9, fgFill = "#C6EFC8", halign = "center")  # green
rep_one  <- createStyle(fontSize = 9, fgFill = "#FFF2CC", halign = "center")  # yellow
rep_none <- createStyle(fontSize = 9, fgFill = "#F4F4F4", halign = "center")  # grey

# Significance highlight for score columns (above MEFF)
sig_style  <- createStyle(fontSize = 9, halign = "right", numFmt = "0.000",
                           fontColour = "#0D3B6E", textDecoration = "bold")

writeData(wb, "Table S1", out, startRow = 1, headerStyle = hdr_style,
          borders = "columns", borderStyle = "thin")

n_rows <- nrow(out) + 1  # +1 for header

# Apply base body styles
addStyle(wb, "Table S1", body_style, rows = 2:n_rows, cols = c(1, 2, 3, 6, 7, 8, 13, 20),
         gridExpand = TRUE)
addStyle(wb, "Table S1", num_style,  rows = 2:n_rows, cols = c(4, 9, 11, 14, 15, 16, 17, 18, 19),
         gridExpand = TRUE)
addStyle(wb, "Table S1", pve_style,  rows = 2:n_rows, cols = 10, gridExpand = TRUE)
addStyle(wb, "Table S1", int_style,  rows = 2:n_rows, cols = c(5, 12), gridExpand = TRUE)

# Alternating row shading for locus groups
locus_ids <- rle(out$Locus)$lengths
locus_rows <- cumsum(c(1, head(locus_ids, -1)))
for (i in seq_along(locus_rows)) {
  if (i %% 2 == 0) {
    r_start <- locus_rows[i] + 1   # +1 for header
    r_end   <- r_start + locus_ids[i] - 1
    addStyle(wb, "Table S1",
             createStyle(fgFill = "#EBF2FA", fontSize = 9),
             rows = r_start:r_end, cols = 1:ncol(out),
             gridExpand = TRUE, stack = TRUE)
  }
}

# Bold significant score cells
for (col_idx in c(9, 14, 15, 16, 17, 18, 19)) {
  col_name <- names(out)[col_idx]
  sig_rows <- which(out[[col_name]] >= MEFF) + 1  # +1 for header
  if (length(sig_rows) > 0) {
    addStyle(wb, "Table S1", sig_style, rows = sig_rows, cols = col_idx,
             gridExpand = FALSE, stack = TRUE)
  }
}

# Replication column color
rep_col <- which(names(out) == "Within-study replication")
for (i in 2:n_rows) {
  val <- out[i-1, `Within-study replication`]
  s <- if (val == "yr.23 + yr.24") rep_both else if (val == "Neither") rep_none else rep_one
  addStyle(wb, "Table S1", s, rows = i, cols = rep_col, stack = TRUE)
}

# Column widths
setColWidths(wb, "Table S1", cols = 1,        widths = 16)  # Locus
setColWidths(wb, "Table S1", cols = 2,        widths = 20)  # Marker
setColWidths(wb, "Table S1", cols = 3,        widths = 8)   # Chr
setColWidths(wb, "Table S1", cols = 4,        widths = 11)  # Position Mb
setColWidths(wb, "Table S1", cols = 5,        widths = 12)  # Position bp
setColWidths(wb, "Table S1", cols = 6,        widths = 10)  # VarType
setColWidths(wb, "Table S1", cols = 7,        widths = 12)  # Primary trait
setColWidths(wb, "Table S1", cols = 8,        widths = 12)  # Model
setColWidths(wb, "Table S1", cols = 9,        widths = 14)  # Score
setColWidths(wb, "Table S1", cols = 10,       widths = 13)  # PVE
setColWidths(wb, "Table S1", cols = 11,       widths = 12)  # AF
setColWidths(wb, "Table S1", cols = 12,       widths = 10)  # N geno
setColWidths(wb, "Table S1", cols = 13,       widths = 22)  # Dosage dist
setColWidths(wb, "Table S1", cols = 14:17,    widths = 14)  # trait scores
setColWidths(wb, "Table S1", cols = 18:19,    widths = 16)  # yr scores
setColWidths(wb, "Table S1", cols = 20,       widths = 18)  # replication

# Freeze top row
freezePane(wb, "Table S1", firstRow = TRUE)

# ---- Legend sheet ----
addWorksheet(wb, "Legend")
legend_text <- data.frame(
  Column = c(
    "Locus",
    "Marker",
    "Chr",
    "Position (Mb)",
    "Position (bp)",
    "Variant type",
    "Primary trait",
    "Dosage model",
    "-log10(p) [BLUE exc.25]",
    "PVE (%) [BLUE exc.25]",
    "Alt allele freq",
    "N genotyped",
    "Dosage dist (D0:D1:D2:D3:D4)",
    "DTFlower/DTFruit/Flow2Fruit/FruitWeight score [BLUE exc.25]",
    "yr.23 score (primary trait)",
    "yr.24 score (primary trait)",
    "Within-study replication"
  ),
  Description = c(
    "1-Mb genomic window used for locus deduplication. Markers in the same window represent a single GWAS locus.",
    "Marker ID in Chr_Position format (Suziblue hap1 linear reference coordinates).",
    "Chromosome.",
    "Marker position in megabases (Mb).",
    "Marker position in base pairs (bp).",
    "SNP = single nucleotide polymorphism; Indel = insertion/deletion <50 bp; Other = complex/multi-allelic.",
    "Trait for which this marker achieved highest -log10(p) in the BLUE_exc.25 primary analysis.",
    "GWASpoly dosage model selected by lowest p-value: additive (0-4 linear), 1-dom-ref (0 vs 1-4), 1-dom-alt (0-3 vs 4), 2-dom-ref (0-1 vs 2-4), 2-dom-alt (0-2 vs 3-4).",
    "Association score (-log10 p-value) from GWASpoly LOCO scan using BLUE_exc.25 phenotype (years 2023+2024 only). Significance threshold: 5.52 (Meff correction, alpha=0.05). Bold values exceed threshold.",
    "Proportion of phenotypic variance explained (PVE) estimated by fit.QTL() in GWASpoly using the best-model scan object.",
    "Alternate allele frequency across all genotyped accessions (0-4 dosage scale).",
    "Number of accessions with valid genotype calls at this marker.",
    "Count of accessions with 0, 1, 2, 3, or 4 copies of the alternate allele (tetraploid dosage 0-4).",
    "Association score for each of the four traits in the BLUE_exc.25 scan. Bold values exceed the Meff significance threshold (5.52).",
    "Association score for the primary trait in the year 2023 individual BLUE scan (within-study replication check).",
    "Association score for the primary trait in the year 2024 individual BLUE scan (within-study replication check).",
    "Green = significant (>=5.52) in both yr.23 and yr.24; Yellow = significant in one year only; Grey = not replicated in either year."
  )
)
writeData(wb, "Legend", legend_text)
setColWidths(wb, "Legend", cols = 1, widths = 45)
setColWidths(wb, "Legend", cols = 2, widths = 90)
addStyle(wb, "Legend", createStyle(fontColour = "white", fgFill = "#0D3B6E",
                                    textDecoration = "bold", fontSize = 10),
         rows = 1, cols = 1:2)

# ── Save ──────────────────────────────────────────────────────────────────────
xlsx_path <- file.path(out_dir, "SupplTable_S1_TopMarkers_LINEAR_BLUEexc25.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)
cat("Excel saved:", xlsx_path, "\n")

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n--- Table S1 summary ---\n")
cat("Total markers:", nrow(out), "\n")
cat("Unique loci (1-Mb windows):", length(unique(out$Locus)), "\n")
cat("Primary trait counts:\n")
print(table(out$`Primary trait`))
cat("\nReplication status:\n")
print(table(out$`Within-study replication`))
cat("\nDosage model counts:\n")
print(table(out$`Dosage model`))
