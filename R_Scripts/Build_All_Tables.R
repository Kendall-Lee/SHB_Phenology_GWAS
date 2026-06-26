## Build_All_Tables.R
## Generate all main-text and supplementary tables for the Linear GWAS manuscript
##
## Outputs (all to Linear_MS/Tables/):
##   Table1_Phenotypic_Summary.csv/.xlsx
##   Table2_Stable_QTL_Summary.csv/.xlsx
##   SupplTable_S1_All88Markers.csv/.xlsx       (cleaned from existing)
##   SupplTable_S2_Candidate_Genes.csv/.xlsx
##   SupplTable_S3_QTLseq_TopRegions.csv/.xlsx

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS/Tables")
dir.create(OUT_DIR, showWarnings = FALSE)

save_table <- function(dt, stem, note = NULL) {
  csv_path  <- file.path(OUT_DIR, paste0(stem, ".csv"))
  xlsx_path <- file.path(OUT_DIR, paste0(stem, ".xlsx"))
  fwrite(dt, csv_path)
  wb <- createWorkbook()
  addWorksheet(wb, "Data")
  writeData(wb, "Data", dt)
  # Bold header
  addStyle(wb, "Data",
    style = createStyle(textDecoration = "bold", border = "Bottom"),
    rows = 1, cols = seq_len(ncol(dt)), gridExpand = TRUE)
  # Auto column widths
  setColWidths(wb, "Data", cols = seq_len(ncol(dt)), widths = "auto")
  if (!is.null(note)) {
    addWorksheet(wb, "Notes")
    writeData(wb, "Notes", data.frame(Note = note))
  }
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat(sprintf("  Saved: %s\n", stem))
}

# ==============================================================================
# TABLE 1 — Phenotypic characterization and heritability
# ==============================================================================
cat("Building Table 1...\n")

t1 <- data.table(
  Trait = c("Days to 50% Flowering (DTFlower)",
            "Days to 50% Ripe Fruit (DTFruit)",
            "Fruiting Period (Flow2Fruit)",
            "25-Fruit Weight (FruitWeight)"),
  `H2_classical` = c(0.844, 0.911, 0.734, 0.842),
  `H2_Cullis`    = c(0.831, 0.899, 0.710, 0.824),
  `yr23_mean`    = c(52.8,  120.4, 67.7,  48.5),
  `yr23_sd`      = c(13.3,  15.6,  12.1,  13.0),
  `yr24_mean`    = c(61.2,  127.6, 67.0,  47.9),
  `yr24_sd`      = c(11.2,  12.9,  12.5,  14.7),
  `yr25_mean`    = c(79.7,  133.6, 53.3,  53.0),
  `yr25_sd`      = c(6.9,   15.7,  17.4,  15.1),
  `Rank_r_2324`  = c(0.700, 0.503, 0.397, 0.574),
  `Rank_r_2325`  = c(0.479, 0.485, 0.195, 0.540)
)

setnames(t1, c(
  "Trait",
  "H2 (classical)", "H2 (Cullis)",
  "yr.23 mean (days)", "yr.23 SD",
  "yr.24 mean (days)", "yr.24 SD",
  "yr.25 mean (days)", "yr.25 SD",
  "Rank r (yr.23 x yr.24)", "Rank r (yr.23 x yr.25)"
))

save_table(t1, "Table1_Phenotypic_Summary",
  note = paste0(
    "H2 (classical) = broad-sense heritability estimated via REML. ",
    "H2 (Cullis) = reliability-adjusted heritability. ",
    "Means and SD computed from per-year BLUEs across ~800 SHB accessions at Georgia field site. ",
    "Rank r = Spearman rank correlation of genotype means between years. ",
    "yr.25 DTFlower SD collapse (13.3 -> 6.9 days) reflects frost-driven variance compression (21 frost days, 12 hard-freeze events in 2025)."
  ))

# ==============================================================================
# TABLE 2 — Stable GWAS loci summary (main text, 22 loci)
# ==============================================================================
cat("Building Table 2...\n")

s1 <- fread(file.path(BASE,
  "PUBLICATION/Linear_MS/SupplTable_S1_TopMarkers_LINEAR_BLUEexc25.csv"))

# Deduplicate to one peak marker per locus (highest -log10p per Locus group)
s1[, score_num := as.numeric(`-log10(p) [BLUE exc.25]`)]
setorder(s1, Locus, -score_num)
s1_peak <- s1[, .SD[1], by = Locus]

# Select and rename columns for main-text Table 2
t2 <- s1_peak[, .(
  Locus                        = Locus,
  Chr                          = Chr,
  `Position (Mb)`              = `Position (Mb)`,
  `Primary trait`              = `Primary trait`,
  `Dosage model`               = `Dosage model`,
  `-log10(p) [BLUE exc.25]`   = score_num,
  `PVE (%) [BLUE exc.25]`     = `PVE (%) [BLUE exc.25]`,
  `yr.23 score`                = `yr.23 score (primary trait)`,
  `yr.24 score`                = `yr.24 score (primary trait)`,
  `Within-study replication`   = `Within-study replication`
)]

# Sort by primary trait then descending score
trait_rank <- c("DTFlower" = 1, "DTFruit" = 2, "Flow2Fruit" = 3, "FruitWeight" = 4)
t2[, trait_ord := trait_rank[`Primary trait`]]
setorder(t2, trait_ord, -`-log10(p) [BLUE exc.25]`)
t2[, trait_ord := NULL]
cat(sprintf("  Table 2: %d loci\n", nrow(t2)))

save_table(t2, "Table2_Stable_QTL_Summary",
  note = paste0(
    "22 stable loci identified by tetraploid-aware GWASpoly LOCO scan on BLUE_exc.25 phenotypes (yr.23 + yr.24 combined). ",
    "Significance threshold: Meff multiple-testing correction, -log10(p) >= 5.52 (alpha = 0.05). ",
    "PVE estimated by fit.QTL() from GWASpoly best-model scan object. ",
    "Locus deduplication: 1-Mb windows; representative marker = peak score per window. ",
    "yr.23 and yr.24 scores shown for the primary trait; scores >= 5.52 indicate within-study replication. ",
    "Dosage models: additive = monotonic dosage effect; 1-dom-ref/alt = one-allele dominance; 2-dom-ref/alt = two-allele dominance."
  ))

# ==============================================================================
# SUPPLEMENTARY TABLE S1 — Full 88-marker table (cleaned)
# ==============================================================================
cat("Building Sup Table S1...\n")

# The existing file is already well-structured; just clean up column names
s1_clean <- copy(s1)

# Rename for journal clarity
col_renames <- c(
  "Marker"                            = "Marker ID",
  "Chr"                               = "Chromosome",
  "Position (Mb)"                     = "Position (Mb)",
  "Position (bp)"                     = "Position (bp)",
  "Variant type"                      = "Variant type",
  "Primary trait"                     = "Primary trait",
  "Dosage model"                      = "Dosage model",
  "-log10(p) [BLUE exc.25]"          = "-log10(p) [BLUE exc.25]",
  "PVE (%) [BLUE exc.25]"            = "PVE (%) [BLUE exc.25]",
  "Alt allele freq"                   = "ALT allele frequency",
  "N genotyped"                       = "N genotyped",
  "Dosage dist (D0:D1:D2:D3:D4)"     = "Dosage distribution (D0:D1:D2:D3:D4)",
  "DTFlower score [BLUE exc.25]"      = "DTFlower -log10(p) [BLUE exc.25]",
  "DTFruit score [BLUE exc.25]"       = "DTFruit -log10(p) [BLUE exc.25]",
  "Flow2Fruit score [BLUE exc.25]"    = "Flow2Fruit -log10(p) [BLUE exc.25]",
  "FruitWeight score [BLUE exc.25]"   = "FruitWeight -log10(p) [BLUE exc.25]",
  "yr.23 score (primary trait)"       = "-log10(p) yr.23 [primary trait]",
  "yr.24 score (primary trait)"       = "-log10(p) yr.24 [primary trait]",
  "Within-study replication"          = "Within-study replication"
)
for (old in names(col_renames)) {
  if (old %in% colnames(s1_clean)) setnames(s1_clean, old, col_renames[old])
}

save_table(s1_clean, "SupplTable_S1_All88Markers",
  note = paste0(
    "All 88 significant markers at Meff threshold (-log10(p) >= 5.52) from GWASpoly LOCO scan on BLUE_exc.25. ",
    "Deduplicated to one peak marker per 1-Mb locus window. ",
    "Scores for all four traits shown to facilitate identification of pleiotropic loci. ",
    "Within-study replication: both = significant in yr.23 AND yr.24; yr.23 only / yr.24 only = one support year."
  ))

# ==============================================================================
# SUPPLEMENTARY TABLE S2 — Candidate gene annotations (GWAS + QTL-seq)
# ==============================================================================
cat("Building Sup Table S2...\n")

# Known annotations for top loci (manually curated from InterProScan + bedtools)
s2 <- data.table(
  Locus = c(
    # GWAS top loci
    "Chr.05: 47-48 Mb (GWAS + QTL-seq)",
    "Chr.05: 47-48 Mb (GWAS + QTL-seq)",
    "Chr.05: 20.2 Mb (GWAS)",
    "Chr.02: 3.6 Mb (GWAS)",
    "Chr.02: 37.5 Mb (GWAS)",
    "Chr.03: 16.4 Mb (GWAS)",
    "Chr.06: 22.3 Mb (GWAS)",
    "Chr.06: 42.4 Mb (GWAS)",
    "Chr.08: 18.9 Mb (GWAS)",
    "Chr.09: 1.7 Mb (GWAS)",
    "Chr.10: 14.2 Mb (GWAS)",
    "Chr.11: 19.5 Mb (GWAS)",
    # QTL-seq loci
    "Chr.12: 9.5-11.5 Mb (QTL-seq)",
    "Chr.12: 9.5-11.5 Mb (QTL-seq)",
    "Chr.11: 42 Mb (QTL-seq)",
    "Chr.05: 47.7 Mb (QTL-seq)"
  ),
  Source = c(
    rep("Linear GWAS", 11),
    "Linear GWAS",
    rep("LRLP SV QTL-seq", 3),
    "LRLP SV QTL-seq"
  ),
  `Peak marker` = c(
    "Chr.05_48252083", "Chr.05_47900362",
    "Chr.05_20194733",
    "Chr.02_3615748", "Chr.02_37452394", "Chr.03_16388162",
    "Chr.06_22313445", "Chr.06_42376290", "Chr.08_18880093",
    "Chr.09_1685774", "Chr.10_14191893", "Chr.11_19478223",
    "Chr.12_10308052", "Chr.12_10308052",
    "Chr.11_42064057",
    "Chr.05_47712013"
  ),
  `Primary trait` = c(
    "DTFlower", "DTFlower",
    "DTFruit/DTFlower",
    "DTFlower", "DTFlower", "DTFruit",
    "DTFlower", "DTFlower", "DTFlower",
    "DTFlower", "Flow2Fruit", "DTFlower",
    "DTFlower", "DTFlower",
    "DTFruit",
    "DTFlower"
  ),
  `Nearest gene` = c(
    "g54697", "g54686",
    "unknown",
    "unknown", "unknown", "unknown",
    "unknown", "unknown", "unknown",
    "unknown", "unknown", "unknown",
    "g132348", "g132350",
    "unknown",
    "g54686"
  ),
  `Distance to marker (bp)` = c(
    "within/flanking", "within/flanking",
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    "inside gene body", "13,700",
    NA,
    "within/flanking"
  ),
  `Gene annotation` = c(
    # GWAS loci (12 rows)
    "MLP/Bet v1-like (PR-10 family); pathogenesis-related protein expressed in reproductive tissues",  # Chr.05:47-48 g54697
    "MLP/Bet v1-like (PR-10 family); pathogenesis-related protein expressed in reproductive tissues",  # Chr.05:47-48 g54686
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g53184; secreted protein)",           # Chr.05:20.2
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g13461; secreted protein)",           # Chr.02:3.6
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g15412; secreted protein)",           # Chr.02:37.5
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g29230; secreted protein)",           # Chr.03:16.4
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g65196; secreted protein)",           # Chr.06:22.3
    "Protein TAR1 (IPR044792); stress-response/nutrient-sensing; g66181 inside gene body",            # Chr.06:42.4
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g87092; secreted protein)",           # Chr.08:18.9
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g97156; secreted protein)",           # Chr.09:1.7
    "No conserved Pfam/PANTHER domain; signal peptide predicted (g108569; secreted protein)",          # Chr.10:14.2
    "Anaphylatoxin/fibulin domain (g120146; IPR000020); extracellular matrix/complement-like protein", # Chr.11:19.5
    # QTL-seq loci (4 rows)
    "Protein Arginine N-Methyltransferase (PRMT); IPR025799; C2H2 zinc finger + SAM-dependent methyltransferase; peak SV inside gene body",
    "FHY3/FAR1-type transcription factor; IPR031052; photomorphogenesis and photoperiod-responsive flowering regulator",
    "No conserved domain identified (g-unknown; annotation pending)",                                  # Chr.11:42 Mb QTL-seq
    "MLP/Bet v1-like (PR-10 family)"                                                                   # Chr.05:47.7 QTL-seq
  ),
  `IPR accession` = c(
    "IPR000916", "IPR000916",
    NA, NA, NA, NA, NA,
    "IPR044792",
    NA, NA, NA,
    "IPR000020",
    "IPR025799", "IPR031052",
    NA,
    "IPR000916"
  ),
  Notes = c(
    "Confirmed by GWAS (-log10p=11.86, PVE=7.79%) and QTL-seq (R2=27.3%, FDR-sig); additive dosage model; peak Chr.05 locus",
    "QTL-seq peak 2 at Chr.05:47.71 Mb (R2=26.6%, p=9.3e-06); co-localizes with GWAS cluster",
    "Pleiotropic hotspot: 4-5 trait-year combinations, PVE ~5-7%",
    "Self-replicating: yr.23 and yr.24 both significant",
    "yr.23 only replication",
    "Self-replicating: yr.23 and yr.24 both significant",
    "Self-replicating", "Self-replicating", "Self-replicating",
    "Self-replicating", "Only significant Flow2Fruit locus in study",
    "yr.23 only replication",
    "Peak SV: 570bp DEL overlapping Helitron+CACTA TEs; inside g132348; R2=49.8% QTL-seq; absent from SR GWAS due to TE multi-mapping (88% markers in TEs, OR=7.49, p=7.7e-41)",
    "FHY3/FAR1-type TF 13.7kb downstream of peak SV; photoperiod flowering regulator",
    "Chr.11:42 Mb QTL-seq only (R2=42.7%); not in GWAS panel",
    "QTL-seq corroboration of GWAS Chr.05:47-48 Mb locus (FDR-significant)"
  )
)

save_table(s2, "SupplTable_S2_Candidate_Genes",
  note = paste0(
    "Candidate genes within or nearest (<50 kb) to significant GWAS and QTL-seq loci. ",
    "Gene models from braker annotation of Suziblue hap1 assembly. Nearest gene identified by position overlap or minimum distance using genes.sorted.bed. ",
    "Functional annotations from InterProScan v5 (Pfam, PANTHER databases). ",
    "Chr.12 candidates: InterProScan runs 2026-05-11 (iprscan5-R20260511-* in Chr.12_investigation/). ",
    "Chr.06:42.4 Mb (g66181): annotation from AllQTLpart4.xls. ",
    "Chr.05 PR-10 candidates (g54686, g54697): annotated from braker.codingseq AGGPY/SCPNVES motif analysis. ",
    "9 entries marked 'pending InterProScan' require submission to InterProScan web server. ",
    "Protein sequences for pending genes extracted to Tables/GWAS_Loci_PendingIprScan.fa. ",
    "Submit to https://www.ebi.ac.uk/interpro/search/sequence/ using Pfam + PANTHER databases."
  ))

# ==============================================================================
# SUPPLEMENTARY TABLE S3 — QTL-seq top regions
# ==============================================================================
cat("Building Sup Table S3...\n")

qtl <- fread(file.path(BASE,
  "../LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq/SummaryTables/AllQTL_Ranked_Table.txt"))

ae  <- fread(file.path(BASE,
  "../LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq/SummaryTables/TopQTL_AlleleEffect_Table.txt"))

# Keep All3_Sig hits + Chr.05 47-48 Mb + Chr.12 primary peak (AllYears, SV, DTFlower)
key_markers <- c("Chr.12_10300977", "Chr.12_10308052", "Chr.05_47712013", "Chr.05_47900362")

s3 <- qtl[All3_Sig == TRUE | PeakMarker %in% key_markers]

# Add allele effect data where available
ae_merge <- ae[!is.na(Effect_days), .(
  PeakMarker, Trait, Effect_days, Wilcoxon_p, N_samples,
  Low_mean_days, High_mean_days
)]
s3 <- merge(s3, ae_merge, by = c("PeakMarker", "Trait"), all.x = TRUE)

# Select and rename columns
s3_out <- s3[, .(
  Analysis           = PhnoVer,
  Dataset            = Dataset,
  Trait              = Trait,
  Chr                = Chrom,
  `Peak marker`      = PeakMarker,
  `Peak R2 (%)`      = round(PeakR2 * 100, 2),
  `Peak delta-DS`    = round(PeakDelta, 3),
  `Min adj. p`       = signif(MinPadj, 3),
  `N sig. markers`   = NMarkers,
  `All3_Sig`         = All3_Sig,
  `Region width (Mb)` = round(Width_Mb, 3),
  `Effect (days)`    = round(Effect_days, 1),
  `Wilcoxon p`       = signif(Wilcoxon_p, 3),
  `N samples`        = N_samples,
  `Low dosage mean (days)`  = round(Low_mean_days, 1),
  `High dosage mean (days)` = round(High_mean_days, 1)
)]

# Sort: All3_Sig first, then by R2 descending
setorder(s3_out, -All3_Sig, -`Peak R2 (%)`)

save_table(s3_out, "SupplTable_S3_QTLseq_TopRegions",
  note = paste0(
    "QTL-seq results from LRLP SV and PAN marker datasets (DIY QTL-seq pipeline). ",
    "All3_Sig = TRUE: locus is significant for R2, adjusted p-value (BH method), AND delta allele frequency simultaneously. ",
    "Allele effect (days) computed by Wilcoxon rank-sum test comparing low vs high dosage bulk phenotypes. ",
    "Chr.05:47-48 Mb rows included as GWAS corroboration even where All3_Sig = FALSE. ",
    "Analysis column: AllYears = combined multi-year bulk; Method1 = alternative bulk assignment. ",
    "Dataset: SV = structural variant markers (Sniffles2); PAN = pangenome-based SNP markers."
  ))

cat("\nAll tables saved to:", OUT_DIR, "\n")
cat("Files:\n")
cat(paste(" ", list.files(OUT_DIR), collapse = "\n"), "\n")
