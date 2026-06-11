# ============================================================================
# QTL_Gene_Analysis.R
#
# PURPOSE:
#   For QTL identified by Linear_QTLSeq.R, this script:
#   1. Filters "promising" QTL by Delta, R2, P-value, and width thresholds
#   2. Finds overlapping and nearby genes (from genes.sorted.bed / braker)
#   3. Extracts SV sequences (REF/ALT) as FASTA for BLASTn
#   4. Extracts protein sequences for overlapping genes as FASTA for BLASTp
#   5. Writes a combined summary table (txt + xlsx)
#   6. Writes a ready-to-run BLAST command sheet
#
# Usage:
#   Rscript QTL_Gene_Analysis.R
#   (edit CONFIGURATION section below to adjust filters)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ─── CONFIGURATION ───────────────────────────────────────────────────────────

LINEAR_DIR <- "/Users/kendalllee/Documents/Blueberry/LRLP_Blueberry/LRLP_QTLVar/DIY_QTLSeq/Linear"
GENE_DIR   <- "/Users/kendalllee/Documents/Blueberry/Gene_Info"
OUT_DIR    <- file.path(LINEAR_DIR, "QTL_Gene_Analysis")

# ── Filters for "promising" QTL ──────────────────────────────────────────────
MIN_ABS_DELTA <- 1.0     # minimum |PeakDelta| to include
MIN_R2        <- 0.15    # minimum PeakR2
MAX_PVAL      <- 0.05    # maximum MinPval (raw, not BH)
MIN_WIDTH_MB  <- 0.0     # set > 0 to exclude single-point hits (e.g. 0.01)

# ── Gene window ──────────────────────────────────────────────────────────────
FLANK_BP      <- 50000   # search ± this distance around QTL region for genes

# ── SV sequence filters (for BLAST FASTA) ────────────────────────────────────
MIN_SV_LEN    <- 50      # only output ALT sequences >= this length for BLASTn

# ─── SETUP ───────────────────────────────────────────────────────────────────

dir.create(OUT_DIR, showWarnings = FALSE)

TRAITS <- c("FruitWT", "DTFruit", "DTFlower", "Flow2Fruit")

# ─── LOAD GENE ANNOTATIONS ───────────────────────────────────────────────────

message("Loading gene annotations ...")

# genes.sorted.bed: Chrom  Start  End  GeneID  (0-based BED coords)
genes_bed <- fread(file.path(GENE_DIR, "genes.sorted.bed"),
                   col.names = c("Chrom", "Start", "End", "GeneID"),
                   header = FALSE)
message(sprintf("  %d genes in genes.sorted.bed", nrow(genes_bed)))

# Parse braker.aa into a named list: GeneID -> protein sequence
message("Parsing braker.aa protein sequences ...")
braker_aa_lines <- readLines(file.path(GENE_DIR, "braker.aa"))
header_idx <- which(startsWith(braker_aa_lines, ">"))
# Extract gene ID from header: ">g1.t1" -> "g1"
aa_headers <- sub("^>([^ ]+).*", "\\1", braker_aa_lines[header_idx])
aa_gene_id <- sub("\\.t.*$", "", aa_headers)   # g1.t1 -> g1

seq_ends <- c(header_idx[-1] - 1, length(braker_aa_lines))
braker_aa <- mapply(function(s, e) {
  paste(braker_aa_lines[(s + 1):e], collapse = "")
}, header_idx, seq_ends, SIMPLIFY = FALSE)
names(braker_aa) <- aa_gene_id
rm(braker_aa_lines); invisible(gc())
message(sprintf("  %d protein sequences loaded", length(braker_aa)))

# ─── LOAD QTL RESULTS ────────────────────────────────────────────────────────

message("\nLoading QTL results ...")

all_qtl   <- list()
all_top25 <- list()

for (trait in TRAITS) {
  qtl_file   <- file.path(LINEAR_DIR, paste0("results_", trait),
                           paste0(trait, "_QTL_regions.txt"))
  top25_file <- file.path(LINEAR_DIR, paste0("results_", trait),
                           paste0(trait, "_Top25_QTL_markers.txt"))

  if (file.exists(qtl_file)) {
    dt <- fread(qtl_file)
    dt[, Trait := trait]
    all_qtl[[trait]] <- dt
    message(sprintf("  %s: %d QTL regions", trait, nrow(dt)))
  } else {
    message(sprintf("  %s: QTL file not found (run still in progress?)", trait))
  }

  if (file.exists(top25_file)) {
    dt2 <- fread(top25_file)
    dt2[, Trait := trait]
    all_top25[[trait]] <- dt2
  }
}

qtl_all   <- rbindlist(all_qtl,   fill = TRUE)
top25_all <- rbindlist(all_top25, fill = TRUE)

# ─── FILTER PROMISING QTL ────────────────────────────────────────────────────

message("\nFiltering promising QTL ...")
message(sprintf("  Criteria: |Delta| >= %.2f, R2 >= %.2f, Pval <= %.3f, Width >= %.3f Mb",
                MIN_ABS_DELTA, MIN_R2, MAX_PVAL, MIN_WIDTH_MB))

promising <- qtl_all[
  abs(PeakDelta) >= MIN_ABS_DELTA &
  PeakR2        >= MIN_R2         &
  MinPval       <= MAX_PVAL       &
  Width_Mb      >= MIN_WIDTH_MB
][order(Trait, -PeakR2)]

message(sprintf("  %d / %d QTL pass filters", nrow(promising), nrow(qtl_all)))
print(promising[, .(Trait, Chrom, Start, End, Width_Mb, PeakPos,
                     PeakDelta, PeakR2, MinPval)])

# ─── FIND OVERLAPPING + NEARBY GENES ─────────────────────────────────────────

message("\nFinding overlapping and nearby genes ...")

find_genes <- function(chrom, start, end, peak_pos, flank = FLANK_BP) {
  # Expand search window
  win_start <- max(0, start - flank)
  win_end   <- end + flank

  hits <- genes_bed[Chrom == chrom & End >= win_start & Start <= win_end]
  if (nrow(hits) == 0) return(NULL)

  hits[, distance_to_peak := pmax(0, pmax(Start, peak_pos) - pmin(End, peak_pos))]
  hits[, relationship := fcase(
    Start <= peak_pos & End >= peak_pos, "overlaps_peak",
    Start >= start    & End <= end,      "within_QTL",
    End   <  start,                      "upstream_of_QTL",
    Start >  end,                        "downstream_of_QTL",
    default =                            "flanks_QTL"
  )]
  hits[order(distance_to_peak)]
}

gene_results <- list()
for (i in seq_len(nrow(promising))) {
  q <- promising[i]
  genes_hit <- find_genes(q$Chrom, q$Start, q$End, q$PeakPos)
  if (!is.null(genes_hit) && nrow(genes_hit) > 0) {
    genes_hit[, `:=`(Trait = q$Trait, QTL_Chrom = q$Chrom,
                      QTL_Start = q$Start, QTL_End = q$End,
                      QTL_PeakPos = q$PeakPos, QTL_PeakR2 = q$PeakR2,
                      QTL_PeakDelta = q$PeakDelta, QTL_MinPval = q$MinPval)]
    gene_results[[i]] <- genes_hit
  }
}

gene_table <- rbindlist(gene_results, fill = TRUE)
gene_table <- gene_table[, .(Trait, QTL_Chrom, QTL_Start, QTL_End,
                               QTL_PeakPos, QTL_PeakDelta, QTL_PeakR2,
                               QTL_MinPval, GeneID,
                               Gene_Start = Start, Gene_End = End,
                               Distance_to_Peak_bp = distance_to_peak,
                               Relationship = relationship)]

message(sprintf("  Found %d gene-QTL associations across %d unique genes",
                nrow(gene_table), uniqueN(gene_table$GeneID)))

# ─── WRITE GENE TABLE ────────────────────────────────────────────────────────

out_gene_txt  <- file.path(OUT_DIR, "QTL_gene_overlap.txt")
out_gene_xlsx <- file.path(OUT_DIR, "QTL_gene_overlap.xlsx")
fwrite(gene_table, out_gene_txt, sep = "\t")
message("  Saved: ", out_gene_txt)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(as.data.frame(gene_table), out_gene_xlsx)
  message("  Saved: ", out_gene_xlsx)
}

# ─── EXTRACT SV SEQUENCES FOR BLASTn ─────────────────────────────────────────
#
# For each promising QTL, grab its peak marker from top25 and write
# REF/ALT sequences as FASTA. Only sequences >= MIN_SV_LEN are written.
# These can be BLASTed against nt (BLASTn) to check conservation / identity.

message("\nExtracting SV sequences for BLASTn ...")

# Get peak marker for each promising QTL
peak_markers <- merge(
  promising[, .(Trait, Chrom, PeakPos)],
  top25_all[, .(Trait, Chrom = Chrom, Position, Marker, REF, ALT,
                DeltaDS_smooth, R2, Pvalue)],
  by.x = c("Trait", "Chrom", "PeakPos"),
  by.y = c("Trait", "Chrom", "Position"),
  all.x = TRUE
)

# If peak not in top25, fall back to any marker within ±500 bp in top25
unmatched <- peak_markers[is.na(Marker)]
if (nrow(unmatched) > 0) {
  matched_fallback <- lapply(seq_len(nrow(unmatched)), function(i) {
    u <- unmatched[i]
    fb <- top25_all[Trait == u$Trait & Chrom == u$Chrom &
                    abs(Position - u$PeakPos) <= 500000][1]
    if (nrow(fb) == 0) return(NULL)
    cbind(u[, .(Trait, Chrom, PeakPos)], fb[, .(Marker, REF, ALT,
                                                   DeltaDS_smooth, R2, Pvalue)])
  })
  fallback_dt <- rbindlist(matched_fallback, fill = TRUE)
  peak_markers <- rbindlist(list(peak_markers[!is.na(Marker)], fallback_dt),
                             fill = TRUE)
}

# Write FASTA per trait
fasta_written <- 0
blast_n_files <- c()

for (trait in unique(peak_markers$Trait)) {
  sub <- peak_markers[Trait == trait & !is.na(Marker)]
  if (nrow(sub) == 0) next

  fasta_lines <- c()
  for (j in seq_len(nrow(sub))) {
    mk  <- sub[j]
    alt <- as.character(mk$ALT)
    ref <- as.character(mk$REF)

    # ALT sequence
    if (!is.na(alt) && alt != "N" && nchar(alt) >= MIN_SV_LEN) {
      fasta_lines <- c(fasta_lines,
        sprintf(">%s|%s:%d|ALT|Delta=%.3f|R2=%.3f",
                mk$Marker, mk$Chrom, mk$PeakPos,
                mk$DeltaDS_smooth, mk$R2),
        alt)
      fasta_written <- fasta_written + 1
    }
    # REF sequence (if it's a real sequence, not just "N")
    if (!is.na(ref) && ref != "N" && nchar(ref) >= MIN_SV_LEN) {
      fasta_lines <- c(fasta_lines,
        sprintf(">%s|%s:%d|REF|Delta=%.3f|R2=%.3f",
                mk$Marker, mk$Chrom, mk$PeakPos,
                mk$DeltaDS_smooth, mk$R2),
        ref)
      fasta_written <- fasta_written + 1
    }
  }

  if (length(fasta_lines) > 0) {
    out_fa <- file.path(OUT_DIR, paste0(trait, "_peak_SV_sequences.fasta"))
    writeLines(fasta_lines, out_fa)
    blast_n_files <- c(blast_n_files, out_fa)
    message(sprintf("  %s: %d sequences -> %s", trait, length(fasta_lines) / 2, out_fa))
  }
}
message(sprintf("  Total SV sequences written: %d", fasta_written))

# ─── EXTRACT PROTEIN SEQUENCES FOR BLASTp ────────────────────────────────────
#
# For each unique gene overlapping a promising QTL, extract its protein
# sequence from braker.aa. Write per-trait and combined FASTA files.

message("\nExtracting protein sequences for BLASTp ...")

blastp_written <- 0
for (trait in unique(gene_table$Trait)) {
  gene_ids <- unique(gene_table[Trait == trait, GeneID])
  seqs     <- braker_aa[names(braker_aa) %in% gene_ids]
  if (length(seqs) == 0) {
    message(sprintf("  %s: no protein sequences found for %d genes", trait, length(gene_ids)))
    next
  }

  # Add functional context to header
  fasta_lines <- unlist(lapply(names(seqs), function(gid) {
    row <- gene_table[Trait == trait & GeneID == gid][1]
    c(sprintf(">%s|%s:%d-%d|%s|R2=%.3f|Delta=%.3f",
              gid, row$QTL_Chrom, row$Gene_Start, row$Gene_End,
              row$Relationship, row$QTL_PeakR2, row$QTL_PeakDelta),
      seqs[[gid]])
  }))

  out_pfa <- file.path(OUT_DIR, paste0(trait, "_QTL_gene_proteins.fasta"))
  writeLines(fasta_lines, out_pfa)
  blastp_written <- blastp_written + length(seqs)
  message(sprintf("  %s: %d proteins -> %s", trait, length(seqs), out_pfa))
}

# Combined protein FASTA (all traits, deduplicated)
all_gene_ids  <- unique(gene_table$GeneID)
combined_seqs <- braker_aa[names(braker_aa) %in% all_gene_ids]
if (length(combined_seqs) > 0) {
  out_combined <- file.path(OUT_DIR, "ALL_traits_QTL_gene_proteins.fasta")
  combined_lines <- unlist(lapply(names(combined_seqs), function(gid) {
    c(paste0(">", gid), combined_seqs[[gid]])
  }))
  writeLines(combined_lines, out_combined)
  message(sprintf("  Combined: %d unique proteins -> %s",
                  length(combined_seqs), out_combined))
}

# ─── WRITE BLAST COMMAND SHEET ───────────────────────────────────────────────

blast_script <- file.path(OUT_DIR, "BLAST_commands.sh")

blast_lines <- c(
  "#!/bin/bash",
  "# ============================================================",
  "# BLAST commands for promising QTL sequences",
  "# Generated by QTL_Gene_Analysis.R",
  "# ============================================================",
  "",
  "# ── SETUP ────────────────────────────────────────────────────",
  "# Load BLAST module or activate conda env first:",
  "# module load blast/2.13.0",
  "# conda activate blast_env",
  "",
  paste0("OUT_DIR=", OUT_DIR),
  "",
  "# ── BLASTn: SV sequences vs NCBI nt ─────────────────────────",
  "# Identifies where SV sequences come from (TEs, known genes, etc.)",
  "# Run on cluster - nt database is large"
)

for (trait in TRAITS) {
  fa <- file.path(OUT_DIR, paste0(trait, "_peak_SV_sequences.fasta"))
  if (file.exists(fa)) {
    blast_lines <- c(blast_lines, "",
      sprintf("# %s", trait),
      sprintf('blastn \\'),
      sprintf('  -query "%s" \\', fa),
      sprintf('  -db nt \\'),
      sprintf('  -out "%s" \\',
              file.path(OUT_DIR, paste0(trait, "_SV_blastn_results.txt"))),
      sprintf('  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\'),
      sprintf('  -evalue 1e-5 -perc_identity 80 -max_target_seqs 5 -num_threads 8')
    )
  }
}

blast_lines <- c(blast_lines, "",
  "# ── BLASTp: Gene proteins vs UniProt/nr ─────────────────────",
  "# Identifies gene function — run against uniprot_sprot for fastest results",
  "",
  sprintf('blastp \\'),
  sprintf('  -query "%s" \\', file.path(OUT_DIR, "ALL_traits_QTL_gene_proteins.fasta")),
  sprintf('  -db uniprot_sprot \\'),
  sprintf('  -out "%s" \\', file.path(OUT_DIR, "ALL_QTL_genes_blastp_results.txt")),
  sprintf('  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \\'),
  sprintf('  -evalue 1e-5 -max_target_seqs 3 -num_threads 8'),
  "",
  "# ── Local BLAST alternative (if cluster not available) ───────",
  "# Download blueberry proteome and BLAST locally:",
  "# makeblastdb -in reference_proteins.fasta -dbtype prot -out ref_prot_db",
  sprintf('# blastp -query "%s" -db ref_prot_db \\',
          file.path(OUT_DIR, "ALL_traits_QTL_gene_proteins.fasta")),
  '#   -out local_blastp_results.txt -outfmt 6 -evalue 1e-5 -num_threads 4',
  "",
  "# ── NCBI web BLAST alternative ───────────────────────────────",
  "# Upload FASTA files at: https://blast.ncbi.nlm.nih.gov/",
  sprintf("# SV sequences (BLASTn):  %s/<trait>_peak_SV_sequences.fasta", OUT_DIR),
  sprintf("# Gene proteins (BLASTp): %s/ALL_traits_QTL_gene_proteins.fasta", OUT_DIR)
)

writeLines(blast_lines, blast_script)
Sys.chmod(blast_script, mode = "0755")
message("\nSaved BLAST command sheet: ", blast_script)

# ─── FINAL SUMMARY ───────────────────────────────────────────────────────────

message("\n", strrep("=", 60))
message("SUMMARY")
message(strrep("=", 60))
message(sprintf("  Promising QTL:          %d across %d traits",
                nrow(promising), uniqueN(promising$Trait)))
message(sprintf("  Gene-QTL associations:  %d (%d unique genes)",
                nrow(gene_table), uniqueN(gene_table$GeneID)))
message(sprintf("  SV FASTA sequences:     %d sequences written", fasta_written))
message(sprintf("  Protein sequences:      %d unique genes", length(combined_seqs)))
message(sprintf("  Output directory:       %s", OUT_DIR))
message("")
message("Files ready for BLAST:")
message(sprintf("  BLASTn (SV seqs):  <trait>_peak_SV_sequences.fasta"))
message(sprintf("  BLASTp (proteins): ALL_traits_QTL_gene_proteins.fasta"))
message(sprintf("  Commands:          BLAST_commands.sh"))
message(sprintf("  Gene overlap:      QTL_gene_overlap.xlsx"))
message("\nNext steps:")
message("  1. Review QTL_gene_overlap.xlsx - check 'overlaps_peak' genes first")
message("  2. Upload *_peak_SV_sequences.fasta to NCBI BLASTn")
message("     -> Insertions from TEs? Known functional elements?")
message("  3. Upload ALL_traits_QTL_gene_proteins.fasta to NCBI BLASTp")
message("     -> What do these genes do? Known flowering/fruit genes?")
message("  4. Cross-reference gene IDs with braker.hap1.gtf for exon structure")
message("  5. Genes with 'overlaps_peak' + high R2 + known function = top candidates")
