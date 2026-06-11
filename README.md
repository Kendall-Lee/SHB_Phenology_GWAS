# SHB Phenology GWAS and QTL-seq Analysis

**"Genomic architecture of phenological timing in Southern Highbush Blueberry (*Vaccinium corymbosum* L.): a chromosome-scale reference genome enables tetraploid GWAS of days to flowering, days to fruit, and fruiting period"**

Kendall Lee et al. — *The Plant Genome* (in preparation)

---

## Overview

This repository contains the analysis scripts and data tables supporting the manuscript. We assembled a chromosome-scale phased reference genome for the SHB cultivar Suziblue (548.4 Mb; 12 chromosomes) and applied tetraploid dosage-aware GWAS (GWASpoly LOCO) to 665 SHB accessions phenotyped over three growing seasons (2023–2025) for four phenology traits. Results were independently validated by bulk segregant analysis (BSA) in a 95-accession long-read low-pass (LRLP) panel sequenced on PacBio Revio.

### Key findings
- **Chr.05:47–48.6 Mb** — major stable locus for DTFlower and DTFruit (−log₁₀p = 11.86, PVE = 7.79%); additive dosage gradient spanning ~21 days; independently confirmed by LRLP BSA (R² = 26–27%, FDR-significant); 20 candidate genes across 5 sub-clusters
- **22 stable GWAS loci** across 8 chromosomes (Meff threshold −log₁₀p ≥ 5.52); Chr.06:22–23 Mb is the only fully within-study-replicated secondary locus
- **Year 2025 frost year** (21 frost days, 12 hard freezes, Jan min −5.3°C) revealed 281 stress-responsive QTL (14× increase over stable analysis), including 3 large-effect loci absent from normal-year scans
- **First SHB-derived phased reference genome** (Suziblue hap1) enabling species-appropriate mapping for GWAS and pangenome construction

---

## Repository structure

```
R_Scripts/       # Publication figure and table generation (R)
BSA_QTLseq/     # In-house LRLP bulk segregant analysis pipeline (R + bash)
Python/          # Utility figure scripts (Python/matplotlib)
Tables/          # CSV tables (main text + supplementary)
docs/            # Figure legends, manuscript notes
```

### R_Scripts

| Script | Output |
|---|---|
| `Fig1_Phenotypic_Context.R` | Fig 1 — trait distributions, frost weather, DTFlower SD collapse |
| `Fig2_Manhattan_TPG.R` | Fig 2 — GWAS Manhattan (BLUE_exc.25 + yr.23/24/25; TPG specs) |
| `Fig2_Manhattan.R` | Fig 2 — earlier draft Manhattan |
| `Fig3_Chr05_Hotspot.R` | Fig 3 — Chr.05:47–48 Mb regional zoom + allele effects + BSA overlay |
| `Fig4_StableQTL_Heatmap.R` | Fig 4 — 22-locus × 12 trait-year heatmap |
| `Fig6_Chr03_StressLocus.R` / `Fig4_Yr25_Stress_QTL.R` | Fig 5 — yr.25 stress QTL PVE bars |
| `FigS2_Yr25_Manhattan.R` | Sup Fig S2 — yr.25 stress GWAS Manhattan (4-trait 2×2 grid) |
| `FigS_Allele_Effects_AllLoci.R` | Sup Fig — dosage boxplots for all stable loci |
| `FigS_Allele_Effects_Secondary.R` | Sup Fig — dosage boxplots for secondary loci |
| `SupFigS7_Suziblue_Assembly.R` | Sup Fig S7 — Suziblue Hi-C contact map + synteny |
| `Build_All_Tables.R` | All main-text and supplementary tables |
| `SupplTable_S1.R` | Supplementary Table S1 generation |
| `Fig4_Chr12_QTLseq.R` | Chr.12 QTL-seq figure (future manuscript) |

### BSA_QTLseq

In-house bulk segregant analysis pipeline for the LRLP long-read panel. Implements PCA-corrected delta-dosage QTL-seq following Takagi et al. adapted for autotetraploid dosage data.

| Script | Purpose |
|---|---|
| `Linear_QTLSeq.R` | Core pipeline: PCA correction, delta-dosage scores, R², FDR via permutation |
| `Run_LinearPCA_FullPipeline.R` | Full pipeline runner — calls all steps in sequence |
| `Enhanced_QTLSeq_Analysis_BLUES.R` | BLUE-based bulk assignment (replaces individual-year extremes) |
| `QTL_seqAlleleffect.R` | Allele effect boxplots at QTL-seq peak markers |
| `QTL_Gene_Analysis.R` | Candidate gene annotation at BSA loci (bedtools intersect + InterProScan) |
| `run_all_traits_BLUES.sh` | Batch runner across all four traits |
| `average_cov.sh` | Samtools breadth-of-coverage summary across LRLP BAMs |

### Python

| Script | Purpose |
|---|---|
| `build_lrlp_coverage_fig.py` | Scatter plot: LRLP mean haploid depth vs. % genome breadth coverage (n=95) |
| `build_methods_clean.py` | Study design flowchart figure (matplotlib, publication style) |

### Tables

| File | Description |
|---|---|
| `Table1_Phenotypic_Summary.csv` | H², phenotypic means by year, rank correlations |
| `Table2_Stable_QTL_Summary.csv` | 22 stable GWAS loci, one peak marker per locus |
| `SupplTable_S1_All88Markers.csv` | All 88 significant markers (full details) |
| `SupplTable_S2_Candidate_Genes.csv` | Candidate gene annotations (GWAS + QTL-seq loci) |
| `SupplTable_S3_QTLseq_TopRegions.csv` | QTL-seq top regions (All3_Sig hits + corroboration) |
| `GWAS_Loci_PendingIprScan.fa` | Protein sequences for 9 GWAS loci awaiting InterProScan annotation |

---

## Data availability

| Data type | Repository | Accession |
|---|---|---|
| Raw phenotyping data | This repository (Tables/) | — |
| Genotyping data (SR SNPs) | NCBI SRA | *[to be deposited]* |
| LRLP long-read sequencing | NCBI SRA | *[to be deposited]* |
| Suziblue hap1 reference | NCBI GenBank | *[to be deposited]* |

---

## Software dependencies

**R scripts**
- R (≥ 4.3)
- GWASpoly (≥ 2.14) — tetraploid dosage GWAS
- data.table, ggplot2, cowplot, scales, openxlsx, readxl
- GWASpoly scan objects (`.RData`) — available from corresponding author on request

**BSA pipeline**
- R (≥ 4.3) with data.table, ggplot2
- samtools (≥ 1.17) — BAM coverage
- bedtools (≥ 2.31) — gene intersection

**Python scripts**
- Python (≥ 3.9) with matplotlib, numpy

---

## Usage

All figure scripts expect GWASpoly scan `.RData` objects at the paths defined in each script's `BASE` variable. Update `BASE` to your local path before running.

```r
# Example
BASE <- "/your/path/to/FINAL_GWAS"
source("R_Scripts/Fig1_Phenotypic_Context.R")
```

Table scripts (`Build_All_Tables.R`) require the 88-marker CSV (`SupplTable_S1_TopMarkers_LINEAR_BLUEexc25.csv`) in `BASE/PUBLICATION/Linear_MS/`.

---

## License

Scripts: MIT License. Data: CC BY 4.0.

## Contact

Kendall Lee — lee.kendall.94@gmail.com
