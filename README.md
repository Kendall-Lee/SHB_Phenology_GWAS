# SHB Phenology GWAS and QTL-seq Analysis

**"Tetraploid-aware GWAS and independent QTL-seq reveal stable loci associated with phenology in Southern Highbush Blueberry (*Vaccinium corymbosum* L.)"**

Kendall Lee et al. — *The Plant Genome* (in preparation)

---

## Overview

This repository contains the R analysis scripts and data tables supporting the manuscript. We applied a tetraploid dosage-aware GWAS (GWASpoly) to ~800 Southern Highbush Blueberry (SHB) accessions phenotyped over three years (2023–2025) for four phenology traits, combined with independent SV-based QTL-seq in an LRLP panel against the Suziblue hap1 reference genome.

### Key findings
- **22 stable GWAS loci** across 9 chromosomes (Meff threshold, −log₁₀p ≥ 5.52); predominantly additive models
- **Chr.05:47–48 Mb** as the anchor locus: highest-confidence signal (−log₁₀p = 11.86, PVE = 7.79%), confirmed by independent QTL-seq (R² = 27%), PR-10/MLP candidate genes (g54686, g54697)
- **Chr.12:9.5–11.5 Mb**: largest-effect QTL-seq locus (R² = 49.8%, +23 days/ALT copy); absent from SR GWAS due to TE-dense region (88% of markers inside TEs); candidate gene g132348 (Protein Arginine N-Methyltransferase)
- Year 2025 frost year (21 frost days, 12 hard freezes) revealed 14× more detectable loci, demonstrating environmental stress amplification of genetic signal

---

## Repository structure

```
R_Scripts/       # Figure and table generation scripts (R)
Tables/          # CSV tables (main text + supplementary)
docs/            # Figure legends, manuscript outline
```

### R_Scripts

| Script | Output |
|---|---|
| `Fig1_Phenotypic_Context.R` | Fig 1 — trait correlations, frost weather, DTFlower SD collapse |
| `Fig2_Manhattan.R` | Fig 2 — linear GWAS Manhattan (BLUE_exc.25 + yr.23/24/25) |
| `Fig3_Chr05_Hotspot.R` | Fig 3 — Chr.05:47–48 Mb regional + allele effects + QTL-seq overlay |
| `Fig4_StableQTL_Heatmap.R` | Fig 4 — 22-locus × 16 trait-year heatmap |
| `Fig4_Chr12_QTLseq.R` | Fig 5 — Chr.12 QTL-seq lollipop + allele effect |
| `Fig4_Yr25_Stress_QTL.R` / `Fig6_Chr03_StressLocus.R` | Sup Fig S2 — year 2025 stress QTL (3-panel) |
| `Build_All_Tables.R` | All main-text and supplementary tables |
| `SupplTable_S1.R` | Supplementary Table S1 generation |

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

- R (≥ 4.3)
- GWASpoly (≥ 2.0) — tetraploid GWAS
- data.table, ggplot2, cowplot, openxlsx, readxl
- GWASpoly scan objects (`.RData`) — available from corresponding author on request

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
