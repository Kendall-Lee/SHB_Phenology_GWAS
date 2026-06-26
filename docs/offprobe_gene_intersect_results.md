# Off-Probe Marker Gene Overlap Analysis

**Script:** `Python/offprobe_gene_intersect.py`  
**Date:** 2026-06-26  
**Output BED:** `Tables/offprobe_genic_markers.bed`

## Background

The GWAS marker panel contains 109,927 SNPs classified as either probe-enriched
(overlapping a VacCap capture probe interval; n = 17,672; 16.1%) or off-probe
(outside all probe intervals; n = 92,255; 83.9%), representing a ~9-fold
enrichment of marker density in probe-targeted regions relative to random
expectation (see `Scripts/classify_markers_by_probe.R`).

The question addressed here: are off-probe markers simply intergenic noise, or
do a meaningful fraction fall in annotated gene bodies — consistent with
homeologous gene copies not covered by the Draper-derived probe design?

## Methods

Gene body intervals were extracted from the BRAKER hap1 annotation
(`Gene_Info/braker.hap1.gtf`; chromosome names Chr.XX.1, stripped to Chr.XX to
match dosage matrix coordinates). Marker positions were intersected with gene
intervals using binary search. Genome-wide genic fraction was calculated from
the Suziblue hap1 assembly (`Enrichment/Suziblue_hap1.fa.fai`; 548,367,372 bp).

## Results

| Class | In gene body | Off gene body | Total | % genic |
|---|---|---|---|---|
| enriched (probe) | 5,214 | 12,458 | 17,672 | 29.5% |
| wgs_only (off-probe) | 29,243 | 63,012 | 92,255 | **31.7%** |

**Genome-wide genic fraction:** 133,451,882 bp / 548,367,372 bp = **24.3%**  
**Off-probe markers in gene bodies:** 29,243 / 92,255 = **31.7%** (1.30× above random expectation)

## Interpretation

Roughly one-third of off-probe markers (31.7%) fall within BRAKER-annotated
gene bodies, exceeding the 24.3% genome-wide genic fraction by 1.30-fold.
This is consistent with partial variant capture in homeologous gene copies or
SHB-specific gene content not represented in the Draper-derived probe design.
The remaining ~68% of off-probe markers are intergenic or in repetitive regions
and represent background whole-genome sequencing coverage from off-target reads.

**Suggested manuscript language:**
> "Although classified as off-probe, roughly one-third of these markers (31.7%)
> fall within annotated gene bodies — above the 24.3% genome-wide genic
> fraction — suggesting that many capture variants in homeologous copies of
> targeted genes not covered by the probe design."
