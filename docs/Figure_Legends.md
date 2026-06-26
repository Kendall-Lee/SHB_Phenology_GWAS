# Figure Legends
## Tetraploid-aware GWAS and independent QTL-seq reveal stable loci associated with phenology in Southern Highbush Blueberry

---

## Main Text Figures

### Figure 1. Phenotypic characterization of four traits across three field years.

**(A)** Pairwise Spearman rank correlations among Days to 50% Flowering (DTFlower), Days to 50% Ripe Fruit (DTFruit), Fruiting Period (Flow2Fruit), and 25-Fruit Weight (FruitWeight) for each year (2023, 2024, 2025). Color scale indicates correlation coefficient (blue = positive, red = negative). The most notable year-to-year shift involves DTFlower×Flow2Fruit, which strengthened from r = −0.36 in yr.23 to r = −0.72 in yr.24 and r = −0.58 in yr.25, reflecting increasing coupling between bloom date and fruiting interval under more stressful conditions. DTFlower×DTFruit was moderately positive across all years (r = 0.56, 0.38, 0.52 in yr.23, yr.24, yr.25, respectively), indicating the two timing traits share genetic control but are not fully correlated. **(B)** Count of frost days (minimum temperature ≤ 32°F) and hard-freeze events (minimum temperature ≤ 28°F) per year at the Alapaha, Georgia, field station. Year 2025 had 3.5× as many frost days (6 → 11 → 21) and 6× as many hard-freeze events (2 → 4 → 12) as 2023. **(C)** Phenotypic standard deviation of DTFlower BLUEs by year (2023: 13.3 days; 2024: 11.2 days; 2025: 6.9 days), illustrating the 48% reduction in phenotypic variance driven by frost-synchronization of flowering in 2025. High H² values across all years (0.73–0.91) confirm that variance compression reflects environmental masking of the genetic signal rather than data-quality loss. BLUEs computed via REML mixed model (lmer: trait ∼ Year + [1|Genotype]) across approximately 800 Southern Highbush Blueberry accessions.

---

### Figure 2. Genome-wide association results for phenology traits across three years.

Manhattan plots from GWASpoly LOCO scan for Days to 50% Flowering (DTFlower), Days to 50% Ripe Fruit (DTFruit), and Fruiting Period (Flow2Fruit). For each trait, the large left panel shows the primary BLUE_exc.25 analysis (combined 2023–2024 BLUEs); three smaller right panels show individual year scans (yr.23, yr.24, yr.25). The horizontal dashed line indicates the Meff multiple-testing threshold (−log₁₀p = 5.52, α = 0.05). Chromosomes are alternately colored. The Chr.05 cluster (47–48 Mb) represents the highest-confidence signal in the study (peak −log₁₀p = 11.86 for DTFlower). Within-study replication is visible as consistent peaks across yr.23 and yr.24 panels. Full marker details are provided in Supplementary Table S1.

---

### Figure 3. Chr.05:47–48 Mb — the anchor locus for blueberry phenology.

**(A)** Regional association plot for the Chr.05:46–49.5 Mb window. Points represent individual marker −log₁₀p scores (BLUE_exc.25 analysis), colored by trait (DTFlower: dark red; DTFruit: blue; Flow2Fruit: teal). The horizontal dashed line indicates the Meff significance threshold. Gene models (braker annotation, Suziblue hap1) are shown as a track below; g54686 and g54697 (MLP/Bet v1-like PR-10 family) are annotated. **(B)** Tetraploid allele dosage effect boxplots for the peak marker Chr.05_48252083. Boxes show per-genotype DTFlower BLUEs for the primary BLUE_exc.25 (dark red border, n = 626), yr.23 (green, n = 669), and yr.24 (purple, n = 648) analyses. Tukey HSD compact letter display indicates pairwise significance between dosage classes. **(C)** Independent QTL-seq signal (R², smoothed) from the LRLP SV panel overlaid on the Chr.05:46–49.5 Mb region. FDR-significant QTL-seq peaks (padj < 0.05) at Chr.05_47712013 (R² = 26.6%, p = 9.3×10⁻⁶) and Chr.05_47900362 (R² = 27.3%, p = 5.0×10⁻⁵) independently corroborate the GWAS signal and localize to the g54686/g54697 gene neighborhood.

---

### Figure 4. Cross-trait, cross-year stability of 22 significant GWAS loci.

Heatmap of −log₁₀p scores for 22 stable loci (rows) across 16 trait × year combinations (columns). Each row represents a unique locus defined by 1-Mb windows; the marker with the highest score per window is shown. Columns are grouped by trait (DTFlower, DTFruit, Flow2Fruit) and within each by analysis year (BLUE_exc.25, yr.23, yr.24, yr.25). Color intensity reflects −log₁₀p; cells below the Meff threshold (5.52) are shown in grey. Row labels include the primary trait designation; loci with signal in both yr.23 and yr.24 individual scans are self-replicating within the study. Chr.05 loci (top rows) dominate the DTFlower signal; Chr.10:14.2 Mb is the sole significant Flow2Fruit locus. Full details are provided in Table 2 and Supplementary Table S1.

---

### Figure 5. Chr.12:9.5–11.5 Mb — a TE-embedded major-effect QTL-seq locus.

**(A)** Regional SV QTL-seq R² profile across Chr.12:9.5–11.5 Mb (LRLP panel, n = 95). Stems show raw R² values for 20 FDR-significant markers (padj < 0.05, allele frequency ≥ 0.21, size ≤ 50 kb); colored by SV type (DEL: red; INS: blue). The grey line is a Gaussian-smoothed R² curve. Gold rug at the bottom indicates EDTA-annotated TE positions; the interval is among the most TE-dense on Chr.12. The dashed red vertical line marks the peak deletion (570 bp DEL, R² = 50%) at Chr.12:10,308,052, which overlaps Helitron and CACTA transposons inside gene g132348 (Protein Arginine N-Methyltransferase, PRMT; IPR025799). The locus was absent from short-read GWAS because 88% of markers in this window fall inside TEs (short reads multi-map at MAPQ=0; OR = 7.49 vs. Chr.12 background, p = 7.7×10⁻⁴¹). **(B)** Allele effect of the Chr.12 QTL. Left: aggregate mean allele dosage (0–1 scale) at the peak SV (Chr.12:10,308,052) by DTFlower phenotype bulk (early, non-extreme, late); 69% of LRLP samples cannot be individually genotyped at this TE-embedded position, so QTL-seq detects the signal via bulk read-level allele frequencies (early bulk AF = 0.28; late bulk AF = 1.00; deltaDS = 2.87). Right: per-sample tetraploid allele dosage (0–4) vs. DTFlower BLUE at the nearest well-covered proxy marker Chr.12:10,399,365 (91 kb from peak, R² = 25.5%, n = 81; Kruskal-Wallis p = 9.5×10⁻⁵). See Supplementary Figure S5 for full TE content and short-read vs. long-read variant density analysis.

---

## Supplementary Figures

### Supplementary Figure S1. Per-locus allele effect plots for all 22 stable GWAS loci.

Tetraploid dosage (0–4) vs. trait BLUE boxplots for the peak marker at each of the 22 stable loci. Each multi-page figure shows three superimposed analyses: BLUE_exc.25 (dark red border), yr.23 (green), and yr.24 (purple) per locus, for the primary associated trait. Wilcoxon rank-sum Tukey HSD compact letters indicate pairwise significance. Source file: `AlleleEffects_LINEAR_BLUEexc25.pdf`.

### Supplementary Figure S2. Year 2025 frost-stress QTL architecture.

**(A)** Count of significant QTL detected at the Meff threshold (−log₁₀p ≥ 5.52) across five analyses: Linear GWAS all-years (358), Pangenome GWAS all-years (345), Linear GWAS BLUE_exc.25 (20), Pangenome GWAS BLUE_exc.25 (23), and year 2025 stress scan (281). The 14-fold increase in the 2025 scan reflects frost-triggered amplification of genetic signal rather than analytical inflation (H² = 0.73–0.91 unchanged across years). **(B)** January minimum temperature recorded at the Georgia field station for 2023–2025 (3.8°C, −1.9°C, −5.3°C), documenting the progressive intensification of winter cold events. **(C)** PVE (%) in year 2025 for exploratory stress-responsive loci (DTFruit Chr.03:38.4 Mb*, DTFlower Chr.11:45.5 Mb, DTFruit Chr.07:25.0 Mb) and two stable phenology loci in yr.25 (DTFlower Chr.05:48.3 Mb, PVE = 3.0%; DTFlower Chr.09:39.6 Mb, PVE = 0.58%, non-significant). Faded bars indicate stable loci; solid bars indicate yr.25-specific candidates. *Chr.03:38.4 Mb signal is Flow2Fruit only; n = 5 rare allele carriers; marker not unique in genome; requires independent replication before biological interpretation.

### Supplementary Figure S3. LD decay curves for the linear GWAS panel.

Pairwise linkage disequilibrium (r²) decay as a function of inter-marker distance (kb) for the linear SNP panel across all chromosomes and per-chromosome subsets. Source file: `Fig_LD_Decay_Combined.pdf`.

### Supplementary Figure S4. Population structure of the linear GWAS panel.

Principal component analysis of the tetraploid SNP genotype matrix (SHB_tet.DP.8.maxmiss.3.minminor.25). PC1 vs. PC2 scores are shown, with accessions colored by geographic origin or breeding program if available. The kinship matrix used for LOCO correction in GWASpoly was derived from this SNP panel.

### Supplementary Figure S5. Mechanistic basis for short-read GWAS failure at Chr.12:9.5–11.5 Mb.

**(A)** SV QTL-seq R² smooth (purple) vs. SNP QTL-seq R² smooth (green) across Chr.12:9.5–11.5 Mb. The SV signal peaks at R² = 50% (Chr.12:10.308 Mb) while the SNP signal is flat (max R² = 11.7%), confirming the causal allele is structural, not a SNP. **(B)** Per-25kb window SNP density (orange = LRLP panel, normalized per 100 samples; blue = SR GWAS panel) and LRLP SV density (purple) with SV QTL-seq R² overlay. SR calls are absent at the QTL peak; the LRLP SV spike precisely co-localizes with R² = 50%. **(C)** Zoom to Chr.12:10.25–10.50 Mb showing TE feature annotations (Helitron = orange; CACTA = purple; LTR = teal) with individual SV positions. The 570 bp peak DEL (red) overlaps both a Helitron and a CACTA transposon within gene g132348. **(D)** SR SNP (blue ticks) and LRLP SNP (amber ticks) positions over the TE rug, confirming SR calls occur only in TE-free gaps while LRLP SNPs penetrate TE-dense zones. Source: `Chr.12_investigation/Fig_Chr12_ManuscriptLocus.pdf`.

### Supplementary Figure S7. Suziblue hap1 chromosome-scale reference genome assembly.

**(A)** Hi-C contact map of the Suziblue phased assembly generated in Juicebox (post-JBAT manual curation). The diagonal blocks represent 24 chromosomes across the two assembled haplotypes (hap1 and hap2; total assembly span 2,147 Mb). Strong intra-chromosomal contact enrichment along the diagonal and the characteristic off-diagonal homeologous cross-pattern confirm chromosome-scale scaffolding consistent with the autotetraploid (*2n = 4x = 48*) genome structure of *V. corymbosum*. **(B)** Dot plot of Suziblue hap1 scaffolds (y-axis) aligned against the diploid *Vaccinium caesariense* W85-20 reference genome (x-axis; Edger et al., 2022) using nucmer. Continuous diagonal synteny blocks confirm collinearity between the assembled Suziblue chromosomes and W85-20 pseudomolecules, validating the chromosome identity assignments (Chr.01–Chr.12) used throughout this study. Source: `Blueberry_Genomes/Suziblue/SuziblueHiCMap.pdf` and `map_Suziblue_Haps_q0_scaffolds_final_to_V_caesariense_W85-20_P0_v2.png`.

### Supplementary Figure S6. Chr.05:47–48 Mb is a novel locus distinct from previously published blueberry chilling requirement QTL.

Synteny-anchored comparison of Chr.05 between the Draper assembly (VaccDscaff7, 41.7 Mb) and the Suziblue hap1 assembly (Chr.05.1, 48.6 Mb). The published chilling-requirement QTL from Lobos et al. (2021; *Frontiers in Plant Science*) maps to bins 113–117 (~0.5–4.3 Mb) and bin 118 (~27.7 Mb) on Draper VaccDscaff7. Through synteny alignment (grey trapezoids; confirmed by three-genome riparian plot), these positions correspond to the proximal ~4 Mb and ~28 Mb regions of Suziblue Chr.05, respectively — approximately 43 Mb distal from the Chr.05:47–48.6 Mb locus identified in this study (dark red). The two loci are therefore genuinely distinct and represent independent genetic contributions to SHB phenology on Chr.05. Source: `Chr.05_compare_lobo/chr05_locus_comparison.pdf`.

---

*Figure file format: PDF (vector) and PNG (300 dpi) at 180 mm width (two-column, per Plant Genome guidelines). All scripts in `R_Scripts/`.*
