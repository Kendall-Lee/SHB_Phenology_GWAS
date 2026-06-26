#!/usr/bin/env python3
"""
offprobe_gene_intersect.py

Classifies GWAS markers as genic or intergenic using the BRAKER hap1 gene
annotation, separately for probe-enriched and off-probe marker classes.

Addresses the question of why 83.9% of GWAS markers fall outside probe
intervals: ~31.7% of off-probe markers still land in annotated gene bodies,
exceeding the genome-wide genic fraction (24.3%), consistent with partial
capture of homeologous gene copies not covered by the Draper-derived probe
design.

INPUTS
------
DOS_FILE   : Dosage matrix (Marker, Chrom, Position, REF, ALT, samples...)
CLASS_FILE : marker_class_stats.csv — output of classify_markers_by_probe.R
             columns: Class (enriched / wgs_only), Missingness, MAF
GTF_FILE   : BRAKER hap1 GTF annotation (chromosome names Chr.XX.1)
FAI_FILE   : Suziblue hap1 .fai for genome size denominator

OUTPUT
------
Prints a summary table to stdout. Also writes:
  - Tables/offprobe_genic_markers.bed  : BED4 of all markers with probe_class
                                         and gene_overlap columns
  - docs/offprobe_gene_intersect_results.md : results summary

USAGE
-----
python3 Python/offprobe_gene_intersect.py
"""

import bisect
import os

REPO_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BED_OUT    = os.path.join(REPO_DIR, "Tables", "offprobe_genic_markers.bed")

BASE       = "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
DOS_FILE   = os.path.join(BASE, "Genotype_Data/SHB_tet.DP.8.maxmiss.3.minminor.25")
CLASS_FILE = os.path.join(BASE, "QTL_Tables/marker_class_stats.csv")
GTF_FILE   = "/Users/kendalllee/Documents/Blueberry/Gene_Info/braker.hap1.gtf"
FAI_FILE   = "/Users/kendalllee/Documents/Blueberry/Enrichment/Suziblue_hap1.fa.fai"

# ── Load marker positions ──────────────────────────────────────────────────────
print("Loading marker positions...")
markers = []
with open(DOS_FILE) as f:
    next(f)
    for line in f:
        p = line.split(',')
        markers.append((p[1], int(p[2])))   # (chrom, pos)
print(f"  {len(markers):,} markers")

classes = []
with open(CLASS_FILE) as f:
    next(f)
    for line in f:
        classes.append(line.strip().split(',')[0])
print(f"  {len(classes):,} class labels")
assert len(markers) == len(classes), "Marker / class count mismatch"

# ── Load BRAKER hap1 gene intervals ───────────────────────────────────────────
# GTF uses Chr.XX.1 naming; strip trailing .1 to match dosage matrix Chr.XX
print("Loading BRAKER hap1 gene intervals...")
gene_by_chr = {}
with open(GTF_FILE) as f:
    for line in f:
        if line.startswith('#'):
            continue
        p = line.strip().split('\t')
        if len(p) < 9 or p[2] != 'gene':
            continue
        raw = p[0]
        parts = raw.split('.')
        if parts[-1] != '1':
            continue                         # hap1 only
        base = '.'.join(parts[:2])           # "Chr.01"
        start = int(p[3]) - 1               # convert to 0-based
        end   = int(p[4])
        gene_by_chr.setdefault(base, []).append((start, end))

for c in gene_by_chr:
    gene_by_chr[c].sort()

n_genes = sum(len(v) for v in gene_by_chr.values())
print(f"  {n_genes:,} gene features across {len(gene_by_chr)} chromosomes")

# ── Genome size from FAI ───────────────────────────────────────────────────────
genome_bp = 0
with open(FAI_FILE) as f:
    for line in f:
        genome_bp += int(line.split('\t')[1])

total_gene_bp = sum(e - s for ivs in gene_by_chr.values() for s, e in ivs)
genic_frac = total_gene_bp / genome_bp
print(f"  Genome size: {genome_bp:,} bp  |  Gene body bp: {total_gene_bp:,}  ({100*genic_frac:.1f}%)")

# ── Interval overlap (binary search) ──────────────────────────────────────────
def in_gene(chrom, pos, intervals):
    if chrom not in intervals:
        return False
    ivs = intervals[chrom]
    idx = bisect.bisect_right([s for s, e in ivs], pos) - 1
    if idx < 0:
        return False
    return ivs[idx][1] >= pos

# ── Load marker IDs for BED output ────────────────────────────────────────────
marker_ids = []
with open(DOS_FILE) as f:
    next(f)
    for line in f:
        marker_ids.append(line.split(',')[0])

# ── Classify markers + write BED ──────────────────────────────────────────────
print("Classifying markers by gene overlap...")
res = {
    'enriched': [0, 0],   # [in_gene, off_gene]
    'wgs_only': [0, 0],
}

os.makedirs(os.path.dirname(BED_OUT), exist_ok=True)
with open(BED_OUT, 'w') as bed:
    bed.write("chrom\tstart\tend\tmarker\tprobe_class\tgene_overlap\n")
    for (chrom, pos), cls, mid in zip(markers, classes, marker_ids):
        hit = in_gene(chrom, pos, gene_by_chr)
        res[cls][0 if hit else 1] += 1
        gene_label = "genic" if hit else "intergenic"
        bed.write(f"{chrom}\t{pos-1}\t{pos}\t{mid}\t{cls}\t{gene_label}\n")

print(f"  BED written: {BED_OUT}")

# ── Report ─────────────────────────────────────────────────────────────────────
print("\n=== Off-probe marker overlap with BRAKER hap1 gene bodies ===")
print(f"{'Class':<12} {'In gene':>8} {'Off gene':>9} {'Total':>8} {'% genic':>9}")
print("-" * 52)
for cls in ['enriched', 'wgs_only']:
    ig, og = res[cls]
    tot = ig + og
    print(f"{cls:<12} {ig:>8,} {og:>9,} {tot:>8,} {100*ig/tot:>8.1f}%")

print()
ig2, og2 = res['wgs_only']
tot2 = ig2 + og2
print(f"Random expectation (genome-wide genic fraction): {100*genic_frac:.1f}%")
print(f"Off-probe markers in gene bodies: {ig2:,} / {tot2:,} ({100*ig2/tot2:.1f}%)")
print(f"Fold enrichment vs. random: {(ig2/tot2)/genic_frac:.2f}x")
