#!/usr/bin/env python3
import json, numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 1.2
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

with open('/tmp/lrlp_depth_breadth.json') as f:
    data = json.load(f)

samples  = sorted(data)
depths   = np.array([data[s]['depth']   for s in samples])
breadths = np.array([data[s]['breadth'] for s in samples])

fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_facecolor('white')

# ── Scatter ───────────────────────────────────────────────────────────────────
ax.scatter(depths, breadths, color='#1F77B4', s=55,
           edgecolors='white', linewidths=0.5, alpha=0.88, zorder=3)

# ── Stats annotation ──────────────────────────────────────────────────────────
# Exclude B825 outlier from summary stats shown
mask = depths <= 35
mn  = depths[mask].mean()
med = np.median(depths[mask])
ax.text(0.97, 0.97,
        f'n = {len(samples)}\nMean depth = {depths.mean():.1f}x\n'
        f'Median depth = {np.median(depths):.1f}x\n'
        f'Mean breadth = {breadths.mean():.1f}%\n'
        f'Median breadth = {np.median(breadths):.1f}%',
        transform=ax.transAxes, fontsize=9,
        va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', fc='#F5F5F5', ec='#BBBBBB', lw=1.0))


# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlabel('Mean Haploid Depth (x)', fontsize=12, fontweight='bold')
ax.set_ylabel('Genome Breadth Coverage (%)', fontsize=12, fontweight='bold')
ax.set_title('LRLP Panel — Sequencing Depth vs. Genome Coverage\n'
             'n = 95 SHB accessions  ·  Suziblue hap1 reference (548.4 Mb)',
             fontsize=12, fontweight='bold', pad=10)

ax.set_xlim(-1, 87)
ax.set_ylim(-2, 102)
ax.tick_params(labelsize=10)
ax.spines[['top','right']].set_visible(False)
ax.grid(True, alpha=0.25, lw=0.8, color='grey')

plt.tight_layout()
OUT = '/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Linear_MS/Figures/LRLP_Depth_vs_Breadth'
plt.savefig(OUT + '.pdf', bbox_inches='tight', facecolor='white', dpi=300)
plt.savefig(OUT + '.png', bbox_inches='tight', facecolor='white', dpi=300)
print(f'Saved: {OUT}.pdf / .png')
