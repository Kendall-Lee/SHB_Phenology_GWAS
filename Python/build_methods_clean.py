#!/usr/bin/env python3
"""
Clean study design flowchart matching KhufuPAN_supp_fig style:
- White boxes, black borders
- Light-coloured panel backgrounds
- Black arrows
- Clean sans-serif, black text
Output: Figures/StudyDesign_Diagram.pdf + .png
"""

import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['pdf.fonttype'] = 42

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

W, H = 13.33, 7.5
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Panel background colours (matching reference) ─────────────────────────────
BLU  = '#9BADC0'   # light steel-blue  (GWAS panel)
GRY  = '#B0B0AE'   # medium grey       (QTL-Seq panel)
AGRY = '#C8C8C6'   # light grey        (Analyzed-by panel)
DRK  = '#525252'   # dark grey         (overlap panel)

# ── Helpers ───────────────────────────────────────────────────────────────────
def panel(l, b, w, h, fc, label='', lw=1.8):
    """Rounded-corner panel background."""
    ax.add_patch(mpatches.FancyBboxPatch(
        (l, b), w, h, boxstyle='round,pad=0.12',
        fc=fc, ec='black', lw=lw, zorder=1))
    if label:
        ax.text(l+0.18, b+h-0.28, label,
                fontsize=13, fontweight='bold', va='top', ha='left',
                color='black' if fc != DRK else 'white', zorder=5)

def rbox(cx, cy, w, h, label, sub='', fs=10, bold=False, fc='white', ec='black', lw=1.4):
    """Rectangle box."""
    ax.add_patch(mpatches.FancyBboxPatch(
        (cx-w/2, cy-h/2), w, h, boxstyle='square,pad=0.0',
        fc=fc, ec=ec, lw=lw, zorder=3))
    fw = 'bold' if bold else 'normal'
    dy = 0.10 if sub else 0
    ax.text(cx, cy+dy, label, ha='center', va='center',
            fontsize=fs, fontweight=fw, color='black', zorder=4)
    if sub:
        ax.text(cx, cy-0.20, sub, ha='center', va='center',
                fontsize=fs-2, color='#444444', zorder=4)

def para(cx, cy, w, h, label, sub='', fs=11, bold=True, skew=0.18):
    """Parallelogram (data node) – matches reference style."""
    s = skew * h
    pts = np.array([
        [cx - w/2 + s, cy + h/2],
        [cx + w/2 + s, cy + h/2],
        [cx + w/2 - s, cy - h/2],
        [cx - w/2 - s, cy - h/2],
    ])
    ax.add_patch(plt.Polygon(pts, closed=True,
                             fc='white', ec='black', lw=1.5, zorder=3))
    fw = 'bold' if bold else 'normal'
    dy = 0.10 if sub else 0
    ax.text(cx, cy+dy, label, ha='center', va='center',
            fontsize=fs, fontweight=fw, color='black', zorder=4)
    if sub:
        ax.text(cx, cy-0.22, sub, ha='center', va='center',
                fontsize=fs-2, color='#333333', zorder=4)

def arr(x1, y1, x2, y2, style='->', lw=1.5, color='black'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color,
                                lw=lw, mutation_scale=14),
                zorder=5)

def txt(x, y, s, fs=10, bold=False, ha='center', va='center', color='black'):
    fw = 'bold' if bold else 'normal'
    ax.text(x, y, s, fontsize=fs, fontweight=fw,
            ha=ha, va=va, color=color, zorder=5)

def hline(y, x0, x1, lw=1.0, color='black', ls='--'):
    ax.plot([x0, x1], [y, y], color=color, lw=lw, ls=ls, zorder=4)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL 1 — GWAS (light blue, left)
# ══════════════════════════════════════════════════════════════════════════════
panel(0.15, 3.80, 6.15, 3.55, BLU, '1.  GWAS')

# Dataset box
rbox(3.22, 6.80, 5.20, 0.56,
     '665 enriched SHB accessions', 'Alapaha, GA  ·  2023, 2024, 2025',
     fs=11, bold=True)

# Two method sub-boxes
rbox(1.55, 5.50, 2.40, 0.52, 'Pangenome', fs=11, bold=True)
rbox(4.90, 5.50, 2.40, 0.52, 'Linear',    fs=11, bold=True)

# Arrows: dataset → methods
arr(3.22, 6.52, 1.55, 5.76)
arr(3.22, 6.52, 4.90, 5.76)

# Output parallelogram
para(3.22, 4.35, 3.80, 0.62, 'Regions / QTL', 'GWAS-identified', fs=11)

# Arrows: methods → output
arr(1.55, 5.24, 3.22, 4.67)
arr(4.90, 5.24, 3.22, 4.67)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL 2 — QTL-Seq (medium grey, right)
# ══════════════════════════════════════════════════════════════════════════════
panel(7.03, 3.80, 6.15, 3.55, GRY, '2.  QTL-Seq')

# Dataset box
rbox(10.10, 6.80, 5.20, 0.56,
     '95 LRLP SHB accessions', 'Phenotype extremes yr.23 + yr.24  ·  PacBio Revio',
     fs=11, bold=True)

# Two method sub-boxes
rbox(8.43, 5.50, 2.40, 0.52, 'Linear',    fs=11, bold=True)
rbox(11.78, 5.50, 2.40, 0.52, 'Pangenome', fs=11, bold=True)

# Arrows
arr(10.10, 6.52, 8.43, 5.76)
arr(10.10, 6.52, 11.78, 5.76)

# Output parallelogram
para(10.10, 4.35, 3.80, 0.62, 'Regions / QTL', 'QTL-Seq-identified', fs=11)

arr(8.43,  5.24, 10.10, 4.67)
arr(11.78, 5.24, 10.10, 4.67)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL 3 — Analyzed by (light grey, full width)
# ══════════════════════════════════════════════════════════════════════════════
panel(0.15, 0.90, 13.03, 2.72, AGRY, 'Analyzed by', lw=1.8)

# Divider lines between columns
hline(3.30, 0.35, 13.00, lw=0.8, color='#666666')
for xd in [4.60, 8.90]:
    ax.plot([xd, xd], [0.98, 3.28], color='#888888', lw=0.8, ls='--', zorder=4)

# Column headers
for cx, hdr in [(2.30, 'Trait'), (6.75, 'Year  /  BLUE'), (11.10, 'Model  &  Threshold')]:
    txt(cx, 3.48, hdr, fs=12, bold=True)

# Column content
trait_items   = ['DTFlower', 'DTFruit', 'Flow2Fruit', 'FruitWeight']
year_items    = ['yr.2023  (individual)', 'yr.2024  (individual)',
                 'yr.2025  (stress)', 'BLUE  (all 3 years)', 'BLUE  excl. yr.2025']
model_items   = ['Additive', '1-dom-ref  /  1-dom-alt',
                 '2-dom-ref  /  2-dom-alt',
                 '─────────────────────',
                 'Bonferroni', 'Meff  (effective # SNPs)']

for j, item in enumerate(trait_items):
    txt(2.30, 3.04 - j*0.47, f'●  {item}', fs=10, ha='center')

for j, item in enumerate(year_items):
    txt(6.75, 3.04 - j*0.40, f'●  {item}', fs=10, ha='center')

for j, item in enumerate(model_items):
    c = '#888888' if '──' in item else 'black'
    txt(11.10, 3.04 - j*0.38, f'●  {item}' if '──' not in item else item,
        fs=10, ha='center', color=c)

# Note
txt(6.67, 0.70,
    '*  The more methods, models, and years that support a QTL  ->  the stronger the evidence',
    fs=10, bold=True, color='#222222')

# ══════════════════════════════════════════════════════════════════════════════
# Arrows: panel outputs → "Analyzed by" and then outputs
# ══════════════════════════════════════════════════════════════════════════════
# GWAS output → Analyzed by
arr(3.22, 4.04, 3.22, 3.62)
# QTL-Seq output → Analyzed by
arr(10.10, 4.04, 10.10, 3.62)

# Outer bracket arrows showing both feed into Analyzed by
ax.annotate('', xy=(6.67, 3.62), xytext=(3.22, 3.45),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.4,
                            connectionstyle='arc3,rad=0.0', mutation_scale=12), zorder=5)
ax.annotate('', xy=(6.67, 3.62), xytext=(10.10, 3.45),
            arrowprops=dict(arrowstyle='->', color='black', lw=1.4,
                            connectionstyle='arc3,rad=0.0', mutation_scale=12), zorder=5)

# ══════════════════════════════════════════════════════════════════════════════
# Bottom connector arrow
# ══════════════════════════════════════════════════════════════════════════════
# Big down arrow from Analyzed-by panel
ax.annotate('', xy=(6.67, 0.82), xytext=(6.67, 0.90),
            arrowprops=dict(arrowstyle='->', color='black', lw=2.0,
                            mutation_scale=16), zorder=5)

plt.tight_layout(pad=0.3)
OUT = '/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Linear_MS/Figures/StudyDesign_Diagram'
plt.savefig(OUT + '.pdf', format='pdf', bbox_inches='tight',
            facecolor='white', dpi=300)
plt.savefig(OUT + '.png', format='png', bbox_inches='tight',
            facecolor='white', dpi=300)
print(f'Saved: {OUT}.pdf / .png')
