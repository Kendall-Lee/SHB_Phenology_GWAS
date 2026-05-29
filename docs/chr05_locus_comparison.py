import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch
from matplotlib.path import Path
import matplotlib.patheffects as pe
import numpy as np

# ── chromosome lengths ──────────────────────────────────────────────────────
SUZI_LEN   = 48.6   # Mb  Suziblue Chr.05.1
VACCD7_LEN = 41.7   # Mb  Draper VaccDscaff7

# Published locus (bins 113-119 on VaccDscaff7)
PUB_CLUSTER1_START = 0.5   # Mb
PUB_CLUSTER1_END   = 4.27  # Mb
PUB_CLUSTER2_START = 27.4  # Mb  (bin 118 ~27.7)
PUB_CLUSTER2_END   = 28.0  # Mb

# User's locus on Suziblue Chr.05
MY_START = 47.0   # Mb
MY_END   = 48.6   # Mb  (end of chromosome)

# ── display geometry ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(13, 5.5))
ax.set_xlim(0, 55)
ax.set_ylim(0, 10)
ax.axis('off')

SCALE   = 55 / max(SUZI_LEN, VACCD7_LEN)   # Mb → display units
Y_TOP   = 7.8   # Draper VaccDscaff7 bar centre
Y_BOT   = 2.2   # Suziblue Chr.05 bar centre
H       = 0.55  # bar half-height

# ── colour palette ──────────────────────────────────────────────────────────
DRAPER_COL  = "#4472C4"
SUZI_COL    = "#548235"
PUB_COL     = "#E07B39"   # orange  = published QTL
MY_COL      = "#C00000"   # red     = this study
RIB_COL     = "#CCCCCC"

# ── helper: draw a filled synteny ribbon ────────────────────────────────────
def ribbon(ax, x1s, x1e, x2s, x2e, y_top, y_bot, colour, alpha=0.35, zorder=1):
    """Draw a trapezoid ribbon between two chromosome bars."""
    xs = [x1s, x1e, x2e, x2s, x1s]
    ys = [y_top, y_top, y_bot, y_bot, y_top]
    ax.fill(xs, ys, color=colour, alpha=alpha, zorder=zorder, linewidth=0)

# ── synteny ribbon (bulk collinear block) ───────────────────────────────────
# VaccDscaff7 full → Suziblue Chr.05 full (roughly collinear)
ribbon(ax,
       x1s=0,                       x1e=VACCD7_LEN * SCALE,
       x2s=0,                       x2e=SUZI_LEN   * SCALE,
       y_top=Y_TOP - H, y_bot=Y_BOT + H,
       colour=RIB_COL, alpha=0.25)

# ── highlight ribbons for the two loci ──────────────────────────────────────
# Equivalent position of user's locus on VaccDscaff7 (for the ribbon endpoint)
my_vaccd7_start = (MY_START / SUZI_LEN) * VACCD7_LEN
my_vaccd7_end   = (MY_END   / SUZI_LEN) * VACCD7_LEN

ribbon(ax,
       x1s=my_vaccd7_start * SCALE, x1e=my_vaccd7_end * SCALE,
       x2s=MY_START        * SCALE, x2e=MY_END         * SCALE,
       y_top=Y_TOP - H, y_bot=Y_BOT + H,
       colour=MY_COL, alpha=0.30, zorder=2)

# Equivalent position of published cluster 1 on Suziblue Chr.05
pub_suzi_start = (PUB_CLUSTER1_START / VACCD7_LEN) * SUZI_LEN
pub_suzi_end   = (PUB_CLUSTER1_END   / VACCD7_LEN) * SUZI_LEN

ribbon(ax,
       x1s=PUB_CLUSTER1_START * SCALE, x1e=PUB_CLUSTER1_END * SCALE,
       x2s=pub_suzi_start      * SCALE, x2e=pub_suzi_end      * SCALE,
       y_top=Y_TOP - H, y_bot=Y_BOT + H,
       colour=PUB_COL, alpha=0.30, zorder=2)

# ── chromosome bars ─────────────────────────────────────────────────────────
def chrom_bar(ax, length, y, colour, label, label_side='left'):
    w = length * SCALE
    rect = Rectangle((0, y - H), w, 2*H,
                      linewidth=1.2, edgecolor='#333333',
                      facecolor=colour, alpha=0.75, zorder=3)
    ax.add_patch(rect)
    # rounded end cap (right)
    cap = mpatches.Ellipse((w, y), 0.18, 2*H,
                            linewidth=1.2, edgecolor='#333333',
                            facecolor=colour, alpha=0.75, zorder=4)
    ax.add_patch(cap)
    if label_side == 'left':
        ax.text(-0.4, y, label, ha='right', va='center', fontsize=10,
                fontweight='bold', color='#222222')
    else:
        ax.text(w + 0.5, y, label, ha='left', va='center', fontsize=10,
                fontweight='bold', color='#222222')

chrom_bar(ax, VACCD7_LEN, Y_TOP, DRAPER_COL,
          "Draper\nVaccDscaff7\n(41.7 Mb)")
chrom_bar(ax, SUZI_LEN,   Y_BOT, SUZI_COL,
          "Suziblue\nChr.05.1\n(48.6 Mb)")

# ── locus highlight boxes ────────────────────────────────────────────────────
def locus_box(ax, x_start, x_end, y, col, zorder=5):
    w = (x_end - x_start) * SCALE
    rect = Rectangle((x_start * SCALE, y - H - 0.04), w, 2*H + 0.08,
                      linewidth=2.0, edgecolor=col, facecolor=col,
                      alpha=0.55, zorder=zorder)
    ax.add_patch(rect)

# Published cluster 1 on VaccDscaff7
locus_box(ax, PUB_CLUSTER1_START, PUB_CLUSTER1_END, Y_TOP, PUB_COL)
# Published cluster 2 on VaccDscaff7 (bin 118)
locus_box(ax, PUB_CLUSTER2_START, PUB_CLUSTER2_END, Y_TOP, PUB_COL)
# User's locus on Suziblue Chr.05
locus_box(ax, MY_START, MY_END, Y_BOT, MY_COL)

# ── annotations ─────────────────────────────────────────────────────────────
# Published QTL label — cluster 1
mid1 = (PUB_CLUSTER1_START + PUB_CLUSTER1_END) / 2 * SCALE
ax.annotate("Published QTL\n(bins 113–117, ~0.5–4.3 Mb)",
            xy=(mid1, Y_TOP + H + 0.08),
            xytext=(mid1, Y_TOP + H + 1.45),
            ha='center', va='bottom', fontsize=8.5, color=PUB_COL, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=PUB_COL, lw=1.5))

# Published QTL label — cluster 2 (bin 118)
mid2 = (PUB_CLUSTER2_START + PUB_CLUSTER2_END) / 2 * SCALE
ax.annotate("Bin 118\n(~27.7 Mb)",
            xy=(mid2, Y_TOP + H + 0.08),
            xytext=(mid2 + 2.5, Y_TOP + H + 1.45),
            ha='center', va='bottom', fontsize=8.0, color=PUB_COL, fontstyle='italic',
            arrowprops=dict(arrowstyle='->', color=PUB_COL, lw=1.2))

# User's locus label
mid_me = (MY_START + MY_END) / 2 * SCALE
ax.annotate("This study\n(47–48.6 Mb)",
            xy=(mid_me, Y_BOT - H - 0.08),
            xytext=(mid_me, Y_BOT - H - 1.5),
            ha='center', va='top', fontsize=8.5, color=MY_COL, fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=MY_COL, lw=1.5))

# Distinct loci callout
ax.text(27.5, (Y_TOP + Y_BOT) / 2,
        "~43 Mb apart\n(distinct loci)",
        ha='center', va='center', fontsize=9.5, color='#555555',
        fontweight='bold', fontstyle='italic',
        bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#AAAAAA', lw=1.2))

# ── scale bars ───────────────────────────────────────────────────────────────
def scale_bar(ax, x0, y, length_mb, label, colour):
    x1 = x0 + length_mb * SCALE
    ax.annotate('', xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle='<->', color=colour, lw=1.2))
    ax.text((x0 + x1) / 2, y - 0.22, label,
            ha='center', va='top', fontsize=7.5, color=colour)

scale_bar(ax, 0.1, Y_TOP - H - 0.7, 10, '10 Mb', DRAPER_COL)
scale_bar(ax, 0.1, Y_BOT - H - 0.7, 10, '10 Mb', SUZI_COL)

# ── tick marks every 10 Mb ──────────────────────────────────────────────────
for mb in range(0, int(VACCD7_LEN) + 1, 10):
    x = mb * SCALE
    ax.plot([x, x], [Y_TOP - H, Y_TOP - H - 0.15], color='#555555', lw=0.8, zorder=5)
    ax.text(x, Y_TOP - H - 0.22, f'{mb}', ha='center', va='top', fontsize=7, color='#555555')

for mb in range(0, int(SUZI_LEN) + 1, 10):
    x = mb * SCALE
    ax.plot([x, x], [Y_BOT - H, Y_BOT - H - 0.15], color='#555555', lw=0.8, zorder=5)
    ax.text(x, Y_BOT - H - 0.22, f'{mb}', ha='center', va='top', fontsize=7, color='#555555')

# ── legend ───────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(facecolor=PUB_COL, alpha=0.7, label='Published QTL (Lobos et al. 2021, FPLS)'),
    mpatches.Patch(facecolor=MY_COL,  alpha=0.7, label='This study — chilling QTL locus'),
    mpatches.Patch(facecolor=RIB_COL, alpha=0.5, label='Synteny block (Suziblue ↔ Draper)'),
]
ax.legend(handles=legend_patches, loc='upper right',
          bbox_to_anchor=(1.0, 0.45), fontsize=8.5,
          framealpha=0.9, edgecolor='#AAAAAA')

# ── supplementary reference note ────────────────────────────────────────────
ax.text(0.99, 0.01,
        "Synteny confirmed by three-genome riparian plot (Supplementary Figure S1)",
        transform=ax.transAxes, ha='right', va='bottom',
        fontsize=7.5, color='#777777', fontstyle='italic')

# ── title ────────────────────────────────────────────────────────────────────
ax.set_title(
    "Chilling requirement loci on Blueberry Chr. 05 / VaccDscaff7 are distinct",
    fontsize=12, fontweight='bold', pad=8, color='#1a1a1a')

plt.tight_layout()
out = "/Users/kendalllee/Documents/chr05_locus_comparison.pdf"
plt.savefig(out, dpi=200, bbox_inches='tight')
out_png = "/Users/kendalllee/Documents/chr05_locus_comparison.png"
plt.savefig(out_png, dpi=200, bbox_inches='tight')
print(f"Saved: {out}")
print(f"Saved: {out_png}")
