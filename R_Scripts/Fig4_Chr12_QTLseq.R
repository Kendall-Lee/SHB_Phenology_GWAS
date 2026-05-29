## Fig4_Chr12_QTLseq.R
## Figure 4 (was Fig 5) -- Chr.12 QTL-seq major-effect locus
## Two-panel main text figure:
##   Panel A: Regional SV QTL-seq lollipop across Chr.12 9.5-11.5 Mb
##   Panel B: Allele effect -- bulk AF at peak SV + per-sample proxy boxplot
##
## Source panels (pre-built, loaded as images):
##   Panel A: Chr.12_investigation/Fig_Chr12_Lollipop.png
##   Panel B: Chr.12_investigation/Chr12_SV_BulkAF.png
##
## The TE mechanistic story (Fig_Chr12_ManuscriptLocus) goes to Sup Fig S5.
##
## Output: Linear_MS/Fig4_Chr12_QTLseq.pdf + .png

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(cowplot)
  library(ggplot2)
})

CHR12DIR <- "/Users/kendalllee/Documents/Blueberry/Chr.12_investigation"
OUT_DIR  <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Linear_MS"

# ── Load pre-built panel images via rasterGrob (no magick needed) -------------
load_png_panel <- function(path) {
  img <- readPNG(path)
  ggdraw() + draw_grob(grid::rasterGrob(img, width = 1, height = 1))
}

panel_A <- load_png_panel(file.path(CHR12DIR, "Fig_Chr12_Lollipop.png"))
panel_B <- load_png_panel(file.path(CHR12DIR, "Chr12_SV_BulkAF.png"))

# ── Title strip ---------------------------------------------------------------
title_strip <- ggdraw() +
  draw_label(
    paste0("Figure 4. Chr.12:9.5-11.5 Mb QTL-seq locus -- ",
           "TE-embedded PRMT (g132348) controlling days to flowering"),
    fontface = "bold", size = 10.5, x = 0.01, hjust = 0
  )

# ── Caption block -------------------------------------------------------------
caption_text <- paste0(
  "A: SV QTL-seq R2 across Chr.12 9.5-11.5 Mb (20 FDR-significant markers; padj<0.05, N>=20 SV carriers). ",
  "Peak DEL 570 bp (red dashed, R2=50%) overlaps Helitron+CACTA TEs inside gene g132348 (PRMT). ",
  "TE rug shown at bottom (gold). ",
  "B: Left -- aggregate mean allele dosage at peak SV (Chr.12:10308052) by DTFlower phenotype bulk; ",
  "early bulk AF=0.28 vs late bulk AF=1.00 (deltaDS=2.87). ",
  "69% of LRLP samples cannot be individually genotyped at this TE-embedded position; ",
  "QTL-seq detects signal via bulk read-level allele frequencies. ",
  "Right -- per-sample tetraploid dosage vs DTFlower BLUE at proxy marker Chr.12:10399365 ",
  "(91 kb from peak, R2=25.5%, n=81; KW p=9.5e-05). ",
  "See Supplementary Figure S5 for full TE content and short-read vs long-read variant density analysis."
)
caption_strip <- ggdraw() +
  draw_label(caption_text, x = 0.01, y = 0.7, hjust = 0, vjust = 1,
             size = 7.5, color = "grey30", lineheight = 1.35)

# ── Assemble ------------------------------------------------------------------
fig4 <- plot_grid(
  title_strip,
  plot_grid(
    panel_A, panel_B,
    ncol = 2, rel_widths = c(1.1, 1),
    labels = c("A", "B"), label_size = 14, label_fontface = "bold",
    label_x = 0.01, label_y = 0.99
  ),
  caption_strip,
  nrow = 3, rel_heights = c(0.04, 1, 0.12)
)

ggsave(file.path(OUT_DIR, "Fig4_Chr12_QTLseq.pdf"),
       fig4, width = 7.09, height = 3.5, units = "in")
ggsave(file.path(OUT_DIR, "Fig4_Chr12_QTLseq.png"),
       fig4, width = 7.09, height = 3.5, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "Fig4_Chr12_QTLseq.pdf")))
