## SupFigS7_Suziblue_Assembly.R
## Supplementary Figure S7 -- Suziblue hap1 reference genome assembly quality
## Panel A: Hi-C contact map (Juicebox; post-JBAT curation)
## Panel B: Suziblue hap1 vs. V. caesariense W85-20 synteny dot plot
##
## Output: Linear_MS/SupFigS7_Suziblue_Assembly.pdf + .png

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(cowplot)
  library(ggplot2)
})

SUZI_DIR <- "/Users/kendalllee/Documents/Blueberry_Genomes/Suziblue"
OUT_DIR  <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Linear_MS"

load_png_panel <- function(path) {
  img <- readPNG(path)
  ggdraw() + draw_grob(grid::rasterGrob(img, width = 1, height = 1))
}

panel_A <- load_png_panel(file.path(SUZI_DIR, "SuziblueHiCMap.png"))
panel_B <- load_png_panel(file.path(SUZI_DIR,
  "map_Suziblue_Haps_q0_scaffolds_final_to_V_caesariense_W85-20_P0_v2.png"))

title_strip <- ggdraw() + draw_label(
  paste0("Supplementary Figure S7. Suziblue hap1 chromosome-scale reference genome assembly"),
  fontface = "bold", size = 10, x = 0.01, hjust = 0
)

fig_s7 <- plot_grid(
  title_strip,
  plot_grid(
    panel_A, panel_B,
    ncol = 2, rel_widths = c(1.05, 1),
    labels = c("A", "B"), label_size = 13, label_fontface = "bold",
    label_x = 0.01, label_y = 0.99
  ),
  nrow = 2, rel_heights = c(0.04, 1)
)

ggsave(file.path(OUT_DIR, "SupFigS7_Suziblue_Assembly.pdf"),
       fig_s7, width = 7.09, height = 3.8, units = "in")
ggsave(file.path(OUT_DIR, "SupFigS7_Suziblue_Assembly.png"),
       fig_s7, width = 7.09, height = 3.8, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "SupFigS7_Suziblue_Assembly.pdf")))
