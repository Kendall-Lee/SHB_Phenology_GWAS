## Fig5_Yr25_Stress_QTL.R  (consolidated from Fig4 + Fig6)
## Single-panel figure: PVE bars for yr.2025 stress-responsive and stable loci.
## Panels A (QTL counts) and B (Jan min temp) removed 2026-06-05:
##   - Panel B redundant with Figure 1B (frost/freeze event counts).
##   - Panel A better stated in results text.
##   - Panel D (Chr.03 allele effect) removed 2026-05-29: n=5 rare allele, off-target marker.
##
## The Plant Genome specs: double-column = 180 mm (7.087 in), 300 DPI, min 8 pt font at print.
## Output: Linear_MS/Figures/Fig5_Yr25_Stress_QTL.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS/Figures")

pub_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text     = element_text(color = "black"),
    axis.title    = element_text(color = "black", face = "bold"),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

trait_colors <- c(
  "DTFlower"    = "#0D3B6E",
  "DTFruit"     = "#1A6FAF",
  "Flow2Fruit"  = "#2FAAB3",
  "FruitWeight" = "#7ECAE0"
)

# ── Panel: PVE bars -----------------------------------------------------------
stress_qtl <- data.table(
  Label = c("Flow2Fruit (2025)\nChr.03:38.4 Mb",
            "DTFlower (2025)\nChr.11:45.5 Mb",
            "DTFruit (2025)\nChr.07:25.0 Mb"),
  PVE   = c(13.26, 6.11, 4.58),
  Trait = c("Flow2Fruit","DTFlower","DTFruit"),
  plab  = c("p = 2.3e-18", "p = 5.8e-9", "p = 5.4e-7"),
  Type  = "Stress-specific (yr.2025 only)"
)
stable_in_yr25 <- data.table(
  Label = c("DTFlower (2025)\nChr.05:48.3 Mb [BLUE=7.79%]",
            "DTFlower (2025)\nChr.09:39.6 Mb [BLUE=7.54%]"),
  PVE   = c(3.00, 0.58),
  Trait = c("DTFlower","DTFlower"),
  plab  = c("p = 5.2e-5", "p = 0.083 (NS)"),
  Type  = "Phenology locus in yr.2025"
)
all_bars <- rbindlist(list(stress_qtl, stable_in_yr25))
all_bars[, Label := factor(Label, levels = rev(c(stress_qtl$Label, stable_in_yr25$Label)))]

p_c <- ggplot(all_bars, aes(y = Label, x = PVE, fill = Trait, alpha = Type)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.55) +
  geom_text(aes(label = sprintf("%.2f%%  %s", PVE, plab)),
            hjust = -0.05, size = 3.0) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Stress-specific (yr.2025 only)" = 1.0,
               "Phenology locus in yr.2025"     = 0.38),
    name = NULL) +
  # x-axis runs to 30 so inline labels clear the plot area
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 5, 10, 15, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "PVE (%) in Year 2025", y = NULL,
       title = "Year 2025 stress-responsive and stable locus performance",
       # manually wrap subtitle so no line exceeds ~90 chars (fits 180 mm at 8 pt)
       subtitle = paste0(
         "Solid bars: loci detected exclusively in yr.2025 (absent from stable-year BLUE_exc.25 scan).\n",
         "Faded bars: stable phenology loci attenuated under frost stress.\n",
         "Chr.03:38.4 Mb (Flow2Fruit) is exploratory: n=5 in rare dosage class, requires replication.")) +
  pub_theme +
  theme(axis.text.y      = element_text(size = 9, lineheight = 1.3),
        plot.margin      = margin(t = 6, r = 10, b = 4, l = 4, unit = "pt"),
        legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        plot.subtitle    = element_text(size = 8, color = "grey40", lineheight = 1.35))

# ── Save — TPG double-column: 180 mm = 7.087 in, 300 DPI -----------------------
W <- 7.087; H <- 4.5
ggsave(file.path(OUT_DIR, "Fig5_Yr25_Stress_QTL.pdf"), p_c, width = W, height = H, units = "in")
ggsave(file.path(OUT_DIR, "Fig5_Yr25_Stress_QTL.png"), p_c, width = W, height = H,
       units = "in", dpi = 300)
cat(sprintf("Saved: %s (%.3f x %.1f in, 300 dpi)\n",
            file.path(OUT_DIR, "Fig5_Yr25_Stress_QTL.pdf"), W, H))
