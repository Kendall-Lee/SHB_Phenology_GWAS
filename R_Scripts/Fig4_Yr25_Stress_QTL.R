## Fig4_Yr25_Stress_QTL.R
## Year 2025 stress-specific QTL figure (3-panel: A=counts, B=temps, C=PVE bars)
## Updated 2026-05-04:
##   - Panel C now includes yr.25 R2 for top two stable loci as comparison bars
##   - Chr.05:48.3 Mb (LINEAR top, BLUE=7.79%) → PVE in yr.25 = 3.00%, p=5.2e-05
##   - Chr.09:39.6 Mb (PAN top, BLUE=7.54%)    → PVE in yr.25 = 0.58%, p=0.083

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(data.table)
})

trait_colors <- c(
  "DTFlower"    = "#0D3B6E",
  "DTFruit"     = "#1A6FAF",
  "Flow2Fruit"  = "#2FAAB3",
  "FruitWeight" = "#7ECAE0"
)
pub_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black", face = "bold"),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )
fig_title <- function(txt, size = 10.5)
  ggdraw() + draw_label(txt, fontface = "bold", size = size, x = 0.01, hjust = 0)

# ── Panel A: QTL counts by analysis ──────────────────────────────────────────
year_counts <- data.table(
  Dataset   = factor(
    c("GWAS Linear\n(2023-2025)", "GWAS Pan.\n(2023-2025)",
      "GWAS Linear\n(redoNo25)",  "GWAS Pan.\n(redoNo25)",
      "GWAS\n(2025 only)"),
    levels = c("GWAS Linear\n(2023-2025)", "GWAS Pan.\n(2023-2025)",
               "GWAS Linear\n(redoNo25)",  "GWAS Pan.\n(redoNo25)",
               "GWAS\n(2025 only)")
  ),
  N_QTL     = c(358, 345, 20, 23, 281),
  Year_type = c("All years","All years","Normal years","Normal years","2025 stress")
)

p4a <- ggplot(year_counts, aes(x = Dataset, y = N_QTL, fill = Year_type)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.65) +
  geom_text(aes(label = N_QTL), vjust = -0.35, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c("All years" = "#1A6FAF", "Normal years" = "#7ECAE0",
               "2025 stress" = "#D84315"),
    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Significant QTL detected", title = "A  QTL counts by analysis") +
  pub_theme +
  theme(legend.position = c(0.72, 0.82), legend.background = element_blank(),
        axis.text.x = element_text(size = 8, lineheight = 1.2))

# ── Panel B: January minimum temperature ─────────────────────────────────────
temp_data <- data.table(
  Year     = factor(c("2023","2024","2025")),
  Min_Temp = c(3.8, -1.9, -5.3)
)

p4b <- ggplot(temp_data, aes(x = Year, y = Min_Temp, fill = Min_Temp)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.55) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = paste0(Min_Temp, "\u00b0C"),
                vjust = ifelse(Min_Temp >= 0, -0.3, 1.3)),
            size = 3.3, fontface = "bold") +
  scale_fill_gradient2(low = "#0D3B6E", mid = "grey92", high = "#D84315",
                       midpoint = 0, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.15))) +
  labs(x = "Year", y = "Min. Jan. Temp (\u00b0C)", title = "B  January minimum temperature") +
  pub_theme

# ── Panel C: PVE bars — stress-specific + phenology loci in yr.25 ────────────
# Stress-specific loci (detected only / predominantly in yr.2025)
stress_qtl <- data.table(
  Label = c("Flow2Fruit (2025)\nChr.03:38.4 Mb",
            "DTFlower (2025)\nChr.11:45.5 Mb",
            "DTFruit (2025)\nChr.07:25.0 Mb"),
  PVE   = c(13.26, 6.11, 4.58),
  Trait = c("Flow2Fruit","DTFlower","DTFruit"),
  plab  = c("p = 2.3e-18", "p = 5.8e-9", "p = 5.4e-7"),
  Type  = "Stress-specific (yr.2025 only)"
)

# Top stable phenology loci — their performance in yr.2025
# Chr.05_48252083: BLUE_exc.25 PVE=7.79%, yr.25 PVE=3.00% (p=5.2e-05)
# Chr.09_39605986: BLUE_exc.25 PVE=7.54%, yr.25 PVE=0.58% (p=0.083, NS)
stable_in_yr25 <- data.table(
  Label = c("DTFlower (2025)\nChr.05:48.3 Mb [LINEAR top; BLUE=7.79%]",
            "DTFlower (2025)\nChr.09:39.6 Mb [PAN top; BLUE=7.54%]"),
  PVE   = c(3.00, 0.58),
  Trait = c("DTFlower","DTFlower"),
  plab  = c("p = 5.2e-5", "p = 0.083 (NS)"),
  Type  = "Phenology locus in yr.2025"
)

all_bars <- rbindlist(list(stress_qtl, stable_in_yr25))
# factor in display order: stress loci on top, phenology at bottom
all_bars[, Label := factor(Label, levels = rev(c(
  stress_qtl$Label, stable_in_yr25$Label
)))]

p4c <- ggplot(all_bars, aes(y = Label, x = PVE, fill = Trait, alpha = Type)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.55) +
  geom_text(aes(label = sprintf("%.2f%%  %s", PVE, plab)),
            hjust = -0.05, size = 2.9) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Stress-specific (yr.2025 only)" = 1.0,
               "Phenology locus in yr.2025"     = 0.38),
    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.65))) +
  labs(x = "PVE (%) in Year 2025",
       y = NULL,
       title = "C  Top stress-responsive QTL and phenology-locus performance in Year 2025",
       subtitle = paste0(
         "Faded bars: phenology loci with high BLUE PVE (7.5-7.8%) that are attenuated in the 2025 stress year.\n",
         "Chr.09:39.6 Mb is not significant in yr.2025 (p>0.05); Chr.05:48.3 Mb retains moderate effect (3.0%, p<0.001).")) +
  pub_theme +
  theme(axis.text.y = element_text(size = 8, lineheight = 1.2),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 7.5, color = "grey40", lineheight = 1.3))

# ── Assemble ──────────────────────────────────────────────────────────────────
fig4 <- plot_grid(
  fig_title("Year 2025 Stress-Specific QTL and Phenology Locus Performance"),
  plot_grid(
    plot_grid(p4a, p4b, ncol = 2, rel_widths = c(1.7, 1)),
    p4c,
    nrow = 2, rel_heights = c(1.0, 1.3)
  ),
  nrow = 2, rel_heights = c(0.04, 1)
)

out_dir <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Figures"
ggsave(file.path(out_dir, "Fig4_Yr25_Stress_QTL.pdf"),
       fig4, width = 7.09, height = 5.5, units = "in")
ggsave(file.path(out_dir, "Fig4_Yr25_Stress_QTL.png"),
       fig4, width = 7.09, height = 5.5, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", file.path(out_dir, "Fig4_Yr25_Stress_QTL.pdf")))
