## Fig6_Yr25_Stress_QTL.R  (was Fig6_Chr03_StressLocus.R)
## Figure 6 -- Year 2025 stress QTL: 3-panel figure
## Panel D (Chr.03:38.4 Mb allele effect) REMOVED 2026-05-29:
##   n=5 in rare dosage class, marker non-unique in genome (6 hits, off-target),
##   absent from stable-year scan (RedoNo25 DTFruit_BLUE score=1.17), not defensible.
##
## Output: Linear_MS/Fig4_Yr25_Stress_QTL.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS")

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

fig_title <- function(txt, size = 10.5)
  ggdraw() + draw_label(txt, fontface = "bold", size = size, x = 0.01, hjust = 0)

# ── Panel A: QTL counts --------------------------------------------------------
year_counts <- data.table(
  Dataset = factor(
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

p_a <- ggplot(year_counts, aes(x = Dataset, y = N_QTL, fill = Year_type)) +
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

# ── Panel B: January minimum temperature ---------------------------------------
temp_data <- data.table(
  Year     = factor(c("2023","2024","2025")),
  Min_Temp = c(3.8, -1.9, -5.3)
)

p_b <- ggplot(temp_data, aes(x = Year, y = Min_Temp, fill = Min_Temp)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.55) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = paste0(Min_Temp, "C"),
                vjust = ifelse(Min_Temp >= 0, -0.3, 1.3)),
            size = 3.3, fontface = "bold") +
  scale_fill_gradient2(low = "#0D3B6E", mid = "grey92", high = "#D84315",
                       midpoint = 0, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.15))) +
  labs(x = "Year", y = "Min. Jan. Temp (C)", title = "B  January minimum temperature") +
  pub_theme

# ── Panel C: PVE bars ---------------------------------------------------------
stress_qtl <- data.table(
  Label = c("DTFruit (2025)\nChr.03:38.4 Mb",
            "DTFlower (2025)\nChr.11:45.5 Mb",
            "DTFruit (2025)\nChr.07:25.0 Mb"),
  PVE   = c(13.26, 6.11, 4.58),
  Trait = c("DTFruit","DTFlower","DTFruit"),
  plab  = c("p = 2.3e-18*", "p = 5.8e-9", "p = 5.4e-7"),
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
            hjust = -0.05, size = 2.9) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Stress-specific (yr.2025 only)" = 1.0,
               "Phenology locus in yr.2025"     = 0.38),
    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.65))) +
  labs(x = "PVE (%) in Year 2025", y = NULL,
       title = "C  Exploratory stress-responsive QTL and stable locus performance in Year 2025",
       subtitle = paste0(
         "Faded bars: stable phenology loci attenuated in 2025 frost year.\n",
         "* Chr.03:38.4 Mb: Flow2Fruit only; exploratory (n=5 rare allele, requires replication).")) +
  pub_theme +
  theme(axis.text.y = element_text(size = 8, lineheight = 1.2),
        legend.position = "bottom", legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 7.5, color = "grey40", lineheight = 1.3))

# ── Assemble ------------------------------------------------------------------
fig6 <- plot_grid(
  fig_title("Figure 6. Year 2025 Frost-Stress QTL Architecture"),
  plot_grid(
    plot_grid(p_a, p_b, ncol = 2, rel_widths = c(1.7, 1)),
    p_c,
    nrow = 2, rel_heights = c(1.0, 1.3)
  ),
  nrow = 2, rel_heights = c(0.04, 1)
)

ggsave(file.path(OUT_DIR, "Fig4_Yr25_Stress_QTL.pdf"), fig6, width = 7.09, height = 5.5, units = "in")
ggsave(file.path(OUT_DIR, "Fig4_Yr25_Stress_QTL.png"), fig6, width = 7.09, height = 5.5, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", file.path(OUT_DIR, "Fig4_Yr25_Stress_QTL.pdf")))
