## Fig1_Phenotypic_Context.R
## Three-panel figure:
##   A — Trait correlation heatmap faceted by year (yr.23, yr.24, yr.25)
##   B — Frost day count and hard-freeze events 2023–2025
##   C — DTFlower phenotypic SD per year (variance collapse)

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(reshape2)
})

pub_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black", face = "bold"),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

yr_colors  <- c("2023" = "#9ECAE1", "2024" = "#3182BD", "2025" = "#08519C")
trait_labels <- c(DTFlower = "Days to\nFlowering",
                  DTFruit  = "Days to\nRipe Fruit",
                  Flow2Fruit = "Fruiting\nPeriod")

# ── Load phenotype data ───────────────────────────────────────────────────────
pheno <- fread("/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/Phenotypes/SHB_BLUE_all_pheno.csv")

traits <- c("DTFlower", "DTFruit", "Flow2Fruit")
years  <- c("23", "24", "25")

# Build 3 correlation matrices (one per year)
cor_list <- lapply(years, function(yr) {
  cols <- paste0(traits, "_yr.", yr)
  mat  <- pheno[, ..cols]
  setnames(mat, cols, traits)
  mat  <- mat[complete.cases(mat)]
  r    <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dt   <- as.data.table(melt(r))
  dt[, Year := paste0("20", yr)]
  dt
})

cor_dt <- rbindlist(cor_list)
setnames(cor_dt, c("Var1", "Var2", "value"), c("Trait1", "Trait2", "r"))
cor_dt[, Year := factor(Year, levels = c("2023", "2024", "2025"))]

# Keep only lower triangle (including diagonal) to avoid redundancy
trait_order <- c("DTFlower", "DTFruit", "Flow2Fruit")
cor_dt[, Trait1 := factor(Trait1, levels = trait_order)]
cor_dt[, Trait2 := factor(Trait2, levels = rev(trait_order))]

# Mask upper triangle
cor_dt[, t1_idx := as.integer(Trait1)]
cor_dt[, t2_idx := as.integer(Trait2)]
cor_dt <- cor_dt[t1_idx >= (4 - t2_idx)]   # lower-left + diagonal

# Correlation labels (hide diagonal 1.00)
cor_dt[, label := ifelse(Trait1 == as.character(Trait2), "",
                          sprintf("%.2f", r))]

# Short axis labels
short_labels <- c(DTFlower = "DTFlw", DTFruit = "DTFrt",
                  Flow2Fruit = "F2Fr")

pA <- ggplot(cor_dt, aes(x = Trait1, y = Trait2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.8, color = "black") +
  facet_wrap(~Year, nrow = 1) +
  scale_fill_gradient2(
    low      = "#4292C6",
    mid      = "white",
    high     = "#08306B",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Pearson r"
  ) +
  scale_x_discrete(labels = short_labels) +
  scale_y_discrete(labels = rev(short_labels)) +
  labs(x = NULL, y = NULL,
       title = "A  Trait correlations by year") +
  pub_theme +
  theme(
    axis.text.x      = element_text(size = 8, angle = 30, hjust = 1),
    axis.text.y      = element_text(size = 8),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.border     = element_rect(fill = NA, color = "grey70"),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width  = unit(0.3, "cm"),
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7)
  )

# ── Panel B: Frost day count + hard-freeze events ────────────────────────────
# Summary from Methods 2.2 (computed from GA weather station daily data)
frost_dt <- data.table(
  Year  = factor(rep(c("2023", "2024", "2025"), 2)),
  Count = c(6, 11, 21,   # frost days (min ≤ 32°F)
            2,  4, 12),  # hard freeze events (min ≤ 28°F)
  Type  = factor(rep(c("Frost days\n(min <=32F)",
                        "Hard freeze events\n(min <=28F)"), each = 3),
                 levels = c("Frost days\n(min <=32F)",
                             "Hard freeze events\n(min <=28F)"))
)

pB <- ggplot(frost_dt, aes(x = Year, y = Count, fill = Year, alpha = Type)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "grey25", linewidth = 0.35) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.7),
            vjust = -0.35, size = 3.5, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = yr_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Frost days\n(min <=32F)"       = 1.0,
               "Hard freeze events\n(min <=28F)" = 0.45),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Year", y = "Event count (Jan-Apr)",
       title = "B  Winter frost and hard-freeze events") +
  pub_theme +
  theme(
    legend.position   = c(0.28, 0.84),
    legend.background = element_blank(),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 7.5, lineheight = 1.2)
  )

# ── Panel C: DTFlower phenotypic SD per year ──────────────────────────────────
sd_dt <- data.table(
  Year = factor(c("2023", "2024", "2025")),
  SD   = c(13.3, 11.2, 6.9)
)

pC <- ggplot(sd_dt, aes(x = Year, y = SD, fill = Year)) +
  geom_col(color = "grey25", linewidth = 0.35, width = 0.55) +
  geom_text(aes(label = sprintf("%.1f days", SD)),
            vjust = -0.35, size = 3.5, fontface = "bold") +
  # annotate the 48% collapse
  annotate("segment",
           x = 1, xend = 3, y = 16.5, yend = 16.5,
           color = "#2171B5", linewidth = 0.6,
           arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed")) +
  annotate("text",
           x = 2, y = 17.6,
           label = "-48% SD collapse (frost year)",
           size = 3.0, color = "#2171B5", fontface = "italic") +
  scale_fill_manual(values = yr_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     limits = c(0, 19)) +
  labs(x = "Year", y = "Phenotypic SD (days)",
       title = "C  DTFlower phenotypic variance by year") +
  pub_theme

# ── Assemble ──────────────────────────────────────────────────────────────────
fig1 <- plot_grid(
  plot_grid(pA, nrow = 1, rel_widths = 1),
  plot_grid(pB, pC, ncol = 2, rel_widths = c(1, 1)),
  nrow = 2,
  rel_heights = c(1.15, 1.0),
  labels = NULL
)

out_dir <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS/PUBLICATION/Figures"
ggsave(file.path(out_dir, "Fig1_Phenotypic_Context.pdf"),
       fig1, width = 7.09, height = 4.5, units = "in")
ggsave(file.path(out_dir, "Fig1_Phenotypic_Context.png"),
       fig1, width = 7.09, height = 4.5, units = "in", dpi = 300)
cat(sprintf("Saved: %s\n", file.path(out_dir, "Fig1_Phenotypic_Context.pdf")))
