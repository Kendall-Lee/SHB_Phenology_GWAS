#!/usr/bin/env Rscript
# Fig4_StableQTL_Heatmap.R
# Figure 4 — Linear GWAS: Stable QTL heatmap
#   Rows:    22 stable loci (peak marker per 1-Mb locus, BLUE_exc.25 analysis)
#            + Chr.03:38.4 Mb stress locus (separator line, yr.25 amplification)
#   Columns: [DTFlower | DTFruit | Flow2Fruit | FruitWeight] × [BLUE_exc.25 | yr.23 | yr.24 | yr.25]
#   Fill:    -log10(p) from GWASpoly scores; cells below Meff threshold (5.52) shown in grey
# Output:   Linear_MS/Fig4_StableQTL_Heatmap.pdf + .png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

BASE    <- "/Users/kendalllee/Documents/Blueberry/FINAL_GWAS"
OUT_DIR <- file.path(BASE, "PUBLICATION/Linear_MS")
THRESH  <- 5.52

# ── Load data ─────────────────────────────────────────────────────────────────

# 88-marker table with all 16 trait×year scores (pre-computed from GWASpoly scans)
raw <- fread(file.path(BASE,
  "All_Markers/LINEAR/BLUE_exc.25/TopMarkers_LINEAR_BLUE_exc25.csv"))

# Fix typo in column name
setnames(raw, "Flow2Fruit_BLUES", "Flow2Fruit_BLUE", skip_absent = TRUE)

# Locus labels from the publication SupplTable
supp <- fread(file.path(BASE,
  "PUBLICATION/Figures/SupplTable_S1_TopMarkers_LINEAR_BLUEexc25.csv"))
setnames(supp, c("Marker", "Locus", "Primary_trait"), c("Marker", "Locus", "Primary_trait"),
         skip_absent = TRUE)
# Standardise column references
locus_map <- supp[, .(Marker, Locus, Primary_trait = `Primary trait`)]

# Join locus labels onto raw data
dat <- merge(raw, locus_map, by = "Marker", all.x = TRUE)

# ── Deduplicate to one peak marker per locus ──────────────────────────────────
dat[, Score_exc25 := as.numeric(Score_exc25)]
setorder(dat, Locus, -Score_exc25)
peak <- dat[!is.na(Locus), .SD[1], by = Locus]

# ── Score columns → long format ───────────────────────────────────────────────
score_cols <- c(
  "DTFlower_BLUE",    "DTFlower_yr.23",    "DTFlower_yr.24",    "DTFlower_yr.25",
  "DTFruit_BLUE",     "DTFruit_yr.23",     "DTFruit_yr.24",     "DTFruit_yr.25",
  "Flow2Fruit_BLUE",  "Flow2Fruit_yr.23",  "Flow2Fruit_yr.24",  "Flow2Fruit_yr.25"
)

long_stable <- melt(
  peak[, c("Locus", "Marker", "Chrom", "Primary_trait", score_cols), with = FALSE],
  id.vars      = c("Locus", "Marker", "Chrom", "Primary_trait"),
  measure.vars = score_cols,
  variable.name = "Col_key",
  value.name    = "Score"
)
long_stable[, Score := as.numeric(Score)]

# ── Add stress locus row (Chr.03:38.4 Mb, yr.25-specific) ────────────────────
# Scores extracted from ALL_TRAITS GWASpoly scan (DTFruit_BLUE = all-year BLUE)
stress_scores <- c(
  DTFlower_BLUE    = 2.80,  DTFlower_yr.23  = 2.68,  DTFlower_yr.24  = 2.55,  DTFlower_yr.25  = 1.08,
  DTFruit_BLUE     = 10.66, DTFruit_yr.23   = 1.79,  DTFruit_yr.24   = 1.86,  DTFruit_yr.25   = 22.69,
  Flow2Fruit_BLUE  = 4.80,  Flow2Fruit_yr.23 = 0.55, Flow2Fruit_yr.24 = 1.33, Flow2Fruit_yr.25 = 7.74
)
long_stress <- data.table(
  Locus        = "Chr.03: 38.4 Mb*",
  Marker       = "Chr.03_38388088",
  Chrom        = "Chr.03",
  Primary_trait = "DTFruit",
  Col_key      = names(stress_scores),
  Score        = as.numeric(stress_scores)
)

long_all <- rbindlist(list(long_stable, long_stress), fill = TRUE)

# ── Parse column key into Trait and Year ──────────────────────────────────────
long_all[, Trait := sub("_(BLUE|yr\\..*)$", "", Col_key)]
long_all[, Year  := fcase(
  grepl("_BLUE$",   Col_key), "BLUE\nexc.25",
  grepl("yr\\.23$", Col_key), "yr.23",
  grepl("yr\\.24$", Col_key), "yr.24",
  grepl("yr\\.25$", Col_key), "yr.25"
)]

# Remove the all-year BLUE cell for the stress locus (different analysis; excluded)
long_all <- long_all[!(Locus == "Chr.03: 38.4 Mb*" & Year == "BLUE\nexc.25")]

# ── Row ordering ──────────────────────────────────────────────────────────────
# Primary trait priority, then descending Score_exc25
trait_order <- c("DTFlower", "DTFruit", "Flow2Fruit")
peak[, sort_key := match(Primary_trait, trait_order)]
setorder(peak, sort_key, -Score_exc25)
locus_order <- c(peak$Locus, "Chr.03: 38.4 Mb*")

# Short display labels for y-axis
locus_labels <- setNames(locus_order, locus_order)
# Add primary trait tag
peak_tag <- peak[, .(Locus, Primary_trait)]
for (i in seq_len(nrow(peak_tag))) {
  lbl <- paste0(peak_tag$Locus[i], "  [", peak_tag$Primary_trait[i], "]")
  locus_labels[peak_tag$Locus[i]] <- lbl
}
locus_labels["Chr.03: 38.4 Mb*"] <- "Chr.03: 38.4 Mb*  [DTFruit]"

long_all[, Locus := factor(Locus, levels = rev(locus_order))]

# ── Column ordering ───────────────────────────────────────────────────────────
year_order  <- c("BLUE\nexc.25", "yr.23", "yr.24", "yr.25")
trait_order2 <- c("DTFlower", "DTFruit", "Flow2Fruit")
long_all[, Trait := factor(Trait, levels = trait_order2)]
long_all[, Year  := factor(Year,  levels = year_order)]

# ── Significance flag & display value ────────────────────────────────────────
long_all[, Sig   := !is.na(Score) & Score >= THRESH]
long_all[, Score_disp := pmin(Score, 23)]  # cap for colour scale

# ── Colour palette ────────────────────────────────────────────────────────────
# Below threshold: light grey; above threshold: blue gradient
# Two-segment fill: grey below threshold, blue scale above

sig_pal <- colorRampPalette(c("#C6DBEF", "#2171B5", "#08306B"))(100)

# ── Plot ──────────────────────────────────────────────────────────────────────
pub_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(color = "black", face = "bold"),
    strip.text      = element_text(face = "bold", size = 9.5),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.spacing.x  = unit(2, "pt"),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8)
  )

# Separate stable and stress rows with a thin visual gap
# Achieved by leaving a gap row — handled via geom_hline
stress_y_pos <- which(rev(locus_order) == "Chr.03: 38.4 Mb*") - 0.5

p <- ggplot(long_all, aes(x = Year, y = Locus, fill = ifelse(Sig, Score_disp, NA_real_))) +
  # Grey tiles for sub-threshold cells
  geom_tile(data = long_all[Sig == FALSE],
            aes(x = Year, y = Locus), fill = "grey88",
            color = "white", linewidth = 0.4, width = 0.92, height = 0.88) +
  # Coloured tiles for significant cells
  geom_tile(data = long_all[Sig == TRUE],
            color = "white", linewidth = 0.4, width = 0.92, height = 0.88) +
  # Score text on significant cells
  geom_text(data = long_all[Sig == TRUE],
            aes(label = ifelse(Score > 22, sprintf("%.1f*", Score),
                               sprintf("%.1f", Score))),
            size = 2.55, color = "white", fontface = "bold") +
  # Separator line above stress locus
  geom_hline(yintercept = stress_y_pos, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  # Facet by trait
  facet_grid(. ~ Trait, scales = "free_x", space = "free_x",
             labeller = labeller(Trait = c(
               DTFlower   = "Days to Flowering",
               DTFruit    = "Days to Ripe Fruit",
               Flow2Fruit = "Fruiting Period"))) +
  scale_fill_gradientn(
    colours   = sig_pal,
    limits    = c(THRESH, 23),
    na.value  = "grey88",
    name      = expression(bold(-log[10](italic(p)))),
    breaks    = c(5.52, 8, 11, 15, 20, 23),
    labels    = c("5.5\n(threshold)", "8", "11", "15", "20", ">=23")
  ) +
  scale_y_discrete(labels = locus_labels) +
  labs(
    title    = "Stable GWAS loci across traits and years - Linear GWAS (Suziblue hap1)",
    subtitle = paste0("Fill = -log10(p) from GWASpoly LOCO scan; grey = below Meff threshold (", THRESH, ").\n",
                      "BLUE exc.25 = yr.23+yr.24 combined BLUE (primary analysis). ",
                      "# Chr.03:38.4 Mb is yr.25 stress-specific (BLUE exc.25 cell excluded). ",
                      "* Score capped at 23 for display (actual = 22.69)."),
    x = NULL, y = NULL
  ) +
  pub_theme +
  theme(
    axis.text.y      = element_text(size = 7.5, family = "mono"),
    axis.text.x      = element_text(size = 7.5),
    legend.position  = "right",
    legend.key.height = unit(1.8, "cm")
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(file.path(OUT_DIR, "Fig4_StableQTL_Heatmap.pdf"), p,
       width = 7.09, height = 5.5, units = "in", device = "pdf")
ggsave(file.path(OUT_DIR, "Fig4_StableQTL_Heatmap.png"), p,
       width = 7.09, height = 5.5, units = "in", dpi = 300, device = "png")

cat(sprintf("Saved Fig4_StableQTL_Heatmap to %s\n", OUT_DIR))
