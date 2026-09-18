FIGURES_WRITTEN <- c("Fig4.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
PUBDIR <- "public_data_tierA/derived/figures_publication"
source("scripts/pub_theme.R"); D <- "public_data_tierA/derived"

## a  553 pathways: local NES vs the decoupling axis
pm <- read.delim(file.path(D, "direction_probe/P1_pathway_decoupling_ranked.tsv"), check.names = FALSE)
pm$name <- sub("__R-HSA.*$", "", pm$pathway)
rna <- "splic|snrnp|intron-containing|3'-end processing|polyadenylation|transport of mature"
trn <- "translation|ribosom|rrna|nonsense|peptide chain|40s|60s"
sts <- "heat shock|perk|hsf1|unfolded protein|nrf2|oxidative stress"
pm$cls <- "other"
pm$cls[grepl(trn, pm$name, ignore.case = TRUE)] <- "translation, rRNA, NMD"
pm$cls[grepl(sts, pm$name, ignore.case = TRUE)] <- "proteostatic stress"
pm$cls[grepl(rna, pm$name, ignore.case = TRUE)] <- "nuclear pre-mRNA processing"
pm$cls <- factor(pm$cls, levels = c("nuclear pre-mRNA processing", "proteostatic stress",
                                    "translation, rRNA, NMD", "other"))
ct <- suppressWarnings(cor.test(pm$local_NES, pm$decoupling_axis_NES, method = "spearman"))
LB <- data.frame(
  name = c("mRNA 3'-end processing", "Processing of Capped Intron-Containing Pre-mRNA",
           "mRNA Splicing - Major Pathway", "Transport of Mature mRNA derived from an Intron-Containing Transcript",
           "PERK regulates gene expression", "Regulation of HSF1-mediated heat shock response",
           "Eukaryotic Translation Elongation", "rRNA processing",
           "Extracellular matrix organization"),
  short = c("3'-end processing", "capped pre-mRNA", "mRNA splicing", "mRNA export",
            "PERK", "HSF1", "translation elongation", "rRNA processing", "ECM"),
  lx = c(-0.55, -0.55, -0.55, -0.55, 3.15, 3.15, 3.15, 3.15, -0.30),
  ly = c(3.78, 3.42, 3.06, 2.70, 3.62, 3.05, 0.42, 2.10, -3.90),
  hj = c(1, 1, 1, 1, 0, 0, 0, 0, 0))
L <- merge(LB, pm[, c("name", "local_NES", "decoupling_axis_NES")], by = "name")
p4a <- ggplot(pm, aes(local_NES, decoupling_axis_NES)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(data = subset(pm, cls == "other"), colour = "#DFE4E9", size = 0.5, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = INK,
              linewidth = 0.3, linetype = "22") +
  geom_segment(data = L, aes(x = local_NES, y = decoupling_axis_NES, xend = lx, yend = ly),
               colour = GREY_L, linewidth = 0.2, inherit.aes = FALSE) +
  geom_point(data = subset(pm, cls != "other"), aes(colour = cls), size = 0.95, stroke = 0) +
  geom_text(data = L, aes(x = lx, y = ly, label = short, hjust = hj),
            size = pt(5.2), family = FONT, colour = GREY_D, inherit.aes = FALSE) +
  annotate("text", x = -2.4, y = -2.6, hjust = 0, vjust = 1, size = pt(5.6), family = FONT,
           colour = INK, lineheight = 1.2,
           label = sprintf("rho = %.3f\nP = %.0e\n553 pathways", unname(ct$estimate), ct$p.value)) +
  annotate("point", x = rep(-2.4, 3), y = c(1.2, 0.75, 0.3), size = 0.95,
           colour = c(BLUE, GOLD, VERM)) +
  annotate("text", x = -2.25, y = c(1.2, 0.75, 0.3), hjust = 0, vjust = 0.42,
           size = pt(5.2), family = FONT, colour = GREY_D,
           label = c("nuclear pre-mRNA processing", "proteostatic stress", "translation, rRNA, NMD")) +
  scale_colour_manual(values = c("nuclear pre-mRNA processing" = BLUE,
                                 "proteostatic stress" = GOLD,
                                 "translation, rRNA, NMD" = VERM, other = "#DFE4E9"),
                      guide = "none") +
  scale_x_continuous(limits = c(-2.5, 4.6), breaks = seq(-2, 3, 1)) +
  scale_y_continuous(limits = c(-4.4, 4.0), breaks = seq(-4, 4, 2)) +
  labs(x = "local anchor NES",
       y = expression(paste("decoupling axis  ", NES[UVA] - NES[senescence]))) +
  theme_pub(7)

## b  module x axis
zm <- read.delim(file.path(D, "direction_probe/P2_module_axis_zmap.tsv"))
ML <- c(splicing = "mRNA splicing", premRNA = "capped pre-mRNA processing",
        mRNA_export = "mRNA export", three_prime = "3'-end processing",
        SF_regulon = "splicing-factor regulon", HSF1_heat_shock = "HSF1 heat shock",
        UPR_PERK = "UPR / PERK", ROS_detox = "ROS detoxification", autophagy = "autophagy",
        NMD = "NMD", translation = "translation", rRNA = "rRNA processing",
        proliferation = "cell cycle", ECM = "ECM organisation", IFN = "interferon")
AL <- c(local = "Repro-CM", UVA_meta = "UVA", senesc_meta = "senescence",
        MPTR_reprogramming = "reprogramming", fibroblast_age = "donor age",
        Tyshkovskiy_ITP_max_lifespan = "lifespan")
zm <- zm[zm$axis %in% names(AL) & zm$module %in% names(ML), ]
zm$m <- factor(ML[zm$module], levels = rev(ML)); zm$a <- factor(AL[zm$axis], levels = AL)
zm$z <- pmax(pmin(zm$matched_z, 8), -8)
p4b <- ggplot(zm, aes(a, m, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", matched_z), colour = abs(z) > 5),
            size = pt(4.8), family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = BLUE, mid = "#F5F7F8", high = VERM, midpoint = 0,
                       limits = c(-8, 8), breaks = c(-8, 0, 8),
                       name = "matched-null z",
                       guide = guide_colourbar(barwidth = unit(28, "pt"), barheight = unit(4, "pt"),
                                               title.position = "left", ticks = FALSE)) +
  labs(x = NULL, y = NULL) +
  theme_pub(7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 6.2),
        legend.position.inside = c(0.52, 1.045), legend.position = "inside",
        legend.direction = "horizontal",
        legend.title = element_text(size = 6, colour = GREY_D, vjust = 1),
        plot.margin = margin(14, 3, 2, 2))

save_fig(compose(p4a, p4b, labels = c("a", "b"), ncol = 2, rel_widths = c(1.15, 1)),
         "Fig4.png", 183, 96)
cat("Fig4 saved\n")
