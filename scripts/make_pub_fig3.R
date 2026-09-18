FIGURES_WRITTEN <- c("Fig3.png")
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

p4 <- read.delim(file.path(D, "direction_probe/P4_decoupling_proliferation_sensitivity.tsv"))
p4$lab <- c("no adjustment", "all axes residualised", "proliferation-neutral genes",
            "proliferation-linked genes")
p4$lab <- factor(p4$lab, levels = rev(p4$lab))
p3a <- ggplot(p4, aes(decoupling_index_D, lab)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = decoupling_index_D, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(size = 1.5, colour = BLUE) +
  geom_text(aes(label = sprintf("%.3f", decoupling_index_D)), hjust = -0.3,
            size = pt(5.4), family = FONT, colour = INK) +
  geom_text(aes(x = 0.012, label = sprintf("n = %s", format(n_genes, big.mark = ","))),
            hjust = 0, vjust = 2.1, size = pt(5), family = FONT, colour = GREY_D) +
  scale_x_continuous(limits = c(0, 0.70), breaks = seq(0, 0.6, 0.2)) +
  labs(x = "decoupling index D", y = NULL) +
  theme_pub(7) + theme_cat_y()

bt <- read.delim(file.path(D, "manuscript_figures_v4/Figure3C_D_by_proliferation_loading_decile.tsv"))
p3b <- ggplot(bt, aes(mean_loading, D)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.25) +
  geom_hline(yintercept = 0.460, colour = VERM, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = GREY, linewidth = 0.35) +
  geom_point(size = 1.3, colour = BLUE) +
  annotate("text", x = max(bt$mean_loading), y = 0.460, label = "primary D", colour = VERM,
           size = pt(5.2), family = FONT, hjust = 1, vjust = -0.6) +
  scale_y_continuous(limits = c(0, 0.86), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "per-gene proliferation loading (decile mean)", y = "decoupling index D") +
  theme_pub(7)

save_fig(compose(p3a, p3b, labels = c("a", "b"), ncol = 2, rel_widths = c(1.12, 1)),
         "Fig3.png", 183, 58)
cat("Fig3 saved\n")
