FIGURES_WRITTEN <- c("Fig8.png")
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
M <- read.delim(file.path(D, "magnitude_vs_alignment/magnitude_vs_alignment.tsv"))
M <- M[M$id != "LOCAL_ReproCM", ]
COL <- c(photoprotection = BLUE, "UV injury" = "#9CC2DE", reprogramming = GREEN,
         "metabolic / culture" = GOLD, "donor age" = PURP, senescence = INK, other = GREY)
M$class <- factor(M$class, levels = names(COL))
ct <- suppressWarnings(cor.test(M$median_abs_logFC, M$rho_anchor, method = "spearman"))
LB <- data.frame(id = c("GSE326951_sauchinone", "GSE109700_deep", "GSE149694_t2iLGoYD13",
                        "GSE179848_Oligomycin"),
                 s = c("sauchinone", "deep senescence", "naive reprogramming", "oligomycin"),
                 lx = c(0.30, 0.62, 0.86, 0.30), ly = c(0.41, -0.185, 0.245, 0.245),
                 hj = c(0, 1, 1, 1))
L <- merge(LB, M[, c("id", "median_abs_logFC", "rho_anchor")], by = "id")
p8a <- ggplot(M, aes(median_abs_logFC, rho_anchor)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK, fill = "#EEF1F4",
              linewidth = 0.3, linetype = "22", alpha = 0.8) +
  geom_segment(data = L, aes(x = median_abs_logFC, y = rho_anchor, xend = lx, yend = ly),
               colour = GREY_L, linewidth = 0.2, inherit.aes = FALSE) +
  geom_point(aes(colour = class), size = 1.15, stroke = 0) +
  geom_text(data = L, aes(x = lx, y = ly, label = s, hjust = hj), size = pt(5),
            family = FONT, colour = INK, inherit.aes = FALSE) +
  annotate("text", x = 0.02, y = 0.46, hjust = 0, size = pt(5.6), family = FONT, colour = INK,
           label = sprintf("rho = %+.3f,  P = %.2f", unname(ct$estimate), ct$p.value)) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(0, 0.95), breaks = seq(0, 0.9, 0.3)) +
  scale_y_continuous(limits = c(-0.22, 0.50), breaks = seq(-0.2, 0.4, 0.2)) +
  labs(x = expression(paste("perturbation magnitude   median |", log[2], " FC|")),
       y = expression(paste("alignment with the anchor   ", rho))) +
  theme_pub(7)

MS <- read.delim(file.path(D, "what_drives_alignment/module_scores_vs_alignment.tsv"))
pr <- c(HSF1 = "HSF1 heat shock", premRNA = "capped pre-mRNA processing",
        splicing = "mRNA splicing", UPR_PERK = "UPR / PERK", prolif = "cell cycle",
        ECM = "ECM organisation", median_abs_logFC = "perturbation magnitude")
B <- do.call(rbind, lapply(names(pr), function(v) {
  k <- is.finite(MS[[v]]); c2 <- suppressWarnings(cor.test(MS[[v]][k], MS$rho_anchor[k], method = "spearman"))
  data.frame(p = pr[[v]], rho = unname(c2$estimate), pv = c2$p.value) }))
B <- B[order(B$rho), ]; B$p <- factor(B$p, levels = B$p)
CARR <- c("HSF1 heat shock", "UPR / PERK", "capped pre-mRNA processing", "mRNA splicing")
B$f <- ifelse(B$p == "perturbation magnitude", "ruled out",
       ifelse(B$p %in% CARR, "carries the decoupling", "other"))
p8b <- ggplot(B, aes(rho, p, fill = f)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_col(width = 0.58) +
  geom_text(aes(label = sprintf("%+.3f  P = %s", rho,
                ifelse(pv < 1e-4, sprintf("%.0e", pv), sprintf("%.3f", pv))),
                hjust = ifelse(rho >= 0, -0.07, 1.07)), size = pt(5.2), family = FONT, colour = INK) +
  scale_fill_manual(values = c(`carries the decoupling` = BLUE, other = "#C3CAD2",
                               `ruled out` = VERM), guide = "none") +
  scale_x_continuous(limits = c(-0.62, 1.04), breaks = seq(-0.5, 1, 0.5)) +
  labs(x = expression(paste("Spearman ", rho, " with alignment to the anchor")), y = NULL) +
  theme_pub(7) + theme_cat_y()
save_fig(compose(p8a, p8b, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.08)),
         "Fig8.png", 183, 68)
cat("Fig8 saved\n")
