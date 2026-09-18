#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "public_data_tierA"
out  <- file.path(root, "derived", "manuscript_figures_v4")
R <- read.delim(file.path(root, "derived", "magnitude_vs_alignment", "magnitude_vs_alignment.tsv"))
## the anchor correlated with itself is 1 by definition and carries no information
R <- R[R$id != "LOCAL_ReproCM", ]
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"; PURPLE <- "#CC79A7"
YELLOW <- "#E69F00"; GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"
COL <- c("Repro-CM" = ORANGE, "photoprotection" = BLUE, "UV injury" = "#8FB8DE",
         "reprogramming" = GREEN, "metabolic / culture" = YELLOW, "donor age" = PURPLE,
         "senescence" = "#444444", "other" = GREY)
R$class <- factor(R$class, levels = names(COL))
ct <- suppressWarnings(cor.test(R$median_abs_logFC, R$rho_anchor, method = "spearman"))
lab <- data.frame(
  id = c("GSE326951_sauchinone", "GSE240226_succinate",
         "GSE109700_deep", "GSE149694_t2iLGoYD13", "GSE222414_UVBinjury",
         "GSE222414_UVB_UIFSP6"),
  short = c("sauchinone", "succinate rescue",
            "deep senescence", "naive reprogramming", "UVB injury", "TAT-UIFSP6"),
  lx = c(0.30, 0.62, 0.62, 0.86, 0.60, 0.22),
  ly = c(0.42, 0.15, -0.185, 0.25, 0.06, 0.02),
  hj = c(0, 0, 1, 1, 0, 1))
L <- merge(lab, R[, c("id", "median_abs_logFC", "rho_anchor")], by = "id")
p <- ggplot(R, aes(median_abs_logFC, rho_anchor)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK,
              fill = "grey85", linewidth = 0.45, linetype = 2, alpha = 0.5) +
  geom_segment(data = L, aes(x = median_abs_logFC, y = rho_anchor, xend = lx, yend = ly),
               colour = "grey70", linewidth = 0.22, inherit.aes = FALSE) +
  geom_point(aes(colour = class), size = 2.2) +
  geom_text(data = L, aes(x = lx, y = ly, label = short, hjust = hj),
            size = 2.4, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = COL, drop = FALSE) +
  scale_x_continuous(limits = c(0, 0.95)) +
  scale_y_continuous(limits = c(-0.22, 0.50)) +
  guides(colour = guide_legend(nrow = 2)) +
  annotate("text", x = 0.02, y = 0.47, hjust = 0, size = 2.9, colour = INK,
           label = sprintf("58 signatures:  rho = %+.3f,  p = %.2f", unname(ct$estimate), ct$p.value)) +
  labs(x = "perturbation magnitude   (median |log2 FC| of the contrast)",
       y = "alignment with the Repro-CM anchor   (Spearman rho)",
       title = "Alignment with the anchor is not a function of how hard the cell was hit",
       subtitle = "Every compendium signature except the anchor itself, plotted against the size of its own transcriptional response",
       caption = paste0("If the anchor simply tracked perturbation strength the cloud would rise to the right. It does not: the relationship is slightly negative\n",
                        "overall (rho = -0.170, p = 0.20) and within photoprotection alone (rho = -0.204, n = 14). The twelve largest perturbations on the map -\n",
                        "four reprogramming media, four senescence contrasts, a UVB injury - align at -0.157 to +0.181. Sauchinone therefore remains an\n",
                        "outlier in alignment that its size does not explain; with no untreated arm in that study, it cannot be resolved further here.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.3),
        plot.caption = element_text(colour = MUTED, size = 7, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8))
ggsave(file.path(out, "Figure8_magnitude_vs_alignment.png"), p,
       width = 172, height = 124, units = "mm", dpi = 300)
cat("Figure 8 written\n")

## ---- 8B : what DOES predict alignment ------------------------------------
M <- read.delim(file.path(root, "derived", "what_drives_alignment", "module_scores_vs_alignment.tsv"))
preds <- c(HSF1 = "HSF1 heat-shock response", premRNA = "capped pre-mRNA processing",
           splicing = "mRNA splicing", UPR_PERK = "UPR / PERK", prolif = "cell cycle",
           ECM = "ECM organization", median_abs_logFC = "perturbation magnitude")
B <- do.call(rbind, lapply(names(preds), function(v) {
  k <- is.finite(M[[v]])
  ct <- suppressWarnings(cor.test(M[[v]][k], M$rho_anchor[k], method = "spearman"))
  data.frame(pred = preds[[v]], rho = unname(ct$estimate), p = ct$p.value, n = sum(k)) }))
B <- B[order(B$rho), ]; B$pred <- factor(B$pred, levels = B$pred)
CARRIERS <- c("HSF1 heat-shock response", "UPR / PERK",
              "capped pre-mRNA processing", "mRNA splicing")
B$fill <- ifelse(B$pred == "perturbation magnitude", "ruled out",
          ifelse(B$pred %in% CARRIERS, "carries the decoupling (Fig. 4B)", "other"))
p8b <- ggplot(B, aes(rho, pred, fill = fill)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%+.3f   p = %s", rho,
                ifelse(p < 1e-4, sprintf("%.0e", p), sprintf("%.3f", p))),
                hjust = ifelse(rho >= 0, -0.06, 1.06)), size = 2.5, colour = INK) +
  scale_fill_manual(values = c(`carries the decoupling (Fig. 4B)` = BLUE, other = GREY,
                               `ruled out` = ORANGE)) +
  scale_x_continuous(limits = c(-0.62, 1.02)) +
  labs(x = "Spearman rho with alignment to the Repro-CM anchor", y = NULL,
       title = "What a perturbation has to do to resemble the Repro-CM response",
       subtitle = sprintf("Module content of %d independent perturbation signatures vs their alignment with the anchor", max(B$n)),
       caption = paste0("The two module families Figure 4B identified as carrying the decoupling also separate\n",
                        "perturbations that resemble the anchor from those that do not. Magnitude does not.\n",
                        "Controlling for magnitude leaves HSF1 at +0.577 and UPR/PERK at +0.364; within the\n",
                        "fifteen agent-alone arms UPR/PERK still predicts alignment (rho = +0.532, p = 0.044).\n",
                        "Part of this is expected - a perturbation sharing the anchor's dominant modules will\n",
                        "correlate with it - but it names which modules, across 59 independent datasets.")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.3),
        plot.caption = element_text(colour = MUTED, size = 7, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8))
ggsave(file.path(out, "Figure8B_what_drives_alignment.png"), p8b,
       width = 168, height = 122, units = "mm", dpi = 300)
cat("Figure 8B written\n")
