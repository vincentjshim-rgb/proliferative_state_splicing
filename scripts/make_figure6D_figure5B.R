#!/usr/bin/env Rscript
# Figure 6D : held-out senescence (HO5, GSE306957)
# Figure 5B : local intron-retention burden (descriptive, n = 1 per condition)
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "public_data_tierA"
out  <- file.path(root, "derived", "manuscript_figures_v4")
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"
GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.4),
        plot.caption = element_text(colour = MUTED, size = 7.2, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8.2))

## ---------------- Figure 6D ----------------
hv <- file.path(root, "derived", "heldout2_senescence")
P <- read.delim(file.path(hv, "HO5_primary_tests.tsv"))
S <- read.delim(file.path(hv, "HO5_sensitivity_tests.tsv"))
D <- rbind(
  data.frame(lab = "Pooled control (primary)",            rho = P$rho[P$id == "HO5a"], adj = "unadjusted",  type = "primary"),
  data.frame(lab = "Pooled control (primary)",            rho = P$rho[P$id == "HO5b"], adj = "proliferation-adjusted", type = "primary"),
  data.frame(lab = "pLenti control",     rho = S$rho[S$id == "HO5s:pLenti"],               adj = "unadjusted", type = "sensitivity"),
  data.frame(lab = "pLenti control",     rho = S$rho[S$id == "HO5s:pLenti:prolifAdj"],     adj = "proliferation-adjusted", type = "sensitivity"),
  data.frame(lab = "siScramble control", rho = S$rho[S$id == "HO5s:siScramble"],           adj = "unadjusted", type = "sensitivity"),
  data.frame(lab = "siScramble control", rho = S$rho[S$id == "HO5s:siScramble:prolifAdj"], adj = "proliferation-adjusted", type = "sensitivity"),
  data.frame(lab = "IFIH1 KO (MDA5)",    rho = S$rho[S$id == "HO5s:IFIH1_KO"],             adj = "unadjusted", type = "sensitivity"),
  data.frame(lab = "IFIH1 KO (MDA5)",    rho = S$rho[S$id == "HO5s:IFIH1_KO:prolifAdj"],   adj = "proliferation-adjusted", type = "sensitivity"),
  data.frame(lab = "DDX58 KO (RIG-I)",   rho = S$rho[S$id == "HO5s:DDX58_KO"],             adj = "unadjusted", type = "sensitivity"),
  data.frame(lab = "DDX58 KO (RIG-I)",   rho = S$rho[S$id == "HO5s:DDX58_KO:prolifAdj"],   adj = "proliferation-adjusted", type = "sensitivity"))
D$lab <- factor(D$lab, levels = rev(unique(D$lab)))
D$adj <- factor(D$adj, levels = c("unadjusted", "proliferation-adjusted"))
disc <- data.frame(rho = c(-0.163, -0.134, -0.100, -0.044))
p6d <- ggplot(D, aes(rho, lab, colour = adj)) +
  annotate("rect", xmin = min(disc$rho), xmax = max(disc$rho), ymin = -Inf, ymax = Inf,
           fill = GREY, alpha = 0.14) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(aes(shape = type), size = 2.6,
             position = position_dodge(width = 0.55)) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 2.3, colour = INK,
            position = position_dodge(width = 0.55), hjust = 1.35) +
  scale_colour_manual(values = c(unadjusted = BLUE, `proliferation-adjusted` = GREEN)) +
  scale_shape_manual(values = c(primary = 16, sensitivity = 17), guide = "none") +
  scale_x_continuous(limits = c(-0.20, 0.02)) +
  labs(x = "Spearman rho (anchor vs senescent-minus-proliferating)", y = NULL,
       title = "D. Held-out senescence: confirmed, and not explained by proliferation",
       subtitle = "GSE306957, MRC5 lung fibroblast, 20 Gy X-ray senescence. Every test p = 1e-4.",
       caption = paste0("Grey band spans the four discovery senescence references (rho = -0.044 to -0.163);\n",
                        "the held-out estimate falls inside it, with no shrinkage. Circles = prespecified primary,\n",
                        "triangles = sensitivity arms. BH family is these two primary tests only.")) +
  base_theme
ggsave(file.path(out, "Figure6D_heldout_senescence.png"), p6d,
       width = 164, height = 92, units = "mm", dpi = 300)

## ---------------- Figure 5B ----------------
b <- read.delim(file.path(root, "derived", "intron_retention", "local_IR_burden.tsv"))
b$label <- factor(b$label, levels = c("UVA 15J, 0 h", "iPSC-CM", "HDF-CM", "Repro-CM"))
b$grp <- ifelse(b$sample == "UVA0", "0 h, no CM (time-confounded)", "24 h conditioned medium")
p5b <- ggplot(b, aes(100 * intronic_read_fraction, label, fill = grp)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * intronic_read_fraction)),
            hjust = -0.25, size = 2.7, colour = INK) +
  scale_fill_manual(values = c(`24 h conditioned medium` = BLUE,
                               `0 h, no CM (time-confounded)` = GREY)) +
  scale_x_continuous(limits = c(0, 7.4), expand = c(0, 0)) +
  labs(x = "Intronic read fraction (%)", y = NULL,
       title = "B. Local intron retention does not support improved splicing",
       subtitle = "6,749 genes covered in all four libraries; unambiguous intronic vs exonic blocks, GENCODE v41",
       caption = paste0("DESCRIPTIVE ONLY: n = 1 per condition, no p-values. Among the three 24 h conditions Repro-CM has the\n",
                        "HIGHEST intron retention, i.e. the opposite of a splicing-efficiency improvement. The 0 h sample differs\n",
                        "about two-fold but is confounded with elapsed time and library batch. Per-gene IR shift correlates with\n",
                        "the expression anchor (rho = -0.25), so the per-gene signal is not independent of expression level.")) +
  base_theme + theme(panel.grid.major.y = element_blank())
ggsave(file.path(out, "Figure5B_local_intron_retention.png"), p5b,
       width = 168, height = 84, units = "mm", dpi = 300)
cat("panels written\n")
