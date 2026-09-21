## DRAFT C  Every gene at once: the age association genome-wide, and each gene's age effect
## against its coupling to the counted division rate.
##
## The forms of the ageing-transcriptome literature (Peters 2015 Manhattan/volcano; the
## effect-size scatters of replication studies):
##   a  volcano of the per-decade age effect across all measured genes in the 107 adult donors
##      (limma, cohort covariates), the 177 splicing genes and the cell-cycle union highlighted
##   b  gene by gene: the age effect (donor cohort) against the correlation with the counted
##      division rate (Cellular Lifespan Study) -- two datasets, one point per gene
## Output: figures_sciadv/draft_C_genomewide.png (not wired into the manuscript)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"
CORE  <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")
DE <- fread(file.path(D, "revision_stats/GSE113957_primary_cohort_per_decade_DE.tsv"))
GV <- fread(file.path(D, "gene_level/gene_vs_counted_rate.tsv"))
## the cell-cycle union used throughout (Reactome Cell Cycle, Mitotic + DNA Replication + M Phase)
GMTD <- tempfile("reactome"); dir.create(GMTD)
unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = GMTD)
GMT <- strsplit(readLines(list.files(GMTD, pattern = "\\.gmt$", full.names = TRUE)[1]), "\t")
CCG <- unique(unlist(lapply(GMT[grepl("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", vapply(GMT, `[`, "", 1))],
                            function(v) v[-(1:2)])))
DE[, gene := toupper(gene)]
DE[, grp := fifelse(gene %in% CORE, "177 splicing genes", fifelse(gene %in% CCG, "cell-cycle union", "other"))]
DE[, grp := factor(grp, levels = c("other", "cell-cycle union", "177 splicing genes"))]
DE[, mlp := -log10(PValue)]
cat(sprintf("a  %d genes; FDR < 0.05: %d (%d down, %d up); splicing set %d, cell-cycle union %d\n",
            nrow(DE), sum(DE$FDR < 0.05), sum(DE$FDR < 0.05 & DE$logFC < 0), sum(DE$FDR < 0.05 & DE$logFC > 0),
            sum(DE$grp == "177 splicing genes"), sum(DE$grp == "cell-cycle union")))
GCOL <- c(other = "#D5DAE0", `cell-cycle union` = RED, `177 splicing genes` = BLUE)
## four labels, placed by hand so that none overlaps (ggrepel is not installed)
SHOW <- c("SRSF1", "SRSF2", "SRSF3", "HNRNPK")
LABS <- DE[gene %in% SHOW]; LABS[, nx := 0.012]
LABS[, ny := c(SRSF1 = 0.25, SRSF2 = -0.25, SRSF3 = 0.0, HNRNPK = 0.0)[gene]]
pA <- ggplot(DE[order(grp)], aes(logFC, mlp, colour = grp)) +
  geom_point(aes(size = grp, alpha = grp)) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_hline(yintercept = -log10(max(DE$PValue[DE$FDR < 0.05])), colour = GREY, linewidth = 0.3, linetype = "22") +
  geom_segment(data = LABS, aes(x = logFC, xend = logFC + nx - 0.003, y = mlp, yend = mlp + ny), colour = GREY, linewidth = 0.25) +
  geom_text(data = LABS, aes(x = logFC + nx, y = mlp + ny, label = gene), family = FONT, size = pt(6.8), fontface = "italic",
            colour = INK, hjust = 0, vjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = GCOL, name = NULL) +
  scale_size_manual(values = c(other = 0.45, `cell-cycle union` = 0.9, `177 splicing genes` = 1.1), guide = "none") +
  scale_alpha_manual(values = c(other = 0.5, `cell-cycle union` = 0.85, `177 splicing genes` = 0.95), guide = "none") +
  annotate("text", x = max(DE$logFC) * 0.98, y = max(DE$mlp) * 0.98, hjust = 1, vjust = 1, family = FONT,
           size = pt(7), colour = INK2, lineheight = 1.05,
           label = sprintf("%s genes, 107 adult donors\ndashed: FDR = 0.05", format(nrow(DE), big.mark = ","))) +
  scale_x_continuous("age effect (log₂ CPM per decade)", labels = num_axis()) +
  scale_y_continuous(expression(-log[10]~italic(P))) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), nrow = 1)) +
  theme_sa(8) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(b = -4),
                      legend.key.size = unit(7, "pt"), legend.text = element_text(size = 7))

## ---- b. each gene's age effect against its coupling to the counted rate ------------------
J <- merge(DE[, .(gene, logFC, FDR, grp)], GV[, .(gene, rho)], by = "gene")
ct <- cor.test(J$rho, J$logFC, method = "spearman", exact = FALSE)
ctS <- cor.test(J$rho[J$grp == "177 splicing genes"], J$logFC[J$grp == "177 splicing genes"], method = "spearman", exact = FALSE)
cat(sprintf("b  %d genes in both datasets: rho(coupling, age effect) = %.3f (P = %.1e); within the 177: %.3f (P = %.2g)\n",
            nrow(J), ct$estimate, ct$p.value, ctS$estimate, ctS$p.value))
LB2 <- J[gene %in% SHOW]; LB2[, ny := c(SRSF1 = 0.055, SRSF2 = 0.085, SRSF3 = -0.075, HNRNPK = -0.105)[gene]]
pB <- ggplot(J[order(grp)], aes(rho, logFC, colour = grp)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) + geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_point(aes(size = grp, alpha = grp)) +
  geom_smooth(data = J, aes(rho, logFC), inherit.aes = FALSE, method = "lm", formula = y ~ x, colour = INK,
              linewidth = 0.5, se = FALSE, linetype = "22") +
  geom_segment(data = LB2, aes(x = rho, xend = rho, y = logFC, yend = logFC + ny * 0.85), colour = GREY, linewidth = 0.25) +
  geom_text(data = LB2, aes(y = logFC + ny, label = gene), family = FONT, size = pt(6.8), fontface = "italic", colour = INK,
            hjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_size_manual(values = c(other = 0.45, `cell-cycle union` = 0.9, `177 splicing genes` = 1.1), guide = "none") +
  scale_alpha_manual(values = c(other = 0.5, `cell-cycle union` = 0.85, `177 splicing genes` = 0.95), guide = "none") +
  annotate("text", x = min(J$rho), y = max(J$logFC), hjust = 0, vjust = 1, family = FONT, size = pt(7), colour = INK,
           lineheight = 1.05, label = sprintf("all genes: %s = %s\n%s\nn = %s genes in both datasets", RHO, num(ct$estimate),
                                              pfmt(ct$p.value), format(nrow(J), big.mark = ","))) +
  scale_x_continuous(sprintf("%s with the counted division rate (328 libraries)", RHO), labels = num_axis()) +
  scale_y_continuous("age effect in the donor cohort (log₂ CPM per decade)", labels = num_axis()) +
  theme_sa(8) + theme(plot.margin = margin(14, 6, 3, 3))
save_fig(lab_grid(pA, pB, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.05)), "draft_C_genomewide.png", 183, 92)
