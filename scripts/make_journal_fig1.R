FIGURES_WRITTEN <- c("Fig1.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/journal_theme.R")
D <- "public_data_tierA/derived"; CS <- file.path(D, "classical_stats")
AR <- arrow(length = unit(1.7, "pt"), type = "closed")

## ============================ A. study design schematic =====================
## coordinate space 0-100 (x) by 0-100 (y); panel is tall and narrow
bx <- function(x, y, w, h, lab, fill = "white", col = INK, fc = INK, tsz = 6.1, lw = 0.32)
  data.frame(x, y, w, h, lab, fill, col, fc, tsz, lw, stringsAsFactors = FALSE)
R <- rbind(
  bx(27, 96,  40, 7.0, "Adult primary HDF"),
  bx(76, 96,  34, 7.0, "UVA 15 J cm⁻²", "#F4EFF9", UV, UV),
  bx(17, 79,  30, 7.0, "HDF-CM"),
  bx(50, 79,  30, 7.0, "Repro-CM", "#E4F1F0", REPRO, REPRO, 6.4),
  bx(83, 79,  30, 7.0, "iPSC-CM"),
  bx(50, 62,  96, 7.4, "RNA-seq of the recipient fibroblast   ·   n = 1 per condition"),
  bx(50, 46,  96, 8.4, ""),
  bx(21, 26,  38, 9.6, "UVA axis\n2 studies", "white", UV, UV),
  bx(79, 26,  38, 9.6, "Senescence axis\n4 studies", "white", SEN, SEN),
  bx(50, 10,  66, 7.4, "D  =  ρ(UVA) − ρ(senescence)", "white", INK, INK, 6.3, 0.6))
R$xmin <- R$x - R$w/2; R$xmax <- R$x + R$w/2
R$ymin <- R$y - R$h/2; R$ymax <- R$y + R$h/2
sg <- function(x, xend, y, yend, a = FALSE, col = GREY, lw = 0.3)
  data.frame(x, xend, y, yend, a, col, lw)
S <- rbind(
  sg(47, 58.5, 96, 96, TRUE, UV, 0.4),                    # HDF -> UVA
  sg(76, 76, 92.5, 86.5), sg(17, 83, 86.5, 86.5),         # UVA -> CM split
  sg(17, 17, 86.5, 83, TRUE), sg(50, 50, 86.5, 83, TRUE), sg(83, 83, 86.5, 83, TRUE),
  sg(17, 17, 75.5, 71), sg(50, 50, 75.5, 71), sg(83, 83, 75.5, 71),
  sg(17, 83, 71, 71), sg(50, 50, 71, 66, TRUE),           # CM -> RNA-seq
  sg(50, 50, 58.3, 51, TRUE),                             # RNA-seq -> anchor
  sg(50, 50, 41.8, 36), sg(21, 79, 36, 36),               # anchor -> two axes
  sg(21, 21, 36, 31.3, TRUE), sg(79, 79, 36, 31.3, TRUE),
  sg(21, 21, 21.2, 17), sg(79, 79, 21.2, 17), sg(21, 79, 17, 17),
  sg(50, 50, 17, 14, TRUE))
p1a <- ggplot() +
  geom_segment(data = S[!S$a, ], aes(x, y, xend = xend, yend = yend),
               colour = S$col[!S$a], linewidth = S$lw[!S$a]) +
  geom_segment(data = S[S$a, ], aes(x, y, xend = xend, yend = yend),
               colour = S$col[S$a], linewidth = S$lw[S$a], arrow = AR) +
  geom_rect(data = R, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = R$fill, colour = R$col, linewidth = R$lw) +
  geom_text(data = R[R$lab != "", ], aes(x, y, label = lab), family = FONT,
            colour = R$fc[R$lab != ""], size = pt(R$tsz[R$lab != ""]), lineheight = 0.94) +
  annotate("text", x = 50, y = 46, parse = TRUE, family = FONT, size = pt(6.2), colour = INK,
    label = "paste(italic(S),' = sign(',italic(d)[H],')  ×  min( |',italic(d)[H],'| , |',italic(d)[I],'| )')") +
  annotate("text", x = 50, y = 68.6, label = "conditioned medium, 24 h", family = FONT,
           size = pt(5.5), colour = GREY, fontface = "italic") +
  annotate("rect", xmin = 8, xmax = 92, ymin = -2.5, ymax = 4.5, fill = NA,
           colour = GREY_L, linewidth = 0.3, linetype = "22") +
  annotate("text", x = 50, y = 1, family = FONT, size = pt(5.9), colour = INK2,
           lineheight = 0.95,
           label = "60-perturbation compendium\n14 preregistered held-out tests") +
  scale_x_continuous(limits = c(-2, 102)) + scale_y_continuous(limits = c(-4, 100.5)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(3, 2, 2, 3))

## ================================ B. PCA ====================================
pca <- read.delim(file.path(CS, "local_pca_scores.tsv"))
ve  <- read.delim(file.path(CS, "local_pca_varexp.tsv"), header = FALSE)
key <- c(UVA0 = UV, HDF = GREY, Repro = REPRO, iPSC = NAVY)
nmk <- c(UVA0 = "UVA 0 h", HDF = "HDF-CM", Repro = "Repro-CM", iPSC = "iPSC-CM")
pca$col <- key[pca$sample]; pca$lab <- nmk[pca$sample]
pca$hj  <- c(UVA0 = -0.16, HDF = -0.16, Repro = 1.14, iPSC = 1.14)[pca$sample]
p1b <- ggplot(pca, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(size = 2.2, colour = pca$col) +
  geom_text(aes(label = lab), hjust = pca$hj, size = pt(6), family = FONT, colour = INK) +
  scale_x_continuous(expand = expansion(mult = 0.30)) +
  scale_y_continuous(expand = expansion(mult = 0.24)) +
  labs(x = sprintf("PC1   (%s%% of variance)", ve$V2[1]),
       y = sprintf("PC2   (%s%%)", ve$V2[2])) +
  theme_j(7)

## ====================== C. sample correlation heatmap =======================
sc <- as.matrix(read.delim(file.path(CS, "local_sample_correlation.tsv"), check.names = FALSE))
dimnames(sc) <- list(nmk[rownames(sc)], nmk[colnames(sc)])
hc  <- hclust(as.dist(1 - sc), method = "average"); ord <- hc$order
lg  <- expand.grid(r = rownames(sc)[ord], c = colnames(sc)[ord], stringsAsFactors = FALSE)
lg$v <- mapply(function(a, b) sc[a, b], lg$r, lg$c)
lg$r <- factor(lg$r, levels = rownames(sc)[ord]); lg$c <- factor(lg$c, levels = colnames(sc)[ord])
dn <- dendro_segments(hc)
dn$x <- 4.56 + 0.70 * dn$x/max(dn$x); dn$xend <- 4.56 + 0.70 * dn$xend/max(dn$xend)
p1c <- ggplot(lg, aes(c, r, fill = v)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sub("^0", ".", sprintf("%.3f", v))), size = pt(5.7), family = FONT,
            colour = ifelse(lg$v > 0.95, "white", INK)) +
  geom_segment(data = dn, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
               colour = INK, linewidth = 0.25) +
  scale_fill_gradientn(colours = c("#F7FBFF","#C6DBEF","#6BAED6","#2166AC"),
                       limits = c(0.86, 1), breaks = c(0.9, 1.0), labels = c(".90","1.0"),
                       name = "Spearman ρ",
                       guide = guide_colourbar(barwidth = unit(30,"pt"), barheight = unit(3.6,"pt"),
                         ticks.colour = "white", frame.colour = NA, title.position = "left",
                         title.vjust = 1, direction = "horizontal")) +
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  coord_cartesian(xlim = c(0.5, 5.42), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_j(7) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x = element_text(angle = 35, hjust = 1),
                     legend.position = "bottom", legend.justification = "left",
                     legend.box.background = element_blank(),
                     legend.margin = margin(t = -3, l = -4),
                     plot.margin = margin(2, 10, 2, 2))

## ================================ D. MA plot ================================
G <- read.delim(file.path(D, "local_repro_anchor/local_Repro_specific_gene_rank.tsv"))
G$cls <- ifelse(abs(G$Repro_specific_score) >= log2(1.25), "|S| ≥ log₂1.25",
         ifelse(G$comparator_consistent %in% c(TRUE,"TRUE"), "direction-consistent", "discordant"))
G$cls <- factor(G$cls, levels = c("discordant","direction-consistent","|S| ≥ log₂1.25"))
G <- G[order(G$cls), ]
COLS <- c(discordant = GREY_L, `direction-consistent` = "#A9C7DF", `|S| ≥ log₂1.25` = REPRO)
p1d <- ggplot(G, aes(mean_log2_expression, logFC_Repro_vs_HDF, colour = cls)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.25) +
  geom_hline(yintercept = c(-1,1)*log2(1.25), colour = GREY, linewidth = 0.22, linetype = "22") +
  geom_point(size = 0.2, alpha = 0.5, stroke = 0) +
  scale_colour_manual(values = COLS, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 3)) +
  coord_cartesian(ylim = c(-3.1, 3.6)) +
  labs(x = expression(paste("mean ", log[2], " expression")),
       y = expression(paste(log[2], "FC  Repro-CM / HDF-CM"))) +
  theme_j(7) + theme(legend.position = "top", legend.justification = "right",
                     legend.box.background = element_blank(),
                     legend.margin = margin(b = -4, r = -2))

## ========================== E. two-comparator quadrants =====================
q  <- G[is.finite(G$logFC_Repro_vs_HDF) & is.finite(G$logFC_Repro_vs_iPSC), ]
nq <- c(pp = sum(q$logFC_Repro_vs_HDF > 0 & q$logFC_Repro_vs_iPSC > 0),
        nn = sum(q$logFC_Repro_vs_HDF < 0 & q$logFC_Repro_vs_iPSC < 0),
        pn = sum(q$logFC_Repro_vs_HDF > 0 & q$logFC_Repro_vs_iPSC < 0),
        np = sum(q$logFC_Repro_vs_HDF < 0 & q$logFC_Repro_vs_iPSC > 0))
L <- 3.6
p1e <- ggplot(q, aes(logFC_Repro_vs_HDF, logFC_Repro_vs_iPSC, colour = cls)) +
  annotate("rect", xmin = 0, xmax = L, ymin = 0, ymax = L, fill = REPRO, alpha = 0.045) +
  annotate("rect", xmin = -L, xmax = 0, ymin = -L, ymax = 0, fill = REPRO, alpha = 0.045) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_point(size = 0.2, alpha = 0.5, stroke = 0) +
  scale_colour_manual(values = COLS, guide = "none") +
  annotate("text", x = c(L-0.15, -L+0.15, L-0.15, -L+0.15), y = c(L-0.2, -L+0.2, -L+0.2, L-0.2),
           hjust = c(1, 0, 1, 0), label = format(nq[c("pp","nn","pn","np")], big.mark = ","),
           family = FONT, size = pt(6.2), fontface = "bold",
           colour = c(REPRO, REPRO, GREY, GREY)) +
  annotate("text", x = 0, y = L - 0.2, family = FONT, size = pt(5.9), colour = INK2,
           label = sprintf("%.1f%% concordant", 100*sum(nq[1:2])/sum(nq))) +
  coord_cartesian(xlim = c(-L, L), ylim = c(-L, L)) +
  labs(x = expression(paste(log[2], "FC vs HDF-CM   (", italic(d)[H], ")")),
       y = expression(paste(log[2], "FC vs iPSC-CM   (", italic(d)[I], ")"))) +
  theme_j(7)

right <- lab_grid(p1b, p1c, p1d, p1e, labels = c("B","C","D","E"), ncol = 2,
                  rel_heights = c(1, 1.06))
fig   <- lab_grid(p1a, right, labels = c("A",""), ncol = 2, rel_widths = c(0.72, 1))
save_fig(fig, "Fig1.png", 183, 108)
