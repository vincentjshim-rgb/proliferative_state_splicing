FIGURES_WRITTEN <- c("Fig3.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

## Fig. 3  Which genes keep an age association once proliferation is accounted for?
## Revised 2026-09-15: rebuilt on the redefined cohort (107 normal donors aged 20-96,
## progeria donors and children removed, repository/instrument/sex as covariates).
## 137 of the 177 genes change with age, 19 keep a decline after adjustment and none
## acquire one; panel e shows how unstable that count is across cohort definitions.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; C <- file.path(D, "cohort_revised")
Z  <- read.delim(file.path(C, "primary/heatmap_matrix.tsv"), check.names = FALSE)
SM <- read.delim(file.path(C, "primary/heatmap_sample_order.tsv"))
G  <- read.delim(file.path(C, "primary/gene_level.tsv"))
KEPT <- readLines(file.path(C, "primary/genes_keeping_decline.txt"))
SENS <- read.delim(file.path(C, "sensitivity_genes.tsv"))
NM <- read.delim(file.path(D, "division_rate/named_splicing_factors.tsv"))
EXPC <- "#B2182B"; OUTC <- "#2166AC"

## ================ a. the 177 genes across 107 donors, ordered by age =======
g <- Z$gene; Mz <- as.matrix(Z[, -1]); rownames(Mz) <- g
ord <- order(G$r_prolif[match(g, G$gene)])          # genes ordered by coupling to proliferation
Mz <- Mz[ord, ]; Mc <- pmax(pmin(Mz, 2), -2)
lg <- expand.grid(gi = seq_len(nrow(Mc)), si = seq_len(ncol(Mc)))
lg$v <- as.vector(Mc)
## label the named splicing factors that belong to the set: eleven of the twenty
## assayed by Holly et al. 2013 (scripts/revision/named_splicing_factors_source.md)
cand <- intersect(NM$gene, rownames(Mc))
sel <- sort(match(cand, rownames(Mc)))
LB <- data.frame(gene = rownames(Mc)[sel], gi = sel)
LB$ly <- seq(3, nrow(Mc) - 2, length.out = nrow(LB))
LB$ret <- LB$gene %in% KEPT
p5a <- ggplot(lg, aes(si, gi, fill = v)) + geom_raster() +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-2, 2), breaks = c(-2, 0, 2),
    labels = num_axis(), name = "expression z",
    guide = guide_colourbar(barwidth = unit(40,"pt"), barheight = unit(3.6,"pt"),
      ticks.colour = "white", frame.colour = NA, title.position = "top")) +
  annotate("rect", xmin = seq_len(nrow(SM)) - 0.5, xmax = seq_len(nrow(SM)) + 0.5,
           ymin = nrow(Mc) + 1.5, ymax = nrow(Mc) + 6.5,
           fill = colorRampPalette(c("#F2E8DC","#8C6D4F"))(100)[cut(SM$age, 100, labels = FALSE)],
           colour = NA) +
  annotate("rect", xmin = seq_len(nrow(SM)) - 0.5, xmax = seq_len(nrow(SM)) + 0.5,
           ymin = nrow(Mc) + 7.5, ymax = nrow(Mc) + 12.5,
           fill = colorRampPalette(c("#FFFFFF","#B2182B"))(100)[cut(SM$prolif, 100, labels = FALSE)],
           colour = NA) +
  annotate("text", x = -2, y = nrow(Mc) + 4, label = "donor age", hjust = 1,
           size = pt(7), family = FONT, colour = INK2) +
  annotate("text", x = -2, y = nrow(Mc) + 10, label = "proliferation", hjust = 1,
           size = pt(7), family = FONT, colour = INK2) +
  annotate("text", x = ncol(Mc)/2, y = -7, label = "107 donors, ordered young to old",
           size = pt(7), family = FONT, colour = INK) +
  annotate("text", x = -20, y = nrow(Mc)/2, hjust = 0.5, angle = 90, size = pt(7),
           family = FONT, colour = INK, lineheight = 0.95,
           label = sprintf("%d genes, ordered by\ncoupling to proliferation", nrow(Mc))) +
  annotate("segment", x = ncol(Mc) + 0.5, xend = ncol(Mc) + 8, y = LB$gi, yend = LB$ly,
           colour = GREY_L, linewidth = 0.2) +
  annotate("text", x = ncol(Mc) + 9.5, y = LB$ly, label = LB$gene, hjust = 0,
           size = pt(7), family = FONT, fontface = "italic",
           colour = ifelse(LB$ret, OUTC, INK2)) +
  annotate("text", x = ncol(Mc)/2, y = -15, size = pt(7), family = FONT, colour = INK2,
           lineheight = 1.05,
           label = sprintf("labelled: the %d named factors in the set;\n%s",
                           nrow(LB),
                           if (any(LB$ret)) "blue, keeps a decline after adjustment"
                           else "none keeps a decline after adjustment")) +
  scale_x_continuous(expand = c(0,0), limits = c(-42, ncol(Mc) + 46)) +
  scale_y_continuous(expand = c(0,0), limits = c(-22, nrow(Mc) + 14)) +
  labs(x = NULL, y = NULL) + coord_cartesian(clip = "off") +
  theme_sa(8) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), legend.position = "bottom",
    legend.justification = "center", legend.margin = margin(t = 2),
    legend.title = element_text(size = 7), legend.text = element_text(size = 7),
    plot.margin = margin(12, 4, 3, 3))

## ============ b. how strongly each gene is coupled to proliferation ========
p5b <- ggplot(G, aes(r_prolif)) +
  geom_histogram(binwidth = 0.06, fill = "#C3CBD4", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = median(G$r_prolif), colour = EXPC, linewidth = 0.5, linetype = "22") +
  annotate("text", x = median(G$r_prolif), y = Inf, vjust = 1.4, hjust = 1.06,
           label = sprintf("median %s", num(median(G$r_prolif))), size = pt(7),
           family = FONT, colour = EXPC) +
  scale_x_continuous(labels = num_axis()) +
  labs(x = "correlation with proliferation", y = "genes") +
  theme_sa(8) + theme(axis.title.y = element_text(margin = margin(r = 6)),
                      plot.margin = margin(3, 9, 3, 3))

## ====== c. genes coupled to proliferation lose their age association ======
G$surv <- ifelse(G$FDR_adj < 0.05, "keeps a decline", "no association after adjustment")
ct <- cor.test(G$r_prolif, G$retained, method = "spearman", exact = FALSE)
p5c <- ggplot(G[is.finite(G$retained), ], aes(r_prolif, retained)) +
  geom_hline(yintercept = c(0, 1), colour = GREY_L, linewidth = 0.3) +
  geom_point(aes(colour = surv), size = 1.5) +
  annotate("text", x = -0.62, y = 1.62, hjust = 0, vjust = 1, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.05,
           label = sprintf("%s\nn = %d genes", rp(ct$estimate, ct$p.value), nrow(G))) +
  scale_colour_manual(values = c(`keeps a decline` = OUTC,
                                 `no association after adjustment` = "#D8A9A4"), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = num_axis()) + scale_y_continuous(labels = num_axis()) +
  coord_cartesian(ylim = c(-0.6, 1.7)) +
  labs(x = "gene's correlation with proliferation", y = "age effect retained") +
  theme_sa(8) + theme(legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7),
                      axis.title.y = element_text(margin = margin(r = 6)),
                      plot.margin = margin(3, 9, 3, 3))

## ============ d. the genes that keep a decline ============================
R <- G[G$gene %in% KEPT, ]
R$b1 <- R$beta * 10; R$b2 <- R$beta_adj * 10        # per decade of donor age
R <- R[order(R$b2), ]; R$gene <- factor(R$gene, levels = R$gene)
p5d <- ggplot(R, aes(y = gene)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(aes(x = b1, xend = b2, yend = gene), colour = GREY, linewidth = 0.4,
               arrow = arrow(length = unit(2.4, "pt"), type = "closed")) +
  geom_point(aes(x = b1), colour = "#E8C4C0", size = 1.5) +
  geom_point(aes(x = b2), colour = OUTC, size = 1.5) +
  annotate("text", x = 0.12, y = nrow(R) + 0.4, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.05,
           label = sprintf("%d of %d genes\nchange with age;\n%d keep a decline\nafter adjustment",
                           sum(G$FDR < 0.05), nrow(G), nrow(R))) +
  annotate("point", x = 0.18, y = 7.2, colour = "#E8C4C0", size = 1.5) +
  annotate("text",  x = 0.30, y = 7.2, hjust = 0, family = FONT,
           size = pt(7), colour = INK2, label = "unadjusted") +
  annotate("point", x = 0.18, y = 5.2, colour = OUTC, size = 1.5) +
  annotate("text",  x = 0.30, y = 5.2, hjust = 0, family = FONT, lineheight = 0.95,
           size = pt(7), colour = OUTC, label = "adjusted for\nproliferation") +
  scale_x_continuous("age effect (z per decade)", breaks = c(-1, -0.5, 0), labels = num_axis(),
                     limits = c(-1.3, 1.55)) +
  theme_sa(8) + theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 7, face = "italic"),
        plot.margin = margin(9, 7, 3, 3))

## ====== e. the count depends on how the cohort is defined =================
SENS$lab <- factor(c("normal adults 20+\n(primary, n = 107)", "all normal donors\n(n = 133)",
                     "normal adults 20-82\n(n = 76)", "published cohort\n(n = 142)")[
                     match(SENS$cohort, c("primary", "normal", "adult2082", "published"))],
                   levels = c("published cohort\n(n = 142)", "all normal donors\n(n = 133)",
                              "normal adults 20+\n(primary, n = 107)", "normal adults 20-82\n(n = 76)"))
SENS$primary <- SENS$cohort == "primary"
p5e <- ggplot(SENS, aes(decline_kept, lab, fill = primary)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = decline_kept), hjust = -0.35, size = pt(7), family = FONT, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = OUTC, `FALSE` = "#C3CBD4"), guide = "none") +
  scale_x_continuous(limits = c(0, 68), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
  annotate("text", x = 0, y = 5.1, hjust = 0, vjust = 0.5, size = pt(7), family = FONT,
           colour = INK2, lineheight = 1.05,
           label = "no gene keeps a decline\nin all four cohorts") +
  coord_cartesian(clip = "off") +
  labs(x = "genes keeping a decline", y = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7, lineheight = 0.95),
                      plot.margin = margin(3, 6, 3, 3))

mid   <- lab_grid(p5b + theme(plot.margin = margin(3, 12, 3, 3)),
                  p5c + theme(plot.margin = margin(3, 12, 3, 3)),
                  labels = c("b", "c"), ncol = 1, rel_heights = c(0.86, 1.14))
right <- lab_grid(p5d, p5e, labels = c("d", "e"), ncol = 1, rel_heights = c(1.25, 1))
save_fig(lab_grid(p5a, mid, right, labels = c("a", "", ""), ncol = 3,
                  rel_widths = c(1.30, 0.86, 1.02)), "Fig3.png", 183, 138)
