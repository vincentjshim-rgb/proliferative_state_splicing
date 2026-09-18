## Supplementary Fig. S2  Which genes keep an age association once proliferation
## is accounted for? (six-figure restructure, 2026-09-17)
## Panels a-d of the former main Fig. 3 (make_sa_fig5gene.R); its panel e, the count
## under four cohort definitions, is now main Fig. 3d (make_sa_fig3merged.R).
## Primary cohort: 107 normal donors aged 20-96 (scripts/revision/cohort_gse113957.R),
## repository, instrument, sex and depth as covariates.
## Two corrections against the former figure:
##   c  18 genes whose retained fraction lies outside the axis were silently not
##      drawn; they are now drawn at the axis limit as open triangles
##   d  gene_level.tsv holds the age coefficient per decade in log2 CPM
##      (run_cohort_core.R fits lc ~ I(age/10) + ...); the former panel multiplied it
##      by ten and called it "z per decade". The table values are now shown as they are.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; C <- file.path(D, "cohort_revised")
Z  <- read.delim(file.path(C, "primary/heatmap_matrix.tsv"), check.names = FALSE)
SM <- read.delim(file.path(C, "primary/heatmap_sample_order.tsv"))
G  <- read.delim(file.path(C, "primary/gene_level.tsv"))
KEPT <- readLines(file.path(C, "primary/genes_keeping_decline.txt"))
NM <- read.delim(file.path(D, "division_rate/named_splicing_factors.tsv"))
EXPC <- "#B2182B"; OUTC <- "#2166AC"; PALE <- "#D8A9A4"
say <- function(...) cat(sprintf(...), "\n")

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
say("a  %d genes x %d donors; named factors in the set: %d (%s); of these keeping a decline: %d",
    nrow(Mc), ncol(Mc), nrow(LB), paste(LB$gene, collapse = ", "), sum(LB$ret))
pa <- ggplot(lg, aes(si, gi, fill = v)) + geom_raster() +
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
  annotate("text", x = ncol(Mc)/2, y = -7, label = sprintf("%d donors, ordered young to old", ncol(Mc)),
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
  scale_x_continuous(expand = c(0,0), limits = c(-42, ncol(Mc) + 58)) +
  scale_y_continuous(expand = c(0,0), limits = c(-22, nrow(Mc) + 14)) +
  labs(x = NULL, y = NULL) + coord_cartesian(clip = "off") +
  theme_sa(8) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), legend.position = "bottom",
    legend.justification = "center", legend.margin = margin(t = 2),
    legend.title = element_text(size = 7), legend.text = element_text(size = 7),
    plot.margin = margin(12, 4, 3, 3))

## ============ b. how strongly each gene is coupled to proliferation ========
say("b  correlation with the proliferation score across %d donors, %d genes: median %.4f, range %.4f to %.4f, negative in %d",
    ncol(Mc), nrow(G), median(G$r_prolif), min(G$r_prolif), max(G$r_prolif), sum(G$r_prolif < 0))
hb <- hist(G$r_prolif, breaks = seq(-0.72, 0.96, 0.06), plot = FALSE)
pb <- ggplot(G, aes(r_prolif)) +
  geom_histogram(breaks = seq(-0.72, 0.96, 0.06), fill = "#C3CBD4", colour = "white", linewidth = 0.2) +
  annotate("segment", x = median(G$r_prolif), xend = median(G$r_prolif), y = 0, yend = max(hb$counts) * 1.12,
           colour = EXPC, linewidth = 0.5, linetype = "22") +
  annotate("text", x = median(G$r_prolif) - 0.03, y = max(hb$counts) * 1.12, vjust = 1, hjust = 1,
           label = sprintf("median %s", num(median(G$r_prolif))), size = pt(7),
           family = FONT, colour = EXPC) +
  scale_x_continuous(labels = num_axis()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = "correlation with proliferation (r)", y = sprintf("genes (n = %d)", nrow(G))) +
  theme_sa(8) + theme(axis.title.y = element_text(margin = margin(r = 6)),
                      plot.margin = margin(8, 12, 3, 3))

## ====== c. genes coupled to proliferation lose their age association ======
YL <- c(-0.6, 1.7)
G$surv <- ifelse(G$FDR_adj < 0.05, "keeps a decline", "no association after adjustment")
ct <- cor.test(G$r_prolif, G$retained, method = "spearman", exact = FALSE)
P <- G[is.finite(G$retained), ]
P$off <- P$retained < YL[1] | P$retained > YL[2]
P$yy <- pmin(pmax(P$retained, YL[1]), YL[2])
say("c  rho = %.4f, P = %.3g, n = %d genes; off the axis: %d (%d above, %d below), of which age-associated before adjustment: %d",
    ct$estimate, ct$p.value, sum(is.finite(G$retained)), sum(P$off), sum(P$retained > YL[2]),
    sum(P$retained < YL[1]), sum(P$off & P$FDR < 0.05))
say("c  coupling of the %d genes keeping a decline: median r %.4f; the other %d: %.4f; Wilcoxon P = %.3g",
    sum(G$gene %in% KEPT), median(G$r_prolif[G$gene %in% KEPT]), sum(!G$gene %in% KEPT),
    median(G$r_prolif[!G$gene %in% KEPT]),
    wilcox.test(G$r_prolif[G$gene %in% KEPT], G$r_prolif[!G$gene %in% KEPT])$p.value)
CLS <- c("keeps a decline", "no association after adjustment", "ratio beyond the axis")
P$cls <- factor(ifelse(P$off, CLS[3], P$surv), levels = CLS)
stopifnot(!any(P$off & P$surv == CLS[1]))
pc <- ggplot(P, aes(r_prolif, yy, colour = cls, shape = cls)) +
  geom_hline(yintercept = c(0, 1), colour = GREY_L, linewidth = 0.3) +
  geom_point(size = 1.4, stroke = 0.4) +
  annotate("text", x = 0.98, y = 1.57, hjust = 1, vjust = 1, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.05,
           label = sprintf("%s = %s, n = %d genes\n%s", RHO, num(ct$estimate),
                           sum(is.finite(G$retained)), pfmt(ct$p.value))) +
  scale_colour_manual(values = setNames(c(OUTC, PALE, PALE), CLS), name = NULL) +
  scale_shape_manual(values = setNames(c(16, 16, 2), CLS), name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  scale_x_continuous(labels = num_axis()) +
  scale_y_continuous(labels = num_axis(), limits = YL, expand = expansion(mult = 0.03)) +
  labs(x = "correlation with proliferation (r)", y = "age effect retained\n(adjusted / unadjusted)") +
  theme_sa(8) + theme(legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7), legend.spacing.y = unit(0, "pt"),
                      axis.title.y = element_text(margin = margin(r = 4), lineheight = 0.95),
                      plot.margin = margin(8, 12, 3, 3))

## ============ d. the genes that keep a decline ============================
## beta and beta_adj are the age coefficients per decade of donor age in log2 CPM
R <- G[G$gene %in% KEPT, ]
R <- R[order(R$beta_adj), ]; R$gene <- factor(R$gene, levels = R$gene)
say("d  %d of %d genes age-associated (FDR < 0.05; all negative: %s); %d keep a decline; new after adjustment: %d",
    sum(G$FDR < 0.05), nrow(G), all(G$beta[G$FDR < 0.05] < 0), nrow(R), sum(G$FDR_adj < 0.05 & G$FDR >= 0.05))
say("d  unadjusted log2 CPM per decade: %.4f to %.4f; adjusted: %.4f to %.4f",
    min(R$beta), max(R$beta), min(R$beta_adj), max(R$beta_adj))
print(data.frame(gene = R$gene, beta = signif(R$beta, 4), beta_adj = signif(R$beta_adj, 4),
                 FDR = signif(R$FDR, 2), FDR_adj = signif(R$FDR_adj, 2), r_prolif = signif(R$r_prolif, 3)),
      row.names = FALSE)
XT <- 0.012
## The arrows ran to the centre of the adjusted symbol, so every 2.4 pt head was hidden
## under it. Each arrow now stops at the edge of that symbol (its radius is about 0.0062
## on this axis at the printed width) and is drawn over the unadjusted symbol. Where the
## change is smaller than the symbol plus a head the two symbols touch or overlap and a
## plain segment is kept: there is no room for a head and no distance to show.
RAD <- 0.0062; HEAD <- 0.0068                       # symbol radius, head length (x units)
R$dx   <- R$beta_adj - R$beta
R$tip  <- R$beta_adj - sign(R$dx) * (RAD + 0.0008)
R$room <- abs(R$dx) >= RAD + HEAD + 0.0008
say("d  arrows drawn for %d of %d genes (change of at least %.4f); plain segment for: %s",
    sum(R$room), nrow(R), RAD + HEAD + 0.0008, paste(R$gene[!R$room], collapse = ", "))
pd <- ggplot(R, aes(y = gene)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(data = R[!R$room, ], aes(x = beta, xend = beta_adj, yend = gene),
               colour = GREY, linewidth = 0.4) +
  geom_point(aes(x = beta), colour = "#E8C4C0", size = 1.5) +
  geom_segment(data = R[R$room, ], aes(x = beta, xend = tip, yend = gene), colour = INK2,
               linewidth = 0.4, linejoin = "mitre",
               arrow = arrow(length = unit(3, "pt"), angle = 25, type = "closed")) +
  geom_point(aes(x = beta_adj), colour = OUTC, size = 1.5) +
  annotate("text", x = XT, y = nrow(R) + 0.4, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.05,
           label = sprintf("%d of %d genes\nchange with age;\n%d keep a decline\nafter adjustment",
                           sum(G$FDR < 0.05), nrow(G), nrow(R))) +
  annotate("point", x = XT + 0.006, y = 7.2, colour = "#E8C4C0", size = 1.5) +
  annotate("text",  x = XT + 0.018, y = 7.2, hjust = 0, family = FONT,
           size = pt(7), colour = INK2, label = "unadjusted") +
  annotate("point", x = XT + 0.006, y = 5.6, colour = OUTC, size = 1.5) +
  annotate("text",  x = XT + 0.018, y = 5.6, hjust = 0, family = FONT, lineheight = 0.95,
           size = pt(7), colour = OUTC, label = "adjusted for\nproliferation") +
  ## the first layer holds a subset of the genes, so the order is fixed explicitly
  scale_y_discrete(limits = levels(R$gene)) +
  scale_x_continuous(sprintf("age effect (log%s CPM per decade)", "₂"),
                     breaks = c(-0.1, -0.05, 0), labels = num_axis(), limits = c(-0.13, 0.155)) +
  theme_sa(8) + theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 7, face = "italic"),
        axis.title.x = element_text(hjust = 0.2),
        plot.margin = margin(9, 7, 3, 3))

mid <- lab_grid(pb, pc, labels = c("b", "c"), ncol = 1, rel_heights = c(0.82, 1.18))
save_fig(lab_grid(pa, mid, pd, labels = c("a", "", "d"), ncol = 3,
                  rel_widths = c(1.28, 0.96, 0.96)), "FigS2.png", 183, 128)
