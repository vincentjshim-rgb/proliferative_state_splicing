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
  annotate("text", x = ncol(Mc)/2, y = -17, size = pt(7), family = FONT, colour = INK2,
           lineheight = 1.05,
           label = sprintf("labelled: the %d named factors in the set;\n%s",
                           nrow(LB),
                           if (any(LB$ret)) "blue, keeps a decline after adjustment"
                           else "none keeps a decline after adjustment")) +
  scale_x_continuous(expand = c(0,0), limits = c(-42, ncol(Mc) + 58)) +
  scale_y_continuous(expand = c(0,0), limits = c(-25, nrow(Mc) + 14)) +
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
## paired bars (2026-09-19) in place of the arrow plot
RL <- rbind(data.frame(gene = R$gene, k = "unadjusted", b = R$beta),
            data.frame(gene = R$gene, k = "adjusted for proliferation", b = R$beta_adj))
RL$k <- factor(RL$k, levels = c("unadjusted", "adjusted for proliferation"))
pd <- ggplot(RL, aes(b, gene, fill = k)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, colour = NA) +
  scale_fill_manual(values = c(unadjusted = "#D9A9A4", `adjusted for proliferation` = OUTC), name = NULL) +
  scale_y_discrete(limits = levels(R$gene)) +
  scale_x_continuous(sprintf("age effect (log%s CPM per decade)", "\u2082"),
                     breaks = c(-0.1, -0.05, 0), labels = num_axis(), limits = c(-0.13, 0.012)) +
  labs(subtitle = sprintf("%d of %d genes change with age;\n%d keep a decline after adjustment",
                          sum(G$FDR < 0.05), nrow(G), nrow(R))) +
  theme_sa(8) + theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 7, face = "italic"),
        legend.position = "bottom", legend.direction = "vertical", legend.margin = margin(t = -4),
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = 7),
        plot.subtitle = element_text(size = 7, colour = INK2, lineheight = 1.05),
        plot.margin = margin(3, 7, 3, 3))

## ============ d-f. moved here from main Fig. 4 on 2026-09-20 =============
## The cohort-definition bars, the share-removed histogram and the count of genes keeping a
## decline under each definition left the main figure, which now shows proliferation alone.
## Every number is unchanged; the code is that of the former Fig. 4 panels b, c (left) and d.
S    <- read.delim(file.path(C, "sensitivity_metrics.tsv"))
SENS <- read.delim(file.path(C, "sensitivity_genes.tsv"))
NL   <- read.delim(file.path(D, "reviewer_sensitivities/adjustment_negative_control.tsv"))
QC   <- read.delim(file.path(D, "reviewer_sensitivities/qc_covariates.tsv"))
NCORE <- length(readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt")))
COHORTS <- c("published", "normal", "primary", "adult2082")
COHLAB  <- c(published = "deposited set\nn = 142",
             normal    = "all normal\nn = 133",
             primary   = "adults 20+\nn = 107, primary",
             adult2082 = "adults 20\u201382\nn = 76")
stopifnot(all(S$n[match(COHORTS, S$cohort)] == c(142, 133, 107, 76)))
K <- S[S$metric == "machinery", ]; K <- K[match(COHORTS, K$cohort), ]
K$row <- rev(seq_along(COHORTS))
FB <- rbind(
  data.frame(cohort = K$cohort, row = K$row, k = "unadjusted", b = K$beta, lo = K$lo, hi = K$hi, p = K$p),
  data.frame(cohort = K$cohort, row = K$row, k = "adjusted for proliferation", b = K$beta_adj,
             lo = K$lo_adj, hi = K$hi_adj, p = K$p_adj))
FB$k <- factor(FB$k, levels = c("unadjusted", "adjusted for proliferation"))
for (i in seq_len(nrow(K)))
  say("d  %-10s n = %3d  unadjusted %.4f [%.4f, %.4f] P = %.3g | adjusted %.4f [%.4f, %.4f] P = %.3g | removed %.1f%%",
      K$cohort[i], K$n[i], K$beta[i], K$lo[i], K$hi[i], K$p[i], K$beta_adj[i], K$lo_adj[i],
      K$hi_adj[i], K$p_adj[i], K$pct_lost[i])
FB$cohort <- factor(FB$cohort, levels = K$cohort[order(K$row, decreasing = TRUE)])
LABK <- K[order(K$row, decreasing = TRUE), ]
LABK$cohort <- factor(LABK$cohort, levels = levels(FB$cohort))
pe <- ggplot(FB, aes(cohort, b, fill = k)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.72), width = 0.22,
                colour = INK2, linewidth = 0.4) +
  geom_text(data = LABK, aes(cohort, y = 0.014, label = sprintf("%.0f%%", pct_lost)),
            inherit.aes = FALSE, family = FONT, size = pt(7), colour = INK2, vjust = 0) +
  scale_fill_manual(values = c(unadjusted = EXPC, `adjusted for proliferation` = GREY), name = NULL) +
  scale_x_discrete(labels = COHLAB) +
  scale_y_continuous(limits = c(-0.235, 0.06), breaks = c(-0.2, -0.1, 0), labels = num_axis(),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "age effect per decade,\npre-mRNA processing score",
       subtitle = "share removed printed above each pair") +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     plot.subtitle = element_text(size = 7, colour = INK2, hjust = 0),
                     legend.margin = margin(b = -6), legend.key.size = unit(6, "pt"),
                     axis.text.x = element_text(size = 7, lineheight = 0.95),
                     plot.margin = margin(3, 4, 3, 6))

OBS <- K[K$cohort == "primary", ]
q <- QC[QC$metric == "machinery" & QC$covariates == "cohort", ]
stopifnot(isTRUE(all.equal(q$beta, OBS$beta)), isTRUE(all.equal(q$pct_lost, OBS$pct_lost)))
sig <- is.finite(NL$p) & NL$p < 0.05
LS <- NL[sig, ]
say("e  %% removed among the %d random sets with an age association: median %.1f; splicing set %.1f%%; at least as much removed: %d of %d",
    sum(sig), median(LS$lost), OBS$pct_lost, sum(LS$lost >= OBS$pct_lost), nrow(LS))
h1 <- hist(LS$lost, breaks = seq(-20, 120, 10), plot = FALSE); top1 <- max(h1$counts)
pf <- ggplot(LS, aes(lost)) +
  geom_histogram(breaks = seq(-20, 120, 10), fill = "#C3CBD4", colour = "white", linewidth = 0.2) +
  annotate("segment", x = median(LS$lost), xend = median(LS$lost), y = 0, yend = top1 * 1.10,
           colour = INK2, linewidth = 0.45, linetype = "22") +
  annotate("segment", x = OBS$pct_lost, xend = OBS$pct_lost, y = 0, yend = top1 * 1.10,
           colour = EXPC, linewidth = 0.6) +
  annotate("text", x = median(LS$lost) - 3, y = top1 * 1.12, hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = INK2, lineheight = 0.95,
           label = sprintf("random sets\nmedian %.0f%%", median(LS$lost))) +
  annotate("text", x = OBS$pct_lost + 3, y = top1 * 1.12, hjust = 0, vjust = 0, family = FONT,
           size = pt(7), colour = EXPC, lineheight = 0.95,
           label = sprintf("splicing set\n%.0f%%", OBS$pct_lost)) +
  scale_x_continuous(breaks = c(0, 50, 100)) +
  coord_cartesian(xlim = c(-20, 120)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, top1 * 1.42)) +
  labs(x = "age effect removed (%)",
       y = sprintf("random gene sets (n = %d of %s\nwith an age association)", nrow(LS),
                   format(nrow(NL), big.mark = ","))) +
  theme_sa() + theme(axis.title.y = element_text(lineheight = 0.95, margin = margin(r = 5)),
                     plot.margin = margin(8, 6, 3, 3))

KL <- setNames(lapply(COHORTS, function(k) readLines(file.path(C, k, "genes_keeping_decline.txt"))), COHORTS)
KL <- lapply(KL, function(x) x[nzchar(x)])
stopifnot(all(lengths(KL) == SENS$decline_kept[match(COHORTS, SENS$cohort)]))
common <- Reduce(intersect, KL)
say("f  genes keeping a decline: %s; common to all four: %d",
    paste(sprintf("%s %d", COHORTS, lengths(KL)), collapse = ", "), length(common))
DB <- data.frame(lab = c(COHLAB[COHORTS], "in all four\ncohort definitions"),
                 n = c(lengths(KL), length(common)),
                 kind = c("other", "other", "primary", "other", "common"))
DB$lab <- factor(DB$lab, levels = rev(DB$lab))
pg <- ggplot(DB, aes(n, lab, fill = kind)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n, x = n + 1.6), hjust = 0, size = pt(7), family = FONT, colour = INK) +
  geom_hline(yintercept = 1.5, colour = GREY_L, linewidth = 0.3) +
  scale_fill_manual(values = c(primary = OUTC, other = "#C3CBD4", common = INK2), guide = "none") +
  scale_x_continuous(limits = c(0, 74), breaks = c(0, 30, 60), expand = expansion(mult = c(0, 0))) +
  labs(x = sprintf("genes keeping a decline\nafter adjustment (of %d)", NCORE), y = NULL) +
  theme_sa() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.y = element_text(size = 7, lineheight = 0.95),
                     axis.title.x = element_text(hjust = 1, lineheight = 0.95),
                     plot.margin = margin(8, 8, 3, 6))

right <- lab_grid(pb, pd, labels = c("b", "c"), ncol = 1, rel_heights = c(0.62, 1.38))
top <- lab_grid(pa, right, labels = c("a", ""), ncol = 2, rel_widths = c(1.2, 1))
bot <- lab_grid(pe, pf, pg, labels = c("d", "e", "f"), ncol = 3, rel_widths = c(1.35, 0.85, 0.8))
save_fig(lab_grid(top, bot, labels = c("", ""), ncol = 1, rel_heights = c(1.36, 0.72)),
         "FigS2.png", 183, 186)
