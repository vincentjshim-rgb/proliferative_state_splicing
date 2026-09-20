## Fig. 3  The machinery's own transcripts: the same genes measured as junction usage.
##
## Every other figure measures these genes as abundance, summarised as a set score. This one
## measures them a different way, from junction counts in the same donors: which junctions the
## transcripts of the machinery genes themselves use, and whether that moves with proliferation.
##   a  SRSF3, drawn from the junction counts: the exon model, and the two junctions that leave
##      the same donor site -- one skipping the short exon at 36,599,820 and one including it --
##      in the donors whose proliferation score is in the top and bottom third
##   b  inclusion of that exon in every donor against the proliferation score
##   c  the class-level answer, which is negative: events in the 177 machinery genes are no more
##      likely to track proliferation than events in other genes
## post hoc; analysis in scripts/revision/run_machinery_transcript_level.R
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"; O <- file.path(D, "machinery_transcript")
X <- readRDS(file.path(O, "junction_detail_SRSF3.rds"))
R <- fread(file.path(O, "psi_vs_proliferation.tsv"))
J <- X$J; ANN <- as.data.table(X$ann); MK <- as.data.frame(X$meta); EX <- as.data.table(X$exons)
GENE <- X$gene
stopifnot(GENE == "SRSF3", ncol(J) == 107)

## the competing pair: same donor site, one junction includes the short exon and one skips it
DON <- 36598984L; INCL_END <- 36599820L; SKIP_END <- 36601151L
i_in <- which(ANN$start == DON & ANN$end == INCL_END)
i_sk <- which(ANN$start == DON & ANN$end == SKIP_END)
stopifnot(length(i_in) == 1, length(i_sk) == 1)
psi <- J[i_in, ] / (J[i_in, ] + J[i_sk, ])
ev  <- R[gene == GENE & start == DON & end == INCL_END]
stopifnot(nrow(ev) == 1)
hi <- MK$prolif >= quantile(MK$prolif, 2/3); lo <- MK$prolif <= quantile(MK$prolif, 1/3)
cat(sprintf("SRSF3 inclusion of the exon at %s: PSI %.3f in the slowest third, %.3f in the fastest\n",
            format(INCL_END, big.mark = ","), mean(psi[lo]), mean(psi[hi])))
cat(sprintf("  junction counts  include %.0f -> %.0f   skip %.0f -> %.0f  (low -> high third)\n",
            mean(J[i_in, lo]), mean(J[i_in, hi]), mean(J[i_sk, lo]), mean(J[i_sk, hi])))
rp <- cor.test(MK$prolif, psi, method = "spearman", exact = FALSE)
cat(sprintf("  PSI vs proliferation score: rho = %.3f, P = %.1e | model dPSI = %.3f, FDR = %.1e\n",
            rp$estimate, rp$p.value, ev$d_prolif, ev$FDR_prolif))

## ====== a. the locus, drawn from the junction counts ==========================
## Exon model from GENCODE v41. The annotation holds one row per transcript, and some
## transcripts read straight through the alternative exon; an exon that strictly contains
## another is that kind of read-through and is not drawn, so the short exon stays visible.
EXU <- unique(EX[, .(s, e)])[order(s, e)]
swallows <- sapply(seq_len(nrow(EXU)), function(i)
  any(EXU$s > EXU$s[i] & EXU$e < EXU$e[i]))
EXU <- EXU[!swallows]
mg <- list(); cur <- c(EXU$s[1], EXU$e[1])
for (k in seq_len(nrow(EXU))[-1]) {
  if (EXU$s[k] <= cur[2]) cur[2] <- max(cur[2], EXU$e[k]) else { mg[[length(mg) + 1]] <- cur
                                                                 cur <- c(EXU$s[k], EXU$e[k]) } }
mg[[length(mg) + 1]] <- cur
EM <- data.table(s = sapply(mg, `[`, 1), e = sapply(mg, `[`, 2))
cat(sprintf("  exon model: %d blocks, %s to %s\n", nrow(EM), format(min(EM$s), big.mark = ","),
            format(max(EM$e), big.mark = ",")))
## the window is the span of the junctions that are drawn, not the whole 3' tail
X0 <- min(EM$s) - 250; X1 <- max(ANN$end) + 350
EM <- EM[s < X1]; EM[e > X1, e := X1]

## every junction covered in half the donors, in each third
JD <- data.table(start = ANN$start, end = ANN$end,
                 low = rowMeans(J[, lo, drop = FALSE]), high = rowMeans(J[, hi, drop = FALSE]))
JL <- rbind(data.table(JD[, .(start, end)], n = JD$low,  grp = "slowest third"),
            data.table(JD[, .(start, end)], n = JD$high, grp = "fastest third"))
JL[, grp := factor(grp, levels = c("slowest third", "fastest third"))]
JL[, base := ifelse(grp == "slowest third", 0, 0)]
JL[, focus := (start == DON & end %in% c(INCL_END, SKIP_END))]
SC <- max(JL$n)
pA <- ggplot() +
  geom_rect(data = EM, aes(xmin = s, xmax = e, ymin = -0.09, ymax = 0.09),
            fill = INK, colour = NA) +
  annotate("segment", x = X0 + 150, xend = X1 - 150, y = 0, yend = 0, colour = INK, linewidth = 0.3) +
  annotate("rect", xmin = INCL_END - 90, xmax = 36599875, ymin = -0.17, ymax = 0.17,
           fill = NA, colour = RED, linewidth = 0.4) +
  geom_curve(data = JL[focus == FALSE], aes(x = start, xend = end, y = 0.1, yend = 0.1,
             linewidth = n), curvature = -0.30, colour = GREY_L, ncp = 12, lineend = "round") +
  geom_curve(data = JL[focus == TRUE & end == SKIP_END],
             aes(x = start, xend = end, y = 0.1, yend = 0.1, linewidth = n),
             curvature = -0.42, colour = BLUE, ncp = 12, lineend = "round") +
  geom_curve(data = JL[focus == TRUE & end == INCL_END],
             aes(x = start, xend = end, y = 0.1, yend = 0.1, linewidth = n),
             curvature = -0.22, colour = RED, ncp = 12, lineend = "round") +
  geom_text(data = JL[focus == TRUE & end == SKIP_END],
            aes(x = (start + end)/2, y = 1.02, label = sprintf("%.0f", n)),
            family = FONT, size = pt(7), colour = BLUE, vjust = 0) +
  ## the skip arc passes above the include arc everywhere they overlap, and at 177 reads it is
  ## thick enough to swallow a count set above the include arc; the count goes inside its own bow
  geom_label(data = JL[focus == TRUE & end == INCL_END],
             aes(x = (start + end)/2, y = 0.235, label = sprintf("%.0f", n)),
             family = FONT, size = pt(7), colour = RED, vjust = 0.5,
             fill = "white", label.size = 0, label.padding = unit(0.6, "pt")) +
  scale_linewidth_continuous(range = c(0.2, 1.9), guide = "none") +
  facet_wrap(~ grp, ncol = 1) +
  geom_text(data = data.table(grp = factor("slowest third", levels = levels(JL$grp))),
            aes(x = INCL_END + 250, y = -0.30), label = "the 26 bp exon", family = FONT,
            size = pt(7), colour = RED, hjust = 0, vjust = 0.5) +
  geom_segment(data = data.table(grp = factor("slowest third", levels = levels(JL$grp))),
               aes(x = INCL_END + 210, xend = INCL_END + 30, y = -0.30, yend = -0.17),
               colour = RED, linewidth = 0.3) +
  scale_x_continuous(sprintf("%s  (kb)", X$locus$chr),
                     limits = c(X0, X1),
                     breaks = seq(36594000, 36602000, by = 1000),
                     labels = function(v) sprintf("%.0f", v/1000)) +
  scale_y_continuous(NULL, limits = c(-0.46, 1.42), breaks = NULL) +
  labs(title = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      panel.spacing = unit(3, "pt"),
                      strip.text = element_text(size = 7.5, hjust = 0, margin = margin(b = 1)),
                      plot.margin = margin(9, 3, 3, 3))

## ====== b. inclusion against the proliferation score ==========================
PB <- data.frame(prolif = MK$prolif, psi = psi, age = MK$age)
pB <- ggplot(PB, aes(prolif, psi)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = INK, fill = GREY, alpha = 0.18,
              linewidth = 0.6) +
  geom_point(size = 1.15, alpha = 0.8, colour = RED) +
  annotate("text", x = min(PB$prolif), y = max(PB$psi), hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.05,
           label = sprintf("%s = %s\n%s\nn = %d donors", RHO, num(rp$estimate),
                           pfmt(rp$p.value), nrow(PB))) +
  scale_y_continuous("exon inclusion (PSI)", labels = num_axis()) +
  scale_x_continuous("proliferation score (20 markers)", labels = num_axis()) +
  theme_sa(8)

## ====== c. the class-level answer, which is negative ==========================
CL <- R[!is.na(gene) & gene != ""]
CL[, cls := ifelse(machinery, "the 177 machinery genes", "all other genes")]
CL[, cls := factor(cls, levels = c("all other genes", "the 177 machinery genes"))]
MMc <- CL[cls == "the 177 machinery genes"]; bgk <- CL[cls == "all other genes"]
wp <- wilcox.test(abs(MMc$t_prolif), abs(bgk$t_prolif))$p.value
pC <- ggplot(CL, aes(abs(t_prolif), colour = cls)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = c(`all other genes` = GREY, `the 177 machinery genes` = BLUE),
                      name = NULL) +
  ## the view stops at 8 but the distribution does not: a scale limit here would drop the 42
  ## events above it before stat_ecdf ran, and the drawn curve would no longer be the ECDF the
  ## printed medians come from. coord_cartesian truncates the view and leaves the statistic whole.
  scale_x_continuous("|t| for the event against the proliferation score", labels = num_axis(0)) +
  coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous("cumulative fraction of events", labels = num_axis(1)) +
  annotate("text", x = 7.9, y = 0.06, hjust = 1, vjust = 0, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.05,
           label = sprintf("median |t| %s against %s\n%s at FDR < 0.05: %.0f%% against %.0f%%\n%s",
                           num(median(abs(MMc$t_prolif))), num(median(abs(bgk$t_prolif))),
                           "events", 100*mean(MMc$FDR_prolif < 0.05), 100*mean(bgk$FDR_prolif < 0.05),
                           pfmt(wp))) +
  annotate("text", x = 0.1, y = 0.99, hjust = 0, vjust = 1, family = FONT, size = pt(7),
           colour = INK2, label = sprintf("%d events beyond the axis, all in other genes",
                                          sum(abs(bgk$t_prolif) > 8))) +
  theme_sa(8) + theme(legend.position = "inside", legend.position.inside = c(0.62, 0.62), legend.key.height = unit(8, "pt"),
                      legend.text = element_text(size = 7))

cat(sprintf("panel c: machinery %d events in %d genes, median |t| %.2f; background %d events, %.2f; Wilcoxon P = %.3f\n",
            nrow(MMc), length(unique(MMc$gene)), median(abs(MMc$t_prolif)), nrow(bgk),
            median(abs(bgk$t_prolif)), wp))

save_fig(lab_grid(lab_grid(pA, labels = "a", ncol = 1),
                  lab_grid(pB, pC, labels = c("b", "c"), ncol = 2, rel_widths = c(1, 1.22)),
                  labels = c("", ""), ncol = 1, rel_heights = c(1.06, 1)),
         "Fig3.png", 183, 128)
