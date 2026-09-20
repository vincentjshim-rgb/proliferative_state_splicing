## Fig. 3  The machinery's own transcripts: the same genes measured as junction usage.
##
## Every other figure measures these genes as abundance, summarised as a set score. This one
## measures them a different way, from junction counts in the same donors: which junctions the
## transcripts of the machinery genes themselves use, and whether that moves with proliferation.
##   a  SRSF3, drawn from the junction counts: the exon model, and the two paths that leave
##      exon 3 -- one skipping the 456-bp exon that GENCODE annotates only in nonsense-mediated-
##      decay transcripts, one including it -- in the donors whose proliferation score is in the
##      top and bottom third
##   b  every event in a machinery gene that moves with proliferation at FDR < 0.05: what the
##      event is, what kind of exon it involves, and how far the measured path shifts
##   c  the events in the named splicing factors, donor by donor
##   d  what kinds of event these are, against the events in every other gene
##   e  the class-level answer, which is negative: events in the 177 machinery genes are no more
##      likely to track proliferation than events in other genes
## post hoc; analysis in scripts/revision/run_machinery_transcript_level.R (tests) and
## scripts/revision/run_machinery_events.R (classification against GENCODE v41)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"; O <- file.path(D, "machinery_transcript")
X  <- readRDS(file.path(O, "junction_detail_SRSF3.rds"))
R  <- fread(file.path(O, "psi_vs_proliferation.tsv"))
EV <- fread(file.path(O, "events_classified_machinery.tsv"))
SM <- fread(file.path(O, "event_type_summary.tsv"))
NE <- readRDS(file.path(O, "named_events_psi.rds"))
EXG <- fread(file.path(O, "gencode_v41_exons_event_genes.tsv"))
J <- X$J; ANN <- as.data.table(X$ann); MK <- as.data.frame(X$meta)
GENE <- X$gene
stopifnot(GENE == "SRSF3", ncol(J) == 107)
PRIME <- function(s) gsub("3'", "3′", gsub("5'", "5′", s, fixed = TRUE), fixed = TRUE)

## ====== a. SRSF3, drawn from the junction counts =============================
## the event: from the donor site of exon 3, one junction skips to exon 5 and one includes the
## exon in between; the classification table gives that exon its bounds and its annotation
CAS <- EV[gene == "SRSF3" & event_type == "cassette exon" & FDR_prolif < 0.05][order(p_prolif)][1]
stopifnot(CAS$element_class == "NMD-only exon", CAS$element_length == 456)
EXS <- CAS$element_start; EXE <- CAS$element_end
DON <- 36598984L; INCL_END <- EXS - 1L; SKIP_END <- 36601151L; DOWN_START <- EXE + 1L
i_in <- which(ANN$start == DON & ANN$end == INCL_END)
i_sk <- which(ANN$start == DON & ANN$end == SKIP_END)
i_dn <- which(ANN$start == DOWN_START & ANN$end == SKIP_END)
stopifnot(length(i_in) == 1, length(i_sk) == 1, length(i_dn) == 1)
psi <- J[i_in, ] / (J[i_in, ] + J[i_sk, ])
hi <- MK$prolif >= quantile(MK$prolif, 2/3); lo <- MK$prolif <= quantile(MK$prolif, 1/3)
cat(sprintf("SRSF3 exon %s-%s (%d bp, %s): PSI %.3f in the slowest third, %.3f in the fastest\n",
            format(EXS, big.mark = ","), format(EXE, big.mark = ","), CAS$element_length,
            CAS$element_class, mean(psi[lo]), mean(psi[hi])))
cat(sprintf("  junction counts  include %.0f -> %.0f (upstream), %.0f -> %.0f (downstream)   skip %.0f -> %.0f\n",
            mean(J[i_in, lo]), mean(J[i_in, hi]), mean(J[i_dn, lo]), mean(J[i_dn, hi]),
            mean(J[i_sk, lo]), mean(J[i_sk, hi])))
rp <- cor.test(MK$prolif, psi, method = "spearman", exact = FALSE)
cat(sprintf("  PSI vs proliferation score: rho = %.3f, P = %.1e | model FDR = %.1e\n",
            rp$estimate, rp$p.value, CAS$FDR_prolif))

## the exon model is the protein-coding transcript; the alternative exon is drawn on its own,
## in red, because no protein-coding transcript carries it
CT <- EXG[gene == GENE & transcript_type == "protein_coding"]
canon <- CT[grepl("Ensembl_canonical", tags)]$transcript[1]
if (is.na(canon)) canon <- CT[, .N, by = transcript][order(-N)]$transcript[1]
EM <- CT[transcript == canon, .(s = start, e = end)][order(s)]
cat(sprintf("  exon model: %s, %d exons\n", canon, nrow(EM)))
X0 <- min(EM$s) - 250; X1 <- max(ANN$end) + 350
EM <- EM[s < X1]; EM[e > X1, e := X1]

JD <- data.table(start = ANN$start, end = ANN$end,
                 low = rowMeans(J[, lo, drop = FALSE]), high = rowMeans(J[, hi, drop = FALSE]))
JL <- rbind(data.table(JD[, .(start, end)], n = JD$low,  grp = "slowest third"),
            data.table(JD[, .(start, end)], n = JD$high, grp = "fastest third"))
JL[, grp := factor(grp, levels = c("slowest third", "fastest third"))]
JL[, kind := ifelse(start == DON & end == SKIP_END, "skip",
             ifelse((start == DON & end == INCL_END) | (start == DOWN_START & end == SKIP_END), "include", "other"))]
NOTE <- data.table(grp = factor("slowest third", levels = levels(JL$grp)))
pA <- ggplot() +
  geom_rect(data = EM, aes(xmin = s, xmax = e, ymin = -0.09, ymax = 0.09), fill = INK, colour = NA) +
  annotate("segment", x = X0 + 150, xend = X1 - 150, y = 0, yend = 0, colour = INK, linewidth = 0.3) +
  annotate("rect", xmin = EXS, xmax = EXE, ymin = -0.09, ymax = 0.09, fill = "#F6D6D2", colour = RED, linewidth = 0.4) +
  geom_curve(data = JL[kind == "other"], aes(x = start, xend = end, y = 0.1, yend = 0.1, linewidth = n),
             curvature = -0.30, colour = GREY_L, ncp = 12, lineend = "round") +
  geom_curve(data = JL[kind == "skip"], aes(x = start, xend = end, y = 0.1, yend = 0.1, linewidth = n),
             curvature = -0.42, colour = BLUE, ncp = 12, lineend = "round") +
  geom_curve(data = JL[kind == "include"], aes(x = start, xend = end, y = 0.1, yend = 0.1, linewidth = n),
             curvature = -0.22, colour = RED, ncp = 12, lineend = "round") +
  geom_text(data = JL[kind == "skip"], aes(x = (start + end)/2, y = 1.02, label = sprintf("%.0f", n)),
            family = FONT, size = pt(7), colour = BLUE, vjust = 0) +
  ## the skip arc passes above the include arcs everywhere they overlap and at 177 reads it is
  ## thick enough to swallow a count set above them, so the counts sit inside their own bows
  geom_label(data = JL[kind == "include"], aes(x = (start + end)/2, y = 0.235, label = sprintf("%.0f", n)),
             family = FONT, size = pt(7), colour = RED, vjust = 0.5, fill = "white", label.size = 0,
             label.padding = unit(0.6, "pt")) +
  scale_linewidth_continuous(range = c(0.2, 1.9), guide = "none") +
  facet_wrap(~ grp, ncol = 1) +
  geom_text(data = NOTE, aes(x = EXE + 120, y = -0.30), family = FONT, size = pt(7), colour = RED,
            hjust = 0, vjust = 0.5, lineheight = 0.95,
            label = sprintf("%d-bp exon, annotated only in\nnonsense-mediated-decay transcripts", CAS$element_length)) +
  geom_segment(data = NOTE, aes(x = EXE + 90, xend = EXE - 60, y = -0.30, yend = -0.12), colour = RED, linewidth = 0.3) +
  scale_x_continuous(sprintf("%s  (kb)", X$locus$chr), limits = c(X0, X1),
                     breaks = seq(36594000, 36602000, by = 1000), labels = function(v) sprintf("%.0f", v/1000)) +
  scale_y_continuous(NULL, limits = c(-0.46, 1.42), breaks = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      panel.spacing = unit(3, "pt"),
                      strip.text = element_text(size = 7.5, hjust = 0, margin = margin(b = 1)),
                      plot.margin = margin(9, 3, 3, 3))

## ====== b. every machinery event that moves with proliferation ===============
SIG <- EV[FDR_prolif < 0.05][order(gene, p_prolif)]
SIG[, what := ifelse(path == "inclusion", "included", "skipped")]
SIG[, kind := fifelse(event_type == "cassette exon", sprintf("%d-bp exon %s", element_length, what),
              fifelse(grepl("splice site", event_type), PRIME(sprintf("%s, %d bp %s", sub("alternative ", "alt. ", event_type), element_length, what)),
                      sprintf("%d-bp segment %s", element_length, what)))]
CLS <- c(`coding exon` = INK2, `NMD-only exon` = RED, `non-coding annotation only` = ORANGE,
         `not an annotated exon` = GREY_L)
CLSLAB <- c("coding exon", "exon only in NMD transcripts", "non-coding annotation only", "not an annotated exon")
SIG[, element_class := factor(element_class, levels = names(CLS))]
SIG[, row := factor(rev(seq_len(.N)), levels = seq_len(.N))]
SIG[event_type == "unannotated splice site", kind := sprintf("unannotated site, %.1f kb %s", element_length / 1000, what)]
XMAX <- max(abs(SIG$dPSI_pp)) * 1.05
pB <- ggplot(SIG, aes(dPSI_pp, row)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(aes(fill = element_class), width = 0.66, orientation = "y") +
  geom_text(aes(x = XMAX * 1.06, label = kind), hjust = 0, family = FONT, size = pt(6.8), colour = INK2) +
  scale_fill_manual(values = CLS, name = NULL, drop = FALSE, labels = CLSLAB) +
  scale_x_continuous("share of the measured path (%),\nfastest third \u2212 slowest third",
                     limits = c(-XMAX, XMAX * 3.2), breaks = seq(-15, 15, 5), labels = num_axis(0)) +
  scale_y_discrete(NULL, labels = setNames(SIG$gene, as.character(SIG$row)), expand = expansion(add = 0.6)) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7, face = "italic"),
                      axis.title.x = element_text(hjust = 0.16, lineheight = 1.0),
                      legend.position = "bottom", legend.text = element_text(size = 7),
                      legend.key.size = unit(7, "pt"), legend.margin = margin(t = -4),
                      ## room on the left for the longest gene symbol, which the grid otherwise cuts
                      plot.margin = margin(4, 3, 3, 9)) +
  guides(fill = guide_legend(nrow = 2))
cat(sprintf("panel b: %d events in %d genes at FDR < 0.05; %d cassette exons, %d in NMD-only exons\n",
            nrow(SIG), uniqueN(SIG$gene), sum(SIG$event_type == "cassette exon"), sum(SIG$element_class == "NMD-only exon")))

## ====== c. the named factors, donor by donor =================================
PD <- copy(NE$psi)
PD[, lab := sprintf("%s\n%d-bp exon %s", gene, element_length, ifelse(path == "inclusion", "included", "skipped"))]
PD[, lab := factor(lab, levels = unique(lab[order(FDR)]))]
PD[, element_class := factor(element_class, levels = names(CLS))]
ST <- PD[, { ct <- cor.test(prolif, psi, method = "spearman", exact = FALSE)
             .(txt = sprintf("%s = %s\n%s", RHO, num(ct$estimate), pfmt(ct$p.value))) }, by = lab]
pC <- ggplot(PD, aes(prolif, 100 * psi)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = INK, fill = GREY, alpha = 0.18, linewidth = 0.5) +
  geom_point(aes(colour = element_class), size = 0.8, alpha = 0.75) +
  geom_text(data = ST, aes(x = Inf, y = Inf, label = txt), hjust = 1.08, vjust = 1.2, family = FONT,
            size = pt(6.8), colour = INK, lineheight = 1.0) +
  facet_wrap(~ lab, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = CLS, guide = "none") +
  scale_x_continuous("proliferation score (20 markers)", labels = num_axis(0)) +
  scale_y_continuous("share of the measured path (%)", labels = num_axis(0), expand = expansion(mult = c(0.05, 0.45))) +
  theme_sa(8) + theme(strip.text = element_text(size = 6.9, hjust = 0, lineheight = 0.95, margin = margin(b = 1.5)),
                      panel.spacing.x = unit(5, "pt"), panel.spacing.y = unit(4, "pt"),
                      axis.text = element_text(size = 6.8), plot.margin = margin(4, 5, 3, 3))

## ====== d. what kinds of event ===============================================
S2 <- SM[genes %in% c("177 machinery genes", "all other genes")]
CP <- melt(S2[, .(genes, `cassette exon` = cassette, `alt. 5′ splice site` = alt5,
                  `alt. 3′ splice site` = alt3, `other / unannotated` = other)],
           id.vars = "genes", variable.name = "type", value.name = "n")
CP[, pct := 100 * n / sum(n), by = genes]
CP[, genes := factor(genes, levels = rev(c("177 machinery genes", "all other genes")))]
GL <- S2[, .(genes, txt = sprintf("%d events, %d (%.0f%%) in exons only in NMD transcripts", n_events, nmd_only, 100 * nmd_only / n_events))]
GL[, genes := factor(genes, levels = levels(CP$genes))]
TYPC <- c(`cassette exon` = BLUE, `alt. 5′ splice site` = TEAL, `alt. 3′ splice site` = PURPLE, `other / unannotated` = GREY_L)
pD <- ggplot(CP, aes(pct, genes, fill = type)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.3) +
  geom_text(data = CP[pct >= 7], aes(label = n), position = position_stack(vjust = 0.5),
            family = FONT, size = pt(7), colour = "white") +
  geom_text(data = GL, aes(x = 0, y = as.integer(genes) + 0.42, label = txt), inherit.aes = FALSE,
            hjust = 0, vjust = 0, family = FONT, size = pt(6.6), colour = INK2) +
  scale_fill_manual(values = TYPC, name = NULL) +
  scale_x_continuous("events at FDR < 0.05 against the proliferation score (%)", expand = c(0, 0),
                     breaks = seq(0, 100, 25)) +
  scale_y_discrete(NULL, expand = expansion(add = c(0.55, 0.95))) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7.2),
                      legend.position = "bottom", legend.text = element_text(size = 7),
                      legend.key.size = unit(7, "pt"), legend.margin = margin(t = -4),
                      plot.margin = margin(4, 8, 3, 3)) +
  guides(fill = guide_legend(nrow = 2))

## ====== e. the class-level answer, which is negative ==========================
CL <- R[!is.na(gene) & gene != ""]
CL[, cls := ifelse(machinery, "the 177 machinery genes", "all other genes")]
CL[, cls := factor(cls, levels = c("all other genes", "the 177 machinery genes"))]
MMc <- CL[cls == "the 177 machinery genes"]; bgk <- CL[cls == "all other genes"]
wp <- wilcox.test(abs(MMc$t_prolif), abs(bgk$t_prolif))$p.value
pE <- ggplot(CL, aes(abs(t_prolif), colour = cls)) +
  stat_ecdf(linewidth = 0.7) +
  scale_colour_manual(values = c(`all other genes` = GREY, `the 177 machinery genes` = BLUE), name = NULL) +
  ## the view stops at 8 but the distribution does not: coord_cartesian leaves the statistic whole
  scale_x_continuous("|t| for the junction against the proliferation score", labels = num_axis(0)) +
  coord_cartesian(xlim = c(0, 8)) +
  scale_y_continuous("cumulative fraction of junctions", labels = num_axis(1)) +
  annotate("text", x = 7.9, y = 0.06, hjust = 1, vjust = 0, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.05,
           label = sprintf("median |t| %s against %s\n%s at FDR < 0.05: %.0f%% against %.0f%%\n%s",
                           num(median(abs(MMc$t_prolif))), num(median(abs(bgk$t_prolif))),
                           "junctions", 100*mean(MMc$FDR_prolif < 0.05), 100*mean(bgk$FDR_prolif < 0.05),
                           pfmt(wp))) +
  annotate("text", x = 0.1, y = 0.99, hjust = 0, vjust = 1, family = FONT, size = pt(7),
           colour = INK2, label = sprintf("%d junctions beyond the axis, all in other genes",
                                          sum(abs(bgk$t_prolif) > 8))) +
  theme_sa(8) + theme(legend.position = "inside", legend.position.inside = c(0.62, 0.62),
                      legend.key.height = unit(8, "pt"), legend.text = element_text(size = 7),
                      plot.margin = margin(4, 3, 3, 3))
cat(sprintf("panel e: machinery %d junctions in %d genes, median |t| %.2f; background %d junctions, %.2f; Wilcoxon P = %.3f\n",
            nrow(MMc), length(unique(MMc$gene)), median(abs(MMc$t_prolif)), nrow(bgk),
            median(abs(bgk$t_prolif)), wp))

save_fig(lab_grid(lab_grid(pA, labels = "a", ncol = 1),
                  lab_grid(pB, pC, labels = c("b", "c"), ncol = 2, rel_widths = c(1.0, 1.0)),
                  lab_grid(pD, pE, labels = c("d", "e"), ncol = 2, rel_widths = c(0.92, 1.08)),
                  labels = c("", "", ""), ncol = 1, rel_heights = c(0.92, 1.12, 0.86)),
         "Fig3.png", 183, 196)
