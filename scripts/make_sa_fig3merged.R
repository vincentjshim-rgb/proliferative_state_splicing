## Fig. 4  What proliferation accounts for (six-figure restructure, 2026-09-17).
## Revised 2026-09-20 (user request): the figure had six plots and proliferation itself was
## hard to see in any of them. It now has three, each about proliferation:
##   a  pre-mRNA processing score against the proliferation score, 107 adult donors
##   b  the age decline itself, with the donors split by proliferation: the pooled slope
##      against the slope within each third, and the model estimate with proliferation held
##   c  the size of the age association against 1,000 expression-matched random gene sets
##   d  every gene: its age effect in the donor cohort against its coupling to the counted rate
## The cohort-definition bars, the share-removed histogram and the count of genes keeping a
## decline moved to Supplementary Fig. S2 (make_sa_figS_gene.R, panels d-f).
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; CR <- file.path(D, "cohort_revised")
M    <- read.delim(file.path(CR, "primary/sample_metrics.tsv"))     # 107 normal adults
S    <- read.delim(file.path(CR, "sensitivity_metrics.tsv"))
NL   <- read.delim(file.path(D, "reviewer_sensitivities/adjustment_negative_control.tsv"))
QC   <- read.delim(file.path(D, "reviewer_sensitivities/qc_covariates.tsv"))
MED  <- read.delim(file.path(D, "chain/mediation.tsv"))
EXPC <- "#B2182B"; BAR <- "#C3CBD4"
say <- function(...) cat(sprintf(...), "\n")
stopifnot(nrow(M) == 107)

## =========== a. pre-mRNA processing score against proliferation ==============
ct <- cor.test(M$prolif, M$machinery)
say("a  r = %.4f [%.4f, %.4f], P = %.3g, n = %d donors", ct$estimate, ct$conf.int[1],
    ct$conf.int[2], ct$p.value, nrow(M))
pa <- ggplot(M, aes(prolif, machinery)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = EXPC, fill = EXPC,
              alpha = 0.16, linewidth = 0.6) +
  geom_point(size = 1.3, colour = EXPC, alpha = 0.85) +
  annotate("text", x = min(M$prolif), y = max(M$machinery), hjust = 0, vjust = 1, family = FONT,
           size = pt(7.5), colour = INK, lineheight = 1.05,
           label = rp(ct$estimate, ct$p.value, lab = "r")) +
  annotate("text", x = max(M$prolif), y = min(M$machinery), hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = INK2, label = sprintf("n = %d donors", nrow(M))) +
  scale_x_continuous(labels = num_axis()) + scale_y_continuous(labels = num_axis()) +
  labs(x = "proliferation score (20 markers)", y = "pre-mRNA processing score") +
  theme_sa() + theme(plot.margin = margin(8, 10, 3, 3))

## =========== b. the age decline, with proliferation held ======================
K <- S[S$metric == "machinery" & S$cohort == "primary", ]
stopifnot(nrow(K) == 1)
q <- QC[QC$metric == "machinery" & QC$covariates == "cohort", ]
stopifnot(isTRUE(all.equal(q$beta, K$beta)), isTRUE(all.equal(q$pct_lost, K$pct_lost)))
M$third <- cut(M$prolif, quantile(M$prolif, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
               labels = c("slowest third", "middle third", "fastest third"))
TS <- do.call(rbind, lapply(levels(M$third), function(k) { s <- M[M$third == k, ]
  f <- lm(machinery ~ I(age/10), s)
  data.frame(third = k, n = nrow(s), slope = coef(f)[2], mean_prolif = mean(s$prolif)) }))
pooled <- lm(machinery ~ I(age/10), M)
say("b  pooled raw slope %.4f per decade; cohort model %.4f -> %.4f with proliferation (%.0f%% removed)",
    coef(pooled)[2], K$beta, K$beta_adj, K$pct_lost)
for (i in seq_len(nrow(TS))) say("b  %-14s n = %2d  slope %.4f per decade  mean proliferation score %+.2f",
                                 TS$third[i], TS$n[i], TS$slope[i], TS$mean_prolif[i])
prop <- MED[MED$quantity == "prop", ]
say("b  mediation: proliferation carries %.0f%% [%.0f, %.0f] of the age effect", 100*prop$estimate, 100*prop$lo, 100*prop$hi)
TCOL <- c(`slowest third` = "#8C6D4F", `middle third` = "#C9A27C", `fastest third` = "#2C6FAF")
## each third's line is labelled at its right-hand end instead of in a key, which did not fit
nd <- data.frame(age = seq(min(M$age), max(M$age), length.out = 50))
nd$fit <- predict(pooled, data.frame(age = nd$age))
END <- do.call(rbind, lapply(levels(M$third), function(k) { s <- M[M$third == k, ]
  f <- lm(machinery ~ age, s)
  sl <- TS$slope[TS$third == k]; if (abs(sl) < 5e-4) sl <- 0          # "-0.000" is not a number
  data.frame(third = k, x = max(M$age) + 1.5, y = predict(f, data.frame(age = max(M$age))),
             lab = sprintf("%s\n%s per decade", k, num(sl, 3))) }))
## the middle and slowest labels would sit two lines apart; spread them a little
END$y <- END$y + c(-0.04, 0.10, 0)[match(END$third, c("slowest third", "middle third", "fastest third"))]
pb <- ggplot(M, aes(age, machinery)) +
  geom_line(data = nd, aes(age, fit), inherit.aes = FALSE, colour = INK, linewidth = 0.55, linetype = "22") +
  geom_smooth(aes(colour = third), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
  geom_point(aes(colour = third), size = 1.3, alpha = 0.8) +
  geom_text(data = END, aes(x, y, label = lab, colour = third), hjust = 0, vjust = 0.5, family = FONT,
            size = pt(7), lineheight = 0.95) +
  annotate("text", x = min(M$age), y = min(M$machinery), hjust = 0, vjust = 0, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.08,
           label = sprintf("all donors (dashed): %s per decade\nproliferation in the model: %s per decade (%.0f%% removed)\nmediation: %.0f%% [%.0f, %.0f] of the age effect",
                           num(K$beta, 3), num(K$beta_adj, 3), K$pct_lost, 100*prop$estimate, 100*prop$lo, 100*prop$hi)) +
  scale_colour_manual(values = TCOL, guide = "none") +
  scale_x_continuous(labels = num_axis(0), breaks = seq(20, 100, 20), limits = c(min(M$age), max(M$age) + 31)) +
  scale_y_continuous(labels = num_axis(), expand = expansion(mult = c(0.30, 0.06))) +
  labs(x = "donor age (years)", y = "pre-mRNA processing score",
       subtitle = "donors split into thirds by proliferation score") +
  theme_sa() + theme(plot.subtitle = element_text(size = 7.5, colour = INK2, hjust = 0, margin = margin(b = 2)),
                     plot.margin = margin(3, 6, 3, 3))

## =========== c. the size of the association against 1,000 random sets =======
OBS <- K
sig <- is.finite(NL$p) & NL$p < 0.05
LS <- NL[sig, ]
say("c  random sets: %d; with an unadjusted age association at P < 0.05: %d (%d negative, %d positive)",
    nrow(NL), sum(sig), sum(LS$beta < 0), sum(LS$beta > 0))
say("c  %% removed among those %d: median %.1f, 2.5%% %.1f, 97.5%% %.1f", sum(sig), median(LS$lost),
    quantile(LS$lost, .025), quantile(LS$lost, .975))
say("c  unadjusted effect, random sets: %.4f to %.4f, largest |effect| %.4f; splicing set %.4f; random sets reaching it: %d",
    min(NL$beta), max(NL$beta), max(abs(NL$beta)), OBS$beta, sum(abs(NL$beta) >= abs(OBS$beta)))
h2 <- hist(NL$beta, breaks = seq(-0.02, 0.0325, 0.0025), plot = FALSE); top2 <- max(h2$counts)
pc <- ggplot(NL, aes(beta)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
  geom_histogram(breaks = seq(-0.02, 0.0325, 0.0025), fill = BAR, colour = NA) +
  annotate("segment", x = OBS$lo, xend = OBS$hi, y = top2 * 0.30, yend = top2 * 0.30,
           colour = EXPC, linewidth = 0.55) +
  annotate("point", x = OBS$beta, y = top2 * 0.30, colour = EXPC, size = 1.7) +
  annotate("text", x = OBS$beta, y = top2 * 0.38, hjust = 0.5, vjust = 0, family = FONT,
           size = pt(7), colour = EXPC, lineheight = 0.95,
           label = sprintf("splicing set\n%s", num(OBS$beta, 3))) +
  annotate("text", x = -0.012, y = top2 * 1.12, hjust = 1, vjust = 1, family = FONT, size = pt(7),
           colour = INK2, lineheight = 0.95,
           label = sprintf("%s random sets\n%s to %s", format(nrow(NL), big.mark = ","),
                           num(min(NL$beta), 3), num(max(NL$beta), 3))) +
  annotate("text", x = -0.19, y = top2 * 1.40, hjust = 0, vjust = 1, family = FONT, size = pt(7),
           colour = INK, lineheight = 0.95,
           label = sprintf("%d of %s random sets\nreach the splicing set's effect",
                           sum(abs(NL$beta) >= abs(OBS$beta)), format(nrow(NL), big.mark = ","))) +
  scale_x_continuous(breaks = c(-0.15, -0.1, -0.05, 0), labels = num_axis()) +
  coord_cartesian(xlim = c(-0.19, 0.04)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, top2 * 1.42)) +
  labs(x = "unadjusted age effect per decade", y = "random gene sets (n = 1,000)") +
  theme_sa() + theme(axis.title.y = element_text(margin = margin(r = 5)),
                     plot.margin = margin(8, 6, 3, 6))

## =========== d. every gene: age effect against coupling to the counted rate ============
## post hoc (2026-09-21). One point per gene measured in both datasets: the per-decade age
## effect in the donor cohort (limma, cohort covariates) against the gene's correlation with
## the counted division rate in the Cellular Lifespan Study.
DE <- read.delim(file.path(D, "revision_stats/GSE113957_primary_cohort_per_decade_DE.tsv"))
GV <- read.delim(file.path(D, "gene_level/gene_vs_counted_rate.tsv"))
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
GMTD <- tempfile("reactome"); dir.create(GMTD)
unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = GMTD)
GMT <- strsplit(readLines(list.files(GMTD, pattern = "\\.gmt$", full.names = TRUE)[1]), "\t")
CCG <- unique(unlist(lapply(GMT[grepl("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", vapply(GMT, `[`, "", 1))],
                            function(v) v[-(1:2)])))
DE$gene <- toupper(DE$gene)
J <- merge(DE[, c("gene", "logFC", "FDR")], GV[, c("gene", "rho")], by = "gene")
J$grp <- ifelse(J$gene %in% CORE, "177 splicing genes", ifelse(J$gene %in% CCG, "cell-cycle union", "other"))
J$grp <- factor(J$grp, levels = c("other", "cell-cycle union", "177 splicing genes"))
J <- J[order(J$grp), ]
ct4 <- cor.test(J$rho, J$logFC, method = "spearman", exact = FALSE)
in177 <- J$grp == "177 splicing genes"
ct4s <- cor.test(J$rho[in177], J$logFC[in177], method = "spearman", exact = FALSE)
say("d  %d genes in both datasets: rho(coupling, age effect) = %.3f (P = %.2g); within the 177: %.3f (P = %.2g)",
    nrow(J), ct4$estimate, ct4$p.value, ct4s$estimate, ct4s$p.value)
say("d  splicing genes: %d of %d fall with age (FDR < 0.05); cell-cycle union %d of %d",
    sum(in177 & J$FDR < 0.05 & J$logFC < 0), sum(in177), sum(J$grp == "cell-cycle union" & J$FDR < 0.05 & J$logFC < 0),
    sum(J$grp == "cell-cycle union"))
GCOL <- c(other = "#D5DAE0", `cell-cycle union` = EXPC, `177 splicing genes` = "#2C6FAF")
SHOW <- c("SRSF1", "SRSF2", "SRSF3", "HNRNPK")
LB2 <- J[J$gene %in% SHOW, ]; LB2$ny <- c(SRSF1 = 0.055, SRSF2 = 0.085, SRSF3 = -0.075, HNRNPK = -0.105)[LB2$gene]
pd <- ggplot(J, aes(rho, logFC, colour = grp)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) + geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_point(aes(size = grp, alpha = grp)) +
  geom_smooth(data = J, aes(rho, logFC), inherit.aes = FALSE, method = "lm", formula = y ~ x, colour = INK,
              linewidth = 0.5, se = FALSE, linetype = "22") +
  geom_segment(data = LB2, aes(x = rho, xend = rho, y = logFC, yend = logFC + ny * 0.85), colour = GREY, linewidth = 0.25) +
  geom_text(data = LB2, aes(y = logFC + ny, label = gene), family = FONT, size = pt(7), fontface = "italic",
            colour = INK, hjust = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = GCOL, name = NULL) +
  scale_size_manual(values = c(other = 0.45, `cell-cycle union` = 0.9, `177 splicing genes` = 1.1), guide = "none") +
  scale_alpha_manual(values = c(other = 0.5, `cell-cycle union` = 0.85, `177 splicing genes` = 0.95), guide = "none") +
  annotate("text", x = min(J$rho), y = max(J$logFC), hjust = 0, vjust = 1, family = FONT, size = pt(7), colour = INK,
           lineheight = 1.05, label = sprintf("all genes: %s = %s\n%s\nn = %s genes in both datasets", RHO, num(ct4$estimate),
                                              pfmt(ct4$p.value), format(nrow(J), big.mark = ","))) +
  scale_x_continuous(sprintf("%s with the counted division rate (328 libraries)", RHO), labels = num_axis()) +
  scale_y_continuous("age effect in the donor cohort (log\u2082 CPM per decade)", labels = num_axis()) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), nrow = 1)) +
  theme_sa() + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(b = -4),
                     legend.key.size = unit(7, "pt"), legend.text = element_text(size = 7),
                     plot.margin = margin(3, 6, 3, 3))

top <- lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(0.82, 1.32))
bot <- lab_grid(pc, pd, labels = c("c", "d"), ncol = 2, rel_widths = c(0.86, 1.28))
save_fig(lab_grid(top, bot, labels = c("", ""), ncol = 1, rel_heights = c(1, 1.02)), "Fig4.png", 183, 150)
