## Fig. 3  What proliferation accounts for (six-figure restructure, 2026-09-17).
## Merges the parts of the former Fig. 3 (make_sa_fig5gene.R) and Fig. 4
## (make_sa_fig4new.R) that carry the main claim, and adds the random-set null
## as a panel of its own:
##   a  pre-mRNA processing score against the proliferation score, 107 adult donors
##   b  age effect per decade before and after adjusting for proliferation, under
##      the four cohort definitions
##   c  the same two models applied to 1,000 expression-matched random gene sets
##      (scripts/revision/run_reviewer_sensitivities.R, part C): share of the age
##      effect removed, and size of the unadjusted age effect
##   d  number of genes keeping a decline after adjustment under the four cohort
##      definitions, and the number common to all four
## The gene-level panels of the former Fig. 3 are now Fig. S2 (make_sa_figS_gene.R),
## the outcome index is Fig. S3 (make_sa_figS_outcome.R) and the age predictor is
## Fig. S5 (make_sa_figS_agepred.R).
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; CR <- file.path(D, "cohort_revised")
M    <- read.delim(file.path(CR, "primary/sample_metrics.tsv"))     # 107 normal adults
S    <- read.delim(file.path(CR, "sensitivity_metrics.tsv"))
SENS <- read.delim(file.path(CR, "sensitivity_genes.tsv"))
NL   <- read.delim(file.path(D, "reviewer_sensitivities/adjustment_negative_control.tsv"))
QC   <- read.delim(file.path(D, "reviewer_sensitivities/qc_covariates.tsv"))
NCORE <- length(readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt")))
EXPC <- "#B2182B"; OUTC <- "#2166AC"; BAR <- "#C3CBD4"
say <- function(...) cat(sprintf(...), "\n")

COHORTS <- c("published", "normal", "primary", "adult2082")
COHLAB  <- c(published = "deposited set\n(n = 142 samples)",
             normal    = "all normal donors\n(n = 133)",
             primary   = "normal adults 20+\n(primary, n = 107)",
             adult2082 = "normal adults 20–82\n(n = 76)")
stopifnot(all(S$n[match(COHORTS, S$cohort)] == c(142, 133, 107, 76)))

## =========== a. pre-mRNA processing score against proliferation ==============
ct <- cor.test(M$prolif, M$machinery)
say("a  r = %.4f [%.4f, %.4f], P = %.3g, n = %d donors", ct$estimate, ct$conf.int[1],
    ct$conf.int[2], ct$p.value, nrow(M))
pa <- ggplot(M, aes(prolif, machinery)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = EXPC, fill = EXPC,
              alpha = 0.16, linewidth = 0.6) +
  geom_point(size = 1.4, colour = EXPC, alpha = 0.85) +
  annotate("text", x = min(M$prolif), y = max(M$machinery), hjust = 0, vjust = 1, family = FONT,
           size = pt(7.5), colour = INK, lineheight = 1.05,
           label = rp(ct$estimate, ct$p.value, lab = "r")) +
  annotate("text", x = max(M$prolif), y = min(M$machinery), hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = INK2, label = sprintf("n = %d donors", nrow(M))) +
  scale_x_continuous(labels = num_axis()) + scale_y_continuous(labels = num_axis()) +
  labs(x = "proliferation score", y = "pre-mRNA processing score") +
  theme_sa() + theme(plot.margin = margin(8, 8, 3, 3))

## =========== b. age effect per decade, before and after adjustment ===========
K <- S[S$metric == "machinery", ]; K <- K[match(COHORTS, K$cohort), ]
K$row <- rev(seq_along(COHORTS))                          # published on top
FB <- rbind(
  data.frame(cohort = K$cohort, row = K$row, k = "unadjusted", b = K$beta, lo = K$lo, hi = K$hi, p = K$p),
  data.frame(cohort = K$cohort, row = K$row, k = "adjusted for proliferation", b = K$beta_adj,
             lo = K$lo_adj, hi = K$hi_adj, p = K$p_adj))
FB$k <- factor(FB$k, levels = c("unadjusted", "adjusted for proliferation"))
FB$y <- FB$row + ifelse(FB$k == "unadjusted", 0.17, -0.17)
FB$plab <- pfmt_v(FB$p)
for (i in seq_len(nrow(K)))
  say("b  %-10s n = %3d  unadjusted %.4f [%.4f, %.4f] P = %.3g | adjusted %.4f [%.4f, %.4f] P = %.3g | removed %.1f%%",
      K$cohort[i], K$n[i], K$beta[i], K$lo[i], K$hi[i], K$p[i], K$beta_adj[i], K$lo_adj[i],
      K$hi_adj[i], K$p_adj[i], K$pct_lost[i])
XP <- 0.012; XR <- 0.098                                  # columns for P and for % removed
pb <- ggplot(FB, aes(b, y, colour = k)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.55) +
  geom_point(size = 1.7) +
  geom_text(aes(x = XP, label = plab), hjust = 0, size = pt(7), family = FONT, show.legend = FALSE) +
  annotate("text", x = XR, y = K$row, hjust = 0, family = FONT, size = pt(7), colour = INK,
           label = sprintf("%.0f%%", K$pct_lost)) +
  annotate("text", x = XR, y = max(K$row) + 0.62, hjust = 0, vjust = 0, family = FONT, size = pt(7),
           colour = INK2, lineheight = 0.95, label = "removed by\nadjustment") +
  scale_colour_manual(values = c(unadjusted = EXPC, `adjusted for proliferation` = GREY), name = NULL) +
  scale_y_continuous(breaks = K$row, labels = COHLAB[K$cohort], limits = c(0.5, max(K$row) + 1.05),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.235, 0.15), breaks = c(-0.2, -0.1, 0), labels = num_axis()) +
  labs(x = "age effect per decade, pre-mRNA processing score", y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -6), legend.key.size = unit(6, "pt"),
                     axis.text.y = element_text(size = 7, lineheight = 0.95),
                     axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = margin(3, 4, 3, 6))

## =========== c. the same two models on 1,000 random gene sets ================
OBS <- K[K$cohort == "primary", ]
## the splicing set's own row was written by the same run; it must agree with Fig. 3b
q <- QC[QC$metric == "machinery" & QC$covariates == "cohort", ]
stopifnot(isTRUE(all.equal(q$beta, OBS$beta)), isTRUE(all.equal(q$pct_lost, OBS$pct_lost)))
sig <- is.finite(NL$p) & NL$p < 0.05
LS <- NL[sig, ]
say("c  random sets: %d; with an unadjusted age association at P < 0.05: %d (%d negative, %d positive)",
    nrow(NL), sum(sig), sum(LS$beta < 0), sum(LS$beta > 0))
say("c  %% removed among those %d: median %.1f, 2.5%% %.1f, 97.5%% %.1f, range %.1f to %.1f",
    sum(sig), median(LS$lost), quantile(LS$lost, .025), quantile(LS$lost, .975), min(LS$lost), max(LS$lost))
say("c  splicing set: %.1f%% removed; random sets with at least as much removed: %d of %d (%.3f)",
    OBS$pct_lost, sum(LS$lost >= OBS$pct_lost), nrow(LS), mean(LS$lost >= OBS$pct_lost))
say("c  unadjusted effect, random sets: %.4f to %.4f, median %.4f, largest |effect| %.4f; splicing set %.4f",
    min(NL$beta), max(NL$beta), median(NL$beta), max(abs(NL$beta)), OBS$beta)
say("c  random sets reaching |effect| >= %.4f: %d of %d; reaching unadjusted P <= %.3g: %d; adjusted P <= %.3g: %d",
    abs(OBS$beta), sum(abs(NL$beta) >= abs(OBS$beta)), nrow(NL), OBS$p, sum(NL$p <= OBS$p),
    OBS$p_adj, sum(NL$p_adj <= OBS$p_adj))
say("c  adjusted effect, random sets: %.4f to %.4f; splicing set %.4f; random sets keeping adjusted P < 0.05: %d of %d",
    min(NL$beta_adj), max(NL$beta_adj), OBS$beta_adj, sum(sig & NL$p_adj < 0.05), sum(sig))

h1 <- hist(LS$lost, breaks = seq(-20, 120, 10), plot = FALSE); top1 <- max(h1$counts)
pc1 <- ggplot(LS, aes(lost)) +
  geom_histogram(breaks = seq(-20, 120, 10), fill = BAR, colour = "white", linewidth = 0.2) +
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
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(-20, 120)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, top1 * 1.42)) +
  labs(x = "age effect removed by adjustment (%)",
       y = sprintf("random gene sets (n = %d of %s\nwith an age association)", nrow(LS),
                   format(nrow(NL), big.mark = ","))) +
  theme_sa() + theme(axis.title.y = element_text(lineheight = 0.95, margin = margin(r = 5)),
                     plot.margin = margin(8, 6, 3, 3))

h2 <- hist(NL$beta, breaks = seq(-0.02, 0.0325, 0.0025), plot = FALSE); top2 <- max(h2$counts)
pc2 <- ggplot(NL, aes(beta)) +
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
           label = sprintf("%d of %s random sets reach\nthe size of the splicing set's effect",
                           sum(abs(NL$beta) >= abs(OBS$beta)), format(nrow(NL), big.mark = ","))) +
  scale_x_continuous(limits = c(-0.19, 0.04), breaks = c(-0.15, -0.1, -0.05, 0), labels = num_axis()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, top2 * 1.42)) +
  labs(x = "unadjusted age effect per decade", y = "random gene sets (n = 1,000)") +
  theme_sa() + theme(axis.title.y = element_text(margin = margin(r = 5)),
                     plot.margin = margin(8, 6, 3, 6))
pc <- plot_grid(pc1, pc2, ncol = 2, rel_widths = c(1, 1.12))

## =========== d. genes keeping a decline, by cohort definition ================
KL <- setNames(lapply(COHORTS, function(k) readLines(file.path(CR, k, "genes_keeping_decline.txt"))), COHORTS)
KL <- lapply(KL, function(x) x[nzchar(x)])
stopifnot(all(lengths(KL) == SENS$decline_kept[match(COHORTS, SENS$cohort)]))
common <- Reduce(intersect, KL)
three <- Reduce(intersect, KL[c("published", "normal", "primary")])
say("d  genes keeping a decline: %s; common to all four: %d",
    paste(sprintf("%s %d", COHORTS, lengths(KL)), collapse = ", "), length(common))
say("d  (for the record) published & normal %d, published & primary %d, normal & primary %d, all three non-empty cohorts %d%s",
    length(intersect(KL$published, KL$normal)), length(intersect(KL$published, KL$primary)),
    length(intersect(KL$normal, KL$primary)), length(three),
    if (length(three)) paste0(": ", paste(sort(three), collapse = ", ")) else "")
say("d  age-associated genes before adjustment: %s",
    paste(sprintf("%s %d", SENS$cohort, SENS$age_assoc), collapse = ", "))
DB <- data.frame(lab = c(COHLAB[COHORTS], "in all four\ncohort definitions"),
                 n = c(lengths(KL), length(common)),
                 kind = c("other", "other", "primary", "other", "common"))
DB$lab <- factor(DB$lab, levels = rev(DB$lab))
pd <- ggplot(DB, aes(n, lab, fill = kind)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n, x = n + 1.6), hjust = 0, size = pt(7), family = FONT, colour = INK) +
  geom_hline(yintercept = 1.5, colour = GREY_L, linewidth = 0.3) +
  scale_fill_manual(values = c(primary = OUTC, other = BAR, common = INK2), guide = "none") +
  scale_x_continuous(limits = c(0, 66), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0))) +
  labs(x = sprintf("genes keeping a decline\nafter adjustment (of %d)", NCORE), y = NULL) +
  theme_sa() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.y = element_text(size = 7, lineheight = 0.95),
                     axis.title.x = element_text(hjust = 1, lineheight = 0.95),
                     plot.margin = margin(8, 8, 3, 6))

top <- lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(0.78, 1.22))
bot <- lab_grid(pc, pd, labels = c("c", "d"), ncol = 2, rel_widths = c(1.34, 0.66))
save_fig(lab_grid(top, bot, labels = c("", ""), ncol = 1, rel_heights = c(1, 0.95)), "Fig3.png", 183, 128)
