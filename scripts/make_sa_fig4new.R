FIGURES_WRITTEN <- c("Fig4.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

## Fig. 4. What adjusting for proliferation removes, and what the outcome measure
## does and does not show, in the redefined GSE113957 cohort.
## Direction (a): the age association of unannotated-junction use is reported as a
## negative result. It is present only in the published 142-sample cohort, which
## contains progeria donors, children and a 31-donor stratum aged 83 and over that
## comes entirely from one cell repository.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; CR <- file.path(D, "cohort_revised")
M  <- read.delim(file.path(CR, "primary/sample_metrics.tsv"))     # 107 normal adults
N  <- read.delim(file.path(CR, "normal/sample_metrics.tsv"))      # 133 normal donors
S  <- read.delim(file.path(CR, "sensitivity_metrics.tsv"))
NCORE <- length(readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt")))
AR <- arrow(length = unit(1.8, "pt"), type = "closed")
EXPC <- "#B2182B"; OUTC <- "#2166AC"
pv <- function(p) vapply(p, pfmt, character(1))      # pfmt() for a column of P values

## ================== a. the two measurements, same libraries =================
bx <- function(x, y, w, h, lab, fill = "white", col = INK, fc = INK, tsz = 7)
  data.frame(x, y, w, h, lab, fill, col, fc, tsz, stringsAsFactors = FALSE)
B <- rbind(
  bx(50, 91, 96, 16, sprintf("%d normal dermal fibroblast donors, 20 to 96 years\nuniformly reprocessed RNA-seq\n(9 progeria donors and 26 donors under 20 excluded)", nrow(M)),
     "#F2F4F6", INK, INK, 7),
  bx(25, 56, 46, 30, "", "#FBECEA", EXPC, EXPC, 7),
  bx(75, 56, 46, 30, "", "#EAF0F8", OUTC, OUTC, 7),
  bx(50, 23, 64, 11, "age effect per decade,\nbefore and after adjusting for proliferation", "white", INK, INK, 7))
B$xmin <- B$x - B$w/2; B$xmax <- B$x + B$w/2; B$ymin <- B$y - B$h/2; B$ymax <- B$y + B$h/2
p4a <- ggplot() +
  geom_segment(data = data.frame(x = c(50, 25, 75, 25, 75), xend = c(50, 25, 75, 25, 75),
                                 y = c(83, 79, 79, 41, 41), yend = c(79, 72.5, 72.5, 29, 29)),
               aes(x, y, xend = xend, yend = yend), colour = GREY, linewidth = 0.32) +
  geom_segment(data = data.frame(x = c(25, 75), xend = c(25, 75), y = c(73.5, 73.5), yend = c(71.5, 71.5)),
               aes(x, y, xend = xend, yend = yend), colour = GREY, linewidth = 0.32, arrow = AR) +
  geom_segment(data = data.frame(x = 25, xend = 75, y = 29, yend = 29),
               aes(x, y, xend = xend, yend = yend), colour = GREY, linewidth = 0.32) +
  geom_segment(data = data.frame(x = 50, xend = 50, y = 29, yend = 28.5),
               aes(x, y, xend = xend, yend = yend), colour = GREY, linewidth = 0.32, arrow = AR) +
  geom_rect(data = B, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = B$fill, colour = B$col, linewidth = 0.4) +
  geom_text(data = B, aes(x, y, label = lab), family = FONT, colour = B$fc, size = pt(B$tsz), lineheight = 1.0) +
  annotate("text", x = 25, y = 66.5, family = FONT, colour = EXPC, size = pt(7), lineheight = 1.05,
           fontface = "bold", label = "machinery") +
  annotate("text", x = 25, y = 61.5, family = FONT, colour = EXPC, size = pt(7), lineheight = 1.05,
           label = "how much spliceosome\nis transcribed") +
  annotate("text", x = 25, y = 53.5, family = FONT, colour = EXPC, size = pt(7), fontface = "italic",
           lineheight = 1.05, label = "SRSF1  SRSF2  HNRNPD\nHNRNPA1  HNRNPK  TRA2B") +
  annotate("text", x = 25, y = 45.5, family = FONT, colour = EXPC, size = pt(7),
           label = sprintf("%d genes, mean z", NCORE)) +
  annotate("text", x = 75, y = 66.5, family = FONT, colour = OUTC, size = pt(7), lineheight = 1.05,
           fontface = "bold", label = "outcome") +
  annotate("text", x = 75, y = 62.0, family = FONT, colour = OUTC, size = pt(7), lineheight = 1.05,
           label = "where the spliceosome\nactually cuts") +
  ## gene model: annotated junctions grey, an unannotated junction blue
  annotate("segment", x = 57, xend = 93, y = 51.0, yend = 51.0, colour = INK2, linewidth = 0.3) +
  annotate("rect", xmin = c(57, 71, 85), xmax = c(65, 79, 93), ymin = 49.8, ymax = 52.2,
           fill = "white", colour = INK2, linewidth = 0.35) +
  annotate("curve", x = 65, xend = 71, y = 52.2, yend = 52.2, curvature = -0.75, ncp = 12,
           colour = GREY, linewidth = 0.4) +
  annotate("curve", x = 79, xend = 85, y = 52.2, yend = 52.2, curvature = -0.75, ncp = 12,
           colour = GREY, linewidth = 0.4) +
  annotate("curve", x = 65, xend = 89, y = 49.8, yend = 49.8, curvature = 0.16, ncp = 16,
           colour = OUTC, linewidth = 0.55) +
  annotate("segment", x = 89, xend = 89, y = 49.8, yend = 52.2, colour = OUTC, linewidth = 0.45) +
  annotate("text", x = 75, y = 56.0, family = FONT, colour = GREY, size = pt(7), label = "annotated") +
  annotate("text", x = 75, y = 44.0, family = FONT, colour = OUTC, size = pt(7), lineheight = 1.05,
           label = "unannotated reads /\nall junction reads") +
  scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(15, 100)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

## =========== b. machinery expression against proliferation =================
ct <- cor.test(M$prolif, M$machinery)
p4b <- ggplot(M, aes(prolif, machinery)) +
  geom_smooth(method = "lm", se = TRUE, colour = EXPC, fill = EXPC, alpha = 0.16, linewidth = 0.6) +
  geom_point(size = 1.5, colour = EXPC, alpha = 0.85) +
  annotate("text", x = min(M$prolif), y = max(M$machinery), hjust = 0, vjust = 1, family = FONT,
           size = pt(7.5), colour = INK, lineheight = 1.05,
           label = rp(ct$estimate, ct$p.value, lab = "r")) +
  annotate("text", x = max(M$prolif), y = min(M$machinery), hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = GREY, label = sprintf("n = %d donors", nrow(M))) +
  scale_x_continuous(labels = num_axis()) + scale_y_continuous(labels = num_axis()) +
  labs(x = "proliferation score", y = "pre-mRNA processing score") + theme_sa()

## =========== c, d. age effect per decade, before and after adjustment ======
COH <- c(published = "published, 142", normal = "all normal, 133",
         primary = "primary: adults 20+, 107", adult2082 = "adults 20 to 82, 76")
forest <- function(metric, scale = 1, xlab, col, note = NULL) {
  K <- S[S$metric == metric, ]
  K <- K[match(names(COH), K$cohort), ]
  d <- rbind(
    data.frame(cohort = COH, k = "unadjusted", b = K$beta * scale, lo = K$lo * scale, hi = K$hi * scale, p = K$p),
    data.frame(cohort = COH, k = "proliferation-adjusted", b = K$beta_adj * scale,
               lo = K$lo_adj * scale, hi = K$hi_adj * scale, p = K$p_adj))
  d$cohort <- factor(d$cohort, levels = rev(COH))
  d$k <- factor(d$k, levels = c("unadjusted", "proliferation-adjusted"))
  rng <- range(c(d$lo, d$hi)); pad <- diff(rng) * 0.42
  ggplot(d, aes(b, cohort, colour = k)) +
    geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.55,
                   position = position_dodge(width = 0.62)) +
    geom_point(size = 1.7, position = position_dodge(width = 0.62)) +
    geom_text(aes(x = hi, label = pv(p)), position = position_dodge(width = 0.62),
              hjust = -0.18, size = pt(7), family = FONT, show.legend = FALSE) +
    scale_colour_manual(values = setNames(c(col, GREY), levels(d$k)), name = NULL) +
    scale_x_continuous(limits = c(rng[1] - pad * 0.15, rng[2] + pad), labels = num_axis()) +
    labs(x = xlab, y = NULL) +
    theme_sa() + theme(legend.position = "top", legend.justification = "left",
                       legend.margin = margin(b = -5), legend.key.size = unit(6, "pt"),
                       axis.text.y = element_text(size = 7.5)) }
p4c <- forest("machinery", 1, "age effect per decade,\npre-mRNA processing score", EXPC)
p4d <- forest("unannot_reads", 1e4,
              sprintf("age effect per decade,\nunannotated read fraction (%s 10%s)", TIMES, sup("-4")), OUTC)

## =========== e. the outcome measure against donor age, 83+ marked ==========
N$grp <- ifelse(N$age >= 83, "aged 83+ (all from one repository)",
                ifelse(N$repo == "AG", "under 83, repository AG", "under 83, other repository"))
ctN <- cor.test(N$age, N$unannot_reads, method = "spearman", exact = FALSE)
ctY <- cor.test(N$age[N$age < 83], N$unannot_reads[N$age < 83], method = "spearman", exact = FALSE)
p4e <- ggplot(N, aes(age, unannot_reads * 1e2)) +
  geom_point(aes(colour = grp, shape = grp), size = 1.5, alpha = 0.9) +
  geom_smooth(data = N[N$age < 83, ], method = "lm", se = FALSE, colour = INK2, linewidth = 0.5,
              linetype = "22") +
  scale_colour_manual(values = setNames(c(ORANGE, INK2, GREY), sort(unique(N$grp))), name = NULL) +
  scale_shape_manual(values = setNames(c(17, 16, 1), sort(unique(N$grp))), name = NULL) +
  annotate("text", x = 0, y = max(N$unannot_reads * 1e2) * 1.045, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.2,
           label = sprintf("all normal donors\n   %s = %s, %s\nunder 83 years\n   %s = %s, %s",
                           RHO, num(ctN$estimate), pfmt(ctN$p.value),
                           RHO, num(ctY$estimate), pfmt(ctY$p.value))) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = "donor age (years)", y = "unannotated junction reads (%)") +
  theme_sa() + theme(legend.position = "bottom", legend.direction = "vertical",
                     legend.key.size = unit(6, "pt"), legend.margin = margin(t = -4),
                     legend.spacing.y = unit(0, "pt"))

## ---- f. a transcriptome-wide age predictor built in the same donors ---------
AP <- read.delim("public_data_tierA/derived/age_predictor/age_predictor_per_donor.tsv")
CVT <- read.delim("public_data_tierA/derived/age_predictor/age_predictor_cv.tsv")
PN <- read.delim("public_data_tierA/derived/age_predictor/age_predictor_permutation_null.tsv")
rfull <- CVT$r[CVT$model == "all genes"]; mfull <- CVT$MAE[CVT$model == "all genes"]
rpro  <- CVT$r[CVT$model == "proliferation score alone"]; mpro <- CVT$MAE[CVT$model == "proliferation score alone"]
rpp   <- cor(AP$pred_full, AP$prolif)
p4f <- ggplot(AP, aes(age, pred_full, colour = prolif)) +
  geom_abline(slope = 1, intercept = 0, colour = GREY, linewidth = 0.3, linetype = "22") +
  geom_point(size = 1.6, alpha = 0.95) +
  scale_colour_gradient2(low = BLUE, mid = GREY_L, high = RED, midpoint = 0,
                         name = "proliferation\nscore", guide = guide_colourbar(barheight = unit(18, "pt"),
                         barwidth = unit(5, "pt"), title.position = "top")) +
  annotate("text", x = 20, y = max(AP$pred_full) * 1.02, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.2,
           label = sprintf("cross-validated %s = %s, MAE %.0f years\nproliferation score alone: %s = %s, MAE %.0f\npredicted age vs proliferation: %s = %s",
                           "r", num(rfull), mfull, "r", num(rpro), mpro, "r", num(rpp))) +
  annotate("text", x = 96, y = min(AP$pred_full), hjust = 1, vjust = 0, family = FONT, size = pt(7),
           colour = GREY, label = sprintf("permuted age: %s = %s to %s", "r", num(min(PN$r)), num(max(PN$r)))) +
  labs(x = "donor age (years)", y = "predicted age (years)") +
  theme_sa() + theme(legend.position = "right", legend.title = element_text(size = 7))

top <- lab_grid(p4a, p4b, p4e, labels = c("a", "b", "c"), ncol = 3, rel_widths = c(1.30, 1, 1.05))
mid <- lab_grid(p4c, p4d, labels = c("d", "e"), ncol = 2, rel_widths = c(1, 1.06))
save_fig(lab_grid(top, mid, lab_grid(p4f, labels = c("f"), ncol = 1), labels = c("", "", ""),
                  ncol = 1, rel_heights = c(1, 0.82, 0.92)),
         "Fig4.png", 183, 178)
