## Supplementary Fig. S5  A transcriptome-wide age predictor built in the same donors
## (six-figure restructure, 2026-09-17). Panel f of the former main Fig. 4
## (make_sa_fig4new.R), extended so that every arm of scripts/revision/run_age_predictor.R
## is shown with its own value:
##   a  cross-validated predicted age against donor age, coloured by proliferation
##   b  the same predictions against the proliferation score
##   c  cross-validated accuracy (r, 95% CI, MAE) of every arm, with the permuted-age
##      null. The arm with proliferation regressed out of every gene gives r = -0.53,
##      not zero: its predictions are almost constant (they span about 1.5 years), and
##      a constant leave-out prediction is anticorrelated with the held-out value.
##      It is drawn with that value and with the spread of its predictions.
## Post hoc analysis; 107 normal adult donors of GSE113957.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
O   <- "public_data_tierA/derived/age_predictor"
AP  <- read.delim(file.path(O, "age_predictor_per_donor.tsv"))
CVT <- read.delim(file.path(O, "age_predictor_cv.tsv"))
PN  <- read.delim(file.path(O, "age_predictor_permutation_null.tsv"))
say <- function(...) cat(sprintf(...), "\n")
cv <- function(m) CVT[CVT$model == m, ]
FULL <- cv("all genes"); PRO <- cv("proliferation score alone")
COVR <- cv("all genes, cohort covariates removed"); NOP <- cv("all genes, proliferation removed")
## the per-donor table and the summary table come from one run and must agree
stopifnot(isTRUE(all.equal(cor(AP$pred_full, AP$age), FULL$r)),
          isTRUE(all.equal(cor(AP$pred_no_prolif, AP$age), NOP$r)),
          isTRUE(all.equal(cor(AP$pred_full, AP$prolif), FULL$r_with_proliferation)))
for (i in seq_len(nrow(CVT)))
  say("cv  %-40s r = %.4f [%.4f, %.4f], P = %.3g, MAE %.2f years, r with proliferation %.4f",
      CVT$model[i], CVT$r[i], CVT$lo[i], CVT$hi[i], CVT$p[i], CVT$MAE[i], CVT$r_with_proliferation[i])
say("null  permuted age, %d draws: r %.4f to %.4f, mean %.4f", nrow(PN), min(PN$r), max(PN$r), mean(PN$r))
MAE0 <- mean(abs(AP$age - mean(AP$age)))
say("ref  predicting the cohort mean age (%.1f years) for every donor: MAE %.2f years", mean(AP$age), MAE0)
say("spread of predictions (years): all genes %.2f to %.2f (SD %.2f); proliferation removed %.2f to %.2f (SD %.2f)",
    min(AP$pred_full), max(AP$pred_full), sd(AP$pred_full),
    min(AP$pred_no_prolif), max(AP$pred_no_prolif), sd(AP$pred_no_prolif))

## ---- a. predicted against chronological age ---------------------------------
pa <- ggplot(AP, aes(age, pred_full, colour = prolif)) +
  geom_abline(slope = 1, intercept = 0, colour = GREY, linewidth = 0.3, linetype = "22") +
  geom_point(size = 1.6, alpha = 0.95) +
  scale_colour_gradient2(low = BLUE, mid = GREY_L, high = RED, midpoint = 0, labels = num_axis(),
                         breaks = c(-3, -2, -1, 0, 1), name = "proliferation\nscore",
                         guide = guide_colourbar(barheight = unit(52, "pt"), barwidth = unit(5, "pt"),
                                                 title.position = "top", ticks.colour = "white")) +
  annotate("text", x = 20, y = 100, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.15,
           label = sprintf("r = %s, %s\nmean absolute error %.0f years\nn = %d donors",
                           num(FULL$r), pfmt(FULL$p), FULL$MAE, nrow(AP))) +
  annotate("text", x = 96, y = 27, hjust = 1, vjust = 0, family = FONT, size = pt(7), colour = GREY,
           label = "dashed line, predicted = actual") +
  scale_x_continuous(limits = c(18, 98), breaks = seq(20, 100, 20)) +
  scale_y_continuous(limits = c(25, 100), breaks = seq(40, 100, 20)) +
  labs(x = "donor age (years)", y = "predicted age (years),\nall genes, cross-validated") +
  theme_sa() + theme(legend.position = "right", legend.title = element_text(size = 7, lineheight = 0.95),
                     legend.text = element_text(size = 7),
                     axis.title.y = element_text(lineheight = 0.95), plot.margin = margin(10, 4, 3, 3))

## ---- b. the predictions against the proliferation score ---------------------
ctb <- cor.test(AP$pred_full, AP$prolif)
r2b <- summary(lm(pred_full ~ prolif, AP))$r.squared
say("b  predicted age vs proliferation score: r = %.4f [%.4f, %.4f], P = %.3g, n = %d donors; R2 = %.4f",
    ctb$estimate, ctb$conf.int[1], ctb$conf.int[2], ctb$p.value, nrow(AP), r2b)
cta <- cor.test(AP$age, AP$prolif)
say("b  (for the record) donor age vs proliferation score: r = %.4f, P = %.3g; Spearman rho = %.4f",
    cta$estimate, cta$p.value, cor(AP$age, AP$prolif, method = "spearman"))
pb <- ggplot(AP, aes(prolif, pred_full)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK2, fill = INK2, alpha = 0.14,
              linewidth = 0.55) +
  geom_point(size = 1.5, colour = INK2, alpha = 0.85) +
  annotate("text", x = max(AP$prolif), y = 100, hjust = 1, vjust = 1, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.15,
           label = sprintf("r = %s, %s\nn = %d donors", num(ctb$estimate), pfmt(ctb$p.value), nrow(AP))) +
  scale_x_continuous(labels = num_axis()) +
  scale_y_continuous(limits = c(25, 100), breaks = seq(40, 100, 20)) +
  labs(x = "proliferation score", y = "predicted age (years)") +
  theme_sa() + theme(plot.margin = margin(10, 8, 3, 6))

## ---- c. every arm with its own value ----------------------------------------
ARM <- data.frame(
  lab = c("all expressed genes",
          "proliferation score alone\n(20 genes)",
          "all genes, cohort covariates\nremoved",
          "all genes, cohort covariates\nand proliferation removed"),
  rbind(FULL, PRO, COVR, NOP)[, c("r", "lo", "hi", "p", "MAE")])
ARM$row <- 6:3
NULLROW <- 2; REFROW <- 1
XM <- 1.04                                             # column for the mean absolute error
pc <- ggplot(ARM) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
  annotate("rect", xmin = min(PN$r), xmax = max(PN$r), ymin = 1.5, ymax = 6.55, fill = GREY_L, alpha = 0.28) +
  geom_segment(aes(x = lo, xend = hi, y = row, yend = row), colour = INK, linewidth = 0.55) +
  geom_point(aes(r, row), colour = INK, size = 1.8) +
  geom_text(aes(x = hi + 0.025, y = row, label = sprintf("r = %s", num(r))), hjust = 0, family = FONT,
            size = pt(7), colour = INK) +
  annotate("point", x = PN$r, y = NULLROW, shape = 1, size = 1.4, stroke = 0.4, colour = GREY) +
  annotate("text", x = max(PN$r) + 0.025, y = NULLROW, hjust = 0, family = FONT, size = pt(7), colour = INK2,
           label = sprintf("r = %s to %s", num(min(PN$r)), num(max(PN$r)))) +
  annotate("text", x = NOP$hi + 0.20, y = 3, hjust = 0, family = FONT, size = pt(7), colour = INK2,
           lineheight = 0.95,
           label = sprintf("predicted ages span only\n%.1f to %.1f years", min(AP$pred_no_prolif),
                           max(AP$pred_no_prolif))) +
  annotate("segment", x = -0.7, xend = 1.17, y = 1.5, yend = 1.5, colour = GREY_L, linewidth = 0.3) +
  annotate("text", x = XM, y = 7.05, hjust = 0.5, vjust = 0.5, family = FONT, size = pt(7), colour = INK2,
           lineheight = 0.95, label = "mean absolute\nerror (years)") +
  annotate("text", x = XM, y = ARM$row, hjust = 0.5, family = FONT, size = pt(7), colour = INK,
           label = sprintf("%.1f", ARM$MAE)) +
  annotate("text", x = XM, y = REFROW, hjust = 0.5, family = FONT, size = pt(7), colour = INK,
           label = sprintf("%.1f", MAE0)) +
  scale_y_continuous(breaks = c(ARM$row, NULLROW, REFROW),
                     labels = c(ARM$lab, sprintf("all genes, permuted ages\n(%d draws)", nrow(PN)),
                                "no model: the cohort's mean\nage for every donor"),
                     limits = c(0.45, 7.55), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.7, 1.17), breaks = seq(-0.5, 0.75, 0.25), labels = num_axis()) +
  labs(x = sprintf("cross-validated correlation between predicted and donor age (r, 95%% CI; n = %d donors)", nrow(AP)),
       y = NULL) +
  theme_sa() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.y = element_text(size = 7, lineheight = 0.95),
                     axis.title.x = element_text(hjust = 0.3),
                     plot.margin = margin(10, 4, 3, 3))

top <- lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(1.22, 0.78))
save_fig(lab_grid(top, pc, labels = c("", "c"), ncol = 1, rel_heights = c(1, 0.95)), "FigS5.png", 183, 132)
