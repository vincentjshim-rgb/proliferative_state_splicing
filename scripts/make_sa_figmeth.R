## Supplementary Fig. S9 (six-figure restructure 2026-09-17; was Fig. S4): an epigenetic
## clock measured in the same cultures does not behave like the transcriptional readouts.
## Preregistered tests M1-M4 (preregistration_methylation_clock_ko.md); outputs from
## run_methylation_clock.R.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
D <- "public_data_tierA/derived/methylation_clock"
M <- read.delim(file.path(D, "methylation_samples.tsv"))
P <- read.delim(file.path(D, "paired_rna_methylation.tsv"))
TB <- read.delim(file.path(D, "preregistered_tests.tsv"))
fisher_ci <- function(r, n) tanh(atanh(r) + c(-1.96, 1.96) / sqrt(n - 3))
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(x) & is.finite(y))) }
## the y-axis titles of a and c ran into the tick labels "60" and "40" (2 pt default gap)
YTITLE <- element_text(margin = margin(r = 6))

## ---- a. the clock advances with cumulative divisions (sanity check) ----------
HC <- M[grepl("^HC", M$cell_line) & M$treatments == "Control" & M$percent_oxygen == 21, ]
sa <- sp(HC$population_doublings_utotal_divisions, HC$horvath)
pa <- ggplot(HC, aes(population_doublings_utotal_divisions, horvath, colour = cell_line)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK2, fill = INK2,
              alpha = 0.12, linewidth = 0.5) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.05, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.1, fill = "white", label.size = 0, label.r = unit(0, "pt"),
           label = sprintf("%s\nn = %d libraries", rp(sa$rho, sa$p), sa$n)) +
  scale_colour_manual(values = c(BLUE, RED, TEAL, ORANGE, PURPLE, GREEN), name = NULL) +
  scale_y_continuous(labels = num_axis(0)) +
  labs(x = "cumulative population doublings", y = "DNAm age (years, uncalibrated)") +
  theme_sa() + theme(legend.position = "top", legend.direction = "horizontal",
                     legend.margin = margin(0, 0, 1, 0), legend.key.size = unit(6, "pt"),
                     legend.key.spacing.x = unit(7, "pt"),
                     legend.text = element_text(size = 7),
                     axis.title.y = YTITLE) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.4)))

## ---- b. but its residual does not track the measured division rate ----------
sb <- sp(M$rate, M$accel)
pb <- ggplot(M, aes(rate, accel)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3, linetype = "22") +
  geom_point(size = 1.1, alpha = 0.55, colour = INK2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = RED, fill = RED,
              alpha = 0.12, linewidth = 0.5) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.05, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.1, fill = "white", label.size = 0, label.r = unit(0, "pt"),
           label = sprintf("%s\nn = %d libraries", rp(sb$rho, sb$p), sb$n)) +
  ## head-room: the white statistic box used to cover about a dozen libraries
  scale_y_continuous(labels = num_axis(0), expand = expansion(mult = c(0.05, 0.26))) +
  labs(x = "measured division rate (divisions per day)",
       y = "clock acceleration (years)") + theme_sa() + theme(axis.title.y = YTITLE)

## ---- c. the two biomarker families are largely independent -------------------
sc <- sp(P$prolif, P$horvath); ci <- fisher_ci(sc$rho, sc$n)
ss <- sp(P$splice, P$horvath)
pc <- ggplot(P, aes(prolif, horvath, colour = cell_line)) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK2, fill = INK2,
              alpha = 0.10, linewidth = 0.5) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.05, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.15, fill = "white", label.size = 0, label.r = unit(0, "pt"),
           label = sprintf("proliferation score: %s = %s [%s, %s]\nsplicing score: %s = %s\nn = %d paired cultures",
                           RHO, num(sc$rho), num(ci[1]), num(ci[2]), RHO, num(ss$rho), sc$n)) +
  scale_colour_manual(values = c(BLUE, RED, ORANGE), name = NULL) +
  scale_x_continuous(labels = num_axis(0)) +   # real minus sign on the negative scores
  scale_y_continuous(labels = num_axis(0), expand = expansion(mult = c(0.05, 0.36))) +
  labs(x = "transcriptomic proliferation score", y = "DNAm age (years, uncalibrated)") +
  theme_sa() + theme(legend.position = "top", legend.direction = "horizontal",
                     legend.margin = margin(0, 0, 1, 0), legend.key.size = unit(6, "pt"),
                     legend.key.spacing.x = unit(7, "pt"),
                     legend.text = element_text(size = 7),
                     axis.title.y = YTITLE) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5)))

## ---- d. four clocks, same question ------------------------------------------
clocks <- c(horvath = "Horvath 2013", skinblood = "skin & blood", phenoage = "PhenoAge", hannum = "Hannum")
Dd <- do.call(rbind, lapply(names(clocks), function(cl) {
  acc <- residuals(lm(M[[cl]] ~ M$days_grown_udays + factor(M$cell_line), na.action = na.exclude))
  s <- sp(M$rate, acc); ci <- fisher_ci(s$rho, s$n)
  data.frame(clock = clocks[[cl]], rho = s$rho, lo = ci[1], hi = ci[2], p = s$p) }))
Dd$clock <- factor(Dd$clock, levels = rev(unname(clocks)))
Dd$lab <- vapply(Dd$p, pfmt, character(1))
pd <- ggplot(Dd, aes(rho, clock)) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.14, colour = INK2, linewidth = 0.4) +
  geom_point(aes(fill = p < 0.05), shape = 21, size = 2.2, colour = INK2, stroke = 0.4) +
  geom_text(aes(x = hi, label = lab), hjust = -0.18, family = FONT, size = pt(7), colour = GREY) +
  scale_fill_manual(values = c(`TRUE` = BLUE, `FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(-0.22, 0.52), labels = num_axis(1)) +
  labs(x = sprintf("%s between clock acceleration and measured division rate", RHO), y = NULL) +
  theme_sa() + theme(axis.text.y = element_text(size = 7.5))

save_fig(lab_grid(lab_grid(pa, pb, labels = c("a", "b"), ncol = 2),
                  lab_grid(pc, pd, labels = c("c", "d"), ncol = 2, rel_widths = c(1, 1.05)),
                  labels = c("", ""), ncol = 1),
         "FigS7.png", 183, 124)
