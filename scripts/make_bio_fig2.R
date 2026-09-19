## Fig. 2  Does splicing-factor expression follow the measured division rate?
## Revised 2026-09-15: the 328 libraries are repeated measures of seven cell lines,
## so the panel reports the mixed-model estimate beside the pooled correlation and
## panel e shows every line separately, including the one that is not significant.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
d  <- read.delim(file.path(D, "division_rate/sample_division_rate.tsv"))
R  <- read.delim(file.path(D, "division_rate/named_splicing_factors.tsv"))
PL <- read.delim(file.path(D, "revision_stats/fig2_per_line_correlations.tsv"))
MM <- read.delim(file.path(D, "revision_stats/fig2_mixed_models.tsv"))
EXPC <- "#B2182B"
LINES <- sort(unique(d$line))
LCOL <- setNames(c(BLUE, RED, TEAL, ORANGE, PURPLE, GREEN, BROWN)[seq_along(LINES)], LINES)

## ====== a. expression against directly measured division rate =============
ct <- cor.test(d$rate, d$splice, method = "spearman", exact = FALSE)
ci <- tanh(atanh(unname(ct$estimate)) + c(-1.96, 1.96) / sqrt(nrow(d) - 3))
mm1 <- MM[MM$model == "random intercept per line", ]
## these two P values lie far below the 2 x 10^-16 floor that pfmt() prints, and the
## manuscript quotes them, so they are written out with the same sci() typography
lab_a  <- sprintf("%s = %s  [%s, %s]\nP = %s",
                  RHO, num(ct$estimate), num(ci[1]), num(ci[2]), sci(ct$p.value, 1))
lab_mm <- sprintf("mixed model, line as random effect\nβ = %s, P = %s",
                  num(mm1$beta), sci(mm1$p, 1))
## The pooled statistic stays in the top-left corner, which holds no data. The
## mixed-model estimate used to sit below it, on top of the points, so it is printed
## in the empty lower-right corner above the n line instead (no white box is drawn,
## because a box would hide cultures).
YSPAN <- diff(range(d$splice))
p2a <- ggplot(d, aes(rate, splice)) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = GREY, alpha = 0.18, linewidth = 0.6) +
  geom_point(aes(colour = line), size = 1.25, alpha = 0.85) +
  annotate("text", x = 0.02, y = max(d$splice) * 1.02, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.05, label = lab_a) +
  annotate("text", x = max(d$rate), y = min(d$splice) + 0.065 * YSPAN, hjust = 1, vjust = 0,
           family = FONT, size = pt(7), colour = INK, lineheight = 1.05, label = lab_mm) +
  annotate("text", x = max(d$rate), y = min(d$splice), hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = GREY,
           label = sprintf("n = %d samples, %d cell lines", nrow(d), length(LINES))) +
  scale_colour_manual(values = LCOL, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
  scale_y_continuous(labels = num_axis()) +
  labs(x = "measured division rate (divisions per day)", y = "pre-mRNA processing score") +
  ## right margin on the legend text keeps each key clear of the preceding label
  theme_sa(8) + theme(legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7, margin = margin(l = 1.5, r = 7)))

## ====== b. individual splicing factors ====================================
## The panel is the twenty factors assayed in senescing human dermal fibroblasts by
## Holly et al. 2013, Table 3 (scripts/revision/named_splicing_factors_source.md).
## Three of them are reported there as HIGHER in late-passage cells, so they are
## drawn open; the four highlighted are those reported as falling with age in more
## than one study cited by this manuscript.
REPEATED <- c("SRSF1", "SRSF2", "SRSF6", "HNRNPK")
UP_IN_AGE <- c("HNRNPUL2", "PNISR", "SF3B1")
SFFILE <- "scripts/revision/named_splicing_factors.txt"
if (file.exists(SFFILE)) { SRC <- readLines(SFFILE); R <- R[R$gene %in% SRC, ] }
R <- R[order(R$rho_rate), ]; R$gene <- factor(R$gene, levels = R$gene)
R$hl  <- ifelse(as.character(R$gene) %in% REPEATED, "falls with age in two or more cited studies", "other")
R$dir <- ifelse(as.character(R$gene) %in% UP_IN_AGE, "reported higher in late-passage cells", "reported lower")
## bars in place of lollipops (2026-09-19): the per-gene form used in the
## splicing-factor literature (Holly 2013, Lee 2016 present per-factor values as bars)
p2b <- ggplot(R, aes(rho_rate, gene, fill = hl)) +
  geom_col(width = 0.7, colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  geom_text(data = R[R$dir == "reported higher in late-passage cells", ],
            aes(x = pmax(rho_rate, 0) + 0.03, label = "\u2020"), hjust = 0, vjust = 0.45,
            family = FONT, size = pt(8), colour = INK2) +
  annotate("text", x = 0.93, y = 1.2, hjust = 1, vjust = 0, family = FONT, size = pt(7),
           colour = GREY, label = sprintf("median %s = %s", RHO, num(median(R$rho_rate)))) +
  scale_fill_manual(values = c(`falls with age in two or more cited studies` = EXPC,
                               other = "#9AA5AE"), name = NULL,
                    breaks = "falls with age in two or more cited studies") +
  scale_x_continuous(limits = c(-0.35, 0.95), breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
                     labels = num_axis()) +
  coord_cartesian(ylim = c(0.5, nrow(R) + 0.4), clip = "off") +
  labs(x = "correlation with measured division rate", y = NULL,
       caption = "\u2020 reported higher, not lower, in late-passage cells") +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7, face = "italic"),
                      legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7),
                      plot.caption = element_text(size = 6.5, colour = INK2, hjust = 0))

## ====== c. untreated cells across the culture lifespan =====================
ex <- d[d$cond == "Control" & d$oxy == 21 & grepl("^HC", d$line), ]
p2d <- ggplot(ex, aes(days, splice, colour = line)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.55) +
  geom_point(aes(size = rate), alpha = 0.85) +
  scale_size_continuous(range = c(0.6, 2.4), name = "divisions\nper day", breaks = c(0.2, 0.4, 0.6)) +
  scale_colour_manual(values = LCOL, guide = "none") +
  ## top-right corner is empty; in the lower-left the note lay under the HC3 fit line
  annotate("text", x = max(ex$days), y = max(ex$splice), hjust = 1, vjust = 1, family = FONT,
           size = pt(7), colour = GREY, lineheight = 1.05,
           label = sprintf("untreated, 21%% O₂\n%d lines, %d samples", length(unique(ex$line)), nrow(ex))) +
  scale_y_continuous(labels = num_axis()) +
  labs(x = "days in culture", y = "pre-mRNA processing score") +
  theme_sa(8) + theme(legend.position = "right", legend.key.size = unit(7, "pt"),
                      legend.title = element_text(size = 7), legend.text = element_text(size = 7))

## ====== e. one correlation per cell line ==================================
PL$line <- factor(PL$line, levels = rev(PL$line))
PL$sig <- ifelse(PL$p < 0.05, "P < 0.05", "not significant")
p2e <- ggplot(PL, aes(rho, line, fill = sig)) +
  geom_col(width = 0.66, colour = NA) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = INK2, linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  scale_fill_manual(values = c(`P < 0.05` = EXPC, `not significant` = "#C9CFD6"), name = NULL) +
  scale_y_discrete(labels = function(x) sprintf("%s (n = %d)", x, PL$n[match(x, as.character(PL$line))])) +
  scale_x_continuous(limits = c(-0.45, 1.05), breaks = c(0, 0.5, 1), labels = num_axis()) +
  labs(x = sprintf("%s with division rate, per cell line", RHO), y = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7),
                      legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7))

## ---- f. how good is the transcriptional proxy that stands in for counting? --------
## Added 2026-09-18.  The counted rate exists only in this resource; every other figure
## uses the 20-marker transcriptional proliferation score instead.  This panel shows how
## well that substitution holds, and which arms it fails on.
LS <- read.delim(file.path(D, "treatment_perturbation/library_scores_vs_counted_rate.tsv"))
rp_pro <- cor(LS$prolif, LS$rate, method = "spearman")
rp_spl <- cor(LS$splice, LS$rate, method = "spearman")
## the arms whose proliferation score sits furthest above what counting supports
LS$resid_arm <- residuals(lm(prolif ~ rate, LS))
armres <- sort(tapply(LS$resid_arm, LS$cond, median), decreasing = TRUE)
NAMED <- names(armres)[1:2]
LS$arm <- ifelse(LS$cond == "Control", "untreated control",
          ifelse(LS$cond %in% NAMED, LS$cond, "other treatment"))
LS$arm <- factor(LS$arm, levels = c("untreated control", NAMED, "other treatment"))
ARMCOL <- setNames(c(GREY, RED, ORANGE, BLUE), levels(LS$arm))
p2f <- ggplot(LS, aes(rate, prolif)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = INK, fill = GREY,
              alpha = 0.18, linewidth = 0.5) +
  geom_point(aes(colour = arm), size = 1.15, alpha = 0.85) +
  scale_colour_manual(values = ARMCOL, name = NULL) +
  scale_x_continuous("measured division rate (divisions per day)", labels = num_axis()) +
  scale_y_continuous("20-marker proliferation score", labels = num_axis()) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.5, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.05,
           label = sprintf("proliferation score %s = %s\nsplicing score %s = %s\nn = %d libraries",
                           RHO, num(rp_pro), RHO, num(rp_spl), nrow(LS))) +
  theme_sa(8) + theme(legend.position = "bottom", legend.margin = margin(t = -5),
                      legend.text = element_text(size = 7),
                      legend.key.height = unit(9, "pt"))
cat(sprintf("\npanel f -- proxy vs counting: proliferation score rho = %.3f, splicing score rho = %.3f\n",
            rp_pro, rp_spl))
cat(sprintf("   arms furthest above the fit: %s\n", paste(NAMED, collapse = ", ")))

top <- lab_grid(p2a, p2b, labels = c("a", "b"), ncol = 2, rel_widths = c(1.18, 1))
bot <- lab_grid(p2d, p2e, labels = c("c", "d"), ncol = 2, rel_widths = c(1.15, 1))
save_fig(lab_grid(top, bot, lab_grid(p2f, labels = "e", ncol = 1), labels = c("", "", ""), ncol = 1,
                  rel_heights = c(1.3, 1, 0.92)), "Fig2.png", 183, 190)
