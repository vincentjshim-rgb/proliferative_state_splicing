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
p2b <- ggplot(R, aes(rho_rate, gene, colour = hl)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho_rate, yend = gene), linewidth = 0.45) +
  geom_point(aes(shape = dir, fill = hl), size = 1.8, stroke = 0.5) +
  annotate("text", x = 0.93, y = 1.2, hjust = 1, vjust = 0, family = FONT, size = pt(7),
           colour = GREY, label = sprintf("median %s = %s", RHO, num(median(R$rho_rate)))) +
  scale_colour_manual(values = c(`falls with age in two or more cited studies` = EXPC,
                                 other = "#9AA5AE"), name = NULL,
                      breaks = "falls with age in two or more cited studies") +
  scale_fill_manual(values = c(`falls with age in two or more cited studies` = EXPC,
                               other = "#9AA5AE"), guide = "none") +
  scale_shape_manual(values = c(`reported lower` = 21, `reported higher in late-passage cells` = 1),
                     name = NULL, breaks = "reported higher in late-passage cells") +
  guides(colour = guide_legend(nrow = 1, order = 1),
         shape  = guide_legend(nrow = 1, order = 2, override.aes = list(colour = INK2))) +
  scale_x_continuous(limits = c(-0.35, 0.95), breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
                     labels = num_axis()) +
  coord_cartesian(ylim = c(0.5, nrow(R) + 0.4), clip = "off") +
  labs(x = "correlation with measured division rate", y = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      axis.text.y = element_text(size = 7, face = "italic"),
                      legend.position = "bottom", legend.box = "vertical",
                      legend.box.just = "left", legend.spacing.y = unit(1, "pt"),
                      legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(5, "pt"),
                      legend.text = element_text(size = 7))

## ====== c. the same factors against time in culture =======================
ctc <- cor.test(R$rho_rate, R$rho_days, method = "spearman", exact = FALSE)
LABC <- data.frame(gene = c("HNRNPD", "HNRNPA1", "HNRNPA0"),
                   hj   = c(-0.15, -0.15, -0.15),
                   vj   = c(0.4, -0.7, 0.4))
LABC <- merge(LABC, R[, c("gene", "rho_rate", "rho_days")], by = "gene")
p2c <- ggplot(R, aes(rho_rate, rho_days)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_point(aes(colour = hl), size = 1.8) +
  scale_colour_manual(values = c(`falls with age in two or more cited studies` = EXPC,
                                 other = "#9AA5AE"), guide = "none") +
  geom_text(data = LABC, aes(label = gene, hjust = hj, vjust = vj),
            size = pt(7), family = FONT, fontface = "italic", colour = INK) +
  annotate("text", x = -0.33, y = -0.86, hjust = 0, vjust = 0, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.05,
           label = sprintf("%s\nn = %d factors\nred: highlighted in b", rp(ctc$estimate, ctc$p.value), nrow(R))) +
  scale_x_continuous(limits = c(-0.35, 1.38), breaks = c(0, 0.4, 0.8),
                     labels = num_axis()) +
  ## upper limit 0.2, not 0.12: IMP3 (+0.14 with days in culture) was being dropped,
  ## so the panel drew 19 of the 20 factors it counts
  scale_y_continuous(limits = c(-0.88, 0.2), labels = num_axis()) +
  labs(x = "correlation with division rate", y = "correlation with days in culture") +
  theme_sa(8)

## ====== d. untreated cells across the culture lifespan ====================
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
p2e <- ggplot(PL, aes(rho, line)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, colour = INK2, linewidth = 0.4) +
  geom_point(aes(fill = sig), shape = 21, size = 2, colour = INK2, stroke = 0.4) +
  geom_text(aes(x = 1.52, label = sprintf("n = %d", n)), hjust = 1, size = pt(7),
            family = FONT, colour = GREY) +
  scale_fill_manual(values = c(`P < 0.05` = EXPC, `not significant` = "white"), name = NULL) +
  scale_x_continuous(limits = c(-0.45, 1.55), breaks = c(0, 0.5, 1), labels = num_axis()) +
  labs(x = sprintf("%s with division rate, per cell line", RHO), y = NULL) +
  theme_sa(8) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
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
bot <- lab_grid(p2c, p2d, p2e, labels = c("c", "d", "e"), ncol = 3, rel_widths = c(1, 1.15, 1.02))
save_fig(lab_grid(top, bot, lab_grid(p2f, labels = "f", ncol = 1), labels = c("", "", ""), ncol = 1,
                  rel_heights = c(1.3, 1, 0.92)), "Fig2.png", 183, 190)
