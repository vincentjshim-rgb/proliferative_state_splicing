## Supplementary Fig. S1 (S5 before the six-figure restructure of 2026-09-17): nine
## perturbations of the same cell lines. Splicing-factor expression stays with the
## proliferation programme; what the drugs pull apart is the programme and the counted
## divisions. Outputs from scripts/revision/run_treatment_perturbation.R, which since
## 2026-09-17 uses the 177-gene pre-mRNA processing score of Fig. 2 (the axis name follows).
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
D <- "public_data_tierA/derived/treatment_perturbation"
N <- read.delim(file.path(D, "matched_deltas.tsv"))
C <- read.delim(file.path(D, "treatment_coefficients.tsv"))

LAB <- c(Contact_Inhibition = "contact inhibition", `Oligomycin+DEX` = "oligomycin + DEX",
         Oligomycin = "oligomycin", `mitoNUITs+DEX` = "mitoNUITs + DEX", mitoNUITs = "mitoNUITs",
         DEX = "DEX", `2-Deoxyglucose` = "2-deoxyglucose",
         betahydroxybutyrate = "β-hydroxybutyrate", Galactose = "galactose")
SHORT <- c(Contact_Inhibition = "contact inh.", `Oligomycin+DEX` = "oligo + DEX",
           Oligomycin = "oligomycin", `mitoNUITs+DEX` = "mitoN + DEX", mitoNUITs = "mitoNUITs",
           DEX = "DEX", `2-Deoxyglucose` = "2-DG", betahydroxybutyrate = "β-HB",
           Galactose = "galactose")
N$lab <- SHORT[N$cond]
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(x) & is.finite(y))) }
MD <- aggregate(cbind(d_rate, d_splice, d_prolif) ~ cond, N, median)
MD$lab <- SHORT[MD$cond]
## Every arm is named, as in panel c. Until 2026-09-17 five were named with nudged text
## and no leader: oligomycin (discussed in the text) was missing, its median lies almost on
## that of mitoNUITs + DEX in panel a, and "2-DG" sat on the fit line away from its point.
## Labels now sit in the sparse parts of each panel, joined to their median by a leader;
## lx/ly is the corner of the label nearest the point, in data units of that panel.
## In panel a the medians of oligomycin and mitoNUITs + DEX overlap, so their leaders
## arrive from different sides (left, and above) and each ends on its own symbol.
LPOS <- read.table(header = TRUE, sep = ",", strip.white = TRUE, check.names = FALSE, text = "
cond,                 ax,    ay,   ahj,  bx,    by,   bhj
Contact_Inhibition,  -0.32, -1.25, 0,   -0.55, -1.52, 0
Oligomycin+DEX,      -0.47, -0.22, 1,    0.20, -1.24, 0
Oligomycin,          -0.33,  0.34, 1,   -0.75,  0.62, 1
mitoNUITs+DEX,       -0.118, 0.62, 0.92, 0.85, -0.42, 0
mitoNUITs,           -0.085,-1.12, 0.5,  0.42, -0.98, 0
DEX,                 -0.075, 1.08, 1,    0.30,  0.95, 1
2-Deoxyglucose,       0.10, -0.30, 0,    1.45, -0.22, 0
betahydroxybutyrate,  0.05, -0.72, 0,    0.62, -0.70, 0
Galactose,            0.09,  1.12, 0,   -0.10,  0.80, 1")
MD <- merge(MD, LPOS, by = "cond")
stopifnot(nrow(MD) == 9)

delta_panel <- function(xv, xlab, lx, ly, lhj, xbreaks = waiver()) {
  s <- sp(N[[xv]], N$d_splice)
  ggplot(N, aes(.data[[xv]], d_splice)) +
    geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
    geom_point(size = 1.1, alpha = 0.35, colour = GREY) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK2, fill = INK2,
                alpha = 0.12, linewidth = 0.5) +
    geom_segment(data = MD, aes(x = .data[[xv]], y = d_splice, xend = .data[[lx]], yend = .data[[ly]]),
                 colour = INK2, linewidth = 0.25) +
    geom_label(data = MD, aes(.data[[lx]], .data[[ly]], label = lab, hjust = .data[[lhj]]),
               vjust = 0.5, size = pt(7), family = FONT, colour = INK, fill = "white",
               label.size = 0, label.r = unit(0, "pt"), label.padding = unit(1.2, "pt")) +
    geom_point(data = MD, aes(.data[[xv]], d_splice), size = 2.3, shape = 21,
               fill = BLUE, colour = "white", stroke = 0.5) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.02, vjust = 1.05, family = FONT,
             size = pt(7), colour = INK, fill = "white", label.size = 0,
             label.r = unit(0, "pt"), lineheight = 1.1,
             label = sprintf("%s\nn = %d treated cultures", rp(s$rho, s$p), s$n)) +
    scale_x_continuous(breaks = xbreaks, labels = num_axis(1)) + scale_y_continuous(labels = num_axis(1)) +
    labs(x = xlab, y = "change in pre-mRNA processing score") + theme_sa() +
    theme(axis.title.y = element_text(margin = margin(r = 4)))   # keeps the title's descenders off the tick labels
}
## breaks that one decimal prints exactly; the default 0.25 steps were being printed as 0.2
pa <- delta_panel("d_rate", "change in measured division rate (divisions per day)", "ax", "ay", "ahj",
                  xbreaks = seq(-0.6, 0.2, 0.2))
pb <- delta_panel("d_prolif", "change in proliferation score", "bx", "by", "bhj")

## ---- c. residual treatment effect under the two adjustments ------------------
keep <- c("splice, measured rate held constant", "splice, proliferation score held constant")
Cc <- C[C$response %in% keep, ]
Cc$model <- factor(ifelse(grepl("measured rate", Cc$response),
                          "measured division rate held constant",
                          "proliferation score held constant"),
                   levels = c("measured division rate held constant",
                              "proliferation score held constant"))
ord <- MD$cond[order(MD$d_rate)]
Cc$term <- factor(Cc$term, levels = rev(ord))
Cc$sig <- Cc$FDR < 0.05
Cc$fillcol <- ifelse(Cc$sig, ifelse(grepl("measured rate", Cc$response), RED, BLUE), "white")
## group = model: without it the fill colour joins the dodge grouping and a filled symbol
## lands on the other model's row, off its own confidence interval
pc <- ggplot(Cc, aes(beta, term, colour = model, group = model)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3, linetype = "22") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45,
                 position = position_dodge(width = 0.62)) +
  geom_point(aes(fill = fillcol), shape = 21, size = 2.1, stroke = 0.5,
             position = position_dodge(width = 0.62)) +
  scale_colour_manual(values = c(RED, BLUE), name = NULL) +
  scale_fill_identity(guide = "none") +
  scale_y_discrete(labels = function(x) LAB[x]) +
  scale_x_continuous(breaks = seq(-0.5, 1, 0.25), labels = num_axis(2)) +
  labs(x = "residual effect of the treatment on the pre-mRNA processing score (z per treatment)",
       y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.direction = "horizontal",
                     legend.margin = margin(0, 0, 1, 0), legend.key.size = unit(7, "pt"),
                     axis.text.y = element_text(size = 7.5),
                     plot.margin = margin(3, 9, 3, 3))   # room for the last tick label

save_fig(lab_grid(lab_grid(pa, pb, labels = c("a", "b"), ncol = 2),
                  lab_grid(pc, labels = c("c"), ncol = 1),
                  labels = c("", ""), ncol = 1, rel_heights = c(1, 1.02)),
         "FigS1.png", 183, 132)
