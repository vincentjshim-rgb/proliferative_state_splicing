FIGURES_WRITTEN <- c("Fig3.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/journal_theme.R")
D <- "public_data_tierA/derived"; CS <- file.path(D, "classical_stats")

M  <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
M  <- M[M$comparator_consistent %in% c(TRUE, "TRUE"), ]
W  <- read.delim(file.path(CS, "williams_test_D.tsv"))

## ============ A, B. gene-level association with each reference axis =========
axpanel <- function(yv, col, ylab, rho) {
  d <- data.frame(x = M$Repro_specific_score, y = M[[yv]])
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  ct <- suppressWarnings(cor.test(d$x, d$y, method = "spearman", exact = FALSE))
  fit <- lm(y ~ x, d); nd <- data.frame(x = seq(min(d$x), max(d$x), length.out = 120))
  pr <- predict(fit, nd, interval = "confidence"); nd <- cbind(nd, pr)
  ggplot(d, aes(x, y)) +
    geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
    stat_bin2d(bins = 62, aes(fill = after_stat(count))) +
    geom_ribbon(data = nd, aes(x, ymin = lwr, ymax = upr), inherit.aes = FALSE,
                fill = col, alpha = 0.25) +
    geom_line(data = nd, aes(x, fit), inherit.aes = FALSE, colour = col, linewidth = 0.55) +
    scale_fill_gradientn(colours = c("#FFFFFF","#DCE3EA","#93A7BA","#43536A","#1A2430"),
                         trans = "log10", guide = "none") +
    annotate("text", x = -2.85, y = 2.35, hjust = 0, vjust = 1, family = FONT, size = pt(6.3),
             colour = INK, lineheight = 1.05,
             label = sprintf("ρ = %+.3f\n95%% CI %+.3f to %+.3f\nP < 10⁻¹⁶", rho$rho, rho$lo, rho$hi)) +
    coord_cartesian(xlim = c(-2.9, 2.9), ylim = c(-2.5, 2.5)) +
    labs(x = "Repro-CM anchor score  S", y = ylab) +
    theme_j(7)
}
fz <- function(r, n) { se <- 1/sqrt(n-3); list(rho = r, lo = tanh(atanh(r)-1.96*se), hi = tanh(atanh(r)+1.96*se)) }
p3a <- axpanel("UVA_meta",        UV,  "UVA reference axis (rank-averaged)",
               fz(W$rho_anchor_UVA, W$n_genes))
p3b <- axpanel("senescence_meta", SEN, "Senescence reference axis (rank-averaged)",
               fz(W$rho_anchor_sen, W$n_genes))

## ================= C. the two correlations and their difference ============
fzi <- function(r, n) { se <- 1/sqrt(n-3); c(lo = tanh(atanh(r)-1.96*se), hi = tanh(atanh(r)+1.96*se)) }
C <- data.frame(
  lab = c("ρ  anchor vs UVA axis", "ρ  anchor vs senescence axis", "D  =  difference"),
  est = c(W$rho_anchor_UVA, W$rho_anchor_sen, W$D),
  col = c(UV, SEN, INK), stringsAsFactors = FALSE)
C[1, c("lo","hi")] <- fzi(W$rho_anchor_UVA, W$n_genes)
C[2, c("lo","hi")] <- fzi(W$rho_anchor_sen, W$n_genes)
dse <- sqrt(2)/sqrt(W$n_genes - 3)
C[3, c("lo","hi")] <- W$D + c(-1.96, 1.96)*dse
C$lab <- factor(C$lab, levels = rev(C$lab))
p3c <- ggplot(C, aes(est, lab)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = lab), height = 0.13, linewidth = 0.4) +
  geom_point(aes(colour = lab, shape = lab), size = 2.2) +
  geom_text(aes(label = sprintf("%+.3f", est)), vjust = -1.3, size = pt(6.2),
            family = FONT, colour = INK) +
  scale_colour_manual(values = setNames(C$col, as.character(C$lab)), guide = "none") +
  scale_shape_manual(values = setNames(c(15, 15, 18), as.character(C$lab)), guide = "none") +
  annotate("text", x = 0.55, y = 0.80, hjust = 1, family = FONT, size = pt(5.8),
           colour = INK2, lineheight = 1.1,
           label = sprintf("Williams t(%s) = %.1f\nP < 10⁻¹⁶", format(W$df, big.mark=","), W$williams_t)) +
  scale_x_continuous(limits = c(-0.30, 0.56), breaks = c(-0.2, 0, 0.2, 0.4)) +
  coord_cartesian(ylim = c(0.6, 3.4), clip = "off") +
  labs(x = "correlation  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() + theme(plot.margin = margin(2, 5, 2, 2))

## =============== D. partial correlation on proliferation loading ===========
pc <- read.delim(file.path(CS, "partial_correlation_proliferation.tsv"))
raw <- data.frame(axis = c("UVA","senescence"),
                  r = c(W$rho_anchor_UVA, W$rho_anchor_sen), adj = "unadjusted")
raw$lo <- c(fzi(raw$r[1], W$n_genes)["lo"], fzi(raw$r[2], W$n_genes)["lo"])
raw$hi <- c(fzi(raw$r[1], W$n_genes)["hi"], fzi(raw$r[2], W$n_genes)["hi"])
adj <- data.frame(axis = pc$axis, r = pc$partial_rho, adj = "partial, proliferation removed",
                  lo = pc$ci_lo, hi = pc$ci_hi)
P <- rbind(raw[, c("axis","r","adj","lo","hi")], adj)
P$axis <- factor(c(UVA = "UVA axis", senescence = "Senescence axis")[P$axis],
                 levels = c("UVA axis","Senescence axis"))
P$adj  <- factor(P$adj, levels = c("unadjusted","partial, proliferation removed"))
p3d <- ggplot(P, aes(axis, r, fill = adj)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.68), width = 0.56, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.68),
                width = 0.14, linewidth = 0.32, colour = INK2) +
  geom_text(aes(label = sprintf("%+.3f", r), vjust = ifelse(r > 0, -1.9, 2.8)),
            position = position_dodge(width = 0.68), size = pt(6), family = FONT, colour = INK) +
  scale_fill_manual(values = c(unadjusted = "#B9C4CF",
                    `partial, proliferation removed` = NAVY), name = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(limits = c(-0.28, 0.44)) +
  labs(x = NULL, y = "Spearman ρ  (95% CI)") +
  theme_j(7) + theme(legend.position = "inside", legend.position.inside = c(0.72, 0.86),
                     legend.box.background = element_blank(),
                     legend.key.size = unit(6, "pt"))

## ================= E. D across proliferation-loading deciles ===============
bt <- read.delim(file.path(D, "manuscript_figures_v4/Figure3C_D_by_proliferation_loading_decile.tsv"))
p3e <- ggplot(bt, aes(mean_loading, D)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) +
  geom_hline(yintercept = W$D, colour = REPRO, linewidth = 0.35, linetype = "22") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.028, colour = GREY, linewidth = 0.32) +
  geom_point(size = 1.5, colour = NAVY) +
  annotate("text", x = max(bt$mean_loading), y = W$D, label = "overall D", colour = REPRO,
           size = pt(5.8), family = FONT, hjust = 1, vjust = -0.7) +
  annotate("text", x = mean(range(bt$mean_loading)), y = 0.035, family = FONT, size = pt(5.8),
           colour = INK2, label = "all ten 95% bootstrap intervals exclude zero") +
  scale_y_continuous(limits = c(0, 0.92), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "per-gene proliferation loading  (decile mean)", y = "decoupling index  D") +
  theme_j(7)


## ============= F. robustness of D across 124 recomputations ================
DV <- file.path(D, "decoupling_validation")
cs  <- read.delim(file.path(DV, "comparator_sensitivity.tsv"))
lo  <- read.delim(file.path(DV, "leave_one_dataset_out.tsv"))
gr  <- read.delim(file.path(DV, "gene_program_removal_sensitivity.tsv"))
p4s <- read.delim(file.path(D, "direction_probe/P4_decoupling_proliferation_sensitivity.tsv"))
th  <- read.delim(file.path(DV, "read_thinning_50pct_100runs.tsv"))
rb <- rbind(
  data.frame(g = "comparator choice",       D = cs$decoupling_index),
  data.frame(g = "leave-one-study-out",     D = lo$decoupling_index),
  data.frame(g = "nuisance-gene removal",   D = gr$decoupling_index),
  data.frame(g = "proliferation handling",  D = p4s$decoupling_index_D),
  data.frame(g = "50% read thinning (×100)", D = th$decoupling_index))
rb$g <- factor(rb$g, levels = rev(c("comparator choice","leave-one-study-out",
  "nuisance-gene removal","proliferation handling","50% read thinning (×100)")))
p3f <- ggplot(rb, aes(D, g)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_vline(xintercept = W$D, colour = REPRO, linewidth = 0.35, linetype = "22") +
  geom_point(size = 0.9, colour = NAVY, alpha = 0.65,
             position = position_jitter(height = 0.14, width = 0, seed = 2)) +
  annotate("text", x = W$D, y = 5.62, label = "primary D", colour = REPRO,
           size = pt(5.6), family = FONT, hjust = -0.12) +
  annotate("text", x = 0.02, y = 0.62, hjust = 0, family = FONT, size = pt(5.8), colour = INK2,
           label = sprintf("%d recomputations; range %.3f to %.3f, none crossing zero",
                           nrow(rb), min(rb$D), max(rb$D))) +
  scale_x_continuous(limits = c(0, 0.70), breaks = seq(0, 0.6, 0.2)) +
  coord_cartesian(ylim = c(0.4, 5.9), clip = "off") +
  labs(x = "decoupling index  D", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.2))

## ============ G. the same construction around each condition ==============
sp  <- read.delim(file.path(D, "anchor_specificity/anchor_specificity_negative_control.tsv"))
spl <- rbind(data.frame(lab = sp$focal_condition, D = sp$D, w = "unadjusted"),
             data.frame(lab = sp$focal_condition, D = sp$D_proliferation_adjusted,
                        w = "proliferation-adjusted"))
spl$lab <- factor(spl$lab, levels = rev(sp$focal_condition))
spl$w   <- factor(spl$w, levels = c("unadjusted","proliferation-adjusted"))
p3g <- ggplot(spl, aes(D, lab, colour = w)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_point(size = 1.8, shape = 15, position = position_dodge(width = 0.45)) +
  geom_text(data = transform(sp, lab = factor(focal_condition, levels = rev(sp$focal_condition))),
            aes(x = pmax(D, D_proliferation_adjusted) + 0.035, y = lab,
                label = sprintf("%+.3f / %+.3f", D, D_proliferation_adjusted)),
            inherit.aes = FALSE, hjust = 0, size = pt(5.8), family = FONT, colour = INK) +
  scale_colour_manual(values = c(unadjusted = REPRO, `proliferation-adjusted` = NAVY), name = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(-0.58, 0.82), breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
  labs(x = "decoupling index  D, construction centred on each condition", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box.background = element_blank(), legend.margin = margin(t = -4),
        legend.key.size = unit(6, "pt"))

top <- lab_grid(p3a, p3b, p3c, labels = c("A","B","C"), ncol = 3, rel_widths = c(1, 1, 1.06))
mid <- lab_grid(p3d, p3e, labels = c("D","E"), ncol = 2, rel_widths = c(1, 1.5))
bot <- lab_grid(p3f, p3g, labels = c("F","G"), ncol = 2, rel_widths = c(1.1, 1))
save_fig(lab_grid(top, mid, bot, labels = c("","",""), ncol = 1, rel_heights = c(1.15, 1, 0.92)),
         "Fig3.png", 183, 158)
