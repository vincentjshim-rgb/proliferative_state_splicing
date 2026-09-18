FIGURES_WRITTEN <- c("Fig2.png")
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
PUBDIR <- "public_data_tierA/derived/figures_benchmark"
D <- "public_data_tierA/derived/benchmark_figs"
Z <- readRDS(file.path(D, "figdata.rds")); R <- Z$R; M <- Z$meta
WB <- read.delim(file.path(D, "within_between_by_class.tsv"))
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", secretome = REPRO,
        `photoprotection / rescue` = AMBER, `UV injury` = UV,
        `metabolic / culture` = GREY, reprogramming = NAVY)
WB <- WB[WB$class %in% names(CC), ]
WB$class <- factor(WB$class, levels = WB$class[order(WB$between)])

## ================ A. between-study reproducibility by class =================
p2a <- ggplot(WB, aes(between, class, colour = class)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = b_lo, xmax = b_hi), height = 0, linewidth = 0.45) +
  geom_point(size = 2.3) +
  geom_text(aes(label = sprintf("%+.3f", between)), hjust = -0.28, vjust = -0.95,
            size = pt(6.1), family = FONT, colour = INK) +
  geom_text(aes(label = sprintf("%d pairs", n_between), x = -0.125), hjust = 0,
            size = pt(5.5), family = FONT, colour = GREY) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(limits = c(-0.13, 0.62), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "median between-study ρ  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.5))

## ============= B. within-study vs between-study, and what survives =========
WB$retain <- 100 * WB$between / WB$within
L <- rbind(transform(WB, v = within,  k = "within a study"),
           transform(WB, v = between, k = "between studies"))
L$k <- factor(L$k, levels = c("within a study","between studies"))
p2b <- ggplot(L, aes(v, class)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_line(aes(group = class), colour = GREY_L, linewidth = 0.5) +
  geom_point(aes(shape = k, colour = class), size = 2.2) +
  geom_text(data = WB[is.finite(WB$retain), ],
            aes(x = 0.92, y = class, label = sprintf("%.0f%%", pmax(retain, 0))),
            hjust = 1, size = pt(5.9), family = FONT, colour = INK) +
  annotate("text", x = 0.92, y = 7.75, label = "retained", hjust = 1, size = pt(5.6),
           family = FONT, colour = INK2, fontface = "bold") +
  scale_shape_manual(values = c(`within a study` = 1, `between studies` = 16), name = NULL) +
  scale_colour_manual(values = CC, guide = "none") +
  guides(shape = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(-0.06, 0.95), breaks = c(0, 0.25, 0.5, 0.75)) +
  coord_cartesian(ylim = c(0.6, 8.0), clip = "off") +
  labs(x = "median Spearman ρ", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom", legend.justification = "center",
        legend.box.background = element_blank(), legend.margin = margin(t = -4))

## ====== C. the pair distributions: senescence vs UV injury ================
uvp <- read.delim(file.path(D, "uv_pairwise.tsv"))
senp <- read.delim(file.path(D, "sen_pairwise.tsv"))
PD <- rbind(data.frame(r = senp$r, k = "senescence"), data.frame(r = uvp$r, k = "UV injury"))
PD$k <- factor(PD$k, levels = c("UV injury","senescence"))
top <- uvp[which.max(uvp$r), ]
p2c <- ggplot(PD, aes(r, k, colour = k)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.42, linewidth = 0.34, fill = NA, show.legend = FALSE) +
  geom_point(size = 1.2, alpha = 0.8, position = position_jitter(height = 0.13, seed = 4)) +
  annotate("segment", x = top$r, xend = top$r, y = 0.68, yend = 0.86, colour = INK, linewidth = 0.3,
           arrow = arrow(length = unit(2, "pt"), type = "closed")) +
  annotate("text", x = top$r, y = 0.62, hjust = 0.5, vjust = 1, size = pt(5.7), family = FONT,
           colour = INK, lineheight = 1.0,
           label = "the only coherent UV pair (+0.348): the two studies\nfrom which a \"UVA reference axis\" would be built") +
  scale_colour_manual(values = c(`UV injury` = UV, senescence = SEN), guide = "none") +
  scale_x_continuous(limits = c(-0.32, 0.86), breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  coord_cartesian(ylim = c(0.18, 2.45), clip = "off") +
  labs(x = "Spearman ρ, individual between-study pairs", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.6))

top_row <- lab_grid(p2a, p2b, labels = c("A","B"), ncol = 2, rel_widths = c(1.2, 1))
save_fig(lab_grid(top_row, p2c, labels = c("","C"), ncol = 1, rel_heights = c(1, 0.72)),
         "Fig2.png", 183, 104)
