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
PUBDIR <- "public_data_tierA/derived/figures_publication"
source("scripts/pub_theme.R")
D <- "public_data_tierA/derived"

## a  alignment QC of the four local libraries
qc <- read.delim(file.path(D, "local_bam_v41/star_qc_four_samples.tsv"))
qc$sample <- factor(qc$sample, levels = rev(qc$sample))
q <- rbind(data.frame(sample = qc$sample, m = "uniquely mapped", v = qc$unique),
           data.frame(sample = qc$sample, m = "assigned to genes", v = qc$assigned))
q$m <- factor(q$m, levels = c("uniquely mapped", "assigned to genes"))
p2a <- ggplot(q, aes(v, sample, colour = m)) +
  geom_segment(aes(x = 80, xend = v, yend = sample), colour = GREY_L, linewidth = 0.4,
               position = position_dodge(width = 0.55)) +
  geom_point(size = 1.4, position = position_dodge(width = 0.55)) +
  geom_text(data = subset(q, sample == levels(q$sample)[1]),
            aes(label = m), position = position_dodge(width = 0.55),
            hjust = -0.16, vjust = 0.4, size = pt(5.2), family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c("uniquely mapped" = BLUE, "assigned to genes" = GREEN),
                      guide = "none") +
  scale_x_continuous(limits = c(80, 108), breaks = c(80, 85, 90, 95)) +
  coord_cartesian(clip = "off") +
  labs(x = "% of input read pairs", y = NULL) +
  theme_pub(7) + theme_cat_y()

## b  anchor construction
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"), check.names = FALSE)
AX$cls <- ifelse(!AX$comparator_consistent, "discordant",
          ifelse(abs(AX$Repro_specific_score) >= log2(1.25), "consistent, strong", "consistent"))
set.seed(1); AXs <- AX[sample(nrow(AX)), ]
p2b <- ggplot(AXs, aes(logFC_Repro_vs_HDF, logFC_Repro_vs_iPSC)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = cls), size = 0.32, stroke = 0, alpha = 0.6) +
  scale_colour_manual(values = c(discordant = "#E4E8EC", consistent = "#A9C6DE",
                                 `consistent, strong` = BLUE), guide = "none") +
  annotate("text", x = -2.9, y = 2.75, hjust = 0, vjust = 1, size = pt(5.4), family = FONT,
           colour = INK, lineheight = 1.15,
           label = "8,921 / 12,319 direction-consistent\n2,204 at |score| ≥ log2(1.25)") +
  annotate("point", x = c(-2.9, -2.9, -2.9), y = c(-2.15, -2.5, -2.85), size = 0.9,
           colour = c("#E4E8EC", "#A9C6DE", BLUE)) +
  annotate("text", x = -2.72, y = c(-2.15, -2.5, -2.85), hjust = 0, vjust = 0.42,
           size = pt(5.2), family = FONT, colour = GREY_D,
           label = c("discordant (score 0)", "direction-consistent", "|score| ≥ log2(1.25)")) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_x_continuous(breaks = seq(-3, 3, 1.5)) + scale_y_continuous(breaks = seq(-3, 3, 1.5)) +
  labs(x = expression(paste(log[2], " FC  Repro-CM vs HDF-CM")),
       y = expression(paste(log[2], " FC  vs iPSC-CM"))) +
  theme_pub(7)

## c  robustness of D
cs <- read.delim(file.path(D, "decoupling_validation/comparator_sensitivity.tsv"))
lo <- read.delim(file.path(D, "decoupling_validation/leave_one_dataset_out.tsv"))
gr <- read.delim(file.path(D, "decoupling_validation/gene_program_removal_sensitivity.tsv"))
th50 <- read.delim(file.path(D, "decoupling_validation/read_thinning_50pct_100runs.tsv"))
p4 <- read.delim(file.path(D, "direction_probe/P4_decoupling_proliferation_sensitivity.tsv"))
rb <- rbind(
  data.frame(g = "comparator", D = cs$decoupling_index),
  data.frame(g = "leave-one-study-out", D = lo$decoupling_index),
  data.frame(g = "gene-programme removal", D = gr$decoupling_index),
  data.frame(g = "proliferation", D = p4$decoupling_index_D),
  data.frame(g = "50% read thinning", D = th50$decoupling_index))
rb$g <- factor(rb$g, levels = rev(c("comparator", "leave-one-study-out", "gene-programme removal",
                                    "proliferation", "50% read thinning")))
p2c <- ggplot(rb, aes(D, g)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_vline(xintercept = 0.460, colour = VERM, linewidth = 0.3, linetype = "22") +
  geom_point(size = 0.85, colour = SLATE, alpha = 0.7,
             position = position_jitter(height = 0.13, width = 0, seed = 2)) +
  annotate("text", x = 0.460, y = 5.62, label = "primary D", colour = VERM,
           size = pt(5.2), family = FONT, hjust = -0.1, vjust = 0.4) +
  scale_x_continuous(limits = c(0, 0.70), breaks = seq(0, 0.6, 0.2)) +
  coord_cartesian(clip = "off") +
  labs(x = "decoupling index D", y = NULL) +
  theme_pub(7) + theme_cat_y()

## d  specificity
sp <- read.delim(file.path(D, "anchor_specificity/anchor_specificity_negative_control.tsv"))
sp$lab <- factor(sp$focal_condition, levels = rev(sp$focal_condition))
spl <- rbind(data.frame(lab = sp$lab, D = sp$D, w = "unadjusted"),
             data.frame(lab = sp$lab, D = sp$D_proliferation_adjusted, w = "proliferation-adjusted"))
p2d <- ggplot(spl, aes(D, lab, colour = w)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_point(size = 1.35, position = position_dodge(width = 0.5)) +
  annotate("point", x = c(0.30, 0.30), y = c(1.34, 1.18), size = 1.2,
           colour = c(SLATE, GREEN)) +
  annotate("text", x = 0.35, y = c(1.34, 1.18), hjust = 0, vjust = 0.42,
           size = pt(5.2), family = FONT, colour = GREY_D,
           label = c("unadjusted", "proliferation-adjusted")) +
  scale_colour_manual(values = c(unadjusted = SLATE, `proliferation-adjusted` = GREEN),
                      guide = "none") +
  scale_x_continuous(limits = c(-0.62, 0.95), breaks = seq(-0.5, 0.5, 0.5)) +
  coord_cartesian(clip = "off") +
  labs(x = "decoupling index D, focal condition swapped", y = NULL) +
  theme_pub(7) + theme_cat_y()

fig2 <- compose(p2a, p2b, p2c, p2d, labels = c("a", "b", "c", "d"), ncol = 2,
                rel_heights = c(1, 1.05))
save_fig(fig2, "Fig2.png", 183, 104)
cat("Fig2 saved\n")
