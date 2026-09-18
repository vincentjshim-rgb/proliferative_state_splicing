#!/usr/bin/env Rscript
# Figure 1 : what the two reference axes do and do not capture
# Figure 2 : the local anchor - quality, construction, robustness, specificity
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "public_data_tierA"; der <- file.path(root, "derived")
out  <- file.path(der, "manuscript_figures_v4")
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"; PURPLE <- "#CC79A7"
GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.3),
        plot.caption = element_text(colour = MUTED, size = 7.1, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8.1))

## ===================== FIGURE 1A =====================
U <- read.delim(file.path(der, "axis_characterisation", "UVA_axis_generalisation.tsv"))
S <- read.delim(file.path(der, "axis_characterisation", "senescence_axis_generalisation.tsv"))
loc <- read.delim(file.path(der, "compendium", "compendium_D_scores.tsv"))
lU <- loc$rho_UVA[loc$id == "LOCAL_ReproCM"]; lS <- loc$rho_senescence[loc$id == "LOCAL_ReproCM"]
mk <- function(x, axis) data.frame(axis = axis, kind = x$kind, rho = x$rho_axis,
                                   lab = paste0(x$dataset, ": ", substr(x$contrast, 1, 34)))
G <- rbind(mk(U, "UVA reference axis"), mk(S, "senescence reference axis"))
G$grp <- ifelse(grepl("\\(axis\\)", G$kind), "defines the axis",
         ifelse(grepl("^injury", G$kind), "independent UV injury",
         ifelse(grepl("^photoprotection", G$kind), "independent photoprotection",
         ifelse(grepl("rescue within", G$kind), "rescue inside the axis study",
         ifelse(grepl("same study", G$kind), "same study as the axis", "independent senescence")))))
G$grp <- factor(G$grp, levels = c("defines the axis", "same study as the axis",
  "independent senescence", "independent photoprotection", "independent UV injury",
  "rescue inside the axis study"))
COL <- c("defines the axis" = "#444444", "same study as the axis" = GREY,
         "independent senescence" = GREEN, "independent photoprotection" = BLUE,
         "independent UV injury" = ORANGE, "rescue inside the axis study" = PURPLE)
G$axis <- factor(G$axis, levels = c("senescence reference axis", "UVA reference axis"))
anno <- data.frame(axis = factor(c("senescence reference axis", "UVA reference axis"),
                                 levels = levels(G$axis)), rho = c(lS, lU))
p1a <- ggplot(G, aes(rho, grp, colour = grp)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_point(size = 2.3, alpha = 0.9) +
  geom_vline(data = anno, aes(xintercept = rho), colour = ORANGE, linetype = 2, linewidth = 0.5) +
  geom_text(data = anno, aes(x = rho, y = 0.62, label = "Repro-CM anchor"),
            inherit.aes = FALSE, colour = ORANGE, size = 2.4, hjust = -0.06) +
  facet_wrap(~ axis, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-0.25, 0.92)) +
  coord_cartesian(clip = "off") +
  labs(x = "Spearman rho with the reference axis", y = NULL,
       title = "A. One axis transfers to independent studies. The other does not.",
       subtitle = "Every compendium contrast in each family, scored on the axis it belongs to",
       caption = paste0("Senescence: axis-defining contrasts median +0.727, four fully independent contrasts median +0.656 - almost no shrinkage.\n",
                        "UVA: axis-defining contrasts median +0.807, but six independent UV-injury contrasts median -0.017.\n",
                        "What loads on it instead is photoprotection (median +0.281; Wilcoxon vs injury p = 0.0087).\n",
                        "The one injury contrast that loads is the earliest timepoint, UVA at 6 h. Rescue contrasts from the axis study are\n",
                        "negative by construction. The Repro-CM anchor sits with photoprotection, not with injury.")) +
  th + theme(strip.background = element_rect(fill = "#F2F5F8", colour = "grey60"),
             strip.text = element_text(face = "bold", size = 8.6, colour = INK),
             panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.3))
ggsave(file.path(out, "Figure1A_axis_generalisation.png"), p1a,
       width = 172, height = 136, units = "mm", dpi = 300)

## ===================== FIGURE 1B =====================
tp <- data.frame(
  lab = c("UVA 6 h\n(single)", "UVA 24 h\n(single)", "UVA 24 h\n(repeated)", "UVA chronic",
          "UVB 1 h", "UVB 4 h"),
  rho = c(0.276, 0.017, -0.030, -0.005, -0.104, -0.089),
  study = c("GSE89005", "GSE89005", "GSE89005", "GSE125429", "GSE116968", "GSE116968"))
tp$lab <- factor(tp$lab, levels = tp$lab)
p1b <- ggplot(tp, aes(lab, rho, fill = study)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%+.3f", rho), vjust = ifelse(rho >= 0, -0.5, 1.4)),
            size = 2.6, colour = INK) +
  scale_fill_manual(values = c(GSE89005 = BLUE, GSE125429 = "#5B9BD5", GSE116968 = ORANGE)) +
  scale_y_continuous(limits = c(-0.17, 0.35)) +
  labs(x = NULL, y = "rho with the UVA reference axis",
       title = "B. Only the earliest injury timepoint loads on the axis",
       subtitle = "Independent UV-injury contrasts, ordered by exposure duration and waveband",
       caption = paste0("Accumulated and chronic UVA, and UVB at either timepoint, sit at or below zero.\n",
                        "The axis behaves as an early photo-adaptive response, not as a damage burden.")) +
  th + theme(panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
             axis.text.x = element_text(size = 7.6))
ggsave(file.path(out, "Figure1B_uv_timepoint_gradient.png"), p1b,
       width = 168, height = 82, units = "mm", dpi = 300)

## ===================== FIGURE 2A - QC =====================
bam <- file.path(der, "local_bam_v41")
qc <- do.call(rbind, lapply(c(HDF = "HDF", REP = "REP", IPS = "IPS", UVA0 = "UVA0"), function(s) {
  L <- readLines(file.path(bam, paste0(s, ".Log.final.out")), warn = FALSE)
  g <- function(p) as.numeric(gsub("%", "", trimws(strsplit(L[grepl(p, L, fixed = TRUE)][1], "|", fixed = TRUE)[[1]][2])))
  rpg <- read.delim(file.path(bam, paste0(s, ".ReadsPerGene.out.tab")), header = FALSE, skip = 4)
  data.frame(sample = c(HDF = "HDF-CM", REP = "Repro-CM", IPS = "iPSC-CM", UVA0 = "UVA 15J 0h")[s],
             input = g("Number of input reads"), unique = g("Uniquely mapped reads %"),
             mismatch = g("Mismatch rate per base"),
             assigned = 100 * sum(rpg[[4]]) / g("Number of input reads")) }))
write.table(qc, file.path(der, "local_bam_v41", "star_qc_four_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ql <- rbind(data.frame(sample = qc$sample, metric = "Uniquely mapped", pct = qc$unique),
            data.frame(sample = qc$sample, metric = "Assigned to genes (reverse strand)", pct = qc$assigned))
ql$sample <- factor(ql$sample, levels = rev(qc$sample))
p2a <- ggplot(ql, aes(pct, sample, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.7),
            hjust = -0.15, size = 2.4, colour = INK) +
  geom_vline(xintercept = 80, colour = GREY, linetype = 2, linewidth = 0.4) +
  scale_fill_manual(values = c("Uniquely mapped" = BLUE, "Assigned to genes (reverse strand)" = GREEN)) +
  scale_x_continuous(limits = c(0, 112), expand = c(0, 0)) +
  labs(x = "Percent of input read pairs", y = NULL,
       title = "A. All four local libraries pass the same alignment standard",
       subtitle = "GENCODE v41, STAR 2.7.3a, identical parameters; 14.9-16.8 M input read pairs per library",
       caption = "Mismatch rate 0.24-0.26% throughout. Technical quality is not the limitation here - biological replication is (n = 1 per condition).") +
  th + theme(panel.grid.major.y = element_blank())
ggsave(file.path(out, "Figure2A_local_qc.png"), p2a, width = 168, height = 78, units = "mm", dpi = 300)

## ===================== FIGURE 2B - construction =====================
AX <- read.delim(file.path(der, "decoupling_validation/primary_gene_axis_matrix.tsv"), check.names = FALSE)
AX$cls <- ifelse(!AX$comparator_consistent, "discordant (score set to 0)",
          ifelse(abs(AX$Repro_specific_score) >= log2(1.25), "consistent, |score| >= log2(1.25)", "consistent"))
AX$cls <- factor(AX$cls, levels = c("discordant (score set to 0)", "consistent", "consistent, |score| >= log2(1.25)"))
set.seed(1); AXs <- AX[sample(nrow(AX)), ]
p2b <- ggplot(AXs, aes(logFC_Repro_vs_HDF, logFC_Repro_vs_iPSC, colour = cls)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = 2, linewidth = 0.35) +
  geom_point(size = 0.7, stroke = 0, alpha = 0.55) +
  scale_colour_manual(values = c("discordant (score set to 0)" = "#D9DEE4",
                                 "consistent" = "#9CC3DE",
                                 "consistent, |score| >= log2(1.25)" = BLUE)) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), nrow = 2)) +
  labs(x = "log2 FC  Repro-CM vs HDF-CM", y = "log2 FC  Repro-CM vs iPSC-CM",
       title = "B. The anchor keeps only what both comparators agree on",
       subtitle = "score = sign(dH) x min(|dH|, |dI|), and 0 whenever the two comparators disagree",
       caption = paste0("8,921 of 12,319 expressed genes (72.4%) agree in direction; 2,204 reach |score| >= log2(1.25).\n",
                        "No p-value is attached to any gene: there is one library per condition.")) +
  th + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
ggsave(file.path(out, "Figure2B_anchor_construction.png"), p2b, width = 130, height = 132, units = "mm", dpi = 300)

## ===================== FIGURE 2C - robustness =====================
cs <- read.delim(file.path(der, "decoupling_validation/comparator_sensitivity.tsv"))
lo <- read.delim(file.path(der, "decoupling_validation/leave_one_dataset_out.tsv"))
gr <- read.delim(file.path(der, "decoupling_validation/gene_program_removal_sensitivity.tsv"))
th50 <- read.delim(file.path(der, "decoupling_validation/read_thinning_50pct_100runs.tsv"))
p4 <- read.delim(file.path(der, "direction_probe/P4_decoupling_proliferation_sensitivity.tsv"))
rb <- rbind(
  data.frame(grp = "comparator", lab = c("Repro vs HDF only", "Repro vs iPSC only", "conservative AND anchor"),
             D = cs$decoupling_index),
  data.frame(grp = "leave-one-study-out", lab = paste("drop", sub("^GSE", "GSE ", lo$dropped_dataset)),
             D = lo$decoupling_index),
  data.frame(grp = "gene-program removal", lab = c("none", "no RPL/RPS/MT", "no translation/rRNA/IFN",
             "no cell cycle/DNA repl.", "no translation/IFN/cell cycle", "no top 250 local", "no top 500 local"),
             D = gr$decoupling_index),
  data.frame(grp = "proliferation", lab = c("no adjustment", "all axes residualised",
             "proliferation-neutral genes", "proliferation-linked genes"), D = p4$decoupling_index_D),
  data.frame(grp = "read thinning", lab = sprintf("50%% thinning, run %d", seq_len(nrow(th50))),
             D = th50$decoupling_index))
rb$grp <- factor(rb$grp, levels = c("comparator", "leave-one-study-out", "gene-program removal",
                                    "proliferation", "read thinning"))
p2c <- ggplot(rb, aes(D, grp)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0.460, colour = ORANGE, linetype = 2, linewidth = 0.45) +
  geom_jitter(height = 0.16, width = 0, size = 1.6, colour = BLUE, alpha = 0.75, stroke = 0) +
  annotate("text", x = 0.460, y = 5.62, label = "primary D = 0.460", colour = ORANGE,
           size = 2.5, hjust = -0.06) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0, 0.72)) +
  coord_cartesian(clip = "off") +
  labs(x = "Decoupling index D", y = NULL,
       title = "C. D under every sensitivity analysis we ran",
       subtitle = "124 recomputations of the same index; none crosses zero",
       caption = paste0("The weakest value, D = 0.241, comes from using iPSC-CM as the only comparator. ",
                        "Read thinning is 100 binomial redraws at 50% depth\nand tests counting stability, not biological replication.")) +
  th
ggsave(file.path(out, "Figure2C_anchor_robustness.png"), p2c, width = 168, height = 84, units = "mm", dpi = 300)

## ===================== FIGURE 2D - specificity =====================
sp <- read.delim(file.path(der, "anchor_specificity", "anchor_specificity_negative_control.tsv"))
sp$lab <- paste0(sp$focal_condition, "\nvs ", sp$comparators)
sp$lab <- factor(sp$lab, levels = rev(sp$lab))
spl <- rbind(data.frame(lab = sp$lab, D = sp$D, w = "unadjusted"),
             data.frame(lab = sp$lab, D = sp$D_proliferation_adjusted, w = "proliferation-adjusted"))
p2d <- ggplot(spl, aes(D, lab, colour = w)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_point(size = 2.7, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sprintf("%+.3f", D)), position = position_dodge(width = 0.5),
            hjust = -0.28, size = 2.4, colour = INK) +
  scale_colour_manual(values = c(unadjusted = BLUE, `proliferation-adjusted` = GREEN)) +
  scale_x_continuous(limits = c(-0.55, 0.68)) +
  labs(x = "Decoupling index D", y = NULL,
       title = "D. The decoupling belongs to Repro-CM, not to the construction",
       subtitle = "The same signed-min consistency rule applied around each condition in turn",
       caption = "If D were an artefact of requiring agreement between two comparators, all three focal conditions would behave alike. They do not.") +
  th
ggsave(file.path(out, "Figure2D_anchor_specificity.png"), p2d, width = 168, height = 74, units = "mm", dpi = 300)
cat("Figures 1A/1B and 2A-2D written\n")
