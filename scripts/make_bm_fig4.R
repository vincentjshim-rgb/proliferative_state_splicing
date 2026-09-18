FIGURES_WRITTEN <- c("Fig4.png")
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
DS <- read.delim(file.path(D, "secretome_within_study.tsv"))
ev <- as.numeric(readLines(file.path(D, "mds_eigen.txt")))
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", secretome = REPRO,
        `photoprotection / rescue` = AMBER, `UV injury` = UV,
        `metabolic / culture` = GREY, reprogramming = NAVY, `Repro-CM` = "#0B3B3A")

## ============== A. secretome: every within-study pair, by study =============
sec <- which(M$class == "secretome")
s <- R[sec, sec]; ds <- M$dataset[sec]; sm <- outer(ds, ds, "==")
W <- data.frame(r = s[upper.tri(s) & sm],  k = "within a study")
B <- data.frame(r = s[upper.tri(s) & !sm], k = "between studies")
PD <- rbind(W, B); PD <- PD[is.finite(PD$r), ]
PD$k <- factor(PD$k, levels = c("between studies","within a study"))
wt <- wilcox.test(W$r, B$r)
p4a <- ggplot(PD, aes(r, k, colour = k)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.34, fill = NA, show.legend = FALSE) +
  geom_point(size = 1.3, alpha = 0.8, position = position_jitter(height = 0.12, seed = 2)) +
  brk(1, 2, 0.90, sprintf("%s   Wilcoxon P = %.4f", sig_stars(wt$p.value), wt$p.value),
      tick = 0.05, size = 5.7) +
  annotate("text", x = c(median(B$r), median(W$r)), y = c(1, 2) - 0.33,
           label = sprintf("%+.3f", c(median(B$r), median(W$r))), size = pt(6), family = FONT,
           colour = c(GREY_D <- INK2, REPRO)) +
  scale_colour_manual(values = c(`within a study` = REPRO, `between studies` = "#9FB6C4"),
                      guide = "none") +
  scale_x_continuous(limits = c(-0.42, 0.92), breaks = c(-0.25, 0, 0.25, 0.5, 0.75)) +
  coord_cartesian(ylim = c(0.55, 2.7)) +
  labs(x = "Spearman ρ between secretome contrasts", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.6))

## ================= B. per-study internal agreement ==========================
DS$lab <- sprintf("%s  (%d arms)", DS$dataset, DS$n)
DS$lab <- factor(DS$lab, levels = rev(DS$lab))
p4b <- ggplot(DS, aes(median_r, lab)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(aes(x = min, xend = max, yend = lab), colour = GREY_L, linewidth = 1.1) +
  geom_point(size = 2.0, colour = REPRO) +
  geom_text(aes(label = sprintf("%+.3f", median_r)), vjust = -1.2, size = pt(6),
            family = FONT, colour = INK) +
  geom_vline(xintercept = median(B$r), colour = "#9FB6C4", linewidth = 0.4, linetype = "22") +
  annotate("text", x = median(B$r), y = 0.62, label = "between-study median", colour = INK2,
           size = pt(5.6), family = FONT, hjust = -0.06) +
  scale_x_continuous(limits = c(0, 0.92), breaks = c(0, 0.25, 0.5, 0.75)) +
  coord_cartesian(ylim = c(0.5, 4.6), clip = "off") +
  labs(x = "Spearman ρ within one experiment (bar = range)", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.3))

## ============ C. MDS of all 76 contrasts ====================================
M$cl <- factor(M$class, levels = names(CC))
lab_ids <- c("LOCAL_ReproCM","GSE109700_deep","GSE179848_ContactInhibition",
             "GSE240226_UVA","SC3","SC8","GSE149694_t2iLGoYD13")
LB <- c(LOCAL_ReproCM="Repro-CM", GSE109700_deep="deep senescence",
        GSE179848_ContactInhibition="contact inhibition", GSE240226_UVA="acute UVA",
        SC3="BM-MSC CM", SC8="endothelial EVs", GSE149694_t2iLGoYD13="naive reprogramming")
M$lab <- LB[M$id]
p4c <- ggplot(M, aes(MDS1, MDS2)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = cl), size = 1.5) +
  geom_point(data = M[M$id %in% lab_ids, ], size = 2.5, shape = 21, colour = INK,
             fill = NA, stroke = 0.45) +
  geom_text(data = M[M$id %in% lab_ids, ],
            aes(label = lab,
                vjust = ifelse(id %in% c("LOCAL_ReproCM","SC3"), 2.0, -1.0),
                hjust = ifelse(id == "SC8", 1.12, ifelse(id == "LOCAL_ReproCM", 1.10, 0.5))),
            size = pt(5.7), family = FONT, colour = INK) +
  scale_colour_manual(values = CC, name = NULL) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.5))) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.7, family = FONT, size = pt(5.8),
           colour = GREY, label = "axis 1 tracks the senescence direction (|r| = 0.74)") +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  labs(x = sprintf("MDS axis 1  (%.1f%% of eigenvalue mass)", ev[1]),
       y = sprintf("MDS axis 2  (%.1f%%)", ev[2])) +
  theme_j(7) + theme(legend.position = "inside", legend.position.inside = c(0.25, 0.86),
                     legend.box.background = element_blank(), legend.key.size = unit(6, "pt"))

top <- lab_grid(p4a, p4b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1))
save_fig(lab_grid(top, p4c, labels = c("","C"), ncol = 1, rel_heights = c(0.72, 1)),
         "Fig4.png", 183, 126)
