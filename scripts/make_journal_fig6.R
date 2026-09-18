FIGURES_WRITTEN <- c("Fig6.png")
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
fzi <- function(r, n) { se <- 1/sqrt(n-3); c(tanh(atanh(r)-1.96*se), tanh(atanh(r)+1.96*se)) }

## ================ A. all preregistered primary tests, by family ============
H <- read.delim(file.path(CS, "heldout_primary_fisherCI.tsv"))
NM <- c(HO1="HO1  OSK +dox, day 4", HO2="HO2  old vs young donor",
        HO3a="HO3a  red light, 1 h", HO3b="HO3b  red light, 4 h",
        HO4a="HO4a  fibroblast day 7", HO4b="HO4b  primed day 13",
        HO5a="HO5a  X-ray senescence", HO5b="HO5b  same, proliferation-adjusted",
        HO6a="HO6a  osmolyte 0.02%", HO6b="HO6b  osmolyte 0.06%",
        HO6c="HO6c  TAT-UIFSP5", HO6d="HO6d  TAT-UIFSP6", HO8="HO8  sauchinone")
H <- H[H$id %in% names(NM), ]; H$nm <- NM[H$id]
FAM <- c("1  photoprotection, age,\n    reprogramming boundary", "2  held-out senescence",
         "3  photoprotection replication", "4  third photoprotection study")
rows <- list(); y <- 0
for (f in sort(unique(H$family))) {
  y <- y - 1
  rows[[length(rows)+1]] <- data.frame(y, kind="header", nm=FAM[f], est=NA, lo=NA, hi=NA,
                                       verdict=NA, fdr=NA)
  s <- H[H$family == f, ]
  for (i in seq_len(nrow(s))) { y <- y - 1
    rows[[length(rows)+1]] <- data.frame(y, kind="test", nm=s$nm[i], est=s$rho[i],
      lo=s$ci_lo[i], hi=s$ci_hi[i], verdict=s$verdict[i],
      fdr=ifelse(is.na(s$BH_FDR[i]), sprintf("P = %.0e", 1e-4),
                 ifelse(s$BH_FDR[i] < 1e-3, sprintf("FDR = %.0e", s$BH_FDR[i]),
                        sprintf("FDR = %.3f", s$BH_FDR[i])))) }
  y <- y - 0.4
}
F <- do.call(rbind, rows)
VC <- c(SUCCESS = REPRO, `BOUNDARY HELD` = NAVY, `NOT SUPPORTED` = GREY)
XL <- c(-0.95, 0.78); AXB <- c(-0.1, 0, 0.1, 0.2, 0.3, 0.4); ybot <- min(F$y) - 1.0
p6a <- ggplot() +
  annotate("rect", xmin = -0.15, xmax = 0.15, ymin = ybot, ymax = 0.4,
           fill = "#F2F4F6", colour = NA) +
  annotate("segment", x = 0, xend = 0, y = ybot, yend = 0.4, colour = INK, linewidth = 0.3) +
  geom_errorbarh(data = F[F$kind=="test",], aes(xmin = lo, xmax = hi, y = y),
                 height = 0.16, colour = GREY, linewidth = 0.33) +
  geom_point(data = F[F$kind=="test",], aes(x = est, y = y, colour = verdict),
             size = 1.7, shape = 15) +
  geom_text(data = F[F$kind=="test",], aes(x = XL[1], y = y, label = nm), hjust = 0,
            size = pt(6), family = FONT, colour = INK2) +
  geom_text(data = F[F$kind=="header",], aes(x = XL[1], y = y, label = paste("Family", nm)),
            hjust = 0, size = pt(6.2), family = FONT, fontface = "bold", colour = INK,
            lineheight = 0.92) +
  geom_text(data = F[F$kind=="test",], aes(x = XL[2], y = y, label = fdr), hjust = 1,
            size = pt(5.6), family = FONT, colour = GREY) +
  scale_colour_manual(values = VC, name = NULL) +
  annotate("text", x = 0, y = 0.75, label = "prespecified null zone  |ρ| < 0.15",
           size = pt(5.5), family = FONT, colour = GREY) +
  annotate("segment", x = min(AXB)-0.02, xend = max(AXB)+0.02, y = ybot, yend = ybot,
           colour = INK, linewidth = 0.3) +
  annotate("segment", x = AXB, xend = AXB, y = ybot, yend = ybot - 0.18, colour = INK, linewidth = 0.3) +
  annotate("text", x = AXB, y = ybot - 0.34, label = sprintf("%.1f", AXB), vjust = 1,
           size = pt(6.5), family = FONT, colour = INK) +
  annotate("text", x = mean(range(AXB)), y = ybot - 1.3, family = FONT, size = pt(7), colour = INK,
           label = "Spearman ρ  (95% CI)") +
  scale_x_continuous(limits = XL) +
  coord_cartesian(ylim = c(ybot - 1.7, 0.95), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_j(7) + theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_blank(), axis.text.y = element_blank(),
                     legend.position = "bottom", legend.justification = "center",
                     legend.box.background = element_blank(),
                     legend.direction = "horizontal", legend.key.size = unit(6, "pt"),
                     legend.spacing.x = unit(9, "pt"), legend.text = element_text(size = 6.2),
                     legend.margin = margin(t = -2))

## ============ B. repair-aligned component, six interventions ==============
HO  <- read.delim(file.path(D, "heldout_validation/heldout_primary_tests.tsv"))
H6  <- read.delim(file.path(D, "heldout3_photoprotection/HO6_primary_tests.tsv"))
R <- data.frame(
  arm   = c("Maifuyin", "succinate", "red light, 1 h", "red light, 4 h", "TAT-UIFSP5", "TAT-UIFSP6"),
  study = c("GSE240226","GSE240226","GSE116968","GSE116968","GSE222414","GSE222414"),
  stage = c("discovery","discovery","held out","held out","held out","held out"),
  rho   = c(0.1271, 0.0935, HO$rho[HO$id=="HO3a"], HO$rho[HO$id=="HO3b"],
            H6$rho[H6$id=="HO6c"], H6$rho[H6$id=="HO6d"]),
  n     = c(2237, 2237, HO$n_shared[HO$id=="HO3a"], HO$n_shared[HO$id=="HO3b"],
            H6$n_shared[H6$id=="HO6c"], H6$n_shared[H6$id=="HO6d"]))
ci <- t(mapply(fzi, R$rho, R$n)); R$lo <- ci[,1]; R$hi <- ci[,2]
R$lab <- factor(sprintf("%s   %s", R$arm, R$study),
                levels = rev(sprintf("%s   %s", R$arm, R$study)))
p6b <- ggplot(R, aes(rho, lab, colour = stage)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16, linewidth = 0.33) +
  geom_point(size = 1.7, shape = 15) +
  geom_text(aes(label = sprintf("%+.3f", rho)), vjust = -1.2, size = pt(6),
            family = FONT, colour = INK) +
  scale_colour_manual(values = c(discovery = GREY, `held out` = REPRO), name = NULL) +
  scale_x_continuous(limits = c(-0.02, 0.20), breaks = c(0, 0.05, 0.1, 0.15)) +
  labs(x = "ρ after removing the injury component  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_text(size = 6.2),
        legend.position = "inside", legend.position.inside = c(0.84, 0.18),
        legend.box.background = element_blank(), legend.key.size = unit(6, "pt"))

## ================= C. held-out senescence, four arms =======================
sp <- read.delim(file.path(D, "heldout2_senescence/HO5_primary_tests.tsv"))
ss <- read.delim(file.path(D, "heldout2_senescence/HO5_sensitivity_tests.tsv"))
gv <- function(i, col = "rho") ss[[col]][ss$id == i]
Cd <- data.frame(
  l = rep(c("pooled control","pLenti","siScramble","IFIH1 KO","DDX58 KO"), each = 2),
  w = rep(c("unadjusted","proliferation-adjusted"), 5),
  r = c(sp$rho[sp$id=="HO5a"], sp$rho[sp$id=="HO5b"],
        gv("HO5s:pLenti"), gv("HO5s:pLenti:prolifAdj"),
        gv("HO5s:siScramble"), gv("HO5s:siScramble:prolifAdj"),
        gv("HO5s:IFIH1_KO"), gv("HO5s:IFIH1_KO:prolifAdj"),
        gv("HO5s:DDX58_KO"), gv("HO5s:DDX58_KO:prolifAdj")))
ci <- t(sapply(Cd$r, fzi, n = 2215)); Cd$lo <- ci[,1]; Cd$hi <- ci[,2]
Cd$l <- factor(Cd$l, levels = rev(unique(Cd$l)))
Cd$w <- factor(Cd$w, levels = c("unadjusted","proliferation-adjusted"))
p6c <- ggplot(Cd, aes(r, l, colour = w)) +
  annotate("rect", xmin = -0.163, xmax = -0.044, ymin = -Inf, ymax = Inf,
           fill = "#F0F2F4", colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.33,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 1.6, shape = 15, position = position_dodge(width = 0.55)) +
  annotate("text", x = -0.104, y = 5.68, label = "range of the four discovery references",
           size = pt(5.5), family = FONT, colour = GREY) +
  scale_colour_manual(values = c(unadjusted = SEN, `proliferation-adjusted` = NAVY), name = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(-0.23, 0.06), breaks = c(-0.2, -0.1, 0)) +
  coord_cartesian(ylim = c(0.5, 6.0)) +
  labs(x = "ρ  anchor vs senescent − proliferating  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box.background = element_blank(), legend.margin = margin(t = -4),
        legend.key.size = unit(6, "pt"))

## =========== D. every UV contrast on the UVA reference axis ================
nw <- read.delim(file.path(D, "heldout3_photoprotection/HO67_axis_loadings.tsv"))
E <- data.frame(
  lab = c("GSE240226  acute UVA","GSE302943  cumulative UVA","GSE89005  UVA 6 h",
          "GSE89005  UVA 24 h","GSE125429  UVA chronic","GSE89005  UVA 24 h repeated",
          "GSE240486  UV","GSE329475  UVA 3 days","GSE116968  UVB 4 h","GSE116968  UVB 1 h",
          "GSE222414  UVB"),
  rho = c(0.816, 0.798, 0.276, 0.017, -0.005, -0.030, nw$rho_UVA[nw$id=="ax1"],
          nw$rho_UVA[nw$id=="HO7"], -0.089, -0.104, nw$rho_UVA[nw$id=="ax2"]),
  g   = c("defines the axis","defines the axis", rep("independent injury", 5),
          "preregistered HO7", rep("independent injury", 3)))
E$g[E$lab == "GSE329475  UVA 3 days"] <- "preregistered HO7"
E$g[E$lab == "GSE240486  UV"] <- "independent injury"
E <- E[order(E$rho), ]; E$lab <- factor(E$lab, levels = E$lab)
p6d <- ggplot(E, aes(rho, lab, colour = g)) +
  annotate("rect", xmin = -0.15, xmax = 0.15, ymin = -Inf, ymax = Inf, fill = "#F2F4F6", colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), linewidth = 0.4) +
  geom_point(size = 1.6, shape = 15) +
  geom_text(aes(label = sprintf("%+.3f", rho), hjust = ifelse(rho >= 0, -0.28, 1.28)),
            size = pt(5.8), family = FONT, colour = INK2) +
  scale_colour_manual(values = c(`defines the axis` = UV, `independent injury` = GREY,
                                 `preregistered HO7` = AMBER), name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_continuous(limits = c(-0.46, 1.06), breaks = c(-0.25, 0, 0.5, 1)) +
  labs(x = "ρ with the UVA reference axis", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_text(size = 6.1),
        legend.position = "bottom", legend.justification = "center",
        legend.box.background = element_blank(), legend.margin = margin(t = -4),
        legend.key.size = unit(6, "pt"))

right <- lab_grid(p6b, p6c, p6d, labels = c("B","C","D"), ncol = 1, rel_heights = c(1, 1.05, 1.3))
save_fig(lab_grid(p6a, right, labels = c("A",""), ncol = 2, rel_widths = c(1.02, 1)),
         "Fig6.png", 183, 150)
