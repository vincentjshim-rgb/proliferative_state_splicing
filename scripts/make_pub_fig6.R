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
PUBDIR <- "public_data_tierA/derived/figures_publication"
source("scripts/pub_theme.R"); D <- "public_data_tierA/derived"

## a  family 1 forest
P <- read.delim(file.path(D, "heldout_validation/heldout_primary_tests.tsv"))
P$lab <- c("HO1  OSK +dox, d4", "HO2  old vs young", "HO3a  red light, 1 h",
           "HO3b  red light, 4 h", "HO4a  fibroblast d7", "HO4b  primed d13")[match(P$id,
           c("HO1","HO2","HO3a","HO3b","HO4a","HO4b"))]
P$lab <- factor(P$lab, levels = rev(P$lab))
P$v <- factor(P$verdict, levels = c("SUCCESS", "BOUNDARY HELD", "NOT SUPPORTED"))
p6a <- ggplot(P, aes(rho, lab)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_linerange(aes(xmin = null_q025, xmax = null_q975), colour = "#D6DCE2", linewidth = 1.1) +
  geom_point(aes(colour = v), size = 1.5) +
  geom_text(aes(label = sprintf("%+.3f", rho)), vjust = -1.15, size = pt(5.2),
            family = FONT, colour = INK) +
  geom_text(aes(x = 0.20, label = sprintf("FDR %.3g", BH_FDR)), hjust = 0,
            size = pt(5), family = FONT, colour = GREY_D) +
  scale_colour_manual(values = c(SUCCESS = GREEN, `BOUNDARY HELD` = BLUE,
                                 `NOT SUPPORTED` = GREY), guide = "none") +
  scale_x_continuous(limits = c(-0.13, 0.33), breaks = c(-0.1, 0, 0.1)) +
  coord_cartesian(clip = "off") +
  labs(x = expression(paste("Spearman ", rho)), y = NULL) +
  theme_pub(7) + theme_cat_y()

## b  rescue replication across three studies
H6 <- read.delim(file.path(D, "heldout3_photoprotection/HO6_primary_tests.tsv"))
R <- data.frame(
  arm = c("Maifuyin", "succinate", "red light 1 h", "red light 4 h", "TAT-UIFSP5", "TAT-UIFSP6"),
  study = c("GSE240226", "GSE240226", "GSE116968", "GSE116968", "GSE222414", "GSE222414"),
  stage = c("discovery", "discovery", "held out", "held out", "held out", "held out"),
  rho = c(0.1271, 0.0935, P$rho[P$id == "HO3a"], P$rho[P$id == "HO3b"],
          H6$rho[H6$id == "HO6c"], H6$rho[H6$id == "HO6d"]))
R$lab <- factor(paste0(R$arm, "  ", R$study), levels = rev(paste0(R$arm, "  ", R$study)))
p6b <- ggplot(R, aes(rho, lab, colour = stage)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_text(aes(label = sprintf("%+.3f", rho)), hjust = -0.35, size = pt(5.2),
            family = FONT, colour = INK) +
  annotate("point", x = c(0.152, 0.152), y = c(5.4, 4.85), size = 1.2,
           colour = c(GREY, GREEN)) +
  annotate("text", x = 0.163, y = c(5.4, 4.85), hjust = 0, vjust = 0.42, size = pt(5),
           family = FONT, colour = GREY_D, label = c("discovery", "held out")) +
  scale_colour_manual(values = c(discovery = GREY, `held out` = GREEN), guide = "none") +
  scale_x_continuous(limits = c(0, 0.23), breaks = c(0, 0.05, 0.10)) +
  labs(x = expression(paste(rho, " after removing the injury component")), y = NULL) +
  theme_pub(7) + theme_cat_y()

## c  held-out senescence
sp <- read.delim(file.path(D, "heldout2_senescence/HO5_primary_tests.tsv"))
ss <- read.delim(file.path(D, "heldout2_senescence/HO5_sensitivity_tests.tsv"))
gv <- function(i) ss$rho[ss$id == i]
Cd <- rbind(
  data.frame(l = "pooled control", w = "unadjusted", r = sp$rho[sp$id == "HO5a"]),
  data.frame(l = "pooled control", w = "adjusted",  r = sp$rho[sp$id == "HO5b"]),
  data.frame(l = "pLenti",     w = "unadjusted", r = gv("HO5s:pLenti")),
  data.frame(l = "pLenti",     w = "adjusted",  r = gv("HO5s:pLenti:prolifAdj")),
  data.frame(l = "siScramble", w = "unadjusted", r = gv("HO5s:siScramble")),
  data.frame(l = "siScramble", w = "adjusted",  r = gv("HO5s:siScramble:prolifAdj")),
  data.frame(l = "IFIH1 KO",   w = "unadjusted", r = gv("HO5s:IFIH1_KO")),
  data.frame(l = "IFIH1 KO",   w = "adjusted",  r = gv("HO5s:IFIH1_KO:prolifAdj")),
  data.frame(l = "DDX58 KO",   w = "unadjusted", r = gv("HO5s:DDX58_KO")),
  data.frame(l = "DDX58 KO",   w = "adjusted",  r = gv("HO5s:DDX58_KO:prolifAdj")))
Cd$l <- factor(Cd$l, levels = rev(unique(Cd$l)))
p6c <- ggplot(Cd, aes(r, l, colour = w)) +
  annotate("rect", xmin = -0.163, xmax = -0.044, ymin = -Inf, ymax = Inf,
           fill = "#EDF1F4", colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.25) +
  geom_point(size = 1.4, position = position_dodge(width = 0.5)) +
  annotate("text", x = -0.104, y = 5.62, label = "discovery range", size = pt(5),
           family = FONT, colour = GREY_D, vjust = 0.3) +
  annotate("point", x = c(-0.052, -0.052), y = c(1.55, 1.22), size = 1.2,
           colour = c(SLATE, GREEN)) +
  annotate("text", x = -0.045, y = c(1.55, 1.22), hjust = 0, vjust = 0.42, size = pt(5),
           family = FONT, colour = GREY_D, label = c("unadjusted", "proliferation-adjusted")) +
  scale_colour_manual(values = c(unadjusted = SLATE, adjusted = GREEN), guide = "none") +
  scale_x_continuous(limits = c(-0.19, 0.075), breaks = c(-0.15, -0.10, -0.05, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression(paste(rho, "  anchor vs senescent − proliferating")), y = NULL) +
  theme_pub(7) + theme_cat_y()

## d  UVA axis does not load independent injury
nw <- read.delim(file.path(D, "heldout3_photoprotection/HO67_axis_loadings.tsv"))
E <- data.frame(
  lab = c("GSE240226 acute UVA", "GSE302943 cumulative UVA", "GSE89005 UVA 6 h",
          "GSE89005 UVA 24 h", "GSE125429 UVA chronic", "GSE89005 UVA 24 h rep.",
          "GSE240486 UV", "GSE329475 UVA 3 d", "GSE116968 UVB 4 h", "GSE116968 UVB 1 h",
          "GSE222414 UVB"),
  rho = c(0.816, 0.798, 0.276, 0.017, -0.005, -0.030, nw$rho_UVA[nw$id == "ax1"],
          nw$rho_UVA[nw$id == "HO7"], -0.089, -0.104, nw$rho_UVA[nw$id == "ax2"]),
  g = c("defines axis", "defines axis", rep("independent injury", 5), "prespecified HO7",
        rep("independent injury", 3)))
E <- E[order(E$rho), ]; E$lab <- factor(E$lab, levels = E$lab)
p6d <- ggplot(E, aes(rho, lab, colour = g)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), colour = GREY_L, linewidth = 0.45) +
  geom_point(size = 1.4) +
  geom_text(aes(label = sprintf("%+.3f", rho), hjust = ifelse(rho >= 0, -0.3, 1.3)),
            size = pt(5), family = FONT, colour = GREY_D) +
  scale_colour_manual(values = c(`defines axis` = INK, `independent injury` = GREY,
                                 `prespecified HO7` = VERM), guide = "none") +
  scale_x_continuous(limits = c(-0.42, 1.02), breaks = c(-0.25, 0, 0.5, 1)) +
  labs(x = expression(paste(rho, " with the UVA reference axis")), y = NULL) +
  theme_pub(7) + theme_cat_y() + theme(axis.text.y = element_text(size = 6))

top <- compose(p6a, p6b, labels = c("a", "b"), ncol = 2, rel_widths = c(1.05, 1))
bot <- compose(p6c, p6d, labels = c("c", "d"), ncol = 2, rel_widths = c(1, 1.12))
save_fig(plot_grid(top, bot, ncol = 1, rel_heights = c(1, 1.15)), "Fig6.png", 183, 112)
cat("Fig6 saved\n")
