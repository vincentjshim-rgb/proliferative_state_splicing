#!/usr/bin/env Rscript
# Figure 6 : prespecified held-out generalization (preregistration SHA256 41e4cf4b...)
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "public_data_tierA"
hv   <- file.path(root, "derived", "heldout_validation")
out  <- file.path(root, "derived", "manuscript_figures_v4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"
GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"
base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.4),
        plot.caption = element_text(colour = MUTED, size = 7.2, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8.2))

P <- read.delim(file.path(hv, "heldout_primary_tests.tsv"))
P$panel_label <- c(
  HO1  = "HO1  GSE297233\nOSK +dox vs -dox, day 4",
  HO2  = "HO2  GSE307377\nold vs young fibroblast",
  HO3a = "HO3a  GSE116968\nred light vs UV, 1 h",
  HO3b = "HO3b  GSE116968\nred light vs UV, 4 h",
  HO4a = "HO4a  GSE149694\nfibroblast day 7 vs day 3",
  HO4b = "HO4b  GSE149694\nprimed day 13 vs day 3")[P$id]
P$pred <- c(HO1 = "null predicted", HO2 = "negative predicted", HO3a = "positive predicted",
            HO3b = "positive predicted", HO4a = "null predicted", HO4b = "null predicted")[P$id]
P$panel_label <- factor(P$panel_label, levels = rev(P$panel_label))
vc <- c(SUCCESS = GREEN, `BOUNDARY HELD` = BLUE, `NOT SUPPORTED` = GREY,
        `REFUTES MODEL` = ORANGE, `REFUTES BOUNDARY` = ORANGE)

pA <- ggplot(P, aes(rho, panel_label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = null_q025, xmax = null_q975), height = 0.15,
                 colour = "grey78", linewidth = 0.7) +
  geom_point(aes(colour = verdict), size = 2.8) +
  geom_text(aes(label = sprintf("%+.3f", rho)), hjust = 0.5, nudge_y = 0.30,
            size = 2.6, colour = INK) +
  geom_text(aes(x = -0.185, label = pred), hjust = 0, size = 2.25, colour = MUTED) +
  geom_text(aes(x = 0.175, label = sprintf("FDR %.3g", BH_FDR)), hjust = 0,
            size = 2.4, colour = MUTED) +
  geom_text(aes(x = 0.255, label = verdict, colour = verdict), hjust = 0, size = 2.4,
            fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = vc) +
  scale_x_continuous(limits = c(-0.19, 0.46),
                     breaks = c(-0.1, 0, 0.1), labels = c("-0.1", "0", "0.1")) +
  labs(x = "Spearman rho (anchor vs held-out contrast)", y = NULL,
       title = "A. Prespecified held-out tests",
       subtitle = "Contrasts, predictions and success criteria fixed before any result was seen",
       caption = paste0("Grey bars: expression-matched permutation 95% NULL intervals, not confidence intervals.\n",
                        "BH correction across these 6 primary tests. HO1/HO4 predicted no positive alignment,\n",
                        "so a null there is consistent with the boundary but is not positive evidence (n = 2-3 per arm).")) +
  base_theme + theme(legend.position = "none",
                     panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3))
ggsave(file.path(out, "Figure6A_heldout_forest.png"), pA, width = 180, height = 104, units = "mm", dpi = 300)

## -- B: rescue-aligned component replicates in an independent study ---------
H6 <- read.delim(file.path(root, "derived", "heldout3_photoprotection", "HO6_primary_tests.tsv"))
R <- data.frame(
  study = c("GSE240226", "GSE240226", "GSE116968", "GSE116968", "GSE222414", "GSE222414"),
  arm   = c("Maifuyin extract", "succinate", "red light, 1 h", "red light, 4 h",
            "TAT-UIFSP5", "TAT-UIFSP6"),
  stage = c("discovery", "discovery", "held-out 1", "held-out 1", "held-out 3", "held-out 3"),
  injury = c("UVA", "UVA", "UVB", "UVB", "UVB", "UVB"),
  rho   = c(0.127086855854731, 0.0935220717798359, P$rho[P$id == "HO3a"], P$rho[P$id == "HO3b"],
            H6$rho[H6$id == "HO6c"], H6$rho[H6$id == "HO6d"]))
R$lab <- paste0(R$arm, "\n", R$study, " (", R$injury, ")")
R$lab <- factor(R$lab, levels = rev(R$lab))
pB <- ggplot(R, aes(rho, lab, colour = stage)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), linewidth = 0.9, colour = "grey82") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%+.3f", rho)), hjust = -0.5, size = 2.7, colour = INK) +
  scale_colour_manual(values = c(discovery = GREY, `held-out 1` = GREEN, `held-out 3` = BLUE)) +
  scale_x_continuous(limits = c(0, 0.178), expand = c(0, 0)) +
  labs(x = "Spearman rho after removing the injury component", y = NULL,
       title = "B. The repair-aligned component replicates independently",
       subtitle = "Six photoprotective interventions, three independent studies, two UV wavebands",
       caption = paste0("Mechanisms are unrelated: a herbal extract, a metabolite, photobiomodulation and a peptide-fusion protein.\n",
                        "In GSE116968 the anchor is NEGATIVELY correlated with the UVB injury axis itself (rho = -0.094, -0.169),\n",
                        "so the positive rescue association is not residual injury. A fourth study (GSE240486) did not replicate;\n",
                        "its own UV injury contrast was inert (rho = -0.030, p = 0.59), and it is shown in the held-out results table.")) +
  base_theme
ggsave(file.path(out, "Figure6B_rescue_replication.png"), pB, width = 162, height = 104, units = "mm", dpi = 300)

## -- C: reprogramming boundary, primary + all sensitivity arms -------------
S <- read.delim(file.path(hv, "heldout_sensitivity_tests.tsv"))
S <- S[grepl("^HO4b:", S$id) | S$id == "HO1s", ]
B <- rbind(
  data.frame(lab = "OSK +dox (primary)",        rho = P$rho[P$id == "HO1"],  p = P$p_matched[P$id == "HO1"],  type = "primary"),
  data.frame(lab = "O4YRSK mutant",             rho = S$rho[S$id == "HO1s"], p = S$p_matched[S$id == "HO1s"], type = "sensitivity"),
  data.frame(lab = "Primed D13 (primary)",      rho = P$rho[P$id == "HO4b"], p = P$p_matched[P$id == "HO4b"], type = "primary"),
  data.frame(lab = sub("^HO4b:", "", S$id[grepl("^HO4b:", S$id)]),
             rho = S$rho[grepl("^HO4b:", S$id)], p = S$p_matched[grepl("^HO4b:", S$id)], type = "sensitivity"))
B$lab <- factor(B$lab, levels = rev(B$lab))
B$flag <- ifelse(B$rho > 0.15 & B$p < 0.05, "exceeds boundary", B$type)
pC <- ggplot(B, aes(rho, lab)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = c(-0.15, 0.15), colour = "grey80", linetype = 2, linewidth = 0.35) +
  geom_point(aes(colour = flag, shape = type), size = 2.6) +
  geom_text(aes(label = sprintf("%+.3f%s", rho, ifelse(p < 0.05, " *", ""))),
            hjust = -0.35, size = 2.5, colour = INK) +
  scale_colour_manual(values = c(primary = BLUE, sensitivity = GREY, `exceeds boundary` = ORANGE),
                      breaks = "exceeds boundary") +
  scale_shape_manual(values = c(primary = 16, sensitivity = 17)) +
  scale_x_continuous(limits = c(-0.17, 0.30)) +
  labs(x = "Spearman rho vs cell-intrinsic reprogramming contrast", y = NULL,
       title = "C. Not a partial-reprogramming mimic",
       subtitle = "Dashed lines mark the prespecified |rho| = 0.15 boundary; * = expression-matched p < 0.05",
       caption = paste0("Both primary contrasts stay inside the boundary. One naive-state sensitivity arm\n",
                        "(t2iLGoY day 13) exceeds it and is reported here rather than omitted; by the\n",
                        "preregistration, sensitivity arms do not change the verdict.")) +
  base_theme
ggsave(file.path(out, "Figure6C_reprogramming_boundary.png"), pC, width = 160, height = 98, units = "mm", dpi = 300)
cat("Figure 6 panels written.\n")

## -- E: HO7, does the UVA axis load independent UV injury? ------------------
AXG <- read.delim(file.path(root, "derived", "axis_characterisation", "UVA_axis_generalisation.tsv"))
new <- read.delim(file.path(root, "derived", "heldout3_photoprotection", "HO67_axis_loadings.tsv"))
E <- rbind(
  data.frame(lab = "GSE240226  acute UVA", rho = 0.816, grp = "defines the axis"),
  data.frame(lab = "GSE302943  cumulative UVA", rho = 0.798, grp = "defines the axis"),
  data.frame(lab = "GSE89005  UVA 6 h", rho = 0.276, grp = "independent injury"),
  data.frame(lab = "GSE89005  UVA 24 h", rho = 0.017, grp = "independent injury"),
  data.frame(lab = "GSE125429  UVA chronic", rho = -0.005, grp = "independent injury"),
  data.frame(lab = "GSE89005  UVA 24 h repeated", rho = -0.030, grp = "independent injury"),
  data.frame(lab = "GSE240486  UV  (held-out 3)", rho = new$rho_UVA[new$id == "ax1"], grp = "independent injury"),
  data.frame(lab = "GSE329475  UVA 3 days  (HO7)", rho = new$rho_UVA[new$id == "HO7"], grp = "prespecified test"),
  data.frame(lab = "GSE116968  UVB 4 h", rho = -0.089, grp = "independent injury"),
  data.frame(lab = "GSE116968  UVB 1 h", rho = -0.104, grp = "independent injury"),
  data.frame(lab = "GSE222414  UVB  (held-out 3)", rho = new$rho_UVA[new$id == "ax2"], grp = "independent injury"))
E <- E[order(E$rho), ]; E$lab <- factor(E$lab, levels = E$lab)
pE <- ggplot(E, aes(rho, lab, colour = grp)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = c(-0.15, 0.15), colour = "grey82", linetype = 3, linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), colour = "grey85", linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%+.3f", rho), hjust = ifelse(rho >= 0, -0.35, 1.35)),
            size = 2.4, colour = INK) +
  scale_colour_manual(values = c(`defines the axis` = "#444444",
                                 `independent injury` = GREY, `prespecified test` = ORANGE)) +
  scale_x_continuous(limits = c(-0.40, 0.98)) +
  labs(x = "rho with the UVA reference axis", y = NULL,
       title = "E. The UVA axis does not load independent UV injury",
       subtitle = "Prespecified test HO7 put the post-hoc reading of the axis at risk of refutation",
       caption = paste0("HO7 fixed the prediction |rho| < 0.15 before the data were opened,\n",
                        "with rho > 0.40 set as the refuting value. A three-day UVA exposure in\n",
                        "dermal fibroblasts returned -0.017. Two further independent injury contrasts\n",
                        "from the same download point the same way, bringing the independent\n",
                        "UV-injury set to eight, none of which loads on the axis.")) +
  base_theme
ggsave(file.path(out, "Figure6E_uva_axis_heldout.png"), pE,
       width = 168, height = 116, units = "mm", dpi = 300)
cat("Figure 6E written\n")
