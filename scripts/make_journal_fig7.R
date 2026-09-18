FIGURES_WRITTEN <- c("Fig7.png")
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
CP <- read.delim(file.path(CS, "compendium_with_CI.tsv"))
cls <- function(id, con) {
  if (grepl("senesc|SIPS|_deep|_early|P27", paste(id, con), ignore.case = TRUE)) "senescence"
  else if (grepl("age|old_vs_young|donor age", con, ignore.case = TRUE)) "donor age"
  else if (grepl("rescue|red_vs|prered|Osmoter|sauchinone|UIFSP|Maifuyin|succinate",
                 paste(id, con), ignore.case = TRUE)) "photoprotection"
  else if (grepl("UV", paste(id, con))) "UV injury"
  else if (grepl("D13|D7|dox|OSK|MPTR|reprogram", paste(id, con), ignore.case = TRUE)) "reprogramming"
  else "other perturbation" }
CP$group <- mapply(cls, CP$id, CP$contrast)
CP$group[CP$id == "LOCAL_ReproCM"] <- "Repro-CM"
CP$indep <- CP$circularity == "independent"
GC <- c(senescence = SEN, `donor age` = "#8C6D4F", `UV injury` = UV,
        photoprotection = AMBER, reprogramming = NAVY, `other perturbation` = GREY,
        `Repro-CM` = REPRO)

## =================== A. D by perturbation class ============================
ord <- names(sort(tapply(CP$D, CP$group, median)))
CP$grp <- factor(CP$group, levels = ord)
## the singleton Repro-CM group is excluded from the omnibus test
kw <- kruskal.test(D ~ droplevels(grp), CP[CP$group != "Repro-CM", ])
wt <- wilcox.test(CP$D[CP$group == "senescence"], CP$D[!CP$group %in% c("senescence","Repro-CM")])
p7a <- ggplot(CP, aes(grp, D)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(aes(colour = grp), outlier.shape = NA, width = 0.55, linewidth = 0.34, fill = NA) +
  geom_point(aes(colour = grp, shape = indep), size = 1.15,
             position = position_jitter(width = 0.16, seed = 3)) +
  geom_hline(yintercept = CP$D[CP$id == "LOCAL_ReproCM"], colour = REPRO,
             linewidth = 0.35, linetype = "22") +
  annotate("text", x = 0.6, y = CP$D[CP$id == "LOCAL_ReproCM"], label = "Repro-CM",
           colour = REPRO, size = pt(5.8), family = FONT, hjust = 0, vjust = -0.6) +
  brk(1, 3.5, 0.60, sprintf("%s   senescence vs rest, Wilcoxon P = %.1e",
      sig_stars(wt$p.value), wt$p.value), tick = 0.05, size = 5.6) +
  annotate("text", x = 6.4, y = 0.80, hjust = 1, family = FONT, size = pt(5.8), colour = INK2,
           label = sprintf("Kruskal-Wallis  χ²(%d) = %.1f,  P = %.1e",
                           kw$parameter, kw$statistic, kw$p.value)) +
  scale_colour_manual(values = GC, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = NULL,
                     labels = c(`TRUE` = "independent of both axes",
                                `FALSE` = "contributes to an axis")) +
  scale_y_continuous(limits = c(-0.88, 0.86), breaks = seq(-0.75, 0.75, 0.25)) +
  labs(x = NULL, y = "decoupling index  D") +
  theme_j(7) + theme(axis.text.x = element_text(angle = 22, hjust = 1, size = 6.3),
                     legend.position = "inside", legend.position.inside = c(0.30, 0.10),
                     legend.box.background = element_blank(), legend.key.size = unit(6, "pt"),
                     legend.direction = "vertical", legend.spacing.y = unit(0, "pt"))

## ================ B. the two axes, with D as distance from identity ========
lab_pts <- c("LOCAL_ReproCM","GSE179848_ContactInhibition","GSE326951_sauchinone",
             "GSE109700_deep","GSE240226_UVA","GSE165177_MPTR","GSE149694_t2iLGoYD13")
CP$show <- CP$id %in% lab_pts
NMS <- c(LOCAL_ReproCM = "Repro-CM", GSE179848_ContactInhibition = "contact inhibition",
         GSE326951_sauchinone = "sauchinone", GSE109700_deep = "deep senescence",
         GSE240226_UVA = "UVA (defines axis)", GSE165177_MPTR = "MPTR",
         GSE149694_t2iLGoYD13 = "naive reprogramming")
CP$nm <- NMS[CP$id]
CP$hj <- ifelse(CP$rho_senescence > 0.3, 1.12, -0.12)
p7b <- ggplot(CP, aes(rho_senescence, rho_UVA)) +
  geom_abline(slope = 1, intercept = 0, colour = GREY_L, linewidth = 0.3, linetype = "22") +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = grp, shape = indep), size = 1.5) +
  geom_point(data = CP[CP$show, ], size = 2.3, shape = 21, colour = INK, fill = NA, stroke = 0.4) +
  geom_text(data = CP[CP$show, ], aes(label = nm, hjust = hj), size = pt(5.7),
            family = FONT, colour = INK, vjust = -0.75) +
  annotate("text", x = 1.0, y = -0.58, hjust = 1, family = FONT, size = pt(5.9), colour = GREY,
           label = "D = vertical distance below the identity line") +
  scale_colour_manual(values = GC, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.4))) +
  coord_cartesian(xlim = c(-0.70, 1.02), ylim = c(-0.62, 1.02)) +
  labs(x = "ρ with the senescence axis", y = "ρ with the UVA axis") +
  theme_j(7) + theme(legend.position = "inside", legend.position.inside = c(0.155, 0.80),
                     legend.box.background = element_blank(), legend.key.size = unit(6, "pt"))

## ==================== C. ranked, independent perturbations =================
IN <- CP[CP$indep, ]; IN <- IN[order(-IN$D), ]
CLEAN <- c(FibroblastD7 = "fibroblast day 7", NHSMD13 = "NHSM medium, day 13",
  O4YRSK = "O4YRSK +dox, day 4", OSK = "OSK +dox, day 4", RSeTD13 = "RSeT medium, day 13",
  `5iLAFD13` = "5iLAF medium, day 13", PrimedD13 = "primed medium, day 13",
  UVBinjury = "UVB injury", age = "donor age, per decade", MPTR = "MPTR reprogramming",
  DDX58_KO = "DDX58 KO senescence", IFIH1_KO = "IFIH1 KO senescence",
  siScramble = "siScramble senescence", pLenti = "pLenti senescence",
  `uv vs control 4h` = "UVB injury, 4 h", `prered vs uv 1h` = "red light rescue, 1 h",
  `red vs control 4h` = "red light alone, 4 h")
raw <- ifelse(is.na(IN$nm), gsub("_", " ", sub("^GSE[0-9]+_", "", IN$id)), IN$nm)
key <- gsub("_", " ", sub("^GSE[0-9]+_", "", IN$id))
names(CLEAN) <- gsub("_", " ", names(CLEAN))
IN$nm2 <- ifelse(key %in% names(CLEAN), CLEAN[key], raw)
IN$nm2 <- substr(IN$nm2, 1, 26)
sel <- rbind(head(IN, 12), tail(IN, 8)); sel$nm2 <- factor(sel$nm2, levels = rev(sel$nm2))
p7c <- ggplot(sel, aes(D, nm2, colour = grp)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = D - 1.96*D_se, xmax = D + 1.96*D_se), height = 0, linewidth = 0.33) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = GC, guide = "none") +
  scale_x_continuous(limits = c(-0.86, 0.66), breaks = seq(-0.75, 0.5, 0.25)) +
  annotate("text", x = -0.84, y = 19.4, hjust = 0, family = FONT, size = pt(5.8), colour = INK2,
           label = "39 perturbations independent of both axes;\ntop 12 and bottom 8 shown") +
  labs(x = "decoupling index  D  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6))

top <- lab_grid(p7a, p7b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1.12))
save_fig(lab_grid(top, p7c, labels = c("","C"), ncol = 1, rel_heights = c(1, 1.05)),
         "Fig7.png", 183, 152)
