#!/usr/bin/env Rscript
# Figure 7 : fibroblast perturbation compendium on the decoupling axis
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
root <- "public_data_tierA"
out  <- file.path(root, "derived", "manuscript_figures_v4")
R <- read.delim(file.path(root, "derived", "compendium", "compendium_D_scores.tsv"))

BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"; PURPLE <- "#CC79A7"
YELLOW <- "#E69F00"; GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"
cls <- function(id, ds, ct) {
  if (id == "LOCAL_ReproCM") return("Repro-CM (this study)")
  if (grepl("senescen|P27|late vs early|SIPS vs", ct, ignore.case = TRUE)) return("senescence")
  if (grepl("donor age|old vs young", ct, ignore.case = TRUE)) return("donor age")
  ## photoprotection first: a rescue contrast also mentions the UV it rescues from
  if (grepl("rescue|Pre-Red|prered|red_vs|Osmoter|UIFSP|Maifuyin|succinate", ct, ignore.case = TRUE))
    return("photoprotection / rescue")
  if (grepl("UVA|UVB|UV |uv_vs|injury|solar", ct, ignore.case = TRUE)) return("UV injury")
  if (ds %in% c("GSE149694", "GSE297233", "GSE165177")) return("reprogramming")
  if (ds == "GSE179848") return("metabolic / culture")
  "other"
}
R$class <- mapply(cls, R$id, R$dataset, R$contrast)
R$class <- factor(R$class, levels = c("Repro-CM (this study)", "photoprotection / rescue",
                  "UV injury", "reprogramming", "metabolic / culture", "donor age",
                  "senescence", "other"))
COL <- c("Repro-CM (this study)" = ORANGE, "photoprotection / rescue" = BLUE,
         "UV injury" = "#8FB8DE", "reprogramming" = GREEN, "metabolic / culture" = YELLOW,
         "donor age" = PURPLE, "senescence" = "#444444", "other" = GREY)
R <- R[order(R$D), ]
R$lab <- paste0(sub("this study", "LOCAL", R$dataset), "  ", substr(R$contrast, 1, 40))
R$lab <- factor(R$lab, levels = R$lab)
R$axis_flag <- ifelse(R$circularity == "IN reference axis", "defines an axis", "not in an axis")

base_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text.y = element_text(colour = MUTED, size = 6.3),
        axis.text.x = element_text(colour = MUTED),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.2),
        plot.caption = element_text(colour = MUTED, size = 7, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(8, "pt"), legend.text = element_text(size = 7.6))

pA <- ggplot(R, aes(D, lab, colour = class)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = D, yend = lab), linewidth = 0.5, colour = "grey85") +
  geom_point(aes(shape = axis_flag), size = 1.9) +
  scale_colour_manual(values = COL) +
  scale_shape_manual(values = c(`defines an axis` = 1, `not in an axis` = 16)) +
  labs(x = "Decoupling index  D", y = NULL,
       title = "A. Fifty-eight fibroblast perturbations on one axis",
       subtitle = "Open circles contributed to a reference axis and are therefore positive controls, not evidence",
       caption = paste0("Senescence contrasts occupy the bottom, including four GSE306957 arms fully independent of the references.\n",
                        "Contact inhibition - arrest WITHOUT senescence - sits on the positive side, so D is not an arrest axis.\n",
                        "All six photoprotection contrasts from studies that did not build the axis have positive D (+0.034 to +0.152);\n",
                        "the two negative rescue points come from the study that defines the UVA axis, where a negative value is structural.\n",
                        "Among the 37 signatures independent of both reference axes, the Repro-CM anchor has the highest D.")) +
  base_theme
ggsave(file.path(out, "Figure7A_compendium_ranked.png"), pA,
       width = 182, height = 216, units = "mm", dpi = 300)

LAB <- data.frame(
  id = c("LOCAL_ReproCM","GSE179848_ContactInhibition","GSE306957_pLenti","GSE109700_deep",
         "GSE240226_UVA","GSE113957_age","GSE165177_MPTR","GSE297233_OSK",
         "GSE179848_2Deoxyglucose","GSE149694_t2iLGoYD13"),
  short = c("Repro-CM","contact inhibition","senescence (held-out)","deep replicative senescence",
            "acute UVA","donor age","MPTR","OSK +dox","2-deoxyglucose","naive reprogramming"),
  lx = c(0.40, 0.30, 0.22, 0.30, 0.62, 0.34, -0.06, 0.16, -0.02, 0.22),
  ly = c(-0.17, 0.04, 0.66, 0.76, -0.05, 0.40, 0.46, -0.15, -0.34, -0.27),
  hj = c(0, 0, 1, 0, 1, 0, 1, 0, 1, 0))
R <- merge(R, LAB, by = "id", all.x = TRUE)
R$lab <- factor(R$lab, levels = levels(R$lab))
lab <- R[!is.na(R$short), ]
pB <- ggplot(R, aes(rho_UVA, rho_senescence)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = 2, linewidth = 0.35) +
  geom_segment(data = lab, aes(x = rho_UVA, y = rho_senescence, xend = lx, yend = ly),
               colour = "grey70", linewidth = 0.22, inherit.aes = FALSE) +
  geom_point(aes(colour = class, shape = axis_flag), size = 2.1) +
  geom_text(data = lab, aes(x = lx, y = ly, label = short, hjust = hj),
            size = 2.3, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = COL) +
  scale_shape_manual(values = c(`defines an axis` = 1, `not in an axis` = 16)) +
  scale_x_continuous(limits = c(-0.25, 0.92)) +
  scale_y_continuous(limits = c(-0.40, 0.86)) +
  guides(colour = guide_legend(nrow = 3, order = 1), shape = guide_legend(nrow = 3, order = 2)) +
  labs(x = expression(rho[UVA]), y = expression(rho[senescence]),
       title = "B. The two axes separate the perturbation classes",
       subtitle = "D is the vertical distance below the dashed identity line; lower-right = decoupled",
       caption = paste0("Senescence and donor-age signatures sit high on the senescence axis.\n",
                        "Contact inhibition - growth arrest without senescence - does not.")) +
  base_theme + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
                     legend.box = "horizontal")
ggsave(file.path(out, "Figure7B_compendium_scatter.png"), pB,
       width = 168, height = 142, units = "mm", dpi = 300)
cat("Figure 7 written\n")
