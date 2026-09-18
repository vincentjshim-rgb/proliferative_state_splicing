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
PUBDIR <- "public_data_tierA/derived/figures_publication"
source("scripts/pub_theme.R"); D <- "public_data_tierA/derived"
R <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
cls <- function(id, ds, ct) {
  if (id == "LOCAL_ReproCM") return("Repro-CM")
  if (grepl("senescen|P27|late vs early|SIPS vs", ct, ignore.case = TRUE)) return("senescence")
  if (grepl("donor age|old vs young", ct, ignore.case = TRUE)) return("donor age")
  if (grepl("rescue|Pre-Red|prered|red_vs|Osmoter|UIFSP|Maifuyin|succinate|sauchinone",
            ct, ignore.case = TRUE)) return("photoprotection")
  if (grepl("UVA|UVB|UV |uv_vs|injury|solar", ct, ignore.case = TRUE)) return("UV injury")
  if (ds %in% c("GSE149694", "GSE297233", "GSE165177")) return("reprogramming")
  if (ds == "GSE179848") return("metabolic / culture")
  "other" }
R$cl <- mapply(cls, R$id, R$dataset, R$contrast)
COL <- c("Repro-CM" = VERM, photoprotection = BLUE, "UV injury" = "#9CC2DE",
         reprogramming = GREEN, "metabolic / culture" = GOLD, "donor age" = PURP,
         senescence = INK, other = GREY)
R$cl <- factor(R$cl, levels = names(COL))
R$axis_flag <- R$circularity == "IN reference axis"
R <- R[order(R$D), ]
R$lab <- factor(seq_len(nrow(R)))
p7a <- ggplot(R, aes(D, lab)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = D, yend = lab), colour = "#EAEEF1", linewidth = 0.35) +
  geom_point(aes(colour = cl, shape = axis_flag), size = 1.15, stroke = 0.45) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16), guide = "none") +
  scale_x_continuous(limits = c(-0.95, 1.16), breaks = seq(-0.75, 0.75, 0.25)) +
  labs(x = "decoupling index D", y = "60 fibroblast perturbations, ranked") +
  theme_pub(7) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_text(size = 6.5, colour = GREY_D))
ann <- data.frame(
  id = c("LOCAL_ReproCM", "GSE179848_ContactInhibition", "GSE109700_deep",
         "GSE240226_UVA", "GSE326951_sauchinone", "GSE113957_age"),
  s  = c("Repro-CM", "contact inhibition", "deep senescence", "acute UVA",
         "sauchinone", "donor age"))
A <- merge(ann, R[, c("id", "D", "lab")], by = "id")
A$y <- as.numeric(as.character(A$lab))
p7a <- p7a +
  geom_text(data = A, aes(x = D, y = y, label = s),
            hjust = ifelse(A$D > 0, -0.18, 1.18), size = pt(5), family = FONT, colour = INK,
            inherit.aes = FALSE) +
  coord_cartesian(clip = "off")
key <- data.frame(cl = factor(names(COL), levels = names(COL)),
                  y = seq(52, 52 - 7 * 3.1, length.out = 8), x = -0.93)
p7a <- p7a +
  geom_point(data = key, aes(x = x, y = y, colour = cl), size = 1.15, inherit.aes = FALSE) +
  geom_text(data = key, aes(x = x + 0.028, y = y, label = cl), hjust = 0, vjust = 0.42,
            size = pt(5), family = FONT, colour = GREY_D, inherit.aes = FALSE) +
  scale_colour_manual(values = COL, guide = "none")

lab2 <- data.frame(
  id = c("LOCAL_ReproCM", "GSE179848_ContactInhibition", "GSE306957_pLenti",
         "GSE109700_deep", "GSE240226_UVA", "GSE113957_age", "GSE165177_MPTR"),
  s = c("Repro-CM", "contact inhibition", "senescence (held out)", "deep senescence",
        "acute UVA", "donor age", "MPTR"),
  lx = c(0.44, 0.30, 0.46, 0.34, 0.62, 0.40, -0.10),
  ly = c(-0.20, 0.055, 0.545, 0.775, -0.055, 0.405, 0.455),
  hj = c(0, 0, 0, 0, 1, 0, 1))
L2 <- merge(lab2, R[, c("id", "rho_UVA", "rho_senescence")], by = "id")
p7b <- ggplot(R, aes(rho_UVA, rho_senescence)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, colour = "#DCE2E7", linetype = "22", linewidth = 0.3) +
  geom_segment(data = L2, aes(x = rho_UVA, y = rho_senescence, xend = lx, yend = ly),
               colour = GREY_L, linewidth = 0.2, inherit.aes = FALSE) +
  geom_point(aes(colour = cl, shape = axis_flag), size = 1.2, stroke = 0.45) +
  geom_text(data = L2, aes(x = lx, y = ly, label = s, hjust = hj), size = pt(5),
            family = FONT, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16), guide = "none") +
  scale_x_continuous(limits = c(-0.30, 0.92), breaks = seq(-0.25, 0.75, 0.25)) +
  scale_y_continuous(limits = c(-0.32, 0.84), breaks = seq(-0.25, 0.75, 0.25)) +
  labs(x = expression(rho[UVA]), y = expression(rho[senescence])) +
  theme_pub(7)
save_fig(compose(p7a, p7b, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.02)),
         "Fig7.png", 183, 92)
cat("Fig7 saved\n")
