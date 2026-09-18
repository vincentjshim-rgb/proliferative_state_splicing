FIGURES_WRITTEN <- c("Fig1.png")
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
dir.create(PUBDIR, showWarnings = FALSE, recursive = TRUE)
source("scripts/pub_theme.R")
D <- "public_data_tierA/derived"

U <- read.delim(file.path(D, "axis_characterisation/UVA_axis_generalisation.tsv"))
S <- read.delim(file.path(D, "axis_characterisation/senescence_axis_generalisation.tsv"))
loc <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
lU <- loc$rho_UVA[loc$id == "LOCAL_ReproCM"]; lS <- loc$rho_senescence[loc$id == "LOCAL_ReproCM"]

grp <- function(k, ax) {
  if (grepl("\\(axis\\)", k)) return("defines axis")
  if (ax == "sen") return(if (grepl("same study", k)) "same study" else "independent")
  if (grepl("^injury", k)) return("independent\nUV injury")
  if (grepl("^photoprotection", k)) return("independent\nphotoprotection")
  "rescue inside\naxis study"
}
G <- rbind(
  data.frame(axis = "UVA reference axis", g = vapply(U$kind, grp, "", ax = "uva"), rho = U$rho_axis),
  data.frame(axis = "Senescence reference axis", g = vapply(S$kind, grp, "", ax = "sen"), rho = S$rho_axis))
G$axis <- factor(G$axis, levels = c("Senescence reference axis", "UVA reference axis"))
G$g <- factor(G$g, levels = rev(c("defines axis", "same study", "independent",
                                  "independent\nphotoprotection", "independent\nUV injury",
                                  "rescue inside\naxis study")))
COL <- c("defines axis" = INK, "same study" = GREY, "independent" = GREEN,
         "independent\nphotoprotection" = BLUE, "independent\nUV injury" = VERM,
         "rescue inside\naxis study" = PURP)
med <- aggregate(rho ~ axis + g, G, median)
anc <- data.frame(axis = factor(c("Senescence reference axis", "UVA reference axis"),
                                levels = levels(G$axis)), rho = c(lS, lU))

p1a <- ggplot(G, aes(rho, g)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(data = anc, aes(xintercept = rho), colour = VERM, linewidth = 0.3, linetype = "22") +
  geom_point(aes(colour = g), size = 1.15, alpha = 0.95,
             position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_point(data = med, shape = "|", size = 2.6, colour = INK, stroke = 0.7) +
  geom_text(data = transform(anc, g = factor(c("independent", "independent\nphotoprotection"),
                                             levels = levels(G$g))),
            aes(x = rho, y = g, label = "Repro-CM"),
            colour = VERM, size = pt(5.2), family = FONT, hjust = -0.12, vjust = -1.1) +
  facet_wrap(~ axis, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(limits = c(-0.28, 0.93), breaks = seq(-0.25, 0.75, 0.25),
                     expand = expansion(mult = 0.01)) +
  coord_cartesian(clip = "off") +
  labs(x = expression(paste("Spearman ", rho, " with the reference axis")), y = NULL) +
  theme_pub(7) + theme_cat_y() +
  theme(axis.text.y = element_text(lineheight = 0.92, hjust = 1),
        panel.spacing.y = unit(9, "pt"),
        strip.text = element_text(face = "plain", size = 7, colour = INK, hjust = 0))

tp <- data.frame(
  lab = c("UVA 6 h", "UVA 24 h", "UVA 24 h,\nrepeated", "UVA, chronic", "UVB 1 h", "UVB 4 h",
          "UVA, 3 d\n(held out)"),
  rho = c(0.276, 0.017, -0.030, -0.005, -0.104, -0.089, -0.0174),
  ho  = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
tp <- tp[order(tp$rho), ]; tp$lab <- factor(tp$lab, levels = tp$lab)
p1b <- ggplot(tp, aes(rho, lab)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_segment(aes(x = 0, xend = rho, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(aes(colour = ho), size = 1.5) +
  geom_text(aes(label = sprintf("%+.3f", rho), hjust = ifelse(rho >= 0, -0.28, 1.28)),
            size = pt(5.4), family = FONT, colour = GREY_D) +
  scale_colour_manual(values = c(`FALSE` = SLATE, `TRUE` = VERM), guide = "none") +
  scale_x_continuous(limits = c(-0.24, 0.42), breaks = c(-0.2, 0, 0.2, 0.4)) +
  labs(x = expression(paste(rho, " with UVA axis")), y = NULL) +
  theme_pub(7) + theme_cat_y() +
  theme(axis.text.y = element_text(lineheight = 0.92, hjust = 1))

ggsave(file.path(PUBDIR, "_p1a.png"), p1a, width = mm(105), height = mm(74), dpi = 300, bg = "white")
cat("p1a ok\n")
ggsave(file.path(PUBDIR, "_p1b.png"), p1b, width = mm(72), height = mm(74), dpi = 300, bg = "white")
cat("p1b ok\n")
fig1 <- compose(p1a, p1b, labels = c("a", "b"), ncol = 2, rel_widths = c(1.34, 1))
save_fig(fig1, "Fig1.png", 183, 74)
cat("Fig1 saved\n")
