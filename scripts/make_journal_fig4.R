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
D <- "public_data_tierA/derived"

## ============== A. 553 Reactome pathways: local NES vs axis difference ======
P1 <- read.delim(file.path(D, "direction_probe/P1_pathway_decoupling_ranked.tsv"))
nmx <- grep("decoupl|diff", names(P1), value = TRUE)
P1$dx <- if (length(nmx)) P1[[nmx[1]]] else P1$UVA_meta_NES - P1$senescence_meta_NES
ct <- suppressWarnings(cor.test(P1$local_NES, P1$dx, method = "spearman", exact = FALSE))
fit <- lm(dx ~ local_NES, P1)
nd <- data.frame(local_NES = seq(min(P1$local_NES), max(P1$local_NES), length.out = 120))
nd <- cbind(nd, predict(fit, nd, interval = "confidence"))
short <- function(x) { x <- sub("__R-HSA.*$", "", x)
  for (k in names(GRK)) x <- gsub(k, GRK[[k]], x, fixed = TRUE); x }
GRK <- list("\u03b5"="e","\u03b1"="a","\u03b2"="b","\u03b3"="g","\u03ba"="k","\u2192"="->","\u2032"="'")
P1$nm <- short(P1$pathway)
key <- rbind(head(P1[order(-P1$dx), ], 3), head(P1[order(P1$dx), ], 2))
key$hj <- ifelse(key$local_NES > 0, 1.07, -0.07)
key$vj <- c(-1.1, 0.5, 2.1, 2.1, -1.1)
p4a <- ggplot(P1, aes(local_NES, dx)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_ribbon(data = nd, aes(local_NES, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = NAVY, alpha = 0.18) +
  geom_point(size = 0.55, colour = GREY, alpha = 0.55, stroke = 0) +
  geom_line(data = nd, aes(local_NES, fit), inherit.aes = FALSE, colour = NAVY, linewidth = 0.55) +
  geom_point(data = key, size = 1.2, colour = INK) +
  geom_text(data = key, aes(label = substr(nm, 1, 30), hjust = hj, vjust = vj),
            size = pt(5.4), family = FONT, colour = INK) +
  annotate("text", x = -2.9, y = 5.2, hjust = 0, vjust = 1, family = FONT, size = pt(6.3),
           colour = INK, lineheight = 1.05,
           label = sprintf("ρ = %+.3f\nP = %.1e\n553 pathways", ct$estimate, ct$p.value)) +
  coord_cartesian(xlim = c(-3.0, 3.4), ylim = c(-5.8, 5.6)) +
  labs(x = "local enrichment (NES) in the Repro-CM anchor",
       y = "NES difference   UVA axis − senescence axis") +
  theme_j(7)

## ================== B. module x axis map, expression-matched z =============
Z <- read.delim(file.path(D, "direction_probe/P2_module_axis_zmap.tsv"))
MODL <- c(premRNA = "Nuclear pre-mRNA processing", splicing = "mRNA splicing",
          three_prime = "mRNA 3'-end processing", mRNA_export = "mRNA export",
          SF_regulon = "Splicing-factor regulon", HSF1_heat_shock = "HSF1 heat-shock response",
          UPR_PERK = "UPR / PERK", autophagy = "Autophagy", ROS_detox = "ROS detoxification",
          IFN = "Interferon signalling", ECM = "Extracellular matrix",
          translation = "Translation", rRNA = "rRNA processing", NMD = "Nonsense-mediated decay",
          proliferation = "Cell cycle")
Z <- Z[Z$module %in% names(MODL), ]
Z$mod <- factor(MODL[Z$module], levels = rev(MODL))
Z$ax  <- factor(c(local = "Repro-CM\nanchor", UVA_meta = "UVA\naxis",
                  senesc_meta = "Senescence\naxis", fibroblast_age = "Donor\nage")[Z$axis],
                levels = c("Repro-CM\nanchor","UVA\naxis","Senescence\naxis","Donor\nage"))
Z <- Z[!is.na(Z$ax), ]
Z$star <- sig_stars(p.adjust(Z$empirical_p, "BH")); Z$star[Z$star == "n.s."] <- ""
Z$zc <- pmax(pmin(Z$matched_z, 8), -8)
p4b <- ggplot(Z, aes(ax, mod, fill = zc)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%+.1f", matched_z)), size = pt(5.6), family = FONT,
            colour = ifelse(abs(Z$zc) > 4.5, "white", INK), vjust = -0.05) +
  geom_text(aes(label = star), size = pt(5.2), family = FONT, vjust = 1.95,
            colour = ifelse(abs(Z$zc) > 4.5, "white", INK2)) +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-8, 8), breaks = c(-8, 0, 8),
                       name = "expression-matched z",
                       guide = guide_colourbar(barwidth = unit(44,"pt"), barheight = unit(3.6,"pt"),
                         ticks.colour = "white", frame.colour = NA, title.position = "top")) +
  scale_x_discrete(expand = c(0,0), position = "top") + scale_y_discrete(expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme_j(7) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.text.x.top = element_text(size = 6.4, lineheight = 0.9,
                                                    margin = margin(b = 3)),
                     legend.position = "bottom", legend.justification = "center",
                     legend.box.background = element_blank(),
                     legend.margin = margin(t = -2))

## ============ C. leading Reactome pathways in the local anchor =============
FG <- read.delim(file.path(D, "local_repro_anchor/local_Repro_specific_Reactome_fgsea.tsv"))
nes <- grep("^NES$", names(FG), value = TRUE); pv <- grep("padj|FDR", names(FG), value = TRUE)[1]
pw  <- grep("pathway", names(FG), value = TRUE)[1]
FG$nm <- ifelse(nchar(short(FG[[pw]])) > 34,
                paste0(substr(short(FG[[pw]]), 1, 32), "\u2026"), short(FG[[pw]]))
FG <- FG[order(-abs(FG[[nes]])), ]
sel <- rbind(head(FG[FG[[nes]] > 0, ], 8), head(FG[FG[[nes]] < 0, ], 6))
sel <- sel[order(sel[[nes]]), ]; sel$nm <- factor(sel$nm, levels = sel$nm)
sel$dir <- ifelse(sel[[nes]] > 0, "up in Repro-CM", "down in Repro-CM")
sel$star <- sig_stars(sel[[pv]])
p4c <- ggplot(sel, aes(.data[[nes]], nm, colour = dir)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = .data[[nes]], yend = nm), linewidth = 0.5) +
  geom_point(aes(size = -log10(.data[[pv]]))) +
  geom_text(aes(label = star, hjust = ifelse(.data[[nes]] > 0, -0.45, 1.45)),
            size = pt(5.6), family = FONT, colour = INK2) +
  scale_colour_manual(values = c(`up in Repro-CM` = REPRO, `down in Repro-CM` = AMBER),
                      name = NULL) +
  scale_size_continuous(range = c(0.8, 2.4), name = "-log10 FDR") +
  guides(colour = guide_legend(order = 1, nrow = 2),
         size   = guide_legend(order = 2, nrow = 2)) +
  scale_x_continuous(limits = c(-3.2, 3.6)) +
  labs(x = "normalised enrichment score", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_text(size = 6.2),
        legend.position = "right", legend.box = "vertical",
        legend.box.background = element_blank(), legend.key.size = unit(6, "pt"),
        legend.spacing.y = unit(1, "pt"), legend.margin = margin(l = -4))

top <- lab_grid(p4a, p4b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1.06))
save_fig(lab_grid(top, p4c, labels = c("","C"), ncol = 1, rel_heights = c(1.28, 1)),
         "Fig4.png", 183, 148)
