## Supplementary Fig. S6 (was main-text Fig. 8, before that Fig. 9): do interventions
## raise splicing-factor expression by more than their effect on the cell cycle predicts?
## Revised 2026-09-15: recomputed contrast table (GSE149694 Fibroblast-D7 is a
## culture-time contrast and leaves the reprogramming class, 8 -> 7), confidence
## intervals for both class correlations, and one name per contrast everywhere.
## Revised 2026-09-17: moved to the supplement (FigS6.png). Panel a no longer plots
## the fitted values of the 63-contrast regression, which contains the interventions
## it was said to predict; each class is predicted from the regression refitted
## without that class (pred_lco from run_interventions.R). Set PRED <- "pred" to get
## the in-sample version back. 2026-09-17: the unreplicated contrast of this study was
## removed from the compendium, so every column here is a public contrast.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
I <- file.path(D, "interventions_revised")
## the shared num_axis() passes the raw break vector to format(), which turns
## ggplot's floating-point zero into 1.1e-16; fix it locally
num_axis <- function(digits = 1) function(x)
  sub("-", MINUS, formatC(round(x, digits + 2), format = "f", digits = digits))

R  <- read.delim(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
S  <- read.delim(file.path(I, "reprog_secretome_splicing.tsv"))
CLS <- read.delim(file.path(I, "class_correlations.tsv"))
NS <- read.delim(file.path(I, "named_SF_reprogramming.tsv"), check.names = FALSE)
ACC <- read.delim(file.path(I, "prediction_accuracy.tsv"))
R$class[is.na(R$class)] <- "metabolic / culture"
PRED <- "pred_lco"                        # "pred" = in-sample fitted value of the 63-contrast fit
REPC <- "#2C4A6E"; SECC <- TEAL

## one name per contrast, used in every panel
NM <- c(GSE149694_t2iLGoYD13 = "naive medium (t2iLGoY), day 13",
        GSE149694_NHSMD13    = "NHSM medium, day 13",
        GSE149694_5iLAFD13   = "5iLAF medium, day 13",
        GSE149694_RSeTD13    = "RSeT medium, day 13",
        GSE149694_PrimedD13  = "primed medium, day 13",
        GSE297233_OSK        = "OSK induction, day 4",
        GSE297233_O4YRSK     = "OCT4YR mutant + SK, day 4")
SHORTNM <- c(GSE149694_t2iLGoYD13 = "naive (t2iLGoY)", GSE149694_NHSMD13 = "NHSM",
             GSE149694_5iLAFD13 = "5iLAF", GSE149694_RSeTD13 = "RSeT",
             GSE149694_PrimedD13 = "primed", GSE297233_OSK = "OSK d4",
             GSE297233_O4YRSK = "OCT4YR + SK")
S$label <- ifelse(S$id %in% names(NM), NM[S$id], sub(" vs .*", "", S$name))

## ====== a. observed against what the cell-cycle effect predicts =============
## the prediction for a class comes from the regression of Fig. 4a refitted without
## that class, so no intervention contributes to its own prediction
S$grp <- ifelse(S$class == "reprogramming", "reprogramming medium", "secretome preparation")
S$grp <- factor(S$grp, levels = c("reprogramming medium", "secretome preparation"))
S$px <- S[[PRED]]
S$resid <- S$spl - S$px
acc <- function(sc) ACC[ACC$scheme == sc & ACC$group == "all interventions", ]
A_OUT <- acc("leave-class-out"); A_IN <- acc("in-sample (63-contrast fit)")
## an unmatched scheme name would give a zero-row data frame and silently empty the
## on-panel note, so both must have matched exactly one row
stopifnot(nrow(A_OUT) == 1, nrow(A_IN) == 1)
A_USE <- if (PRED == "pred_lco") A_OUT else A_IN
stopifnot(abs(median(abs(S$resid)) - A_USE$median_abs_residual) < 1e-9, nrow(S) == A_USE$n_contrasts)
xlim <- c(-0.30, 0.55); ylim <- c(-0.34, 0.82)
## named: every contrast more than 0.10 from its prediction. Label and leader
## positions are set by hand so that nothing overlaps.
LABS <- data.frame(
  id = c("GSE149694_t2iLGoYD13", "GSE149694_NHSMD13", "SC11", "SC10", "SC7",
         "GSE297233_OSK", "GSE297233_O4YRSK"),
  lx = c( 0.235, -0.018, -0.018, 0.285, 0.20,  0.20,  0.012),
  ly = c( 0.55,   0.43,   0.26,  0.07, -0.05, -0.20, -0.29),
  hj = c( 0,      1,      1,     0,     0,     0,     0),
  sx = c( NA,    -0.012, -0.012, 0.28,  0.19,  0.19,  0.006),
  sy = c( NA,     0.41,   0.25,  0.08, -0.06, -0.19, -0.275))
stopifnot(setequal(LABS$id, S$id[abs(S$resid) > 0.10]))
LABS <- merge(LABS, S[, c("id", "px", "spl", "label")], by = "id")
LABS$label[LABS$id == "GSE149694_t2iLGoYD13"] <- "naive medium\n(t2iLGoY), day 13"
LABS$label[LABS$id == "SC10"] <- "primary MSC secretome"
LABS$label[LABS$id == "SC11"] <- "PRF serum"
LABS$label[LABS$id == "SC7"]  <- "hTSC secretome"
NOTE <- sprintf(paste0("dashed line: observed = predicted\n",
                       "median |residual| %s against\nmedian |change| %s (n = %d contrasts);\n",
                       "%s with the in-sample 63-contrast fit"),
                num(A_OUT$median_abs_residual, 3), num(A_OUT$median_abs_change, 3),
                A_OUT$n_contrasts, num(A_IN$median_abs_residual, 3))
GV <- c("reprogramming medium", "secretome preparation")
pA <- ggplot(S, aes(px, spl)) +
  geom_abline(slope = 1, intercept = 0, colour = GREY, linewidth = 0.4, linetype = "22") +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_segment(aes(xend = px, yend = px), colour = GREY_L, linewidth = 0.35) +
  geom_segment(data = LABS[!is.na(LABS$sx), ], aes(x = px, y = spl, xend = sx, yend = sy),
               colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = grp, shape = grp, fill = grp), size = 1.9, stroke = 0.6) +
  geom_text(data = LABS, aes(lx, ly, label = label, hjust = hj),
            size = pt(7), family = FONT, colour = INK, lineheight = 0.95, vjust = 0.5) +
  annotate("label", x = xlim[1] + 0.008, y = ylim[2] + 0.004, hjust = 0, vjust = 1, size = pt(7), family = FONT,
           colour = INK2, lineheight = 1.1, label = NOTE, fill = "white", label.size = 0,
           label.r = unit(0, "pt"), label.padding = unit(1.5, "pt")) +
  scale_colour_manual(values = c(`reprogramming medium` = REPC, `secretome preparation` = SECC),
                      breaks = GV, name = NULL) +
  scale_fill_manual(values = c(`reprogramming medium` = REPC, `secretome preparation` = SECC),
                    breaks = GV, name = NULL) +
  scale_shape_manual(values = c(`reprogramming medium` = 21, `secretome preparation` = 21),
                     breaks = GV, name = NULL) +
  scale_x_continuous(breaks = seq(-0.25, 0.50, 0.25), labels = num_axis(2)) +
  scale_y_continuous(breaks = seq(-0.25, 0.75, 0.25), labels = num_axis(2)) +
  guides(colour = guide_legend(nrow = 1)) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(x = "change predicted from the change in cell-cycle genes\n(63-contrast regression refitted without the class being predicted)",
       y = "observed change in\npre-mRNA processing genes") +
  theme_sa() + theme(legend.position = "bottom", legend.margin = margin(t = -5))

## ====== b. the coupling within each intervention class =====================
## solid line: regression within the class. Dashed grey: the regression across the
## contrasts outside the class, which is what panel a predicts from.
sub <- R[R$class %in% c("reprogramming", "secretome"), ]
NSER <- c(reprogramming = length(unique(S$series[S$class == "reprogramming"])),
          secretome = length(unique(S$series[S$class == "secretome"])))
st <- do.call(rbind, lapply(CLS$class, function(k) { r <- CLS[CLS$class == k, ]
  data.frame(class = k, lab = sprintf("%s\nn = %d contrasts from %d series", k, r$n, NSER[[k]]),
             txt = sprintf("%s = %s [%s, %s]\n%s", RHO, num(r$rho), num(r$lo), num(r$hi),
                           pfmt(r$p))) }))
OUT <- do.call(rbind, lapply(CLS$class, function(k) { m <- lm(spl ~ cc, R[R$class != k, ])
  data.frame(class = k, a = coef(m)[[1]], b = coef(m)[[2]]) }))
OUT$lab <- st$lab[match(OUT$class, st$class)]
sub$lab <- st$lab[match(sub$class, st$class)]
note <- data.frame(lab = st$lab[st$class == "secretome"], x = 0.50, y = -0.26,
                   txt = "dashed: regression across the\ncontrasts outside the class")
pB <- ggplot(sub, aes(cc, spl)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_abline(data = OUT, aes(intercept = a, slope = b), colour = GREY, linewidth = 0.4, linetype = "22") +
  geom_smooth(aes(colour = class, fill = class), method = "lm", formula = y ~ x, se = TRUE,
              alpha = 0.15, linewidth = 0.55, show.legend = FALSE) +
  geom_point(aes(colour = class), size = 1.6, show.legend = FALSE) +
  ## starts just right of x = 0: from the panel edge the zero line ran through the text
  geom_text(data = st, aes(x = 0.012, y = Inf, label = txt), hjust = 0, vjust = 1.25,
            size = pt(7), family = FONT, colour = INK, lineheight = 1.05, inherit.aes = FALSE) +
  geom_text(data = note, aes(x, y, label = txt), hjust = 1, vjust = 0.5, size = pt(7), family = FONT,
            colour = INK2, lineheight = 1.0, inherit.aes = FALSE) +
  facet_wrap(~lab, nrow = 1) +
  scale_colour_manual(values = c(reprogramming = REPC, secretome = SECC)) +
  scale_fill_manual(values = c(reprogramming = REPC, secretome = SECC)) +
  scale_x_continuous(labels = num_axis(1)) +
  scale_y_continuous(labels = num_axis(1), expand = expansion(mult = c(0.05, 0.28))) +
  labs(x = "change in cell-cycle genes", y = "change in pre-mRNA\nprocessing genes") +
  theme_sa() + theme(strip.text = element_text(lineheight = 1.0))

## ====== c. named splicing factors across the reprogramming contrasts =======
pick <- intersect(names(SHORTNM), names(NS))
L <- do.call(rbind, lapply(pick, function(k)
  data.frame(gene = NS$gene, v = NS[[k]], set = SHORTNM[k], stringsAsFactors = FALSE)))
ordg <- names(sort(sapply(split(L$v, L$gene), mean, na.rm = TRUE)))
ordc <- names(sort(sapply(split(L$v, L$set), mean, na.rm = TRUE), decreasing = TRUE))
L$gene <- factor(L$gene, levels = ordg)
L$set  <- factor(L$set,  levels = ordc)
## Genes absent from a contrast's gene table were filled with #EDEFF1, which cannot be
## told from the near-zero fill (#F7F7F7). Absent cells are now a mid grey that is not
## on the colour scale, and the grey has its own key.
ABSENT <- "gene absent from the table"; ABSC <- "#AEB4BB"
cat(sprintf("  panel c: absent cells per contrast: %s\n",
            paste(sprintf("%s %d", levels(L$set), tapply(is.na(L$v), L$set, sum)), collapse = "; ")))
pC <- ggplot(L, aes(set, gene, fill = pmax(pmin(v, 1.2), -1.2))) +
  geom_tile(colour = "white", linewidth = 0.35) +
  ## the key for absent cells is drawn only when some cell is actually absent
  (if (any(is.na(L$v))) list(
     geom_tile(data = L[is.na(L$v), ], aes(alpha = ABSENT), fill = ABSC, colour = "white",
               linewidth = 0.35),
     scale_alpha_manual(values = setNames(1, ABSENT), name = NULL,
                        guide = guide_legend(order = 2, override.aes = list(fill = ABSC)))) else NULL) +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-1.2, 1.2), breaks = c(-1, 0, 1),
    labels = c(paste0(MINUS, "1"), "0", "1"), na.value = ABSC,
    name = sprintf("log₂ FC"), guide = guide_colourbar(order = 1, barwidth = unit(40, "pt"),
      barheight = unit(4, "pt"), ticks.colour = "white", frame.colour = NA,
      title.position = "top")) +
  scale_x_discrete(expand = c(0, 0), position = "top") + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_sa() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(angle = 40, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7, face = "italic"),
    legend.position = "bottom", legend.margin = margin(t = 1),
    legend.box = "horizontal", legend.box.just = "bottom", legend.spacing.x = unit(10, "pt"),
    legend.text = element_text(size = 7), legend.key.size = unit(8, "pt"),
    plot.margin = margin(3, 34, 3, 3))

left <- lab_grid(pA, pB, labels = c("a", "b"), ncol = 1, rel_heights = c(1.30, 1))
save_fig(lab_grid(left, pC, labels = c("", "c"), ncol = 2, rel_widths = c(1.38, 1)),
         "FigS4.png", 183, 136)
