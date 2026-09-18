## Fig. 4 (six-figure restructure 2026-09-17; was Fig. 5, before that Fig. 6): does
## the coupling between splicing-factor expression and the cell cycle generalise
## across the published literature?
## Revised 2026-09-15: contrast table recomputed (donor age in the redefined
## cohort; GSE149694 Fibroblast-D7 moved to metabolic / culture), raw P values
## added to the class test, and a specificity panel against expression-matched
## random gene sets.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
## the shared num_axis() passes the raw break vector to format(), which turns
## ggplot's floating-point zero into 1.1e-16; fix it locally
num_axis <- function(digits = 1) function(x)
  sub("-", MINUS, formatC(round(x, digits + 2), format = "f", digits = digits))
R <- read.delim(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
R <- R[is.finite(R$spl) & is.finite(R$cc), ]
CL <- read.delim(file.path(D, "revision_stats/fig6_class_residuals.tsv"))
NL <- read.delim(file.path(D, "revision_stats/fig6_null_sets.tsv"))
CC <- c(senescence = RED, `donor age` = BROWN, secretome = TEAL,
        `photoprotection / rescue` = ORANGE, `UV injury` = PURPLE,
        `metabolic / culture` = GREY, reprogramming = "#2C4A6E")
R$class[is.na(R$class)] <- "metabolic / culture"
ST <- read.delim(file.path(D, "figures_sciadv/supp_table1.tsv"), check.names = FALSE)
NSTUDY <- length(unique(ST[[1]][grepl("^GSE", ST[[1]])]))

XLAB <- "change in cell-cycle genes (mean log₂ FC vs background)"
YLAB <- "change in pre-mRNA processing genes"

## ====== a. the 63 contrasts ================================================
## landmark contrasts, placed so that no leader line crosses another
NAMES <- c(GSE109700_deep = "deep replicative\nsenescence",
           GSE93535_SIPS = "stress-induced\nsenescence",
           GSE179848_ContactInhibition = "contact inhibition",
           GSE113957_age = "donor age",
           GSE297233_OSK = "OSK induction, day 4",
           GSE149694_t2iLGoYD13 = "naive medium, day 13",
           GSE329475_UVA3d = "UVA, 3 days")
## lx/ly: where the text sits; sx/sy: where its leader line stops, just short of it
LPOS <- data.frame(
  id = names(NAMES),
  lx = c(-1.70, -1.35, -1.28, -0.98,  0.78, -0.02,  0.86),
  ly = c(-1.26, -0.25,  0.03,  0.40, -0.62,  0.70,  0.68),
  hj = c( 0,     0,     0,     0,     1,     1,     1),
  sx = c(-1.72, -0.86, -0.90, -0.62,  0.33,  0.06,  0.80),
  sy = c(-1.30, -0.25,  0.03,  0.40, -0.62,  0.66,  0.60))
## 2026-09-17: "UVA, 3 days" sat on the regression line and was moved into clear
## space (position only, no value changed). The unreplicated contrast of this study
## was removed from the compendium on the same day.
LB <- merge(LPOS, R[, c("id", "cc", "spl")], by = "id"); LB$lab <- NAMES[LB$id]
ct <- cor.test(R$cc, R$spl, method = "spearman", exact = FALSE)
f  <- lm(spl ~ cc, R)
## sensitivity values printed in the grey note, computed here from the contrast table
## exactly as scripts/revision/run_stats_supplements.R (C2, C3) computes them
RHO_NOSHARE <- cor(R$cc, R$spl_noshare, method = "spearman")      # shared genes removed
NOSEN       <- R$class != "senescence"
RHO_NOSEN   <- cor(R$cc[NOSEN], R$spl[NOSEN], method = "spearman")
PERSTUDY    <- aggregate(cbind(spl, cc) ~ study, FUN = mean,
                         data = transform(R, study = sub("_.*$", "", id)))
RHO_STUDY   <- cor(PERSTUDY$cc, PERSTUDY$spl, method = "spearman")
## number of genes the 177-gene set shares with the cell-cycle set (same definitions)
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
GMTD <- tempfile("reactome"); dir.create(GMTD)
unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = GMTD)
GMT  <- strsplit(readLines(list.files(GMTD, pattern = "\\.gmt$", full.names = TRUE)[1]), "\t")
CCG  <- unique(unlist(lapply(GMT[grepl("^Cell Cycle, Mitotic|^DNA Replication|^M Phase",
                                       vapply(GMT, `[`, "", 1))], function(v) v[-(1:2)])))
N_SHARED <- length(intersect(CORE, CCG))
cat(sprintf("  panel a note: shared genes %d, rho without them %.3f; without senescence %.3f (n = %d); per study %.3f (n = %d)\n",
            N_SHARED, RHO_NOSHARE, RHO_NOSEN, sum(NOSEN), RHO_STUDY, nrow(PERSTUDY)))
nd <- data.frame(cc = seq(min(R$cc), max(R$cc), length.out = 80))
nd <- cbind(nd, predict(f, nd, interval = "confidence"))
pA <- ggplot(R, aes(cc, spl)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_ribbon(data = nd, aes(cc, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = GREY, alpha = 0.22) +
  geom_line(data = nd, aes(cc, fit), inherit.aes = FALSE, colour = INK, linewidth = 0.6) +
  geom_segment(data = LB, aes(x = cc, y = spl, xend = sx, yend = sy), inherit.aes = FALSE,
               colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = class), size = 1.7) +
  geom_text(data = LB, aes(lx, ly, label = lab, hjust = hj), inherit.aes = FALSE,
            family = FONT, size = pt(7), colour = INK2, vjust = 0.5, lineheight = 0.95) +
  annotate("text", x = -1.93, y = 0.72, hjust = 0, vjust = 1, family = FONT,
           size = pt(8), colour = INK, lineheight = 1.05,
           label = sprintf("%s = %s\n%s = %s\n%s", RHO, num(ct$estimate), R2,
                           num(summary(f)$r.squared), pfmt(ct$p.value))) +
  ## four lines now (the shared-gene result was cited to this panel but not drawn),
  ## anchored at the bottom edge so that the note grows upwards into empty space
  annotate("label", x = 0.85, y = -1.585, hjust = 1, vjust = 0, family = FONT,
           size = pt(7), colour = GREY, lineheight = 1.05, fill = "white", label.size = 0,
           label.r = unit(0, "pt"), label.padding = unit(1.5, "pt"),
           label = sprintf("%d contrasts, %d studies\nwithout the %d shared genes %s = %s\nwithout senescence %s = %s\none value per study %s = %s",
                           nrow(R), NSTUDY, N_SHARED, RHO, num(RHO_NOSHARE),
                           RHO, num(RHO_NOSEN), RHO, num(RHO_STUDY))) +
  scale_colour_manual(values = CC, name = NULL) +
  scale_x_continuous(XLAB, labels = num_axis()) +
  scale_y_continuous(YLAB, labels = num_axis()) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.6))) +
  coord_cartesian(xlim = c(-1.95, 0.85), ylim = c(-1.50, 0.75)) +
  theme_sa() + theme(legend.position = "bottom", legend.margin = margin(t = -4),
                     legend.text = element_text(size = 7))

## ====== b. the same within each class ======================================
keep <- names(which(table(R$class) >= 5))
Rk <- R[R$class %in% keep, ]
SHORT <- c(`metabolic / culture` = "metabolic,\nculture", `photoprotection / rescue` = "photo-\nprotection",
           reprogramming = "reprogramming", secretome = "secretome", senescence = "senescence",
           `UV injury` = "UV injury")
st <- do.call(rbind, lapply(sort(unique(Rk$class)), function(k) { s <- Rk[Rk$class == k, ]
  ctk <- cor.test(s$cc, s$spl, method = "spearman", exact = FALSE)
  data.frame(class = k, lab = SHORT[k],
             txt = sprintf("%s = %s\nn = %d", RHO, num(ctk$estimate), nrow(s))) }))
st$lab <- factor(st$lab, levels = SHORT[sort(unique(Rk$class))])
Rk$lab <- factor(SHORT[Rk$class], levels = levels(st$lab))
pB <- ggplot(Rk, aes(cc, spl)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_abline(slope = coef(f)[2], intercept = coef(f)[1], colour = GREY, linewidth = 0.4,
              linetype = "22") +
  geom_point(aes(colour = class), size = 1.4) +
  geom_text(data = st, aes(x = Inf, y = -Inf, label = txt), hjust = 1.1, vjust = -0.35,
            size = pt(7), family = FONT, colour = INK, lineheight = 1.0, inherit.aes = FALSE) +
  facet_wrap(~lab, nrow = 1) +
  scale_colour_manual(values = CC, guide = "none") +
  ## whole-number ticks: "−2.0", "−1.0", "0.0" ran together in the narrow facets
  scale_x_continuous(XLAB, breaks = c(-2, -1, 0), labels = num_axis(0)) +
  scale_y_continuous("change in pre-mRNA\nprocessing genes", labels = num_axis()) +
  theme_sa() + theme(strip.text = element_text(size = 7.2, lineheight = 0.95))

## ====== c. what is left once the cell-cycle change is accounted for ========
CL <- CL[order(CL$mean), ]
CL$lab <- factor(sprintf("%s (n = %d)", CL$class, CL$n),
                 levels = sprintf("%s (n = %d)", CL$class, CL$n))
## pfmt() takes one value at a time
CL$ptxt <- vapply(CL$p_raw, pfmt, character(1))
CL$clab <- sprintf("%s (n = %d)", CL$class, CL$n)
pC <- ggplot(CL, aes(mean, lab)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = class), height = 0, linewidth = 0.5) +
  geom_point(aes(colour = class), size = 1.9) +
  geom_text(aes(x = hi + 0.015, label = ptxt),
            hjust = 0, family = FONT, size = pt(7), colour = INK2) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous("residual (class mean, 95% CI)",
                     limits = c(-0.32, 0.45), breaks = seq(-0.2, 0.2, 0.2), labels = num_axis()) +
  scale_y_discrete(expand = expansion(add = c(0.75, 0.75))) +
  ## the note used to sit inside the panel, where the zero line ran through it
  ## BH-corrected P is identical for all six classes; read it from the file rather
  ## than hard-coding it, which left a stale 0.94 behind when the set changed
  labs(subtitle = sprintf("no detectable residual; BH P = %s for all six",
                          num(max(CL$p_BH), 2))) +
  theme_sa() + theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 7.2),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 7, colour = GREY, face = "italic", hjust = 1,
                                     margin = margin(b = 2)),
        plot.margin = margin(3, 5, 8, 3))

## ====== d. is the splicing set special? ====================================
obs <- unname(ct$estimate)
pD <- ggplot(NL, aes(rho)) +
  geom_histogram(bins = 34, fill = GREY_L, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs, colour = RED, linewidth = 0.6) +
  geom_vline(xintercept = median(NL$rho), colour = INK, linewidth = 0.4, linetype = "22") +
  annotate("text", x = obs - 0.04, y = Inf, hjust = 1, vjust = 1.4, family = FONT,
           size = pt(7.5), colour = RED, lineheight = 1.0,
           label = sprintf("observed\n%s = %s", RHO, num(obs))) +
  annotate("text", x = median(NL$rho) - 0.03, y = Inf, hjust = 1, vjust = 1.4, family = FONT,
           size = pt(7), colour = INK,
           label = sprintf("median\n%s", num(median(NL$rho)))) +
  ## label held clear of the bars, with a leader down to the largest random set
  annotate("segment", x = max(NL$rho), xend = max(NL$rho), y = 80, yend = 7,
           colour = GREY, linewidth = 0.25) +
  annotate("text", x = 0.62, y = 96, hjust = 0.5, vjust = 0.5, family = FONT,
           size = pt(7), colour = INK2, lineheight = 1.05,
           label = sprintf("maximum\n%s", num(max(NL$rho)))) +
  scale_x_continuous(sprintf("%s with the cell-cycle change", RHO),
                     limits = c(-0.5, 1), breaks = seq(-0.4, 0.8, 0.4), labels = num_axis()) +
  scale_y_continuous("random gene sets", expand = expansion(mult = c(0, 0.30))) +
  theme_sa()

top <- lab_grid(pA, pC, labels = c("a", "c"), ncol = 2, rel_widths = c(1.12, 1))
bot <- lab_grid(pB, pD, labels = c("b", "d"), ncol = 2, rel_widths = c(2.3, 1))
save_fig(plot_grid(top, bot, ncol = 1, rel_heights = c(1.42, 1)), "Fig4.png", 183, 138)
