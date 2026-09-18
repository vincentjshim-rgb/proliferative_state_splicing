## Fig. 6. Is the coupling a property of culture?  Preregistered GTEx boundary test.
## Main figure (Fig6.png), 96-gene splicing set throughout:
##   a  samples from the GTEx donor pool
##   b  the 96 splicing genes against the proliferation score in culture, skin, muscle
##   c  partial rho with proliferation for every programme in every tissue
##   d  age effect per decade before and after proliferation adjustment in two
##      fibroblast cohorts, and (set apart, post hoc) the two adjusted estimates combined
##   e  the six preregistered hypotheses
## Supplementary figure (FigS10.png), written by the same script:
##   a  regression slope by tissue and variance-matched subsampling        (post hoc)
##   b  coupling before and after cell-composition adjustment              (post hoc)
##   c  donor concordance between culture and skin (preregistered secondary analysis)
##      with the genotype-driven positive control                          (post hoc)
## The composition-adjustment values are read from composition_adjusted.tsv when the
## figure is rendered; nothing from that table is typed into this script.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
O  <- "public_data_tierA/derived/gtex_boundary"
PH <- "public_data_tierA/derived/gtex_posthoc"
CR <- "public_data_tierA/derived/cohort_revised"
TS <- readRDS(file.path(O, "tissue_scores.rds"))
S1 <- fread(file.path(O, "programme_by_tissue.tsv"))
H  <- fread(file.path(O, "primary_hypotheses.tsv"))
S3 <- fread(file.path(O, "donor_concordance_culture_vs_legskin.tsv"))
CNT <- fread(file.path(PH, "sample_counts.tsv"))
PC  <- fread(file.path(PH, "donor_identity_positive_control.tsv"))
RR  <- fread(file.path(PH, "range_restriction.tsv"))
VM  <- fread(file.path(PH, "variance_matched_summary.tsv"))
VD  <- fread(file.path(PH, "variance_matched_draws.tsv"))
CA  <- fread(file.path(PH, "composition_adjusted.tsv"))
SENS <- fread(file.path(CR, "sensitivity_metrics.tsv"))
PRIM <- fread(file.path(CR, "primary/sample_metrics.tsv"))
TK   <- c("culture","legskin","pubskin","muscle")
TCOL <- c(culture = RED, legskin = ORANGE, pubskin = TEAL, muscle = PURPLE)
nice <- function(x) { x <- sub("cell cycle \\(positive control\\)", "cell cycle, mitotic", x)
                      sub("^splicing 96$", "96 splicing genes", x) }
G <- function(tk, pg, col) S1[tissue == tk & programme == pg][[col]]
NF <- setNames(CNT$final_donors, CNT$tissue)          # donors in each tissue
NC <- setNames(CNT$complete_covariates, CNT$tissue)   # donors with complete covariates
NPAIR <- CNT[tissue == "culture"]$paired_with_legskin

## ---------------- a. design -------------------------------------------------------
## the culture box is taller than the others because it carries a third line
BX <- data.table(tk = TK, ymin = c(3.57, 2.62, 1.67, 0.72), ymax = c(4.67, 3.42, 2.47, 1.52),
  lab = c("cultured fibroblasts", "skin, sun-exposed", "skin, not sun-exposed", "skeletal muscle"),
  site = c("", "lower leg, ", "suprapubic, ", ""))
BX[, ymid := (ymin + ymax) / 2]
BX[, n := sprintf("%sn = %d donors", site, NF[tk])]
BX[, extra := ifelse(tk == "culture", sprintf("%d with complete covariates", NC["culture"]), "")]
BX[, `:=`(y1 = ifelse(tk == "culture", ymid + 0.31, ymid + 0.16),
          y2 = ifelse(tk == "culture", ymid,        ymid - 0.16),
          y3 = ymid - 0.31)]
pA <- ggplot(BX) +
  geom_rect(aes(xmin = 0.02, xmax = 0.98, ymin = ymin, ymax = ymax, colour = tk),
            fill = "white", linewidth = 0.45) +
  geom_rect(aes(xmin = 0.02, xmax = 0.06, ymin = ymin, ymax = ymax, fill = tk), colour = NA) +
  geom_text(aes(x = 0.085, y = y1, label = lab), hjust = 0, family = FONT, size = pt(7.5), colour = INK) +
  geom_text(aes(x = 0.085, y = y2, label = n), hjust = 0, family = FONT, size = pt(7), colour = INK2) +
  geom_text(aes(x = 0.085, y = y3, label = extra), hjust = 0, family = FONT, size = pt(7), colour = INK2) +
  annotate("segment", x = 1.02, xend = 1.02, y = BX$ymid[1], yend = BX$ymid[2], colour = GREY, linewidth = 0.35) +
  annotate("segment", x = 0.99, xend = 1.02, y = BX$ymid[1:2], yend = BX$ymid[1:2], colour = GREY, linewidth = 0.35) +
  annotate("text", x = 1.05, y = mean(BX$ymid[1:2]), hjust = 0, family = FONT, size = pt(7), colour = INK2,
           lineheight = 1.0, label = sprintf("%d donors\nin both", NPAIR)) +
  annotate("text", x = 0.02, y = 4.99, hjust = 0, family = FONT, size = pt(7.5), colour = INK2,
           label = "GTEx v8, one pipeline (recount3)") +
  scale_colour_manual(values = TCOL, guide = "none") + scale_fill_manual(values = TCOL, guide = "none") +
  coord_cartesian(xlim = c(0, 1.34), ylim = c(0.62, 5.15), clip = "off") + theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(12, 3, 3, 3))

## ---------------- b. the splicing genes against proliferation --------------------
## the y axis is given headroom so that the statistic never sits on the points
SC <- rbindlist(lapply(c("culture","legskin","muscle"), function(tk) { T <- TS[[tk]]
  data.table(tk = tk, x = T$d$prolif, y = T$S[, "splicing 96"]) }))
SC[, tk := factor(tk, levels = c("culture","legskin","muscle"))]
LB <- data.table(tk = factor(c("culture","legskin","muscle"), levels = levels(SC$tk)),
  txt = sapply(c("culture","legskin","muscle"), function(tk)
    sprintf("partial %s = %s\n%s\nn = %d donors", RHO, num(G(tk, "splicing 96", "rho_prolif")),
            pfmt(G(tk, "splicing 96", "p_prolif")), G(tk, "splicing 96", "n"))))
STRIP <- c(culture = "cultured fibroblasts", legskin = "skin, sun-exposed", muscle = "skeletal muscle")
pB <- ggplot(SC, aes(x, y)) +
  geom_point(aes(colour = tk), size = 0.5, alpha = 0.4) +
  geom_smooth(aes(colour = tk, fill = tk), method = "lm", formula = y ~ x, linewidth = 0.55, alpha = 0.15) +
  geom_text(data = LB, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.06,
            vjust = 1.12, family = FONT, size = pt(7), colour = INK, lineheight = 1.05) +
  facet_wrap(~tk, nrow = 1, scales = "free", labeller = as_labeller(STRIP)) +
  scale_colour_manual(values = TCOL, guide = "none") + scale_fill_manual(values = TCOL, guide = "none") +
  scale_x_continuous(labels = num_axis()) +
  scale_y_continuous(labels = num_axis(), expand = expansion(mult = c(0.04, 0.40))) +
  labs(x = "proliferation score (20 cell-cycle markers)", y = "96 splicing genes") +
  theme_sa() + theme(strip.text = element_text(size = 8))

## ---------------- c. every programme in every tissue -----------------------------
## transposed relative to the published version so that no label falls below 7 pt
M <- copy(S1)[, tk := factor(tissue, levels = TK)]
ord <- M[tissue == "culture"][order(-rho_prolif)]$programme
M[, pg := factor(nice(programme), levels = nice(ord))]
HL <- data.table(tk = factor(c("culture","legskin","muscle","culture"), levels = TK),
                 pg = factor(nice(c("splicing 96","splicing 96","splicing 96","collagen formation")),
                             levels = nice(ord)))
NPROG <- sapply(TK, function(tk) { n <- unique(S1[tissue == tk]$n); stopifnot(length(n) == 1); n })
YLAB <- setNames(sprintf("%s, n = %d donors", c("cultured fibroblasts", "skin, sun-exposed",
                                                "skin, not exposed", "skeletal muscle"), NPROG), TK)
pC0 <- ggplot(M, aes(pg, tk)) +
  geom_tile(aes(fill = pmax(pmin(rho_prolif, 1), -1)), colour = "white", linewidth = 0.4) +
  geom_tile(data = HL, fill = NA, colour = INK, linewidth = 0.6) +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-1, 1), breaks = c(-1, 0, 1),
    labels = c(paste0(MINUS, "1"), "0", "1"),
    name = sprintf("partial %s with\nproliferation", RHO),
    guide = guide_colourbar(barwidth = unit(5, "pt"), barheight = unit(42, "pt"),
                            title.position = "top", ticks.colour = "white", frame.colour = NA)) +
  scale_y_discrete(limits = rev(TK), labels = YLAB, expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme_sa() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, colour = INK),
    axis.text.y = element_text(size = 7.5, colour = INK),
    legend.title = element_text(size = 7.5, lineheight = 1.05), legend.text = element_text(size = 7),
    legend.margin = margin(0, 0, 0, 0), plot.margin = margin(3, 6, 3, 3))
## the colour bar goes into the empty corner under the tissue names, left of the
## slanted programme names, instead of adding a row of its own below them
LEGC <- get_plot_component(pC0 + theme(legend.position = "right"), "guide-box-right")
pC <- ggdraw(pC0 + theme(legend.position = "none")) +
  draw_grob(LEGC, x = 0.035, y = 0.08, width = 0.13, height = 0.44, hjust = 0, vjust = 0)

## ---------------- d. the raw age effect depends on how proliferation moves with age
## The fibroblast cohort is the redefined one (normal donors aged 20+, technical
## covariates). That redefinition was made after the outcomes were known, so nothing
## on this panel is labelled preregistered. The combined row is post hoc and is drawn
## apart from the two cohort rows, below a thin rule and as a diamond.
K <- SENS[metric == "machinery96" & cohort == "primary"]
rho_ap <- cor(PRIM$age, PRIM$prolif, method = "spearman")
d <- copy(TS$culture$d); d[, s := TS$culture$S[, "splicing 96"]]
tech <- c("rin","isch","loglib")[c(mean(is.finite(d$rin)) > 0.9, mean(is.finite(d$isch)) > 0.9, TRUE)]
cov <- c("sex", if (nlevels(droplevels(d$hardy)) > 1) "hardy", tech)
ga0 <- lm(as.formula(paste("s ~ I(age/10) +", paste(cov, collapse = "+"))), d); ga1 <- update(ga0, . ~ . + prolif)
ci <- function(m, term) { s <- summary(m)$coef; k <- confint(m)
  data.table(b = s[term, 1], lo = k[term, 1], hi = k[term, 2], p = s[term, 4]) }
dm <- d[as.integer(rownames(model.frame(ga0)))]           # the donors the model used
rho_gt <- cor(dm$prolif, dm$age, method = "spearman")
## GTEx releases age as ten-year bands; the models use the band midpoints
gt_age <- sprintf("%d to %d years", min(dm$age) - 5L, max(dm$age) + 4L)
## post hoc: the two independent cohorts combined by inverse variance. Only the
## adjusted estimates are pooled; unadjusted effects depend on how proliferation
## happens to move with age in each collection and are not the same quantity.
MT <- fread("public_data_tierA/derived/residual_meta/residual_age_effect_meta.tsv")
MT <- MT[grepl("96", gene_set) & grepl("^combined", cohort)]; stopifnot(nrow(MT) == 1)
E <- rbind(
  data.table(block = 1, model = c("unadjusted", "proliferation-adjusted"),
             b = c(K$beta, K$beta_adj), lo = c(K$lo, K$lo_adj), hi = c(K$hi, K$hi_adj), p = c(K$p, K$p_adj)),
  data.table(block = 2, model = c("unadjusted", "proliferation-adjusted"),
             rbind(ci(ga0, "I(age/10)"), ci(ga1, "I(age/10)"))),
  ## 2026-09-18: the inverse-variance row was withdrawn.  With two studies the
  ## between-study variance is not estimable, and the two estimates are not the same
  ## quantity -- GSE113957 is an attenuation, GTEx a suppression effect.
  NULL)
YB <- c(3.10, 1.60); YRULE <- NA_real_; YHEAD <- 4.02; YLIM <- c(0.75, 4.25)
E[, y := YB[block] + ifelse(model == "unadjusted", 0.26, -0.26)]
E[, model := factor(model, levels = c("unadjusted", "proliferation-adjusted"))]
E[, est := sprintf("%s [%s, %s]", num(b, 3), num(lo, 3), num(hi, 3))]
E[, ptx := pfmt_v(p)]
print(E)
CO <- data.table(y = YB, lab = c(
  sprintf("GSE113957 fibroblasts\n%d donors, %d to %d years\nproliferation vs age %s = %s",
          K$n, as.integer(min(PRIM$age)), as.integer(max(PRIM$age)), RHO, num(rho_ap)),
  sprintf("GTEx fibroblasts\n%d donors, %s\nproliferation vs age %s = %s", nobs(ga0), gt_age, RHO, num(rho_gt)),
  NULL))
void_d <- theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
ysc <- scale_y_continuous(limits = YLIM, expand = c(0, 0))
rule <- NULL  ## the rule separated the pooled row, which was withdrawn
head_d <- function(x, lab) annotate("text", x = x, y = YHEAD, hjust = 0, family = FONT, size = pt(7),
                                    colour = GREY, label = lab)
pDl <- ggplot() + rule +
  geom_text(data = CO, aes(0.04, y, label = lab), hjust = 0, family = FONT, size = pt(7), colour = INK,
            lineheight = 1.05) +
  geom_text(data = E, aes(0.63, y, label = model), hjust = 0, family = FONT, size = pt(7), colour = INK2) +
  head_d(0.04, "fibroblast cohort") + head_d(0.63, "model") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ysc + void_d +
  theme(plot.margin = margin(3, 0, 3, 3))
pDm <- ggplot(E, aes(b, y, colour = model)) + rule +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.55) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = c(unadjusted = GREY, `proliferation-adjusted` = BLUE), guide = "none") +
  scale_x_continuous("age effect per decade, 96 splicing genes",
                     limits = c(-0.20, 0.05), breaks = seq(-0.20, 0.05, 0.05), labels = num_axis()) +
  ysc + labs(y = NULL) + theme_sa() +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"), plot.margin = margin(3, 0, 3, 0))
pDr <- ggplot(E) + rule +
  geom_text(aes(0.06, y, label = est), hjust = 0, family = FONT, size = pt(7), colour = INK) +
  geom_text(aes(0.70, y, label = ptx), hjust = 0, family = FONT, size = pt(7), colour = INK) +
  head_d(0.06, "estimate [95% CI]") + head_d(0.70, "two-sided P") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ysc + void_d +
  theme(plot.margin = margin(3, 6, 3, 0))
pD <- plot_grid(pDl, pDm, pDr, nrow = 1, align = "h", axis = "tb", rel_widths = c(0.385, 0.345, 0.27))

## ---------------- e. preregistered scorecard --------------------------------------
sgn <- function(x, digits) paste0(ifelse(x > 0, "+", ""), num(x, digits))
V <- data.table(
  id = paste0("H", 1:6),
  what = c("culture: splicing genes track proliferation",
           "skin: coupling weaker than in culture",
           "skeletal muscle: no coupling",
           "culture: age effect shrinks with adjustment",
           "culture: collagen formation independent of proliferation",
           "skin: collagen formation falls with age"),
  crit = c(sprintf("%s ≥ 0.5", RHO), "one-sided P < 0.05", sprintf("%s < 0.3", RHO),
           "≥ 50% smaller", sprintf("|%s| < 0.3", RHO), "negative, P < 0.05"),
  obs = c(num(G("culture","splicing 96","rho_prolif")),
          sprintf("%s vs %s", num(G("legskin","splicing 96","rho_prolif")),
                  num(G("culture","splicing 96","rho_prolif"))),
          num(G("muscle","splicing 96","rho_prolif")),
          sprintf("no unadjusted effect (%s)", pfmt(G("culture","splicing 96","p_age"))),
          num(G("culture","collagen formation","rho_prolif")),
          sprintf("%s per decade (%s)", sgn(G("legskin","collagen formation","beta_age"), 3),
                  pfmt(G("legskin","collagen formation","p_age")))),
  verdict = H$verdict)
V[, vshort := fifelse(verdict == "supported", "supported",
              fifelse(grepl("^not evaluable", verdict), "not evaluable", "not supported"))]
V[, vcol := c(supported = PASS, `not evaluable` = NEUTRAL, `not supported` = FAIL)[vshort]]
V[, y := 6:1]
pE <- ggplot(V) +
  geom_hline(yintercept = seq(0.5, 6.5, 1), colour = GREY_L, linewidth = 0.25) +
  geom_text(aes(0.00, y, label = id), hjust = 0, family = FONT, size = pt(7.5), fontface = "bold", colour = INK) +
  geom_text(aes(0.035, y, label = what), hjust = 0, family = FONT, size = pt(7.5), colour = INK) +
  geom_text(aes(0.45, y, label = crit), hjust = 0, family = FONT, size = pt(7), colour = INK2) +
  geom_text(aes(0.62, y, label = obs), hjust = 0, family = FONT, size = pt(7), colour = INK) +
  geom_text(aes(0.90, y, label = vshort), hjust = 0, family = FONT, size = pt(7.5), fontface = "bold",
            colour = V$vcol) +
  annotate("text", x = c(0.035, 0.45, 0.62, 0.90), y = 6.95, hjust = 0, vjust = 0.5, family = FONT,
           size = pt(7), colour = GREY,
           label = c("preregistered hypothesis", "criterion", "observed", "verdict")) +
  annotate("text", x = 0, y = 0.15, hjust = 0, family = FONT, size = pt(7), colour = GREY,
           label = sprintf("H2 was preregistered as a one-sided test (P = %s); every other P is two-sided", sci(H[id == "H2"]$p))) +
  coord_cartesian(xlim = c(0, 1.02), ylim = c(0.05, 7.2), clip = "off") + theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(6, 6, 3, 6))

## ---------------- f. does cell composition explain the residual coupling in skin? ----
## Post hoc, added 2026-09-18.  The composition adjustment originally ran on the splicing
## sets only, so it could not distinguish "composition removes splicing's coupling" from
## "composition flattens the whole proliferation axis".  Run across programmes it
## separates them: the positive control is untouched and only splicing falls.
## NB: `CA` is already the splicing-only composition table loaded at the top and used by
## Supplementary Fig. S10b -- this one needs its own name.
CAP <- fread(file.path(PH, "composition_adjusted_all_programmes.tsv"))
CAP <- CAP[tissue %in% c("legskin", "pubskin")]
PGF <- c("cell cycle (positive control)", "splicing 96", "splicing 177",
         "mRNA splicing", "pre-mRNA processing", "collagen formation")
stopifnot(all(PGF %in% CAP$set))
CF <- CAP[set %in% PGF]
CF <- melt(CF, id.vars = c("tissue", "set"),
           measure.vars = c("rho_tech", "rho_tech_comp"),
           variable.name = "adj", value.name = "rho")
CF[, adj := factor(adj, levels = c("rho_tech", "rho_tech_comp"),
                   labels = c("technical covariates", "+ cell-composition markers"))]
CF[, pgf := factor(nice(set), levels = rev(nice(PGF)))]
CF[, tkf := factor(tissue, levels = c("legskin", "pubskin"),
                   labels = c("skin, sun-exposed (n = 739)", "skin, not exposed (n = 626)"))]
## the two splicing definitions that were never selected in these data are marked, so a
## reader can see the drop is not a property of the selected set
UNSEL <- data.table(pgf = factor(nice(c("mRNA splicing", "pre-mRNA processing")),
                                 levels = rev(nice(PGF))))
pF <- ggplot(CF, aes(rho, pgf, colour = adj)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_line(aes(group = interaction(tkf, pgf)), colour = GREY, linewidth = 0.4) +
  geom_point(size = 2) +
  facet_wrap(~ tkf, nrow = 1) +
  scale_colour_manual(values = c("technical covariates" = GREY,
                                 "+ cell-composition markers" = BLUE), name = NULL) +
  scale_x_continuous(sprintf("partial %s with proliferation score", RHO),
                     limits = c(-0.15, 0.78), breaks = seq(0, 0.75, 0.25), labels = num_axis()) +
  labs(y = NULL) + theme_sa() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 7.5, colour = INK),
        strip.text = element_text(size = 7.5, colour = INK),
        legend.position = "top", legend.text = element_text(size = 7),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(3, 6, 3, 3))
cat("\npanel f -- composition adjustment by programme:\n")
print(dcast(CF, set ~ tkf + adj, value.var = "rho"), digits = 3)

## 2026-09-18: the preregistered scorecard was a table drawn as a figure panel, and the
## figure had grown past a printable page.  It is written out as Supplementary Table 10
## instead; every verdict, including the failure and the unevaluable test, is unchanged.
fwrite(V[, .(hypothesis = id, test = what, criterion = crit, observed = obs, verdict = verdict)],
       file.path(O, "preregistered_scorecard.tsv"), sep = "\t")
cat("\nSupplementary Table 10 (preregistered scorecard) written\n")

top <- lab_grid(pA, pB, labels = c("a", "b"), ncol = 2, rel_widths = c(0.94, 2.0))
save_fig(plot_grid(top, lab_grid(pC, labels = "c", ncol = 1), lab_grid(pD, labels = "d", ncol = 1),
                   lab_grid(pF, labels = "e", ncol = 1), ncol = 1,
                   rel_heights = c(53, 53, 44, 46)), "Fig6.png", 183, 196)

## ================= Supplementary Fig. S10: post hoc tissue analyses and donor concordance
TN <- c(culture = "cultured\nfibroblasts", legskin = "skin,\nsun-exposed",
        pubskin = "skin,\nnot exposed", muscle = "skeletal\nmuscle")
## ---- S10a. slope by tissue, with variance-matched subsampling (post hoc)
R9 <- RR[set == "splicing 96"][, tissue := factor(tissue, levels = TK)]
XN <- setNames(sprintf("%s\nn = %d", TN[as.character(R9$tissue)], R9$n), as.character(R9$tissue))
vmw <- VM[scheme == "window"]; vmk <- VM[scheme == "kernel"]
nsub <- unique(VD$n); stopifnot(length(nsub) == 1, vmw$draws == vmk$draws)
vm_txt <- sprintf(paste0("post hoc. Cultured fibroblasts subsampled\nto the leg-skin SD of proliferation (%s);\n",
                         "%d draws of %d donors, median [95%% range]\n",
                         "hard window: partial %s = %s [%s, %s]\n",
                         "kernel weights: partial %s = %s [%s, %s]"),
                  num(vmw$target_sd), vmw$draws, nsub,
                  RHO, num(vmw$rho_median), num(vmw$rho_lo), num(vmw$rho_hi),
                  RHO, num(vmk$rho_median), num(vmk$rho_lo), num(vmk$rho_hi))
pSa <- ggplot(R9, aes(as.integer(tissue), slope, colour = tissue)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.4) +
  geom_errorbar(aes(ymin = slope_lo, ymax = slope_hi), width = 0.12, linewidth = 0.55) +
  geom_point(size = 2.1) +
  ## the last label goes to the left of its point so that the axis need not be stretched for it
  geom_text(aes(x = as.integer(tissue) + ifelse(tissue == "muscle", -0.14, 0.14),
                hjust = ifelse(tissue == "muscle", 1, 0),
                label = sprintf("slope %s\nSD %s", num(slope), num(sd_prolif))),
            vjust = 0.5, size = pt(7), family = FONT, lineheight = 1.05, show.legend = FALSE) +
  scale_colour_manual(values = TCOL, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = XN[TK], limits = c(0.6, 4.4), expand = c(0, 0)) +
  scale_y_continuous(labels = num_axis(), expand = expansion(mult = c(0.10, 0.62))) +
  labs(x = NULL, y = "slope: 96 splicing genes per unit\nproliferation score (95% CI)") +
  annotate("text", x = 1.55, y = Inf, hjust = 0, vjust = 1.08, family = FONT, size = pt(7),
           colour = INK2, lineheight = 1.15, label = vm_txt) +
  theme_sa()

## ---- S10b. coupling before and after composition adjustment (post hoc)
## every value, n and marker list on this panel comes from composition_adjusted.tsv
CC <- CA[set == "splicing 96"][match(TK, tissue)][, tissue := factor(tissue, levels = TK)]
stopifnot(nrow(CC) == 4, !anyNA(CC$rho_tech), !anyNA(CC$rho_tech_comp))
CL <- melt(CC[, .(tissue, before = rho_tech, after = rho_tech_comp)], id.vars = "tissue",
           variable.name = "k", value.name = "rho")
CL[, k := factor(k, levels = c("before","after"),
                 labels = c("technical covariates", "+ cell-composition markers"))]
CC[, ntxt := ifelse(n_tech == n_comp, sprintf("n = %d donors", n_comp),
                    sprintf("n = %d and %d donors", n_tech, n_comp))]
YN <- setNames(sprintf("%s\n%s", gsub("\n", " ", TN[as.character(CC$tissue)]), CC$ntxt), as.character(CC$tissue))
mk_full <- c(kerat = "keratinocyte", fibro = "fibroblast", immune = "immune",
             myofibre = "myofibre", satellite = "satellite cell")
mk_txt <- function(s) { p <- strsplit(s, "+", fixed = TRUE)[[1]]
  paste("+", paste(ifelse(p %in% names(mk_full), mk_full[p], p), collapse = ", ")) }
xr <- range(c(0, CL$rho)); XL <- xr + c(-0.10, 0.14) * diff(xr)
## the marker list of each tissue is written under its own row, below the value
## labels and on the side away from the points, so that it follows the table
## without running into a value
MK <- if ("markers" %in% names(CC)) CC[, .(tissue, txt = vapply(markers, mk_txt, ""),
        right = pmax(rho_tech, rho_tech_comp) < mean(XL))] else NULL
pSb <- ggplot(CL, aes(rho, tissue, colour = k)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
  geom_line(aes(group = tissue), colour = GREY, linewidth = 0.4) +
  geom_point(size = 2.1) +
  geom_text(aes(label = num(rho), vjust = ifelse(k == "technical covariates", -1.1, 1.9)),
            size = pt(7), family = FONT, show.legend = FALSE) +
  { if (!is.null(MK)) geom_text(data = MK, aes(x = ifelse(right, XL[2], 0.03), y = tissue, label = txt,
                                               hjust = ifelse(right, 1, 0)),
                                inherit.aes = FALSE, vjust = 3.7, size = pt(7), family = FONT, colour = INK2) } +
  annotate("text", x = 0.03, y = 4.5, hjust = 0, family = FONT, size = pt(7), colour = INK2, label = "post hoc") +
  scale_colour_manual(values = c(`technical covariates` = INK2, `+ cell-composition markers` = ORANGE),
                      name = NULL) +
  scale_y_discrete(limits = rev(TK), labels = YN, expand = expansion(add = c(0.75, 0.55))) +
  scale_x_continuous(labels = num_axis(), limits = XL, expand = c(0, 0)) +
  labs(x = sprintf("partial %s with proliferation, 96 splicing genes", RHO), y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.location = "plot", legend.justification = "right",
                     legend.margin = margin(b = -4), legend.key.size = unit(6, "pt"),
                     axis.text.y = element_text(size = 7.5, lineheight = 1.05),
                     plot.margin = margin(3, 10, 3, 3))

## ---- S10c. culture erases the donor, and the control that shows the design can
##            detect a donor (the old main panel d). The concordance of programme
##            scores is the preregistered secondary analysis S3; the genotype-driven
##            positive control was added post hoc and is labelled so.
band <- 1.96 / sqrt(NPAIR - 3)
S3 <- S3[order(rho)]; stopifnot(all(S3$n_donors == NPAIR))
## dots are stacked greedily: each goes into the lowest row in which it keeps its distance
SEP <- 0.024; rows <- list(); S3[, k := NA_integer_]
for (i in seq_len(nrow(S3))) { k <- 1L
  while (k <= length(rows) && any(abs(rows[[k]] - S3$rho[i]) < SEP)) k <- k + 1L
  rows[[k]] <- c(if (k <= length(rows)) rows[[k]], S3$rho[i]); S3$k[i] <- k }
S3[, yy := 0.66 + (k - 1) * 0.085]
PCg <- PC[in_culture == TRUE]; stopifnot(all(PCg$n == NPAIR))
sexrng <- range(PCg[group == "sex"]$rho); eq <- PCg[group == "eqtl"][order(-rho)]
PTS <- rbind(
  data.table(x = S3$rho, y = S3$yy, what = "programme score"),
  data.table(x = PCg[group == "sex"]$rho, y = 0.14, what = "sex-chromosome gene"),
  data.table(x = eq$rho, y = 0.36, what = "structural-variant gene"))
pSc <- ggplot() +
  annotate("rect", xmin = -band, xmax = band, ymin = -Inf, ymax = Inf, fill = GREY_L, alpha = 0.55) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_point(data = PTS, aes(x, y, colour = what, shape = what), size = 1.6, alpha = 0.9) +
  annotate("text", x = -0.52, y = 1.18, hjust = 0, vjust = 1, family = FONT, size = pt(7), colour = INK,
           lineheight = 1.1,
           label = sprintf("%d programme scores\n%s = %s to %s\n%s FDR < 0.05", nrow(S3), RHO,
                           num(min(S3$rho)), num(max(S3$rho)),
                           ifelse(any(S3$FDR < 0.05), sprintf("%d with", sum(S3$FDR < 0.05)), "none with"))) +
  annotate("text", x = max(eq$rho) - 0.025, y = 0.36, hjust = 1, family = FONT, size = pt(7), colour = BLUE,
           parse = TRUE, label = sprintf("italic('%s')~'(%s = %s)'", paste(eq$gene, collapse = ", "), RHO,
                                         paste(num(eq$rho), collapse = ", "))) +
  annotate("text", x = sexrng[1] - 0.025, y = 0.14, hjust = 1, family = FONT, size = pt(7), colour = TEAL,
           label = sprintf("%d sex-chromosome genes (%s = %s to %s)", sum(PCg$group == "sex"), RHO,
                           num(sexrng[1]), num(sexrng[2]))) +
  annotate("text", x = -0.52, y = 0.25, hjust = 0, family = FONT, size = pt(7), colour = INK2,
           lineheight = 1.1, label = "positive control, post hoc:\ngenotype-driven genes") +
  annotate("text", x = 0.86, y = 1.18, hjust = 1, vjust = 1, family = FONT, size = pt(7), colour = INK2,
           lineheight = 1.1,
           label = sprintf("n = %d donors who gave both samples\ngrey band, 95%% range under no association", NPAIR)) +
  scale_colour_manual(values = c(`programme score` = INK2, `sex-chromosome gene` = TEAL,
                                 `structural-variant gene` = BLUE), guide = "none") +
  scale_shape_manual(values = c(`programme score` = 16, `sex-chromosome gene` = 17,
                                `structural-variant gene` = 15), guide = "none") +
  scale_x_continuous(sprintf("Spearman %s between a donor's cultured fibroblasts and the same donor's sun-exposed skin", RHO),
                     limits = c(-0.55, 0.9), breaks = seq(-0.4, 0.8, 0.2), labels = num_axis()) +
  scale_y_continuous(NULL, limits = c(-0.02, 1.22), breaks = NULL) +
  theme_sa() + theme(axis.line.y = element_blank())

save_fig(plot_grid(lab_grid(pSa, pSb, labels = c("a", "b"), ncol = 2, rel_widths = c(0.95, 1.1)),
                   lab_grid(pSc, labels = "c", ncol = 1), ncol = 1, rel_heights = c(84, 58)),
         "FigS10.png", 183, 142)
