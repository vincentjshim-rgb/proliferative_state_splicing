## Fig. 7. What splicing-factor abundance reports: summary of the study (2026-09-19).
##   a  schematic of the reading the paper arrives at, carrying the paper's own numbers
##   b  the coupling between the splicing set and proliferation in every dataset, in the
##      order of the Results (Spearman rho with a Fisher-z 95% CI). Every value in panel b is
##      computed here from the deposited tables; the numbers printed in panel a are the
##      manuscript's and are checked against the same tables before the figure is written.
## Skeletal muscle, where the proliferation score recovers no programme, is reported in
## Fig. 6b and Supplementary Table 9 and is not drawn here.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"
## The schematic is laid out from text extents, which ggplot measures on a null pdf device
## before the PNG device exists; register the family there with the Helvetica metrics it is
## built on (Nimbus Sans is metrically Helvetica), so the measurement is exact and silent.
if (!FONT %in% names(grDevices::pdfFonts()))
  do.call(grDevices::pdfFonts, setNames(list(grDevices::Type1Font(FONT,
    c("Helvetica.afm", "Helvetica-Bold.afm", "Helvetica-Oblique.afm",
      "Helvetica-BoldOblique.afm", "Symbol.afm"))), FONT))

## ---------------- b. the same coupling in every dataset ----------------------------
fz <- function(r, n) { z <- atanh(r); se <- 1/sqrt(n - 3)
  c(lo = tanh(z - qnorm(0.975) * se), hi = tanh(z + qnorm(0.975) * se)) }
DR <- fread(file.path(D, "division_rate/sample_division_rate.tsv"))                    # Fig. 2a
CC <- fread(file.path(D, "hallmark_revised/division_rate_cellcycle_removed.tsv"))       # Fig. 5a
DS <- fread(file.path(D, "machinery_specificity/donor_scores.tsv"))[primary == TRUE]    # Fig. 3a
LO <- fread(file.path(D, "revision_stats/leave_out_defining_series.tsv"))[set == "all contrasts"]  # Fig. 4a
GT <- fread(file.path(D, "gtex_boundary/programme_by_tissue.tsv"))[programme == "splicing 96"]     # Fig. 6a
CA <- fread(file.path(D, "gtex_posthoc/composition_adjusted.tsv"))[set == "splicing 96"]           # Fig. 6d
CR <- fread(file.path(D, "cohort_revised/sensitivity_metrics.tsv"))
stopifnot(nrow(DR) == 328, nrow(DS) == 107, LO$n_contrasts == 63, nrow(GT) == 4)

g <- function(tk, col) GT[tissue == tk][[col]]
row <- function(dataset, unit, r, n, lo = NA, hi = NA, tissue = "culture", adj = NA_real_) {
  if (is.na(lo)) { k <- fz(r, n); lo <- k[["lo"]]; hi <- k[["hi"]] }
  data.table(dataset = dataset, unit = unit, r = r, n = n, lo = lo, hi = hi, tissue = tissue, adj = adj) }
B <- rbind(
  row("Counted division rate, 177-gene set",
      sprintf("%d libraries, 7 cell lines (Fig. 2a)", nrow(DR)),
      cor(DR$rate, DR$splice, method = "spearman"), nrow(DR)),
  row("Same libraries, unselected Reactome definition",
      sprintf("mRNA splicing, %d libraries (Fig. 5a)", nrow(DR)),
      CC[programme == "mRNA splicing", rho_rate][1], nrow(DR)),
  row("Donor cohort, proliferation score",
      sprintf("%d adult donors (Fig. 3a)", nrow(DS)),
      cor(DS$prolif, DS$splicing, method = "spearman"), nrow(DS)),
  row("63 contrasts, Δ splicing set vs Δ cell-cycle genes",
      sprintf("%d GEO series (Fig. 4a)", LO$n_series),
      LO$rho, LO$n_contrasts, LO$lo, LO$hi),
  row("GTEx cultured fibroblasts, 96-gene set",
      sprintf("%d donors, partial %s (Fig. 6a)", g("culture", "n"), RHO),
      g("culture", "rho_prolif"), g("culture", "n")),
  row("GTEx skin, sun-exposed",
      sprintf("%d donors, partial %s (Fig. 6a)", g("legskin", "n"), RHO),
      g("legskin", "rho_prolif"), g("legskin", "n"), tissue = "skin",
      adj = CA[tissue == "legskin", rho_tech_comp]),
  row("GTEx skin, not exposed",
      sprintf("%d donors, partial %s (Fig. 6a)", g("pubskin", "n"), RHO),
      g("pubskin", "rho_prolif"), g("pubskin", "n"), tissue = "skin",
      adj = CA[tissue == "pubskin", rho_tech_comp]))
B[, lab := paste0(dataset, "\n", unit)]
B[, lab := factor(lab, levels = rev(lab))]
B[, xlab := pmax(hi, r) + 0.025]
## study-level value for the contrast row: one averaged value per study, secretome
## preparations counted separately, exactly as run_stats_supplements.R computes it
CT <- fread(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
CT[, study := sub("_.*$", "", id)]
ST <- CT[, .(spl = mean(spl), cc = mean(cc)), by = study]
r_study <- cor(ST$cc, ST$spl, method = "spearman"); n_study <- nrow(ST)
stopifnot(n_study == 29, abs(r_study - 0.79) < 0.006)
B[, study := NA_real_]; B[4, study := r_study]
fwrite(B[, .(dataset = gsub("\n", " ", dataset), unit, n, rho = r, lo, hi, tissue,
             rho_after_composition_markers = adj, rho_one_value_per_study = study,
             n_studies = ifelse(is.na(study), NA_integer_, n_study))],
       file.path(PUBDIR, "fig7_coupling_strip.tsv"), sep = "\t")
cat("Fig. 7b values:\n"); print(B[, .(dataset = gsub("\n", " ", dataset), n, r = round(r, 3),
                                      lo = round(lo, 3), hi = round(hi, 3), adj = round(adj, 3))])

pB <- ggplot(B, aes(x = r, y = lab)) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_col(aes(fill = tissue), width = 0.6, orientation = "y") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.22, linewidth = 0.35, colour = INK,
                orientation = "y") +
  geom_point(data = B[!is.na(adj)], aes(x = adj), shape = 21, fill = "white", colour = INK,
             size = 1.9, stroke = 0.45) +
  geom_point(data = B[!is.na(study)], aes(x = study), shape = 23, fill = INK, colour = INK, size = 1.9) +
  geom_text(aes(x = xlab, label = num(r)), hjust = 0, size = pt(7), family = FONT, colour = INK) +
  annotate("point", x = 0.50, y = 2.5, shape = 21, fill = "white", colour = INK, size = 1.9, stroke = 0.45) +
  annotate("text", x = 0.525, y = 2.5, hjust = 0, size = pt(7), family = FONT, colour = INK2,
           label = "after cell-composition markers (post hoc)") +
  annotate("point", x = 0.50, y = 3.45, shape = 23, fill = INK, colour = INK, size = 1.9) +
  annotate("text", x = 0.525, y = 3.45, hjust = 0, size = pt(7), family = FONT, colour = INK2,
           label = sprintf("one value per study (%s = %s, n = %d)", RHO, num(r_study), n_study)) +
  scale_fill_manual(values = c(culture = BLUE, skin = ORANGE),
                    labels = c(culture = "cultured fibroblasts", skin = "bulk skin"), name = NULL) +
  scale_x_continuous(sprintf("Spearman %s with the proliferation measure named in each row (partial %s for the GTEx rows)\n95%% confidence interval from Fisher's z", RHO, RHO),
                     limits = c(-0.12, 1.06), breaks = seq(-0.2, 1, 0.2), labels = num_axis(1),
                     expand = c(0, 0)) +
  labs(y = NULL) +
  theme_sa() +
  theme(axis.text.y = element_text(hjust = 0, lineheight = 0.92), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.80, 0.11), legend.key.height = unit(6, "pt"),
        legend.spacing.y = unit(1, "pt"))

## ---------------- a. schematic, carrying the paper's numbers -----------------------
## the numbers printed below are checked against the tables that produced them
chk <- function(x, y, tol = 0.006) stopifnot(abs(x - y) < tol)
chk(0.71, B$r[1]); chk(0.86, B$r[4], 0.01); chk(0.74, B$r[5]); chk(0.23, B$r[6]); chk(0.06, B$r[7])
chk(0.10, B$adj[6]); chk(-0.06, B$adj[7])
CCY <- fread(file.path(D, "gtex_boundary/programme_by_tissue.tsv"))[programme == "cell cycle (positive control)"]
chk(0.67, CCY[tissue == "legskin", rho_prolif]); chk(0.60, CCY[tissue == "pubskin", rho_prolif])
prim <- CR[metric == "machinery" & cohort == "primary" & covariates == TRUE]
stopifnot(nrow(prim) == 1); chk(-0.133, prim$beta, 0.002); chk(-0.040, prim$beta_adj, 0.002)
chk(70, prim$pct_lost, 1)

box <- function(x1, y1, x2, y2, fill = "white", col = INK, lwd = 0.35, lty = "solid")
  annotate("rect", xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = fill, colour = col,
           linewidth = lwd, linetype = lty)
tx <- function(x, y, s, size = 7, hjust = 0.5, vjust = 0.5, face = "plain", col = INK, lh = 0.95)
  annotate("text", x = x, y = y, label = s, size = pt(size), hjust = hjust, vjust = vjust,
           fontface = face, colour = col, family = FONT, lineheight = lh)
arr <- function(x1, y1, x2, y2, lwd = 0.6, col = INK, lty = "solid")
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2, linewidth = lwd, colour = col,
           linetype = lty, arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
           linejoin = "mitre")
CULT <- "#E9F0F8"; SKIN <- "#FBEFE0"; BAND <- "#EEF1F4"
c1 <- 30; c2 <- 90; c3 <- 150            # column centres (mm); columns are 54 mm wide
L1 <- c1 - 27; R1 <- c1 + 27; L2 <- c2 - 28; R2 <- c2 + 28; L3 <- c3 - 27; R3 <- c3 + 27

pA <- ggplot() +
  ## column headers
  tx(c1, 83, "Cultured human fibroblasts", 8, face = "bold") +
  tx(c2, 83, "The age effect in culture", 8, face = "bold") +
  tx(c3, 83, "Bulk skin, same donor pool (GTEx)", 8, face = "bold") +
  ## ---- column 1: what moves proliferative state, and what follows it
  box(L1, 66.5, R1, 79.5, fill = CULT) +
  tx(L1 + 2, 76.6, "slower", 7, hjust = 0, face = "bold") +
  tx(L1 + 17, 76.6, "donor age, passage,\ncontact inhibition", 7, hjust = 0) +
  tx(L1 + 2, 71.6, "faster", 7, hjust = 0, face = "bold") +
  tx(L1 + 17, 71.6, "most secretome preparations", 7, hjust = 0) +
  tx(L1 + 2, 68.4, "either way", 7, hjust = 0, face = "bold") +
  tx(L1 + 17, 68.4, "reprogramming media", 7, hjust = 0) +
  arr(c1, 66.5, c1, 62.8) +
  box(L1, 49.5, R1, 62.5, fill = CULT) +
  tx(c1, 60.3, "Proliferative state", 8, face = "bold") +
  tx(c1, 54.4, "counted division rate, or a 20-marker\ntranscriptional score\n(checked against counting, Fig. 2e)", 7, lh = 0.92) +
  arr(L1 + 11, 49.5, L1 + 11, 41.8) +
  tx(L1 + 14, 45.6, "growth programme:\ntranscribed with the cell cycle", 7, hjust = 0) +
  box(L1, 30, R1, 41.5, fill = CULT) +
  tx(c1, 39, "Splicing-machinery transcripts", 8, face = "bold") +
  tx(c1, 34.5, "177 pre-mRNA processing genes:\nspliceosome, capping, 3′-end processing", 7) +
  tx(c1, 22.8, sprintf("%s = 0.71 with the counted division rate\n63 contrasts: %s = 0.86\n8 of 23 ageing programmes follow the rate too\ncollagen formation indifferent (H5, %s = 0.05)", RHO, RHO, RHO), 7, col = INK2) +
  ## ---- column 2: what the age effect in culture contains
  tx(c2, 77.5, "107 adult donors: −0.133 SD per decade", 7) +
  annotate("rect", xmin = L2, xmax = L2 + 0.70 * 56, ymin = 68, ymax = 74, fill = BLUE, colour = NA) +
  annotate("rect", xmin = L2 + 0.70 * 56, xmax = R2, ymin = 68, ymax = 74, fill = GREY_L, colour = NA) +
  tx(L2 + 0.35 * 56, 71, "via proliferative state, 70%", 7, col = "white") +
  tx(L2 + 0.85 * 56, 71, "residual", 7) +
  tx(L2 + 0.85 * 56, 65.8, "−0.040 per decade", 7, col = INK2) +
  tx(c2, 57.5, "The residual is not robust:\nit halves with a second programme in the model,\nand in GTEx cultures it appears only by suppression.\nNo age effect independent of proliferation is reported.", 7) +
  tx(c2, 44, "The same adjustment removes a median 56% of the\nage effect of the 205 expression-matched random sets\nthat have one; none of 1,000 reaches this size.", 7, col = INK2) +
  tx(c2, 32.5, "Splicing outcome (unannotated junctions) in the\nsame libraries carries no age association.", 7, col = INK2) +
  ## ---- column 3: bulk skin
  box(L3, 67, R3, 79.5, fill = SKIN) +
  tx(c3, 73.3, "most cells are not cycling;\nproliferation varies little between donors", 7) +
  arr(c3, 67, c3, 62.3, lwd = 0.45, col = GREY, lty = "22") +
  box(L3, 52, R3, 62, fill = SKIN) +
  tx(c3, 59, "Proliferative state", 8, face = "bold") +
  tx(c3, 55, "the same score still recovers the mitotic\ncell cycle (0.67 sun-exposed, 0.60 not exposed)", 7) +
  arr(L3 + 11, 52, L3 + 11, 42.3, lwd = 0.45, col = GREY, lty = "22") +
  tx(L3 + 14, 47, "coupling weakens", 7, hjust = 0, col = INK2) +
  box(L3, 30, R3, 42, fill = SKIN) +
  tx(c3, 39.6, "Splicing-machinery transcripts", 8, face = "bold") +
  tx(c3, 34.2, sprintf("partial %s = 0.23 and 0.06 at the two sites;\n0.10 and −0.06 after cell-composition\nmarkers (post hoc)", RHO), 7, lh = 0.92) +
  tx(c3, 22.8, "collagen formation stays small (0.03 and 0.09;\nthe preregistered test, H5, was in culture);\nage effects within GTEx tissues are not read\n(the collagen positive control, H6, failed)", 7, col = INK2) +
  ## ---- reading
  box(3, 3, 177, 13.5, fill = BAND, col = NA) +
  tx(6, 8.3, "Reading", 8, hjust = 0, face = "bold") +
  tx(22, 8.3, "In a cultured fibroblast, a change in splicing-factor transcripts is first of all a change in proliferative state;\nin skin the coupling weakens, and the claim stops at the culture dish until it is shown in the tissue of interest.", 7, hjust = 0) +
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 86), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(1, 1, 1, 1))

save_fig(plot_grid(lab_grid(pA, labels = "a", ncol = 1),
                   lab_grid(pB, labels = "b", ncol = 1),
                   ncol = 1, rel_heights = c(86, 66)), "Fig7.png", 180, 152)
