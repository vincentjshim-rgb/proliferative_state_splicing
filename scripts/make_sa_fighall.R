## Fig. 5 (six-figure restructure 2026-09-17; was Fig. 6, before that Fig. 7): which
## ageing-associated programmes in cultured fibroblasts are readouts of division
## rate, and which are not?
## Revised 2026-09-15:
##   - every programme is also scored after removing the Reactome cell-cycle genes
##     it contains (reviewer 1, major 7a), which is now panel a;
##   - the donor-age half uses the redefined cohort (107 normal adults, repository,
##     instrument and sex covariates), in which collagen formation is no longer
##     age-associated, so the old "age-associated yet division-independent" reading
##     is gone (panels b and e);
##   - programmes with FDR >= 0.05 are drawn as open circles instead of alpha = 0,
##     which made them invisible.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
O  <- "public_data_tierA/derived/hallmark_audit"
OR <- "public_data_tierA/derived/hallmark_revised"
CCR <- read.delim(file.path(OR, "division_rate_cellcycle_removed.tsv"))
AGE <- read.delim(file.path(OR, "donor_age_primary_cohort.tsv"))
SGR <- read.delim(file.path(OR, "senescence_dissection_cellcycle_union.tsv"))
SG  <- read.delim(file.path(O, "senescence_geneset_dissected.tsv"))
IN <- readRDS(file.path(O, "audit_inputs.rds"))
z <- IN$z; d0 <- IN$meta; P <- IN$pathways
NPRIM <- 107L
## the shared num_axis() passes the raw break vector to format(), which turns
## ggplot's floating-point zero into 1.1e-16; fix it locally
num_axis <- function(digits = 1) function(x)
  sub("-", MINUS, formatC(round(x, digits + 2), format = "f", digits = digits))

CLASS <- c(
 `chromatin organisation`="RNA / chromatin", `cellular senescence`="composite",
 `pre-mRNA processing`="RNA / chromatin", `mRNA splicing`="RNA / chromatin",
 `cell cycle (positive control)`="composite", `interferon signalling`="immune",
 `DNA double-strand repair`="genome upkeep", `telomere maintenance`="genome upkeep",
 `base excision repair`="genome upkeep", `interleukin signalling`="immune",
 `rRNA processing`="RNA / chromatin", `mitochondrial translation`="metabolic",
 `translation`="RNA / chromatin", `TCA cycle`="metabolic",
 `respiratory electron transport`="metabolic", `proteasome`="proteostasis",
 `chaperone / HSF1`="proteostasis", `autophagy`="proteostasis",
 `mTOR signalling`="proteostasis", `unfolded protein response`="proteostasis",
 `collagen formation`="matrix", `NF-kB / TNF`="immune",
 `insulin / IGF signalling`="metabolic", `extracellular matrix`="matrix")
PAL <- c(`RNA / chromatin`=BLUE, `genome upkeep`=PURPLE, proteostasis=GREEN,
         matrix=BROWN, immune=RED, metabolic=TEAL, composite=GREY)

J <- merge(CCR, AGE[, c("programme", "rho_age", "beta_decade", "p_age", "FDR_age",
                        "beta_adj", "FDR_adj")], by = "programme", all.x = TRUE)
J$class <- factor(CLASS[J$programme], levels = names(PAL))
J$short <- sub("cell cycle \\(positive control\\)", "cell cycle, mitotic", J$programme)
CCrho <- J$rho_rate[J$programme == "cell cycle (positive control)"]
ord <- order(J$rho_rate)
J$short <- factor(J$short, levels = J$short[ord])
N <- nrow(J)

## ---------------- a. division rate, before and after removing cell-cycle genes ----
## the two programmes whose gene sets are entirely inside the cell-cycle union
## cannot be scored after removal; they are marked rather than dropped
J$noCC_missing <- !is.finite(J$rho_rate_noCC)
## two-column heatmap (2026-09-19) in place of the dumbbell plot: the programme x
## condition matrix with printed values is the form used for programme-level
## summaries in the ageing transcriptome literature
if (!exists("DIVERGE")) DIVERGE <- c(BLUE, "#F7F7F7", RED)
JH <- rbind(data.frame(short = J$short, class = J$class, col = "all genes",
                       rho = J$rho_rate, miss = FALSE),
            data.frame(short = J$short, class = J$class, col = "without\ncell-cycle genes",
                       rho = J$rho_rate_noCC, miss = J$noCC_missing))
JH$col <- factor(JH$col, levels = c("all genes", "without\ncell-cycle genes"))
JH$lab <- ifelse(JH$miss, "\u2020", ifelse(is.finite(JH$rho), sub("-", MINUS, sprintf("%.2f", JH$rho)), ""))
JH$rho[JH$miss] <- NA
LABCOL <- unname(PAL[as.character(J$class[match(levels(droplevels(J$short)), as.character(J$short))])])
pA <- ggplot(JH, aes(col, short)) +
  geom_tile(aes(fill = pmax(pmin(rho, 1), -1)), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), family = FONT, size = pt(6.5), colour = INK) +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-1, 1), breaks = c(-1, 0, 1),
                       labels = c(paste0(MINUS, "1"), "0", "1"), na.value = "#E6E8EB",
                       name = sprintf("%s with measured division rate", RHO),
                       guide = guide_colourbar(title.position = "top", barwidth = unit(48, "pt"),
                                               barheight = unit(4, "pt"), ticks.colour = "white")) +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = "\u2020 all measured genes are cell-cycle genes\nrow colour: programme class (key in b)") +
  theme_sa() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 7.2, colour = LABCOL),
        axis.text.x = element_text(size = 7, lineheight = 0.95),
        legend.position = "bottom", legend.key.height = unit(5, "pt"),
        legend.key.width = unit(28, "pt"), legend.title = element_text(size = 7),
        legend.text = element_text(size = 7), legend.margin = margin(t = 2),
        plot.caption = element_text(size = 6.5, colour = INK2, hjust = 0, lineheight = 1.05),
        plot.margin = margin(3, 4, 3, 3))

## ---------------- b. division readout against donor-age readout -------------
J$sig <- J$FDR_age < 0.05
LP <- data.frame(
 prog = c("extracellular matrix", "insulin / IGF signalling", "NF-kB / TNF", "proteasome",
          "unfolded protein response", "collagen formation",
          "mRNA splicing", "cellular senescence", "cell cycle (positive control)"),
 lab = c("extracellular matrix", "insulin / IGF", "NF-kB / TNF", "proteasome",
         "unfolded protein response", "collagen formation",
         "mRNA splicing", "cellular\nsenescence", "cell cycle"),
 ## 2026-09-17: every label now sits beside its own point, on one side of the zero
 ## lines, with the leader running up to it. Before, "cell cycle" was printed level
 ## with the cellular-senescence point and three labels crossed the vertical zero line.
 lx = c(-0.45, -0.45, -0.03,  0.25,  0.03, -0.45,  0.94,  0.71,  0.75),
 ly = c(-0.11,  0.15,  0.22,  0.05, -0.56, -0.52, -0.61, -0.371, -0.43),
 hj = c( 0,     0,     1,     0,     0,     0,     1,     0,     0),
 sx = c(-0.34, -0.30, -0.14,  0.245, 0.07, -0.05,  0.70,  0.70,  0.74),
 sy = c(-0.13,  0.13,  0.20,  0.035,-0.54, -0.50, -0.59, -0.371, -0.435))
LP <- merge(LP, J[, c("programme", "rho_rate", "rho_age", "class")],
            by.x = "prog", by.y = "programme")
pB <- ggplot(J, aes(rho_rate, rho_age)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_segment(data = LP, aes(x = rho_rate, y = rho_age, xend = sx, yend = sy),
               colour = GREY_L, linewidth = 0.22) +
  geom_point(aes(colour = class, fill = ifelse(sig, "sig", "ns")), size = 1.9,
             shape = 21, stroke = 0.6) +
  geom_text(data = LP, aes(lx, ly, label = lab, colour = class, hjust = hj),
            family = FONT, size = pt(7), vjust = 0.5, lineheight = 0.95) +
  scale_colour_manual(values = PAL, name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.7))) +
  scale_fill_manual(values = c(sig = INK2, ns = "white"), guide = "none") +
  scale_x_continuous(sprintf("%s with measured division rate (328 libraries)", RHO),
      limits = c(-0.47, 0.95), breaks = seq(-0.4, 0.8, 0.2), labels = num_axis()) +
  scale_y_continuous(sprintf("%s with donor age (107 donors; open, FDR > 0.05)", RHO),
      limits = c(-0.64, 0.32), breaks = seq(-0.6, 0.2, 0.2), labels = num_axis()) +
  theme_sa() +
  theme(legend.position = "bottom", legend.justification = "left", legend.margin = margin(t = -4),
        legend.text = element_text(size = 7), legend.key.size = unit(7, "pt"))

## ---------------- c. the senescence set is largely a cell-cycle set ---------
SGR$group <- factor(SGR$group, levels = c("also in cell cycle", "senescence set only"))
nlab <- sprintf("%s\nn = %d", sub(" ", "\n", levels(SGR$group)), as.integer(table(SGR$group)))
AGG   <- cor(colMeans(z[as.character(SGR$gene), ]), d0$rate, method = "spearman")
AGGnc <- cor(colMeans(z[as.character(SGR$gene[SGR$group == "senescence set only"]), ]),
             d0$rate, method = "spearman")
wt <- wilcox.test(rho ~ group, data = SGR)
## SenMayo (Saul et al. 2022), the panel used most often in intervention work, scored
## the same way: it contains almost no cell-cycle genes and shows no positive coupling
SMT <- read.delim("public_data_tierA/derived/reviewer_sensitivities/senmayo_division_rate.tsv")
SMrho <- SMT$rho_set[SMT$set == "SenMayo (Saul 2022)"]
SMn   <- SMT$n_measured[SMT$set == "SenMayo (Saul 2022)"]
SMcc  <- SMT$n_cellcycle[SMT$set == "SenMayo (Saul 2022)"]
pC <- ggplot(SGR, aes(group, rho)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3, linetype = "21") +
  geom_hline(yintercept = AGG, colour = RED, linewidth = 0.45) +
  geom_hline(yintercept = AGGnc, colour = ORANGE, linewidth = 0.45, linetype = "21") +
  geom_violin(aes(fill = group), colour = NA, alpha = 0.20, width = 0.78, scale = "width") +
  geom_point(aes(colour = group), size = 0.85, alpha = 0.7,
             position = position_jitter(width = 0.13, height = 0, seed = 1)) +
  geom_boxplot(width = 0.17, outlier.shape = NA, fill = NA, colour = INK, linewidth = 0.32) +
  ## 2026-09-17: the three set-level values are labelled in a gutter to the right of
  ## the violins; they used to be printed across the points and the box
  annotate("text", x = 2.45, y = AGG, vjust = -0.25, hjust = 0, family = FONT, size = pt(7),
      colour = RED, lineheight = 0.95,
      label = sprintf("all %d genes\n%s = %s", nrow(SGR), RHO, num(AGG))) +
  annotate("text", x = 2.45, y = AGGnc, vjust = -0.25, hjust = 0, family = FONT, size = pt(7),
      colour = ORANGE, lineheight = 0.95,
      label = sprintf("%d non-cell-cycle\ngenes, %s = %s", sum(SGR$group == "senescence set only"),
                      RHO, num(AGGnc))) +
  geom_hline(yintercept = SMrho, colour = BLUE, linewidth = 0.45, linetype = "12") +
  annotate("text", x = 2.45, y = SMrho, vjust = 1.25, hjust = 0, family = FONT, size = pt(7),
      colour = BLUE, lineheight = 0.95,
      label = sprintf("SenMayo panel\n%s = %s", RHO, num(SMrho))) +
  scale_fill_manual(values = c(GREY, GREY_L), guide = "none") +
  scale_colour_manual(values = c(INK2, GREY), guide = "none") +
  scale_x_discrete(NULL, labels = nlab, expand = expansion(add = c(0.50, 2.00))) +
  scale_y_continuous(sprintf("gene-wise %s with division rate", RHO), limits = c(-0.72, 1.02),
      breaks = seq(-0.6, 0.8, 0.2), labels = num_axis()) +
  sigbar(1, 2, 0.88, pfmt(wt$p.value), 0.04, size = 7) +
  theme_sa() + theme(axis.text.x = element_text(size = 7.2, lineheight = 0.95),
                     ## the reference lines stop short of the gene names of panel d
                     plot.margin = margin(3, 11, 3, 3))

## ---------------- d. canonical senescence markers behave correctly ----------
EXP <- c(LMNB1 = "lost in senescence", E2F1 = "lost in senescence",
         CDKN1A = "induced (p21)", CDKN2A = "induced (p16)", CDKN2B = "induced (p15)")
M <- SG[SG$marker, ]; M <- M[order(M$rho), ]
M$gene <- factor(M$gene, levels = M$gene)
M$exp  <- EXP[as.character(M$gene)]
M$col  <- ifelse(is.na(M$exp), GREY, ifelse(M$rho > 0, BLUE, RED))
pD <- ggplot(M, aes(rho, gene)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(fill = M$col, width = 0.66) +
  geom_text(aes(x = ifelse(rho > 0, -0.05, 0.05), label = ifelse(is.na(exp), "", exp),
                hjust = ifelse(rho > 0, 1, 0)), family = FONT, size = pt(7), colour = INK2) +
  scale_x_continuous(sprintf("%s with division rate", RHO), limits = c(-0.95, 0.95),
      breaks = seq(-0.8, 0.8, 0.4), labels = num_axis()) +
  annotate("text", x = -0.92, y = 13.2, hjust = 0, vjust = 0.5, family = FONT, size = pt(7),
      colour = RED, label = "higher in slow") +
  annotate("text", x = 0.92, y = 13.2, hjust = 1, vjust = 0.5, family = FONT, size = pt(7),
      colour = BLUE, label = "higher in fast") +
  theme_sa() + theme(axis.title.y = element_blank(), axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), axis.text.y = element_text(size = 7.2, face = "italic"),
      plot.margin = margin(12, 3, 3, 3)) + coord_cartesian(clip = "off")

## ---------------- e. the two behaviours, shown as data ----------------------
## division rate from the 328 counted libraries; donor age in the redefined cohort,
## scored from the same recount3 matrix used everywhere else (cached on first run)
CACHE <- file.path(PUBDIR, "cache_fig6_primary_scores.rds")
SPL <- "mRNA Splicing - Major Pathway"; COL <- "Collagen formation"
if (file.exists(CACHE)) { PS <- readRDS(CACHE) } else {
  suppressPackageStartupMessages({library(data.table); library(edgeR)})
  source("scripts/revision/cohort_gse113957.R")
  DON <- gse113957_donors(); DON <- DON[gse113957_keep(DON, "primary"), ]
  gs <- fread(cmd = "zcat public_data_tierA/recount3/sra.gene_sums.SRP144355.G026.gz | tail -n +2")
  GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
  mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
  key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
  rownames(GM) <- key[rownames(GM)]
  GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
  GM <- GM[, intersect(DON$srr, colnames(GM)), drop = FALSE]
  DON <- DON[match(colnames(GM), DON$srr), ]
  y <- DGEList(GM); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
  zz <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
  PS <- list(age = DON$age,
             spl = colMeans(zz[intersect(rownames(zz), P[[SPL]]), , drop = FALSE]),
             col = colMeans(zz[intersect(rownames(zz), P[[COL]]), , drop = FALSE]))
  saveRDS(PS, CACHE)
}
mk <- function(sc, xv, pmodel = NULL) { ct <- cor.test(sc, xv, method = "spearman", exact = FALSE)
  ## for donor age the P printed is the one from the cohort model (log depth,
  ## repository, instrument, sex), not the raw correlation, so that this panel
  ## agrees with panel b, where collagen formation is not significant
  data.frame(x = xv, y = sc, lab = rp(ct$estimate, if (is.null(pmodel)) ct$p.value else pmodel)) }
pmod <- function(nm) AGE$p_age[match(nm, AGE$programme)]
sp <- function(Dt, xlab, ylab, col, sub) ggplot(Dt, aes(x, y)) +
  geom_point(colour = col, size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, colour = col, fill = col, alpha = 0.14, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = Dt$lab[1], hjust = -0.12, vjust = 1.25,
           family = FONT, size = pt(7), colour = INK, lineheight = 0.95) +
  ## head-room so that the statistic is not printed over the points
  scale_y_continuous(labels = num_axis(), expand = expansion(mult = c(0.05, 0.45))) +
  labs(x = xlab, y = ylab, subtitle = sub) + theme_sa() +
  theme(plot.subtitle = element_text(size = 7.2, colour = INK, hjust = 0.5, margin = margin(b = 1)),
        plot.margin = margin(2, 4, 2, 2), axis.title = element_text(size = 7.5))
sc_spl <- colMeans(z[intersect(rownames(z), P[[SPL]]), , drop = FALSE])
sc_col <- colMeans(z[intersect(rownames(z), P[[COL]]), , drop = FALSE])
pE <- plot_grid(
  sp(mk(sc_spl, d0$rate), NULL, "mRNA splicing", BLUE, sprintf("n = %d libraries", ncol(z))),
  sp(mk(PS$spl, PS$age, pmod("mRNA splicing")),  NULL, NULL,  BLUE, sprintf("n = %d donors", NPRIM)),
  sp(mk(sc_col, d0$rate), "divisions per day", "collagen formation", BROWN, NULL),
  sp(mk(PS$col, PS$age, pmod("collagen formation")),  "donor age (years)",  NULL,   BROWN, NULL),
  ncol = 2, align = "hv", axis = "tblr", rel_heights = c(1.16, 1))

top <- lab_grid(pA, pB, labels = c("a", "b"), ncol = 2, rel_widths = c(0.9, 1.1))
## panel c is wider than before so that its gutter holds the three set-level labels
bot <- lab_grid(pC, pD, pE, labels = c("c", "d", "e"), ncol = 3, rel_widths = c(1.07, 1.00, 1.28))
save_fig(plot_grid(top, bot, ncol = 1, rel_heights = c(1.34, 1)), "Fig6.png", 183, 164)
