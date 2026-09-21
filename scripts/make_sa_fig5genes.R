## Fig. 5  Six named splicing factors, gene by gene, across the ageing samples (2026-09-21).
##
## The form used for ageing transcriptomes (Fleischer 2018; the GTEx age analyses; Holly 2013):
## the expression of individual genes along time in culture and across age groups. One row per
## gene, one column per dataset; the top row is the proliferation score in the same samples, so
## that a reader sees whether the genes fall where proliferation falls.
##   a  the counted resource: untreated cultures of HC1 and HC2 along days in culture
##   b  the 107-donor cohort by decade of age (points coloured by proliferation third)
##   c  GTEx cultured fibroblasts by age band
##   d  GTEx sun-exposed skin by age band
## Expression is standardised within each dataset (z of log2 CPM). Every number printed below
## is quoted in the Results paragraph that cites this figure.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table); library(edgeR)})
D <- "public_data_tierA/derived"
GENES <- c("HNRNPD", "HNRNPA1", "HNRNPK", "SRSF1", "SRSF2", "SRSF3")   # named factors inside the 177
say <- function(...) cat(sprintf(...), "\n")
sp <- function(x, y) { ct <- cor.test(x, y, method = "spearman", exact = FALSE); c(rho = unname(ct$estimate), p = ct$p.value) }

## ---- column a: the counted resource ---------------------------------------------------
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds")); z <- IN$z; m <- IN$meta
ctl <- m$cond == "Control" & m$oxy == 21 & m$line %in% c("HC1", "HC2")
DR <- fread(file.path(D, "treatment_perturbation/library_scores_vs_counted_rate.tsv"))
PS <- DR[match(m$sample[ctl], sample)]; stopifnot(!anyNA(PS$prolif))
C1 <- rbindlist(c(list(data.table(gene = "proliferation score", x = m$days[ctl], y = as.numeric(scale(PS$prolif)))),
                  lapply(GENES, function(g) data.table(gene = g, x = m$days[ctl], y = z[g, ctl]))))
C1[, line := rep(m$line[ctl], length(GENES) + 1)]
say("a  %d untreated HC1/HC2 cultures, days %d to %d", sum(ctl), round(min(m$days[ctl])), round(max(m$days[ctl])))
for (g in c("proliferation score", GENES)) { r <- sp(C1[gene == g]$x, C1[gene == g]$y)
  say("   %-20s rho(days) = %+.2f  P = %.1e", g, r["rho"], r["p"]) }

## ---- column b: the donor cohort -------------------------------------------------------
Zc <- fread(file.path(D, "cohort_revised/primary/heatmap_matrix.tsv"))
SM <- fread(file.path(D, "cohort_revised/primary/heatmap_sample_order.tsv"))
G13 <- fread("supplementary_data/SupplementaryData13_gene_age_effects.tsv")
Zm <- as.matrix(Zc[, -1]); rownames(Zm) <- Zc$gene
stopifnot(all(GENES %in% rownames(Zm)), all(colnames(Zm) == SM$srr), nrow(SM) == 107)
dec <- function(a) { b <- pmin(floor(a / 10) * 10, 80); factor(ifelse(b >= 80, "80+", sprintf("%ds", b)),
                                                                levels = c("20s", "30s", "40s", "50s", "60s", "70s", "80+")) }
SM[, third := cut(prolif, quantile(prolif, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                  labels = c("slowest", "middle", "fastest"))]
C2 <- rbindlist(c(list(data.table(gene = "proliferation score", x = SM$age, y = as.numeric(scale(SM$prolif)))),
                  lapply(GENES, function(g) data.table(gene = g, x = SM$age, y = Zm[g, ]))))
C2[, band := dec(x)]; C2[, third := rep(SM$third, length(GENES) + 1)]
say("b  %d donors; per decade %s", nrow(SM), paste(sprintf("%s %d", names(table(dec(SM$age))), table(dec(SM$age))), collapse = ", "))
for (g in c("proliferation score", GENES)) { r <- sp(C2[gene == g]$x, C2[gene == g]$y)
  k <- G13[gene == g]
  say("   %-20s rho(age) = %+.2f  P = %.1e%s", g, r["rho"], r["p"],
      if (nrow(k)) sprintf("  | model: %+.3f log2CPM/decade FDR %.1e; adjusted FDR %.2f", k$log2FC_per_decade, k$FDR, k$FDR_adj) else "") }
say("   oldest donors (80+) in the slowest third: %d of %d", sum(SM$third == "slowest" & SM$age >= 80), sum(SM$age >= 80))

## ---- columns c and d: GTEx cultured fibroblasts and sun-exposed skin ------------------
## normalised as run_gtex_boundary.R does (filterByExpr, TMM, log-CPM), within each tissue
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
G <- "public_data_tierA/gtex"
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
hn <- names(fread(cmd = sprintf("zcat %s/gtex.gene_sums.SKIN.G026.gz | grep -v '^##' | head -2", G)))
BANDS <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
gt <- function(tk) {
  d <- as.data.table(TS[[tk]]$d)
  idcol <- if (sum(as.character(d$rail_id) %in% hn) > sum(d$external_id %in% hn)) "rail_id" else "external_id"
  ids <- as.character(d[[idcol]]); stopifnot(all(ids %in% hn))
  gx <- fread(cmd = sprintf("zcat %s/gtex.gene_sums.SKIN.G026.gz | grep -v '^##'", G), select = c(hn[1], ids))
  sym <- key[sub("\\..*$", "", gx[[1]])]
  X <- as.matrix(gx[, -1, with = FALSE])[, ids]; rm(gx)
  ok <- !is.na(sym) & sym != ""; X <- rowsum(X[ok, , drop = FALSE], sym[ok])
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  lg <- cpm(y, log = TRUE); rm(X, y)
  zg <- t(scale(t(lg[GENES, , drop = FALSE])))
  mid <- as.numeric(substr(as.character(d$age_bin), 1, 2)) + 4.5
  out <- rbindlist(c(list(data.table(gene = "proliferation score", band = as.character(d$age_bin), mid = mid,
                                     y = as.numeric(scale(d$prolif)))),
                     lapply(GENES, function(g) data.table(gene = g, band = as.character(d$age_bin), mid = mid, y = zg[g, ]))))
  out[, band := factor(sub("-", "–", band), levels = sub("-", "–", BANDS))]
  say("%s  %d donors", tk, nrow(d))
  for (g in c("proliferation score", GENES)) { r <- sp(out[gene == g]$mid, out[gene == g]$y)
    say("   %-20s rho(age band) = %+.2f  P = %.1e", g, r["rho"], r["p"]) }
  out }
C3 <- gt("culture"); C4 <- gt("legskin")

## ---- drawing ----------------------------------------------------------------------------
LV <- c("proliferation score", GENES)
for (x in list(C1, C2, C3, C4)) x[, gene := factor(gene, levels = LV)]
TCOL <- c(slowest = "#8C6D4F", middle = "#C9A27C", fastest = BLUE)
lab <- function(v) ifelse(v == "proliferation score", "atop('proliferation', 'score')", sprintf("italic('%s')", v))
strip_it <- function(p, strips = FALSE) {
  p <- p + facet_wrap(~ gene, ncol = 1, scales = "free_y", strip.position = "right",
                      labeller = as_labeller(lab, label_parsed)) +
    theme_sa(8) + theme(strip.placement = "outside", panel.spacing.y = unit(3, "pt"),
                        axis.text = element_text(size = 7), plot.margin = margin(3, 2, 3, 3),
                        plot.subtitle = element_text(size = 7.2, colour = INK2, lineheight = 1.0))
  if (strips) p + theme(strip.text.y = element_text(size = 7, angle = 0, hjust = 0)) else p + theme(strip.text.y = element_blank()) }
p1 <- ggplot(C1, aes(x, y)) +
  geom_point(aes(colour = line), size = 0.8, alpha = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = INK, linewidth = 0.5, span = 0.9) +
  scale_colour_manual(values = c(HC1 = BLUE, HC2 = RED), name = NULL) +
  scale_x_continuous("days in culture", breaks = c(50, 100, 150, 200)) +
  scale_y_continuous("expression (z within dataset)", labels = num_axis(0)) +
  labs(subtitle = sprintf("untreated HC1 and HC2\nn = %d cultures", sum(ctl)))
p1 <- strip_it(p1) + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(b = -4),
                          legend.key.size = unit(6, "pt"), legend.text = element_text(size = 7))
p2 <- ggplot(C2, aes(band, y)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, fill = "#EEF1F4", colour = INK2, linewidth = 0.3) +
  geom_jitter(aes(colour = third), width = 0.16, height = 0, size = 0.7, alpha = 0.8) +
  scale_colour_manual(values = TCOL, name = NULL) +
  scale_x_discrete("donor age") + scale_y_continuous(NULL, labels = num_axis(0)) +
  labs(subtitle = sprintf("donor cohort, by decade; n = %d\npoints: proliferation third", nrow(SM)))
p2 <- strip_it(p2) + guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 1.5))) +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(b = -4),
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = 7, margin = margin(r = 2)),
        legend.spacing.x = unit(2, "pt"), legend.key.spacing.y = unit(0, "pt"))
bx <- function(dt, sub, strips = FALSE) { p <- ggplot(dt, aes(band, y)) +
  geom_boxplot(outlier.size = 0.4, outlier.colour = GREY, width = 0.62, fill = "#EEF1F4", colour = INK2, linewidth = 0.3) +
  stat_summary(fun = median, geom = "line", aes(group = 1), colour = INK, linewidth = 0.45) +
  scale_x_discrete("age band (years)") + scale_y_continuous(NULL, labels = num_axis(0)) + labs(subtitle = sub)
  strip_it(p, strips) + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) }
p3 <- bx(C3, sprintf("GTEx fibroblast cultures\nn = %d donors", nrow(TS$culture$d)))
p4 <- bx(C4, sprintf("GTEx sun-exposed skin\nn = %d donors", nrow(TS$legskin$d)), strips = TRUE)
save_fig(lab_grid(p1, p2, p3, p4, labels = c("a", "b", "c", "d"), ncol = 4, rel_widths = c(1.15, 1.15, 0.95, 1.12),
                  align = "h", axis = "tb"),
         "Fig5.png", 183, 190)
