## Statistical supplements requested by the pre-submission review (steps 3-4).
##   A. Fig. 2: the 328 libraries are repeated measures of seven cell lines.
##      Per-line correlations with exact P where n is small, mixed models with
##      the line as a random effect, and confidence intervals for the headline r.
##   B. Fig. 1D: exact P for the four culture-time series (n = 6 to 25).
##   C. Fig. 6: the contrast-level coupling with (i) the donor-age contrast
##      recomputed in the redefined cohort, (ii) expression-matched random gene
##      sets as a null, (iii) the senescence class removed, (iv) one value per
##      study, and (v) class residual tests reported with raw as well as BH P.
##   D. Sensitivity: would the 177-gene set change if its GSE113957 input were
##      computed in the redefined cohort?
## Outputs: public_data_tierA/derived/revision_stats/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(nlme); library(limma); library(fgsea)})
D <- "public_data_tierA/derived"
O <- file.path(D, "revision_stats"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
fisher_ci <- function(r, n) { z <- atanh(r); s <- 1/sqrt(n - 3); tanh(z + c(-1.96, 1.96) * s) }
sp <- function(x, y, exact = NULL) {
  n <- sum(is.finite(x) & is.finite(y))
  if (is.null(exact)) exact <- n < 10
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = exact))
  ci <- fisher_ci(unname(ct$estimate), n)
  list(rho = unname(ct$estimate), p = ct$p.value, n = n, lo = ci[1], hi = ci[2], exact = exact)
}

## ============================ A. division rate ==============================
T <- read.delim(file.path(D, "division_rate/sample_division_rate.tsv"))
a <- sp(T$rate, T$splice, exact = FALSE)
say("A. pooled across %d libraries from %d lines: rho = %.3f [%.2f, %.2f], P = %.1e",
    a$n, length(unique(T$line)), a$rho, a$lo, a$hi, a$p)

PL <- do.call(rbind, lapply(sort(unique(T$line)), function(l) {
  s <- T[T$line == l, ]; r <- sp(s$rate, s$splice)
  data.frame(line = l, n = r$n, rho = r$rho, lo = r$lo, hi = r$hi, p = r$p,
             exact = r$exact, days_from = min(s$days), days_to = max(s$days)) }))
print(PL, digits = 3, row.names = FALSE)
write.table(PL, file.path(O, "fig2_per_line_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

T$z_rate <- as.numeric(scale(T$rate)); T$z_spl <- as.numeric(scale(T$splice))
mm <- function(m, lab) { if (is.null(m)) return(data.frame(model = lab, beta = NA, se = NA, df = NA, p = NA))
  tt <- summary(m)$tTable["z_rate", ]
  data.frame(model = lab, beta = tt[1], se = tt[2], df = tt[3], p = tt[5], row.names = NULL) }
MM <- rbind(
  mm(lme(z_spl ~ z_rate, random = ~ 1 | line, data = T, method = "REML"), "random intercept per line"),
  mm(tryCatch(lme(z_spl ~ z_rate, random = ~ z_rate | line, data = T, method = "REML",
                  control = lmeControl(opt = "optim")), error = function(e) NULL), "random intercept and slope"),
  mm(lme(z_spl ~ z_rate + days, random = ~ 1 | line, data = T, method = "REML"), "random intercept, plus days in culture"))
print(MM, digits = 3, row.names = FALSE)
write.table(MM, file.path(O, "fig2_mixed_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

A <- aggregate(cbind(rate, splice) ~ line + cond + oxy, data = T, FUN = mean)
ag <- sp(A$rate, A$splice, exact = FALSE)
H <- T[T$cond == "Control" & T$oxy == 21 & grepl("^HC", T$line), ]
hh <- sp(H$rate, H$splice, exact = FALSE)
say("   line x condition x oxygen means: rho = %.2f [%.2f, %.2f], P = %.1e, n = %d series", ag$rho, ag$lo, ag$hi, ag$p, ag$n)
say("   untreated healthy lines at 21%% O2: rho = %.2f [%.2f, %.2f], P = %.1e, n = %d", hh$rho, hh$lo, hh$hi, hh$p, hh$n)

## ============================ B. culture time ===============================
L <- read.delim(file.path(D, "core_depth/GSE179848_longitudinal_core.tsv"))
CT <- do.call(rbind, lapply(sort(unique(L$line)), function(l) {
  s <- L[L$line == l, ]; r <- sp(s$days, s$score)
  data.frame(line = l, n = r$n, days = round(max(s$days)), rho = r$rho, p = r$p, exact = r$exact) }))
say("\nB. culture-time series (exact P where n < 10):")
print(CT, digits = 3, row.names = FALSE)
write.table(CT, file.path(O, "fig1d_culture_time.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ====================== C. contrast-level coupling ==========================
source("scripts/revision/cohort_gse113957.R")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
CC <- unique(unlist(P[grep("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(P))]))

## C1. donor-age contrast recomputed in the redefined cohort ------------------
DON <- gse113957_donors()
fp  <- read.delim(gzfile("public_data_tierA/aging/GSE113957_fpkm.txt.gz"), check.names = FALSE)
sym <- toupper(sub("\\|.*$", "", fp[["Annotation/Divergence"]]))
scols <- names(fp)[grepl("^[0-9]+_", names(fp))]
X <- as.matrix(fp[, scols]); storage.mode(X) <- "numeric"
X <- rowsum(X, sym); X <- X[rownames(X) != "" & rownames(X) != "-", ]
id_fp <- sub("^([0-9]+)_.*$", "\\1", scols); id_dn <- sub("^([0-9]+)_.*$", "\\1", DON$title)
dn <- DON[match(id_fp, id_dn), ]
keep <- !is.na(dn$srr) & gse113957_keep(dn, "primary")
say("\nC1. donor-age contrast: %d of %d FPKM columns are primary-cohort donors", sum(keep), length(scols))
d <- dn[keep, ]
des <- model.matrix(~ I(age/10) + sex + instr + repo, data = d)
fit <- eBayes(lmFit(log2(X[, keep, drop = FALSE] + 0.1), des), trend = TRUE, robust = TRUE)
tab <- topTable(fit, coef = "I(age/10)", number = Inf, sort.by = "none")
newDE <- data.frame(dataset = "GSE113957",
                    contrast = "per_10_years_normal_adults_20plus_adjusted_for_sex_platform_and_repository",
                    gene = rownames(tab), logFC = tab$logFC, AveExpr = tab$AveExpr,
                    statistic = tab$t, PValue = tab$P.Value, FDR = tab$adj.P.Val)
write.table(newDE[order(newDE$PValue), ], file.path(O, "GSE113957_primary_cohort_per_decade_DE.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
oldDE <- read.delim(file.path(D, "figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv"))
cmn <- intersect(oldDE$gene, newDE$gene)
say("    agreement with the published contrast: r = %.3f over %d genes; genes down in both: %d",
    cor(oldDE$logFC[match(cmn, oldDE$gene)], newDE$logFC[match(cmn, newDE$gene)]), length(cmn),
    sum(oldDE$logFC[match(cmn, oldDE$gene)] < 0 & newDE$logFC[match(cmn, newDE$gene)] < 0))

## C2. the 64 contrasts, with the new donor-age vector ------------------------
Z <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); META <- Z$meta
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
vec <- lapply(split(N, N$id), function(x) setNames(x$logFC, x$gene))
for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
  x <- read.delim(f); vec[[sub("\\.tsv$", "", basename(f))]] <- setNames(x$logFC, toupper(x$gene)) }
lgf <- list(GSE109700_deep = "figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
            GSE179848_late = "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
            GSE93535_SIPS  = "figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
            GSE191055_P27  = "figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
            GSE109700_early = "figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv")
for (id in names(lgf)) { f <- file.path(D, lgf[[id]]); if (file.exists(f)) {
  x <- read.delim(f); vec[[id]] <- setNames(x$logFC, toupper(x$gene)) } }
vec[["GSE113957_age"]] <- setNames(newDE$logFC, toupper(newDE$gene))     # redefined cohort
## 2026-09-17: the unreplicated contrast of this study is no longer part of the
## compendium. The study is now a reanalysis of public data only; the earlier work it
## follows is cited rather than reanalysed. Set INCLUDE_LOCAL to TRUE to restore it.
INCLUDE_LOCAL <- FALSE
if (INCLUDE_LOCAL) {
  loc <- AX[AX$comparator_consistent %in% c(TRUE, "TRUE") & AX$Repro_specific_score != 0, ]
  vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))
}

dlt <- function(v, set) { g <- intersect(names(v), set); if (length(g) < 15) return(NA)
  bg <- v[is.finite(v)]; mean(bg[names(bg) %in% set]) - mean(bg) }
R <- do.call(rbind, lapply(names(vec), function(i)
  data.frame(id = i, spl = dlt(vec[[i]], CORE), cc = dlt(vec[[i]], CC))))
R$class <- META$class[match(R$id, META$id)]
FIX <- c(GSE116968_uv_vs_control_1h = "UV injury", GSE116968_uv_vs_control_4h = "UV injury",
         GSE179848_late = "senescence",
         ## culture-time control in the reprogramming series, not a reprogramming condition
         GSE149694_FibroblastD7 = "metabolic / culture")
hit <- R$id %in% names(FIX); R$class[hit] <- FIX[R$id[hit]]
R <- R[is.finite(R$spl) & is.finite(R$cc), ]
obs <- sp(R$cc, R$spl, exact = FALSE); r2 <- summary(lm(spl ~ cc, R))$r.squared
say("\nC2. %d contrasts: rho = %.3f [%.2f, %.2f], P = %.1e, R2 = %.3f", nrow(R), obs$rho, obs$lo, obs$hi, obs$p, r2)
noshare <- setdiff(CORE, CC)
Rn <- do.call(rbind, lapply(names(vec), function(i) data.frame(id = i, spl = dlt(vec[[i]], noshare))))
R$spl_noshare <- Rn$spl[match(R$id, Rn$id)]
say("    with the %d genes shared with the cell-cycle set removed: rho = %.3f, R2 = %.3f",
    length(intersect(CORE, CC)), cor(R$cc, R$spl_noshare, method = "spearman"),
    summary(lm(spl_noshare ~ cc, R))$r.squared)
write.table(R, file.path(O, "fig6_contrasts_revised.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## class residuals, raw and BH-corrected
R$resid <- resid(lm(spl ~ cc, R))
cls <- names(which(table(R$class) >= 3))
CL <- do.call(rbind, lapply(cls, function(k) { v <- R$resid[R$class == k]; tt <- t.test(v)
  data.frame(class = k, n = length(v), mean = mean(v), lo = tt$conf.int[1], hi = tt$conf.int[2], p_raw = tt$p.value) }))
CL$p_BH <- p.adjust(CL$p_raw, "BH")
print(CL, digits = 3, row.names = FALSE)
write.table(CL, file.path(O, "fig6_class_residuals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## C3. leverage and null sets ------------------------------------------------
nos <- R$class != "senescence"
say("\nC3. without the %d senescence contrasts: rho = %.2f, R2 = %.2f, n = %d",
    sum(!nos), cor(R$cc[nos], R$spl[nos], method = "spearman"), summary(lm(spl ~ cc, R[nos, ]))$r.squared, sum(nos))
R$study <- sub("_.*$", "", R$id)
St <- aggregate(cbind(spl, cc) ~ study, data = R, FUN = mean)
say("    one averaged value per study: rho = %.2f, R2 = %.2f, n = %d",
    cor(St$cc, St$spl, method = "spearman"), summary(lm(spl ~ cc, St))$r.squared, nrow(St))

set.seed(20260915)
lev <- tapply(rowMeans(log2(X + 1)), rownames(X), max)
meas <- table(unlist(lapply(vec[R$id], function(v) names(v)[is.finite(v)])))
univ <- intersect(names(meas)[meas >= 0.8 * nrow(R)], names(lev))
lv <- lev[univ]; bins <- cut(lv, quantile(lv, seq(0, 1, .05)), include.lowest = TRUE, labels = FALSE)
names(bins) <- univ
core_b <- table(bins[intersect(CORE, univ)]); pool <- setdiff(univ, c(CORE, CC))
draw <- function() unlist(lapply(names(core_b), function(b) sample(pool[bins[pool] == as.integer(b)], core_b[[b]])))
NL <- t(replicate(1000, { g <- draw(); s <- sapply(R$id, function(i) dlt(vec[[i]], g)); ok <- is.finite(s)
  c(rho = cor(R$cc[ok], s[ok], method = "spearman"), r2 = summary(lm(s[ok] ~ R$cc[ok]))$r.squared) }))
say("    expression-matched random sets (%d genes, 1000 draws): rho median %.2f [%.2f, %.2f], max %.2f; empirical P = %.3f",
    length(intersect(CORE, univ)), median(NL[, 1]), quantile(NL[, 1], .025), quantile(NL[, 1], .975),
    max(NL[, 1]), mean(NL[, 1] >= obs$rho))
write.table(data.frame(NL), file.path(O, "fig6_null_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

CTX <- do.call(rbind, lapply(c("Translation", "Metabolism of RNA", "Collagen formation",
                               "Extracellular matrix organization", "Unfolded Protein Response (UPR)",
                               "Cellular Senescence"), function(nm) {
  g <- setdiff(P[[nm]], CC); s <- sapply(R$id, function(i) dlt(vec[[i]], g)); ok <- is.finite(s)
  data.frame(programme = nm, n_genes = length(g), rho = cor(R$cc[ok], s[ok], method = "spearman"),
             r2 = summary(lm(s[ok] ~ R$cc[ok]))$r.squared) }))
print(CTX, digits = 3, row.names = FALSE)
write.table(CTX, file.path(O, "fig6_other_programmes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ================== D. would the gene set itself change? ====================
## The 177-gene set was defined from three series, one of which is GSE113957 in
## its published (adults 22-89, HGPS already excluded) form. Redefine it with the
## donor-age contrast recomputed in the primary cohort and see what changes.
say("\nD. set-definition sensitivity")
SER <- list(GSE179848_late = file.path(D, "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv"),
            GSE307377_age  = file.path(D, "compendium/signatures/GSE307377_age.tsv"))
vecs <- c(list(GSE113957_age = setNames(newDE$logFC, toupper(newDE$gene))),
          lapply(SER, function(f) { x <- read.delim(f); setNames(x$logFC, toupper(x$gene)) }))
gg  <- Reduce(union, lapply(vecs, names))
mat <- sapply(vecs, function(v) v[gg]); rownames(mat) <- gg
mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
DOWN <- rownames(mat)[rowSums(mat < 0) == ncol(mat)]
pre <- unique(unlist(P[grep("mRNA Splicing|Processing of Capped Intron|3'-end processing", names(P))]))
v3 <- intersect(DOWN, pre)
say("    genes measured in all three series: %d | falling in all three: %d | redefined set: %d genes",
    nrow(mat), length(DOWN), length(v3))
say("    overlap with the 177-gene set: %d (%.0f%% of the 177, %.0f%% of the redefined set)",
    length(intersect(v3, CORE)), 100 * length(intersect(v3, CORE)) / length(CORE),
    100 * length(intersect(v3, CORE)) / max(1, length(v3)))
writeLines(v3, file.path(O, "splicing_set_primary_cohort_sensitivity.txt"))
s3 <- sapply(R$id, function(i) dlt(vec[[i]], v3)); ok3 <- is.finite(s3)
say("    contrast-level coupling with the redefined set: rho = %.2f, R2 = %.2f (published set %.2f, %.2f)",
    cor(R$cc[ok3], s3[ok3], method = "spearman"), summary(lm(s3[ok3] ~ R$cc[ok3]))$r.squared, obs$rho, r2)
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()

## ================================================================================
## Leave-out check for the set-defining series (2026-09-18, referee request).
## Fourteen of the 63 contrasts come from the three series that defined the splicing
## set.  This asks whether the compendium relationship depends on them.
## NB: join contrasts to series through the id, NOT by row position against
## supp_table1.tsv -- the two files are sorted differently and a positional join
## silently mislabels every contrast (see CLAUDE.md section 20).
## ================================================================================
{
  CO <- read.delim(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
  SG <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
  CO$series <- sub("^(GSE[0-9]+).*$", "\\1", CO$id)
  CO$series[!grepl("^GSE", CO$series)] <- NA
  CO$series[is.na(CO$series)] <- SG$dataset[match(CO$id[is.na(CO$series)], SG$id)]
  stopifnot(!anyNA(CO$series))
  DEFINING <- c("GSE179848", "GSE113957", "GSE307377")
  fit1 <- function(x) {
    ct <- suppressWarnings(cor.test(x$cc, x$spl, method = "spearman", exact = FALSE))
    z <- atanh(ct$estimate); se <- 1 / sqrt(nrow(x) - 3)
    data.frame(n_contrasts = nrow(x), n_series = length(unique(x$series)),
               rho = unname(ct$estimate), lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
               p = ct$p.value, r2 = summary(lm(spl ~ cc, x))$r.squared)
  }
  LO <- rbind(data.frame(set = "all contrasts", fit1(CO)),
              data.frame(set = "set-defining series dropped",
                         fit1(CO[!(CO$series %in% DEFINING), ])))
  write.table(LO, file.path(D, "revision_stats/leave_out_defining_series.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("\nleave-out check (set-defining series):\n"); print(LO, digits = 3)
}
