## Re-express the study's inferential statistics in the idiom conventional for
## ageing / molecular-biology journals:
##   - GSEA (fgsea) with NES and BH-FDR as the primary signature-association test
##   - Spearman rho reported with a Fisher-z 95% CI
##   - DerSimonian-Laird random-effects meta-analysis across studies, with Q and I^2
##   - Williams'/Steiger test for the difference between two dependent correlations
##   - partial correlation and linear models for the proliferation adjustment
##   - Wilcoxon rank-sum / Kruskal-Wallis for group comparisons
## The expression-matched permutation null is retained as a calibration check only.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(metafor); library(fgsea)})
D   <- "public_data_tierA/derived"
OUT <- file.path(D, "classical_stats"); dir.create(OUT, showWarnings = FALSE)
set.seed(20260912)
say <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------- 1. local data
gc_file <- function(s) file.path(D, "local_bam_v41", paste0(s, ".ReadsPerGene.out.tab"))
samples <- c(UVA0 = "UVA0", HDF = "HDF", Repro = "REP", iPSC = "IPS")
rd <- lapply(samples, function(s) {
  x <- read.delim(gc_file(s), header = FALSE, stringsAsFactors = FALSE)
  x <- x[!grepl("^N_", x$V1), ]
  setNames(x$V4, sub("\\..*$", "", x$V1))          # column 4 = reverse-stranded
})
g  <- Reduce(intersect, lapply(rd, names))
cnt <- do.call(cbind, lapply(rd, function(v) v[g]))
rownames(cnt) <- g

## map ENSG -> symbol using the locked anchor table
rank_tab <- read.delim(file.path(D, "local_repro_anchor/local_Repro_specific_gene_rank.tsv"))
dge <- DGEList(cnt); dge <- calcNormFactors(dge, method = "TMM")
cpm_all <- cpm(dge)
keep <- rowSums(cpm_all >= 1) >= 2
lcpm <- log2(cpm_all[keep, ] + 0.5)
say("local libraries: %s | genes after CPM>=1 in >=2: %d", paste(colnames(lcpm), collapse=", "), nrow(lcpm))
write.table(data.frame(gene = rownames(lcpm), lcpm), file.path(OUT, "local_logCPM_4samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## sample-sample Spearman + PCA on the 2000 most variable genes
sc <- cor(lcpm, method = "spearman")
write.table(round(sc, 4), file.path(OUT, "local_sample_correlation.tsv"), sep = "\t", quote = FALSE)
v  <- apply(lcpm, 1, var); top <- names(sort(v, decreasing = TRUE))[1:2000]
pc <- prcomp(t(lcpm[top, ]), scale. = FALSE)
pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
write.table(data.frame(sample = rownames(pc$x), pc$x[, 1:3], check.names = FALSE),
            file.path(OUT, "local_pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(paste0("PC", seq_along(pv), "\t", pv), file.path(OUT, "local_pca_varexp.tsv"))
say("PCA variance explained: PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%", pv[1], pv[2], pv[3])

## ------------------------------------------------- 2. gene-axis matrix and rho CIs
M <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
MA <- M                                   # full shared universe (all measured genes)
M  <- M[M$comparator_consistent %in% c(TRUE, "TRUE"), ]
say("shared universe across all six studies: %d genes", nrow(MA))
say("comparator-consistent genes shared with both axes: %d", nrow(M))

fisher_ci <- function(r, n, conf = 0.95) {
  z <- atanh(r); se <- 1/sqrt(n - 3); q <- qnorm(1 - (1 - conf)/2)
  c(lo = tanh(z - q*se), hi = tanh(z + q*se), z = z, se = se)
}
rho_row <- function(x, y, label) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  ci <- fisher_ci(unname(ct$estimate), n)
  data.frame(contrast = label, n = n, rho = unname(ct$estimate),
             ci_lo = ci["lo"], ci_hi = ci["hi"], z = ci["z"], se = ci["se"],
             S = unname(ct$statistic), p = ct$p.value, row.names = NULL)
}

DS <- list(GSE240226 = "UVA", GSE302943 = "UVA", GSE179848 = "senescence",
           GSE109700 = "senescence", GSE191055 = "senescence", GSE93535 = "senescence")
lab <- c(GSE240226 = "GSE240226  acute UVA", GSE302943 = "GSE302943  cumulative UVA",
         GSE179848 = "GSE179848  longitudinal RS", GSE109700 = "GSE109700  deep RS",
         GSE191055 = "GSE191055  late passage", GSE93535  = "GSE93535  SIPS")

## PRIMARY meta-analysis input: the locked per-study estimates, each computed on that
## study's own shared gene set with the full anchor vector (discordant genes scored 0).
PM <- read.delim(file.path(D, "decoupling_validation/primary_dataset_metrics.tsv"))
PM$acc <- sub("_.*$", "", PM$dataset_contrast)
per <- data.frame(contrast = PM$acc, axis = unlist(DS[PM$acc]),
                  n = PM$n_shared, rho = PM$rho_rank_stat, row.names = NULL)
ci <- t(mapply(fisher_ci, per$rho, per$n))
per <- cbind(per, ci_lo = ci[,"lo"], ci_hi = ci[,"hi"], z = ci[,"z"], se = ci[,"se"])
per$label <- lab[per$contrast]
write.table(per, file.path(OUT, "per_study_rho_fisherCI.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(per[, c("label","axis","n","rho","ci_lo","ci_hi")], digits = 3)

## SENSITIVITY: the same six correlations restricted to the 2,237-gene common set
sens <- do.call(rbind, lapply(names(DS), function(d) {
  r <- rho_row(M$Repro_specific_score, M[[d]], d); r$axis <- DS[[d]]; r }))
sens$label <- lab[sens$contrast]
write.table(sens, file.path(OUT, "per_study_rho_commonset_sensitivity.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
say("common-set sensitivity: UVA %+.3f/%+.3f  senescence %+.3f/%+.3f/%+.3f/%+.3f",
    sens$rho[1], sens$rho[2], sens$rho[3], sens$rho[4], sens$rho[5], sens$rho[6])

## ---------------------------------------- 3. random-effects meta-analysis per axis
meta_out <- do.call(rbind, lapply(c("UVA", "senescence"), function(a) {
  s <- per[per$axis == a, ]
  m <- rma(yi = s$z, sei = s$se, method = "DL")
  pooled <- tanh(as.numeric(m$b)); ci <- tanh(c(m$ci.lb, m$ci.ub))
  data.frame(axis = a, k = m$k, pooled_rho = pooled, ci_lo = ci[1], ci_hi = ci[2],
             z = as.numeric(m$b), se = as.numeric(m$se), p = m$pval,
             Q = m$QE, Q_df = m$k - 1, Q_p = m$QEp, I2 = m$I2, tau2 = m$tau2)
}))
write.table(meta_out, file.path(OUT, "meta_analysis_axes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
for (i in seq_len(nrow(meta_out))) with(meta_out[i,],
  say("%s axis: pooled rho = %+.3f [%.3f, %.3f], k=%d, Q=%.1f (df=%d, P=%.3g), I2=%.1f%%",
      axis, pooled_rho, ci_lo, ci_hi, k, Q, Q_df, Q_p, I2))

## difference between the two pooled estimates (subgroup contrast in Fisher-z space)
du <- meta_out[meta_out$axis == "UVA", ]; ds <- meta_out[meta_out$axis == "senescence", ]
dz <- du$z - ds$z; dse <- sqrt(du$se^2 + ds$se^2)
dci <- dz + c(-1.96, 1.96) * dse
say("pooled z difference = %.3f [%.3f, %.3f], Z = %.2f, P = %.3g",
    dz, dci[1], dci[2], dz/dse, 2*pnorm(-abs(dz/dse)))
write.table(data.frame(delta_z = dz, se = dse, ci_lo = dci[1], ci_hi = dci[2],
            Z = dz/dse, p = 2*pnorm(-abs(dz/dse))),
            file.path(OUT, "meta_axis_difference.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ------------------- 4. Williams' test for two dependent correlations (meta axes)
r_au <- cor(M$Repro_specific_score, M$UVA_meta,        method = "spearman")
r_as <- cor(M$Repro_specific_score, M$senescence_meta, method = "spearman")
r_us <- cor(M$UVA_meta,             M$senescence_meta, method = "spearman")
n    <- nrow(M)
williams <- function(r12, r13, r23, n) {
  Rdet <- 1 - r12^2 - r13^2 - r23^2 + 2*r12*r13*r23
  rbar <- (r12 + r13)/2
  t <- (r12 - r13) * sqrt(((n - 1)*(1 + r23)) /
       (2*((n - 1)/(n - 3))*Rdet + rbar^2*(1 - r23)^3))
  c(t = t, df = n - 3, p = 2*pt(-abs(t), n - 3))
}
w <- williams(r_au, r_as, r_us, n)
say("rho(anchor,UVA)=%+.3f  rho(anchor,sen)=%+.3f  rho(UVA,sen)=%+.3f | D=%+.3f",
    r_au, r_as, r_us, r_au - r_as)
say("Williams t(%d) = %.1f, P = %.3g   [gene-level n; anti-conservative, see Methods]",
    w["df"], w["t"], w["p"])
write.table(data.frame(rho_anchor_UVA=r_au, rho_anchor_sen=r_as, rho_UVA_sen=r_us,
            D=r_au-r_as, n_genes=n, williams_t=w["t"], df=w["df"], p=w["p"]),
            file.path(OUT, "williams_test_D.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ------------------------------------- 5. fgsea: anchor gene sets vs each reference
up <- MA$gene[MA$Repro_specific_score >=  log2(1.25)]
dn <- MA$gene[MA$Repro_specific_score <= -log2(1.25)]
say("anchor gene sets in the shared universe: up %d, down %d", length(up), length(dn))
sets <- list(`Repro-CM up` = up, `Repro-CM down` = dn)
tie_frac <- sapply(names(DS), function(d) {
  v <- MA[[d]]; 1 - length(unique(v))/length(v) })
gsea_ok <- names(DS)[tie_frac < 0.20]
say("GSEA references: %s | excluded for tied ranks: %s",
    paste(gsea_ok, collapse=", "),
    paste(sprintf("%s (%.0f%% tied)", names(DS)[tie_frac >= 0.20],
                  100*tie_frac[tie_frac >= 0.20]), collapse=", "))
write.table(data.frame(contrast = names(DS), tie_fraction = round(tie_frac, 4),
            used_for_gsea = tie_frac < 0.20),
            file.path(OUT, "gsea_rank_tie_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
gs <- do.call(rbind, lapply(gsea_ok, function(d) {
  st <- setNames(MA[[d]], MA$gene); st <- st[is.finite(st)]
  set.seed(1); r <- fgsea(sets, st, minSize = 15, nproc = 1)
  data.frame(contrast = d, axis = DS[[d]], pathway = r$pathway, NES = r$NES,
             pval = r$pval, size = r$size)
}))
gs$FDR <- p.adjust(gs$pval, "BH")
gs$label <- lab[gs$contrast]
write.table(gs, file.path(OUT, "fgsea_anchor_vs_references.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(gs[, c("contrast","pathway","NES","pval","FDR")], digits = 3)

## running enrichment curves against the two meta axes, for plotting
for (ax in c("UVA_meta", "senescence_meta")) {
  st <- setNames(MA[[ax]], MA$gene); st <- sort(st[is.finite(st)], decreasing = TRUE)
  for (nm in names(sets)) {
    es <- plotEnrichmentData(sets[[nm]], st)
    write.table(es$curve, file.path(OUT, sprintf("gsea_curve_%s_%s.tsv", ax,
                gsub("[^A-Za-z]", "", nm))), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(es$ticks, file.path(OUT, sprintf("gsea_ticks_%s_%s.tsv", ax,
                gsub("[^A-Za-z]", "", nm))), sep="\t", quote=FALSE, row.names=FALSE)
  }
}

## ---------------------------------- 6. proliferation: partial correlation + models
pl <- read.delim(file.path(D, "direction_probe/P3_proliferation_loading_GSE179848.tsv"))
M$prolif <- pl$proliferation_loading[match(M$gene, pl$gene)]
sub <- M[is.finite(M$prolif), ]
pcor <- function(x, y, z) {
  rx <- resid(lm(rank(x) ~ rank(z))); ry <- resid(lm(rank(y) ~ rank(z)))
  ct <- cor.test(rx, ry); ci <- fisher_ci(unname(ct$estimate), length(rx))
  c(r = unname(ct$estimate), lo = ci["lo"], hi = ci["hi"], p = ct$p.value)
}
pu <- pcor(sub$Repro_specific_score, sub$UVA_meta,        sub$prolif)
ps <- pcor(sub$Repro_specific_score, sub$senescence_meta, sub$prolif)
say("partial rho | proliferation: UVA %+.3f [%.3f,%.3f]  senescence %+.3f [%.3f,%.3f]  D=%+.3f",
    pu["r"], pu["lo.lo"], pu["hi.hi"], ps["r"], ps["lo.lo"], ps["hi.hi"], pu["r"]-ps["r"])
write.table(data.frame(axis=c("UVA","senescence"),
            partial_rho=c(pu["r"],ps["r"]), ci_lo=c(pu["lo.lo"],ps["lo.lo"]),
            ci_hi=c(pu["hi.hi"],ps["hi.hi"]), p=c(pu["p.p"],ps["p.p"]), n=nrow(sub)),
            file.path(OUT, "partial_correlation_proliferation.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
say("anchor vs proliferation loading: rho = %+.3f",
    cor(sub$Repro_specific_score, sub$prolif, method="spearman"))

## ------------------------------- 7. compendium: Fisher-z CIs and group comparisons
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
CP$rho_UVA_lo <- NA; CP$rho_UVA_hi <- NA; CP$rho_sen_lo <- NA; CP$rho_sen_hi <- NA
for (i in seq_len(nrow(CP))) {
  a <- fisher_ci(CP$rho_UVA[i],        CP$n_genes[i])
  b <- fisher_ci(CP$rho_senescence[i], CP$n_genes[i])
  CP$rho_UVA_lo[i] <- a["lo"]; CP$rho_UVA_hi[i] <- a["hi"]
  CP$rho_sen_lo[i] <- b["lo"]; CP$rho_sen_hi[i] <- b["hi"]
  CP$D_se[i] <- sqrt(2)/sqrt(CP$n_genes[i] - 3)
}
write.table(CP, file.path(OUT, "compendium_with_CI.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cls <- function(id, con) {
  if (grepl("senesc|SIPS|_deep|_early|P27", paste(id, con), ignore.case=TRUE)) "senescence"
  else if (grepl("age|old_vs_young|donor age", con, ignore.case=TRUE)) "donor age"
  else if (grepl("rescue|red_vs|prered|Osmoter|sauchinone|UIFSP|Maifuyin|succinate", paste(id,con), ignore.case=TRUE)) "photoprotection"
  else if (grepl("UV", paste(id, con))) "UV injury"
  else if (grepl("D13|D7|dox|OSK|MPTR|reprogram", paste(id, con), ignore.case=TRUE)) "reprogramming"
  else "other perturbation" }
CP$group <- mapply(cls, CP$id, CP$contrast)
kw <- kruskal.test(D ~ factor(group), data = CP)
say("Kruskal-Wallis across %d groups: chi2 = %.1f (df = %d), P = %.3g",
    length(unique(CP$group)), kw$statistic, kw$parameter, kw$p.value)
gsum <- aggregate(D ~ group, CP, function(x) c(n=length(x), med=median(x), q1=quantile(x,.25), q3=quantile(x,.75)))
gsum <- data.frame(group=gsum$group, as.data.frame(gsum$D))
write.table(gsum, file.path(OUT, "compendium_group_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(gsum, digits = 3)
writeLines(capture.output(kw), file.path(OUT, "compendium_kruskal.txt"))

sen <- CP$D[CP$group == "senescence"]; oth <- CP$D[CP$group != "senescence"]
wt <- wilcox.test(sen, oth); say("senescence vs rest (Wilcoxon): W = %.0f, P = %.3g", wt$statistic, wt$p.value)

## ------------------------------------------ 8. held-out tests with Fisher-z CIs
KEEP <- c("id","contrast","n_shared","rho","p_matched","BH_FDR","verdict")
rd_ho <- function(f, fam) {
  if (!file.exists(f)) return(NULL)
  x <- read.delim(f); for (k in KEEP) if (!k %in% names(x)) x[[k]] <- NA
  cbind(x[, KEEP], family = fam) }
ho <- do.call(rbind, list(
  rd_ho(file.path(D, "heldout_validation/heldout_primary_tests.tsv"), 1),
  rd_ho(file.path(D, "heldout2_senescence/HO5_primary_tests.tsv"),    2),
  rd_ho(file.path(D, "heldout3_photoprotection/HO6_primary_tests.tsv"), 3),
  rd_ho(file.path(D, "heldout3_photoprotection/HO67_axis_loadings.tsv"), 3),
  rd_ho(file.path(D, "heldout4_sauchinone/HO8_primary_test.tsv"),    4)))
ho <- ho[!is.na(ho$rho) & !is.na(ho$n_shared), ]
ci <- t(mapply(fisher_ci, ho$rho, ho$n_shared))
ho$ci_lo <- ci[,"lo"]; ho$ci_hi <- ci[,"hi"]
write.table(ho, file.path(OUT, "heldout_primary_fisherCI.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("held-out primary tests with CIs: %d", nrow(ho))

## ---------------------------------------------------------------- 9. session
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done -> %s", OUT)
