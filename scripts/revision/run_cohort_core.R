## Recompute the GSE113957 donor-cohort results under the cohort definition fixed
## in scripts/revision/cohort_gse113957.R (primary: normal donors aged 20+).
##
## Junction-derived quantities (unannotated-junction use, entropy) are per-sample
## statistics and are taken from the existing metrics file unchanged; everything
## that depends on normalisation across samples (machinery score, proliferation
## score, per-gene models) is recomputed within each cohort.
##
## Outputs: public_data_tierA/derived/cohort_revised/<cohort>/...
##          public_data_tierA/derived/cohort_revised/sensitivity_*.tsv
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
source("scripts/revision/cohort_gse113957.R")
D  <- "public_data_tierA/derived"
OUT <- file.path(D, "cohort_revised"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
CORE  <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))     # 177 genes
CORE96 <- readLines(file.path(D, "conserved_core/age_down_splicing_core.txt"))       # preregistered 96
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")

## ---- donors, and the per-sample junction metrics computed earlier ------------
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))
MET <- MET[, .(srr, depth, unannot_reads, unannot_junc, entropy, log_depth)]
DON <- merge(DON, MET, by = "srr")                       # 142 with junction metrics
say("donors with junction metrics: %d", nrow(DON))

## ---- gene counts ------------------------------------------------------------
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
say("gene matrix: %d genes x %d samples", nrow(GM), ncol(GM))

## ---- scores and models within one cohort ------------------------------------
score_cohort <- function(keep) {
  d <- DON[keep, ]
  X <- GM[, d$srr, drop = FALSE]
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
  z <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
  d$machinery  <- colMeans(z[intersect(rownames(z), CORE), , drop = FALSE])
  d$machinery96 <- colMeans(z[intersect(rownames(z), CORE96), , drop = FALSE])
  d$prolif     <- colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE])
  list(d = d, lc = lc)
}

## age effect per decade before and after adjustment, with technical covariates
age_models <- function(d, metric, covar = TRUE) {
  base <- "I(age/10) + log_depth"
  if (covar) {
    if (length(unique(d$repo))  > 1) base <- paste(base, "+ repo")
    if (length(unique(d$instr)) > 1) base <- paste(base, "+ instr")
    if (length(unique(d$sex))   > 1) base <- paste(base, "+ sex")
  }
  f0 <- as.formula(paste(metric, "~", base))
  f1 <- as.formula(paste(metric, "~", base, "+ prolif"))
  c0 <- summary(lm(f0, d))$coef["I(age/10)", ]; c1 <- summary(lm(f1, d))$coef["I(age/10)", ]
  ci <- function(fit) { s <- summary(fit)$coef["I(age/10)", ]; s[1] + c(-1, 1) * qt(.975, fit$df.residual) * s[2] }
  data.frame(metric = metric, n = nrow(d),
             beta = c0[1], lo = ci(lm(f0, d))[1], hi = ci(lm(f0, d))[2], p = c0[4],
             beta_adj = c1[1], lo_adj = ci(lm(f1, d))[1], hi_adj = ci(lm(f1, d))[2], p_adj = c1[4],
             pct_lost = 100 * (1 - c1[1] / c0[1]), row.names = NULL)
}

## per-gene age models within the 177-gene set
gene_models <- function(d, lc, covar = TRUE) {
  g <- intersect(rownames(lc), CORE)
  base <- "v ~ I(age/10) + log_depth"
  if (covar) {
    if (length(unique(d$repo))  > 1) base <- paste(base, "+ repo")
    if (length(unique(d$instr)) > 1) base <- paste(base, "+ instr")
    if (length(unique(d$sex))   > 1) base <- paste(base, "+ sex")
  }
  f0 <- as.formula(base); f1 <- as.formula(paste(base, "+ prolif"))
  res <- do.call(rbind, lapply(g, function(gn) {
    d$v <- lc[gn, ]
    a <- summary(lm(f0, d))$coef["I(age/10)", ]; b <- summary(lm(f1, d))$coef["I(age/10)", ]
    data.frame(gene = gn, r_prolif = cor(d$v, d$prolif),
               beta = a[1], p = a[4], beta_adj = b[1], p_adj = b[4],
               retained = b[1] / a[1], row.names = NULL) }))
  res$FDR <- p.adjust(res$p, "BH"); res$FDR_adj <- p.adjust(res$p_adj, "BH")
  res$retained[!is.finite(res$retained)] <- NA
  res
}

COHORTS <- c(primary = "normal donors aged 20+ (primary)",
             normal  = "all normal donors (children included)",
             adult2082 = "normal adults aged 20-82",
             published = "the published 142-donor cohort")
SENS <- list(); GSENS <- list(); SC <- list()
for (k in names(COHORTS)) {
  keep <- gse113957_keep(DON, k)
  s <- score_cohort(keep); d <- s$d
  dir.create(file.path(OUT, k), showWarnings = FALSE)
  ## covariates only where the cohort definition does not already control them
  cov <- k != "published"
  M <- do.call(rbind, lapply(c("machinery", "machinery96", "unannot_reads", "unannot_junc", "entropy"),
                             function(m) age_models(d, m, covar = cov)))
  M$cohort <- k; M$covariates <- cov
  G <- gene_models(d, s$lc, covar = cov)
  kept <- G$gene[G$FDR < .05 & G$FDR_adj < .05 & G$beta < 0 & G$beta_adj < 0]
  gsum <- data.frame(cohort = k, n = nrow(d),
                     rho_age_prolif = cor(d$age, d$prolif, method = "spearman"),
                     r_machinery_prolif = cor(d$machinery, d$prolif),
                     age_assoc = sum(G$FDR < .05), adj_assoc = sum(G$FDR_adj < .05),
                     decline_kept = length(kept), adj_only = sum(G$FDR_adj < .05 & G$FDR >= .05),
                     rho_coupling_retained = cor(G$r_prolif, G$retained, method = "spearman", use = "pair"))
  fwrite(d, file.path(OUT, k, "sample_metrics.tsv"), sep = "\t")
  fwrite(G, file.path(OUT, k, "gene_level.tsv"), sep = "\t")
  ## heatmap input: the 177 genes as gene-wise z-scores, donors ordered young to old
  gset <- intersect(rownames(s$lc), CORE); ordm <- order(d$age)
  Zh <- t(scale(t(s$lc[gset, d$srr[ordm], drop = FALSE])))
  write.table(data.frame(gene = rownames(Zh), Zh, check.names = FALSE),
              file.path(OUT, k, "heatmap_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  fwrite(d[ordm, c("srr", "age", "sex", "repo", "instr", "prolif", "machinery", "unannot_reads")],
         file.path(OUT, k, "heatmap_sample_order.tsv"), sep = "\t")
  writeLines(kept, file.path(OUT, k, "genes_keeping_decline.txt"))
  SENS[[k]] <- M; GSENS[[k]] <- gsum; SC[[k]] <- list(d = d, G = G, kept = kept)
  say("== %-10s n=%3d | rho(age,prolif)=%+.2f | machinery: P %.2g -> %.2g (%.0f%% lost) | genes: %d age, %d kept",
      k, nrow(d), gsum$rho_age_prolif, M$p[M$metric == "machinery"], M$p_adj[M$metric == "machinery"],
      M$pct_lost[M$metric == "machinery"], gsum$age_assoc, gsum$decline_kept)
}
S <- rbindlist(SENS); fwrite(S, file.path(OUT, "sensitivity_metrics.tsv"), sep = "\t")
GS <- rbindlist(GSENS); fwrite(GS, file.path(OUT, "sensitivity_genes.tsv"), sep = "\t")

say("\n=== age effect per decade, machinery (177 genes) ===")
print(as.data.frame(S[S$metric == "machinery", .(cohort, n, beta, p, beta_adj, p_adj, pct_lost)]), digits = 3, row.names = FALSE)
say("\n=== age effect per decade, use of unannotated junctions (outcome) ===")
print(as.data.frame(S[S$metric == "unannot_reads", .(cohort, n, beta, p, beta_adj, p_adj)]), digits = 3, row.names = FALSE)
say("\n=== gene level ===")
print(as.data.frame(GS), digits = 3, row.names = FALSE)

core_kept <- Reduce(intersect, lapply(SC, `[[`, "kept"))
say("\ngenes keeping a decline in all four cohort definitions: %d%s", length(core_kept),
    if (length(core_kept)) paste0(" -> ", paste(core_kept, collapse = ", ")) else "")
writeLines(core_kept, file.path(OUT, "genes_keeping_decline_all_cohorts.txt"))
saveRDS(SC, file.path(OUT, "cohort_results.rds"))
