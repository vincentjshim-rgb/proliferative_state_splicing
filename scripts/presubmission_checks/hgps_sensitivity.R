.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
O <- "public_data_tierA/derived/outcome_vs_expression"
S <- "scripts/presubmission_checks"
M  <- read.delim(file.path(O, "GSE113957_sample_metrics.tsv"))
st <- read.delim(file.path(S, "gse113957_status.tsv"))
M$disease <- st$disease[match(M$srr, st$srr)]
cat("HGPS in cohort:", sum(M$disease == "HGPS"), "| ages:", paste(sort(M$age[M$disease == "HGPS"]), collapse = ","), "\n")
cat("HGPS vs normal (same age < 10): prolif", round(mean(M$prolif[M$disease=="HGPS"]),3), "vs", round(mean(M$prolif[M$disease=="Normal" & M$age < 10]),3),
    "| unannot_reads", signif(mean(M$unannot_reads[M$disease=="HGPS"]),3), "vs", signif(mean(M$unannot_reads[M$disease=="Normal" & M$age < 10]),3), "\n\n")
fit <- function(d, metric) {
  f0 <- lm(d[[metric]] ~ I(d$age/10) + d$log_depth); f1 <- lm(d[[metric]] ~ I(d$age/10) + d$prolif + d$log_depth)
  c0 <- summary(f0)$coef[2, ]; c1 <- summary(f1)$coef[2, ]
  sprintf("beta %+.4g P %.2g -> adj %+.4g P %.2g (%.0f%% lost)", c0[1], c0[4], c1[1], c1[4], 100*(1 - c1[1]/c0[1])) }
sets <- list(`all 142 (as published)` = M, `133 normal (HGPS removed)` = M[M$disease == "Normal", ],
             `adults >= 20, normal` = M[M$disease == "Normal" & M$age >= 20, ])
for (nm in names(sets)) { d <- sets[[nm]]
  cat(sprintf("== %s  n=%d | rho(age,prolif)=%.2f | rho(machinery,prolif) r=%.2f | rho(age,unannot)=%.2f P=%.1g\n", nm, nrow(d),
      cor(d$age, d$prolif, method="spearman"), cor(d$machinery, d$prolif),
      cor(d$age, d$unannot_reads, method="spearman"), cor.test(d$age, d$unannot_reads, method="spearman", exact=FALSE)$p.value))
  for (m in c("machinery", "unannot_reads", "entropy")) cat("   ", sprintf("%-14s", m), fit(d, m), "\n") }
## gene level: z-scores are a per-gene linear transform of logCPM, so model P values are unchanged
Z  <- read.delim(file.path(O, "core_gene_heatmap_matrix.tsv"), check.names = FALSE)
SM <- read.delim(file.path(O, "heatmap_sample_order.tsv"))
X <- as.matrix(Z[, -1]); rownames(X) <- Z$gene
stopifnot(ncol(X) == nrow(SM)); SM$disease <- st$disease[match(SM$srr, st$srr)]
SM$log_depth <- M$log_depth[match(SM$srr, M$srr)]
gl <- function(keep) {
  d <- SM[keep, ]; Xk <- X[, keep, drop = FALSE]
  r <- t(apply(Xk, 1, function(v) { f0 <- lm(v ~ I(d$age/10) + d$log_depth); f1 <- lm(v ~ I(d$age/10) + d$prolif + d$log_depth)
        c(b = coef(f0)[2], p = summary(f0)$coef[2, 4], ba = coef(f1)[2], pa = summary(f1)$coef[2, 4]) }))
  fdr <- p.adjust(r[, 2], "BH"); fdra <- p.adjust(r[, 4], "BH")
  sprintf("n=%d | age-assoc %d | adj-assoc %d | decline kept (raw & adj, both negative) %d | adj-only %d",
          nrow(d), sum(fdr < .05), sum(fdra < .05), sum(fdr < .05 & fdra < .05 & r[,1] < 0 & r[,3] < 0), sum(fdra < .05 & fdr >= .05)) }
cat("\n== gene level (177 genes)\n")
cat("   all 142         :", gl(rep(TRUE, nrow(SM))), "\n")
cat("   133 normal      :", gl(SM$disease == "Normal"), "\n")
cat("   adults>=20 norm :", gl(SM$disease == "Normal" & SM$age >= 20), "\n")
