## Two small tables that give a source to numbers quoted in the text and had none.
##
##   A. Which technical covariate separates the residual age effect from zero?
##      The proliferation-adjusted age model of the 177-gene score is refitted in the
##      107-donor primary cohort of GSE113957 under five covariate sets: sequencing
##      depth only; depth + cell repository; depth + sequencing instrument; depth + sex;
##      and the full cohort model (depth + repository + instrument + sex). Each row is
##      compared with the unadjusted model that carries the same covariates.
##      -> reviewer_sensitivities/covariate_attribution.tsv
##
##   B. The correlation of the whole 177-gene set and of its two subsets (127 core
##      spliceosome genes; 50 genes of capping, 3'-end processing and export) with the
##      measured division rate in the 328 counted libraries of GSE179848. These were
##      printed by run_reviewer_sensitivities.R but never written to a table. All three
##      are scored in the expression matrix of the programme audit, as there.
##      -> reviewer_sensitivities/subset_division_rate.tsv
##
## The scores are not recomputed: A reads the per-donor table written by
## run_cohort_core.R, and the cohort is checked against scripts/revision/cohort_gse113957.R.
## Run from the project root:  Rscript scripts/revision/run_covariate_attribution.R
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "reviewer_sensitivities")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
LOG <- file(file.path(O, "covariate_attribution.log"), open = "wt")
say <- function(...) { x <- sprintf(...); cat(x, "\n"); cat(x, "\n", file = LOG) }

## ================= A. covariate attribution ===================================
P <- fread(file.path(D, "cohort_revised/primary/sample_metrics.tsv"))
DON <- gse113957_donors()
stopifnot(nrow(P) == 107,
          setequal(P$srr, DON$srr[gse113957_keep(DON, "primary")]),
          all(P$disease == "Normal"), all(P$age >= 20),
          all(c("machinery", "prolif", "log_depth", "repo", "instr", "sex") %in% names(P)))
say("primary cohort: %d donors, ages %.0f-%.0f | repository %s | instrument %s | sex %s",
    nrow(P), min(P$age), max(P$age),
    paste(names(table(P$repo)), table(P$repo), collapse = ", "),
    paste(names(table(P$instr)), table(P$instr), collapse = ", "),
    paste(names(table(P$sex)), table(P$sex), collapse = ", "))
xt <- table(P$instr, P$repo)
say("instrument by repository: %s",
    paste(sprintf("%s: %s", rownames(xt),
                  apply(xt, 1, function(r) paste(colnames(xt), r, collapse = ", "))), collapse = "; "))

SETS <- list("depth"                               = "log_depth",
             "depth + repository"                  = "log_depth + repo",
             "depth + instrument"                  = "log_depth + instr",
             "depth + sex"                         = "log_depth + sex",
             "depth + repository + instrument + sex (cohort model)" =
                                                     "log_depth + repo + instr + sex")
age_row <- function(fit) {
  s <- summary(fit)$coef["I(age/10)", ]
  ci <- s[1] + c(-1, 1) * qt(.975, fit$df.residual) * s[2]
  c(beta = s[[1]], lo = ci[1], hi = ci[2], p = s[[4]])
}
A <- rbindlist(lapply(names(SETS), function(nm) {
  f0 <- lm(as.formula(paste("machinery ~ I(age/10) +", SETS[[nm]])), P)
  f1 <- lm(as.formula(paste("machinery ~ I(age/10) +", SETS[[nm]], "+ prolif")), P)
  a0 <- age_row(f0); a1 <- age_row(f1)
  data.table(covariates = nm, n = nrow(P),
             beta_unadj = a0[["beta"]], p_unadj = a0[["p"]],
             beta_adj = a1[["beta"]], lo = a1[["lo"]], hi = a1[["hi"]], p = a1[["p"]],
             pct_removed = 100 * (1 - a1[["beta"]] / a0[["beta"]]))
}))
print(A, digits = 3)
fwrite(A, file.path(O, "covariate_attribution.tsv"), sep = "\t")
for (i in seq_len(nrow(A)))
  say("   %-55s unadjusted %+.4f (P = %.2g) | adjusted %+.4f [%+.4f, %+.4f] P = %.3g | %.1f%% removed",
      A$covariates[i], A$beta_unadj[i], A$p_unadj[i], A$beta_adj[i], A$lo[i], A$hi[i], A$p[i],
      A$pct_removed[i])

## the covariate terms themselves in the full proliferation-adjusted model
full <- summary(lm(machinery ~ I(age/10) + log_depth + repo + instr + sex + prolif, P))$coef
for (term in setdiff(rownames(full), "(Intercept)"))
  say("   full adjusted model, %-22s %+.4f (P = %.2g)", term, full[term, 1], full[term, 4])

## agreement with the cohort model already reported (sensitivity_metrics.tsv)
ref <- fread(file.path(D, "cohort_revised/sensitivity_metrics.tsv"))[metric == "machinery" & cohort == "primary"]
stopifnot(abs(ref$beta_adj - A$beta_adj[5]) < 1e-10, abs(ref$pct_lost - A$pct_removed[5]) < 1e-8)
say("   full-model row reproduces cohort_revised/sensitivity_metrics.tsv (%.4f, %.1f%% removed)",
    ref$beta_adj, ref$pct_lost)

## ================= B. subsets against the measured division rate ==============
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
spliceo <- readLines(file.path(O, "core_spliceosome_genes.txt"))
stopifnot(length(CORE) == 177, all(spliceo %in% CORE))
rest <- setdiff(CORE, spliceo)
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z; rate <- IN$meta$rate
stopifnot(ncol(z) == length(rate))
B <- rbindlist(lapply(list(list("whole set (177 genes)", CORE),
                           list("core spliceosome (127 genes)", spliceo),
                           list("capping, 3'-end, export (50 genes)", rest)), function(x) {
  g <- intersect(rownames(z), x[[2]]); s <- colMeans(z[g, , drop = FALSE])
  ok <- is.finite(s) & is.finite(rate)
  ct <- cor.test(s[ok], rate[ok], method = "spearman", exact = FALSE)
  ci <- tanh(atanh(ct$estimate) + c(-1, 1) * qnorm(.975) / sqrt(sum(ok) - 3))
  data.table(gene_set = x[[1]], genes_in_set = length(x[[2]]), genes_measured = length(g),
             n_libraries = sum(ok), rho_rate = unname(ct$estimate), lo = ci[1], hi = ci[2],
             p = ct$p.value)
}))
print(B, digits = 3)
fwrite(B, file.path(O, "subset_division_rate.tsv"), sep = "\t")
for (i in seq_len(nrow(B)))
  say("   %-36s division-rate rho = %+.4f [%.3f, %.3f] (P = %.1g), %d of %d genes measured, n = %d",
      B$gene_set[i], B$rho_rate[i], B$lo[i], B$hi[i], B$p[i], B$genes_measured[i],
      B$genes_in_set[i], B$n_libraries[i])
close(LOG)
## the directory's sessionInfo.txt belongs to run_reviewer_sensitivities.R
writeLines(capture.output(sessionInfo()), file.path(O, "covariate_attribution_sessionInfo.txt"))
