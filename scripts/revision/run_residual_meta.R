## One number for the residual age effect.
##
## The only positive claim that survives proliferation adjustment is a small decline
## in splicing-factor expression per decade, and the paper currently shows it as two
## cohort estimates side by side. This combines them by inverse variance and reports
## the heterogeneity, for both versions of the gene set. Post hoc.
##
## Only the proliferation-adjusted estimates are pooled: unadjusted effects depend on
## how proliferation happens to move with age in each collection, so pooling them
## would combine quantities that are not the same thing (Discussion, and the rule in
## the handover notes). The two cohorts are independent; the four definitions of the
## fibroblast cohort are not, and are reported as sensitivity rather than pooled.
## Output: public_data_tierA/derived/residual_meta/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table)})
D <- "public_data_tierA/derived"
O <- file.path(D, "residual_meta"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

SENS <- fread(file.path(D, "cohort_revised/sensitivity_metrics.tsv"))
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
d <- as.data.table(TS$culture$d)

gtex_effect <- function(setname) {
  d[, s := TS$culture$S[, setname]]
  tech <- c("rin", "isch", "loglib")[c(mean(is.finite(d$rin)) > 0.9,
                                       mean(is.finite(d$isch)) > 0.9, TRUE)]
  cov <- c("sex", if (nlevels(droplevels(d$hardy)) > 1) "hardy", tech)
  m <- lm(as.formula(paste("s ~ I(age/10) + prolif +", paste(cov, collapse = "+"))), d)
  s <- summary(m)$coef; k <- confint(m)
  data.table(cohort = "GTEx cultured fibroblasts", n = nobs(m), beta = s["I(age/10)", 1],
             se = s["I(age/10)", 2], lo = k["I(age/10)", 1], hi = k["I(age/10)", 2],
             p = s["I(age/10)", 4])
}

meta <- function(E) {
  w <- 1 / E$se^2
  b <- sum(w * E$beta) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (E$beta - b)^2); df <- nrow(E) - 1
  data.table(cohort = "combined (inverse variance)", n = sum(E$n), beta = b, se = se,
             lo = b - 1.96 * se, hi = b + 1.96 * se, p = 2 * pnorm(-abs(b / se)),
             Q = Q, Q_p = pchisq(Q, df, lower.tail = FALSE),
             I2 = max(0, 100 * (Q - df) / Q))
}

run <- function(which_metric, setname, label) {
  K <- SENS[metric == which_metric & cohort == "primary"]
  se_k <- (K$hi_adj - K$lo_adj) / (2 * 1.96)
  E <- rbind(
    data.table(cohort = "GSE113957 fibroblast donors", n = K$n, beta = K$beta_adj, se = se_k,
               lo = K$lo_adj, hi = K$hi_adj, p = K$p_adj),
    gtex_effect(setname))
  M <- meta(E)
  out <- rbind(E, M, fill = TRUE)
  out[, gene_set := label]
  say("\n-- %s --", label)
  print(out[, .(cohort, n, beta = round(beta, 4), lo = round(lo, 4), hi = round(hi, 4),
                p = signif(p, 2), I2 = round(I2))])
  out
}

## second pair: the Reactome programme the set was drawn from, which is the only other
## definition scored identically in both cohorts (the 177-gene set was not scored in GTEx)
DS <- fread(file.path(D, "machinery_specificity/donor_scores.tsv"))
m <- lm(premrna_reactome ~ I(age/10) + prolif + log_depth + repo + instr + sex, DS)
sm <- summary(m)$coef; km <- confint(m)
K2 <- data.table(cohort = "GSE113957 fibroblast donors", n = nobs(m), beta = sm["I(age/10)", 1],
                 se = sm["I(age/10)", 2], lo = km["I(age/10)", 1], hi = km["I(age/10)", 2],
                 p = sm["I(age/10)", 4])
E2 <- rbind(K2, gtex_effect("pre-mRNA processing"))
T2 <- rbind(E2, meta(E2), fill = TRUE)[, gene_set := "Reactome pre-mRNA processing"]
say("\n-- Reactome pre-mRNA processing, scored identically in both cohorts --")
print(T2[, .(cohort, n, beta = round(beta, 4), lo = round(lo, 4), hi = round(hi, 4),
             p = signif(p, 2), I2 = round(I2))])
TAB <- rbind(run("machinery96", "splicing 96", "preregistered 96-gene set"), T2)
fwrite(TAB, file.path(O, "residual_age_effect_meta.tsv"), sep = "\t")

## the four definitions of the fibroblast cohort, as sensitivity and not pooled
S4 <- SENS[metric %in% c("machinery", "machinery96"),
           .(metric, cohort, n, beta_adj, lo_adj, hi_adj, p_adj)][order(metric, cohort)]
print(S4, digits = 3)
fwrite(S4, file.path(O, "cohort_definition_sensitivity.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
