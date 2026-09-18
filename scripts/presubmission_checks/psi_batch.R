.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(limma); library(data.table)})
S <- "scripts/presubmission_checks"
X <- readRDS("public_data_tierA/derived/psi_age/psi_matrix.rds"); PSI <- X$PSI; ANN <- X$ANN; M <- as.data.frame(X$M)
md <- read.delim(file.path(S, "gse113957_meta2.tsv"))
idcol <- intersect(c("srr", "external_id"), names(M))[1]
if (is.na(idcol)) { cat("M columns:", paste(names(M), collapse = ", "), "\n"); quit() }
M <- cbind(M, md[match(M[[idcol]], md$srr), c("disease", "prefix", "instr")]); M$repo <- ifelse(M$prefix == "AG", "AG", "other")
stopifnot(ncol(PSI) == nrow(M))
lgt <- function(x) { x <- pmin(pmax(x, 0.005), 0.995); log(x/(1-x)) }
lg <- lgt(PSI)
run <- function(keep, form) { fit <- eBayes(lmFit(lg[, keep], model.matrix(form, data = M[keep, ]))); topTable(fit, coef = 2, number = Inf, sort.by = "none") }
col <- res <- fread("public_data_tierA/derived/psi_age/psi_age_events.tsv")
ev <- which(res$gene == "COL1A2" & res$FDR < 0.05)
show <- function(lab, keep, extra = "") {
  a0 <- run(keep, as.formula(paste("~ I(age/10) + log_depth", extra)))
  a1 <- run(keep, as.formula(paste("~ I(age/10) + prolif + log_depth", extra)))
  cat(sprintf("%-40s n=%3d | age-assoc %3d | after proliferation %3d | COL1A2 age events P: %s\n", lab, sum(keep),
      sum(a0$adj.P.Val < .05), sum(a1$adj.P.Val < .05 & a0$adj.P.Val < .05),
      paste(sprintf("%.1g", a0$P.Value[ev]), collapse = ", "))) }
N <- M$disease == "Normal"
show("published, all 142", rep(TRUE, nrow(M)))
show("normal 133", N)
show("normal 133 + repository + instrument", N, "+ repo + instr")
show("normal, age < 83", N & M$age < 83)
show("normal, adults 20-82", N & M$age >= 20 & M$age < 83)
cat("\nCOL1A2 event PSI (strongest), median by age band x repository (normal):\n")
top <- ev[which.min(res$p[ev])]; M$psi <- PSI[top, ]
M$band <- cut(M$age, c(0, 20, 40, 60, 83, 100), right = FALSE)
print(round(with(M[N, ], tapply(psi, list(band, repo), median)), 3))
