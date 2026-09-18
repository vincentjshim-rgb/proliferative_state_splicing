.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
S <- "scripts/presubmission_checks"
O <- "public_data_tierA/derived/outcome_vs_expression"
Z  <- read.delim(file.path(O, "core_gene_heatmap_matrix.tsv"), check.names = FALSE)
SM <- read.delim(file.path(O, "heatmap_sample_order.tsv")); M <- read.delim(file.path(O, "GSE113957_sample_metrics.tsv"))
md <- read.delim(file.path(S, "gse113957_meta2.tsv"))
SM <- cbind(SM, md[match(SM$srr, md$srr), c("disease","prefix","instr")]); SM$log_depth <- M$log_depth[match(SM$srr, M$srr)]
SM$repo <- ifelse(SM$prefix == "AG", "AG", "other")
X <- as.matrix(Z[, -1]); rownames(X) <- Z$gene
gl <- function(keep, extra = "") {
  d <- SM[keep, ]; Xk <- X[, keep, drop = FALSE]
  f0 <- as.formula(paste("v ~ I(age/10) + log_depth", extra)); f1 <- as.formula(paste("v ~ I(age/10) + prolif + log_depth", extra))
  r <- t(apply(Xk, 1, function(v) { d$v <- v; a <- summary(lm(f0, d))$coef[2, ]; b <- summary(lm(f1, d))$coef[2, ]; c(a[1], a[4], b[1], b[4]) }))
  fdr <- p.adjust(r[, 2], "BH"); fdra <- p.adjust(r[, 4], "BH")
  kept <- rownames(r)[fdr < .05 & fdra < .05 & r[, 1] < 0 & r[, 3] < 0]
  list(txt = sprintf("n=%3d | age-assoc %3d | adj-assoc %3d | decline kept %3d | adj-only %d", nrow(d), sum(fdr < .05), sum(fdra < .05), length(kept), sum(fdra < .05 & fdr >= .05)), kept = kept) }
base <- gl(rep(TRUE, nrow(SM)))
res <- list(
  "published, all 142"                     = base,
  "normal 133"                             = gl(SM$disease == "Normal"),
  "normal 133 + repository + instrument"   = gl(SM$disease == "Normal", "+ repo + instr"),
  "normal, age < 83"                       = gl(SM$disease == "Normal" & SM$age < 83),
  "normal, adults 20-82"                   = gl(SM$disease == "Normal" & SM$age >= 20 & SM$age < 83))
for (nm in names(res)) cat(sprintf("%-40s %s\n", nm, res[[nm]]$txt))
core <- Reduce(intersect, lapply(res, `[[`, "kept"))
cat("\ngenes keeping a decline in every version:", length(core), "->", paste(core, collapse = ", "), "\n")
