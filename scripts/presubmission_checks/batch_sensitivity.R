.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
S <- "scripts/presubmission_checks"
M  <- read.delim("public_data_tierA/derived/outcome_vs_expression/GSE113957_sample_metrics.tsv")
md <- read.delim(file.path(S, "gse113957_meta2.tsv")); M <- cbind(M, md[match(M$srr, md$srr), c("disease","prefix","instr","sex")])
M$repo <- ifelse(M$prefix == "AG", "AG", "other")
cat("HGPS by prefix:", paste(M$prefix[M$disease=="HGPS"], collapse=","), "\n\n")
fit <- function(d, metric, extra = "") {
  f0 <- as.formula(paste(metric, "~ I(age/10) + log_depth", extra)); f1 <- as.formula(paste(metric, "~ I(age/10) + prolif + log_depth", extra))
  c0 <- summary(lm(f0, d))$coef[2, ]; c1 <- summary(lm(f1, d))$coef[2, ]
  sprintf("P %.2g -> adj P %.2g (%3.0f%% of effect lost)", c0[4], c1[4], 100*(1 - c1[1]/c0[1])) }
N <- M[M$disease == "Normal", ]
sets <- list(
  "published: all 142"                                    = list(M, ""),
  "normal 133"                                            = list(N, ""),
  "normal 133 + repository + instrument covariates"       = list(N, "+ repo + instr"),
  "normal, age < 83"                                      = list(N[N$age < 83, ], ""),
  "normal, AG lines only"                                 = list(N[N$repo == "AG", ], ""),
  "normal, AG lines, age < 83"                            = list(N[N$repo == "AG" & N$age < 83, ], ""),
  "normal, adults 20-82"                                  = list(N[N$age >= 20 & N$age < 83, ], ""))
for (nm in names(sets)) { d <- sets[[nm]][[1]]; ex <- sets[[nm]][[2]]
  cat(sprintf("== %-50s n=%3d  rho(age,unannot)=%+.2f  rho(age,prolif)=%+.2f\n", nm, nrow(d),
      cor(d$age, d$unannot_reads, method="spearman"), cor(d$age, d$prolif, method="spearman")))
  for (m in c("machinery", "unannot_reads")) cat("     ", sprintf("%-14s", m), fit(d, m, ex), "\n") }
cat("\nmedian unannotated fraction x 1e4, normal donors, by age band and repository:\n")
N$band <- cut(N$age, c(0, 20, 40, 60, 83, 100), right = FALSE)
print(round(with(N, tapply(unannot_reads, list(band, repo), median)) * 1e4, 2))
print(with(N, table(band, repo)))
