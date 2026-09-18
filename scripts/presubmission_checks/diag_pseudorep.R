.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(nlme))
T <- read.delim("public_data_tierA/derived/division_rate/sample_division_rate.tsv")
ct <- cor.test(T$rate, T$splice, method = "spearman", exact = FALSE)
cat(sprintf("pooled (as published): rho = %.3f, P = %.1g, n = %d libraries from %d lines\n", ct$estimate, ct$p.value, nrow(T), length(unique(T$line))))
cat("\nwithin-line Spearman:\n")
for (l in sort(unique(T$line))) { s <- T[T$line == l, ]; c1 <- cor.test(s$rate, s$splice, method = "spearman", exact = FALSE)
  cat(sprintf("  %-8s n=%3d  rho=%+.2f  P=%.1g   (days %.0f-%.0f)\n", l, nrow(s), c1$estimate, c1$p.value, min(s$days), max(s$days))) }
T$z_rate <- as.numeric(scale(T$rate)); T$z_spl <- as.numeric(scale(T$splice))
m1 <- lme(z_spl ~ z_rate, random = ~ 1 | line, data = T, method = "REML")
m2 <- tryCatch(lme(z_spl ~ z_rate, random = ~ z_rate | line, data = T, method = "REML", control = lmeControl(opt = "optim")), error = function(e) NULL)
m3 <- lme(z_spl ~ z_rate + days, random = ~ 1 | line, data = T, method = "REML")
pr <- function(m, lab) { if (is.null(m)) return(cat(lab, ": did not converge\n")); tt <- summary(m)$tTable["z_rate", ]
  cat(sprintf("  %-46s beta(z) = %.3f  SE %.3f  df %d  P = %.1g\n", lab, tt[1], tt[2], tt[3], tt[5])) }
cat("\nmixed models (standardised):\n")
pr(m1, "random intercept per line"); pr(m2, "random intercept + slope per line"); pr(m3, "random intercept, + days in culture")
## line x condition means: one value per culture series, a conservative unit
A <- aggregate(cbind(rate, splice) ~ line + cond + oxy, data = T, FUN = mean)
ca <- cor.test(A$rate, A$splice, method = "spearman", exact = FALSE)
cat(sprintf("\nline x condition x oxygen means: rho = %.2f, P = %.1g, n = %d series\n", ca$estimate, ca$p.value, nrow(A)))
## healthy controls only (untreated, 21% O2) - the cleanest culture-lifespan subset
H <- T[T$cond == "Control" & T$oxy == 21 & grepl("^HC", T$line), ]
ch <- cor.test(H$rate, H$splice, method = "spearman", exact = FALSE)
cat(sprintf("untreated healthy lines at 21%% O2: rho = %.2f, P = %.1g, n = %d libraries, %d lines\n", ch$estimate, ch$p.value, nrow(H), length(unique(H$line))))
mh <- lme(scale(splice) ~ scale(rate), random = ~ 1 | line, data = H)
cat(sprintf("   same, mixed model beta = %.2f, P = %.1g\n", fixef(mh)[2], summary(mh)$tTable[2, 5]))
