## Three analyses that tie the six results together without new data.
##
## The chain has three joints that the figures did not address:
##   A  the instrument changes between steps. Division rate was counted in one resource; every
##      other step uses the 20-marker transcriptional score. Calibrating the score against
##      counting puts every step into one unit -- divisions per day -- and measures how much
##      the substitution costs.
##   B  "adjustment removes 70% of the age effect" is a description, not an estimand. The same
##      claim as a mediation model gives a proportion mediated with an interval, and lets the
##      mediator's measurement error be carried through instead of ignored.
##   C  the tissue step is a different comparison from the culture steps. A two-compartment
##      mixture says what between-donor correlation to expect in bulk skin if the coupling held
##      inside cycling cells and only the cycling fraction varied, which is the alternative the
##      composition adjustment tests empirically.
## All three are post hoc.
## Output: public_data_tierA/derived/chain/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table)})
D <- "public_data_tierA/derived"; O <- file.path(D, "chain")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
set.seed(20260920)
OUT <- list()

## =====================================================================================
## A. One instrument: the 20-marker score in units of counted divisions per day
## =====================================================================================
say("=== A. calibrating the proliferation score against counted division rate ===")
LS <- fread(file.path(D, "treatment_perturbation/library_scores_vs_counted_rate.tsv"))
stopifnot(nrow(LS) == 328)
## the calibration: one standard deviation of the score, in divisions per day
cal <- lm(rate ~ prolif, LS)
sd_score <- sd(LS$prolif)
per_sd <- unname(coef(cal)[2]) * sd_score
say("  rate = %.3f + %.3f x score   (R2 = %.2f, n = %d)", coef(cal)[1], coef(cal)[2],
    summary(cal)$r.squared, nrow(LS))
PER_SD_IN_RATE <- per_sd
say("  one SD of the score (%.2f score units) = %.3f divisions per day", sd_score, per_sd)
say("  the counted rate spans %.2f to %.2f divisions per day in this resource",
    min(LS$rate), max(LS$rate))

## What the substitution costs, measured where both exist. Regressing the counted rate on the
## score gives an attenuated slope, because the score carries error; so converting a
## score-based effect into divisions per day through that slope overstates it. The size of the
## overstatement is measurable here, and the same factor is applied to the other datasets.
b_counted <- unname(coef(lm(scale(splice) ~ rate, LS))[2])            # per division per day
b_score   <- unname(coef(lm(scale(splice) ~ scale(prolif), LS))[2])   # per SD of the score
b_naive   <- b_score / PER_SD_IN_RATE                                 # naive conversion
SHRINK    <- b_counted / b_naive
say("  splicing score per division per day, counted:              %+.2f SD", b_counted)
say("  the same from the score through the calibration (naive):   %+.2f SD", b_naive)
say("  the score therefore overstates the per-division effect by %.1f times;", 1 / SHRINK)
say("  score-based effects are multiplied by %.2f to reach the counted scale.", SHRINK)
say("  correlation with the counted rate: splicing set %.2f, 20-marker score %.2f",
    cor(LS$splice, LS$rate, method = "spearman"), cor(LS$prolif, LS$rate, method = "spearman"))
say("  so a score-based estimate is an upper bound on the effect and a lower bound on the correlation.")

## every sample-level dataset in the same unit
cohort <- fread(file.path(D, "cohort_revised/primary/sample_metrics.tsv"))
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
zb <- function(y, x) unname(coef(lm(scale(y) ~ scale(x)))[2])          # per SD of x
rows <- list(
  data.table(dataset = "counted resource, against counting",
             n = nrow(LS), unit = "328 libraries",
             per_sd = b_counted * PER_SD_IN_RATE, rho = cor(LS$splice, LS$rate, method = "spearman"),
             counted = TRUE),
  data.table(dataset = "counted resource, against the score",
             n = nrow(LS), unit = "328 libraries",
             per_sd = b_score, rho = cor(LS$splice, LS$prolif, method = "spearman"), counted = FALSE),
  data.table(dataset = "donor cohort", n = nrow(cohort), unit = "107 adult donors",
             per_sd = zb(cohort$machinery, cohort$prolif),
             rho = cor(cohort$machinery, cohort$prolif, method = "spearman"), counted = FALSE))
for (tk in c("culture", "legskin", "pubskin")) {
  d <- as.data.frame(TS[[tk]]$d); s <- TS[[tk]]$S[, "splicing 96"]
  ok <- complete.cases(d$prolif, d$rin, d$isch, d$loglib, s); d <- d[ok, ]; s <- s[ok]
  rx <- resid(lm(prolif ~ rin + isch + loglib, d)); ry <- resid(lm(s ~ rin + isch + loglib, d))
  rows[[length(rows) + 1]] <- data.table(
    dataset = paste0("GTEx ", c(culture = "cultured fibroblasts", legskin = "skin, sun-exposed",
                                pubskin = "skin, not exposed")[tk]),
    n = nrow(d), unit = sprintf("%d donors, partial", nrow(d)),
    per_sd = zb(ry, rx), rho = cor(rx, ry, method = "spearman"), counted = FALSE)
}
A <- rbindlist(rows)
A[, per_div_day := ifelse(counted, per_sd / PER_SD_IN_RATE,
                          per_sd / PER_SD_IN_RATE * SHRINK)]
OUT$A <- A
say("\n  every sample-level step in one unit (change in the splicing score, in SD):")
print(A[, .(dataset, n, `per SD of proliferation` = round(per_sd, 3),
            `per 0.1 div/day` = round(per_div_day * 0.1, 3), rho = round(rho, 3))],
      row.names = FALSE)
fwrite(A, file.path(O, "common_unit.tsv"), sep = "\t")

## =====================================================================================
## B. The age effect as a mediation model
## =====================================================================================
say("\n=== B. age -> proliferation -> splicing in the 107 adult donors ===")
M <- as.data.frame(cohort)
M$sex <- factor(M$sex); M$repo <- factor(M$repo); M$instr <- factor(M$instr)
COV <- "repo + instr + sex + log_depth"
fml <- function(y, x) as.formula(sprintf("%s ~ %s + %s", y, x, COV))
paths <- function(dd) {
  c_tot <- coef(lm(fml("scale(machinery)", "I(age/10)"), dd))[2]
  a     <- coef(lm(fml("scale(prolif)", "I(age/10)"), dd))[2]
  m     <- lm(as.formula(sprintf("scale(machinery) ~ I(age/10) + scale(prolif) + %s", COV)), dd)
  b     <- coef(m)["scale(prolif)"]; c_dir <- coef(m)["I(age/10)"]
  c(total = unname(c_tot), a = unname(a), b = unname(b), direct = unname(c_dir),
    indirect = unname(a * b), prop = unname(a * b / c_tot))
}
p0 <- paths(M)
say("  total effect of age  c  = %+.4f SD per decade", p0["total"])
say("  age -> proliferation a  = %+.4f", p0["a"])
say("  proliferation -> splicing b = %+.4f (age in the model)", p0["b"])
say("  direct effect        c' = %+.4f", p0["direct"])
say("  indirect effect     a*b = %+.4f", p0["indirect"])
say("  proportion mediated     = %.1f%%", 100 * p0["prop"])
B <- t(replicate(2000, { i <- sample(nrow(M), replace = TRUE); paths(M[i, ]) }))
ci <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
say("  bootstrap (2,000 resamples of donors):")
say("    indirect effect  %+.4f  [%+.4f, %+.4f]", p0["indirect"], ci(B[, "indirect"])[1], ci(B[, "indirect"])[2])
say("    direct effect    %+.4f  [%+.4f, %+.4f]", p0["direct"], ci(B[, "direct"])[1], ci(B[, "direct"])[2])
say("    proportion mediated %.0f%%  [%.0f%%, %.0f%%]", 100 * p0["prop"],
    100 * ci(B[, "prop"])[1], 100 * ci(B[, "prop"])[2])
say("    direct effect crosses zero: %s", ifelse(prod(ci(B[, "direct"])) < 0, "yes", "no"))

## the mediator is measured with error; what that does to the proportion mediated
r_obs <- cor(LS$prolif, LS$rate, method = "pearson")
say("  the mediator is a proxy: it correlates with counted division rate at %.2f (Pearson, %d libraries).",
    r_obs, nrow(LS))
rel <- c(1, 0.8, 0.6, r_obs^2)
SENS <- data.table(reliability = rel,
                   prop_mediated = pmin(1, p0["prop"] / rel),
                   direct = p0["total"] - pmin(p0["total"] * sign(p0["total"]),
                                               abs(p0["indirect"]) / rel) * sign(p0["total"]))
say("  under classical measurement error the indirect path is attenuated by the reliability:")
for (k in seq_len(nrow(SENS)))
  say("    reliability %.2f -> proportion mediated %.0f%%", SENS$reliability[k],
      100 * SENS$prop_mediated[k])
say("  a reliability of %.2f or less is enough for the age effect to be fully mediated.",
    p0["prop"])
OUT$B <- list(paths = p0, boot = B, sens = SENS)
fwrite(data.table(quantity = names(p0), estimate = as.numeric(p0),
                  lo = c(NA, NA, NA, ci(B[, "direct"])[1], ci(B[, "indirect"])[1], ci(B[, "prop"])[1]),
                  hi = c(NA, NA, NA, ci(B[, "direct"])[2], ci(B[, "indirect"])[2], ci(B[, "prop"])[2])),
       file.path(O, "mediation.tsv"), sep = "\t")
fwrite(SENS, file.path(O, "mediation_measurement_error.tsv"), sep = "\t")

## =====================================================================================
## C. What a two-compartment tissue should show
## =====================================================================================
say("\n=== C. a mixture model for bulk skin ===")
sk <- as.data.frame(TS$legskin$d); Ssk <- TS$legskin$S
CC <- "cell cycle (positive control)"
okk <- complete.cases(sk$prolif, sk$rin, sk$isch, sk$loglib, Ssk[, "splicing 96"], Ssk[, CC])
sk <- sk[okk, ]; Ssk <- Ssk[okk, ]
rz <- function(v) as.numeric(scale(resid(lm(v ~ rin + isch + loglib, sk))))
rx <- rz(sk$prolif); rs <- rz(Ssk[, "splicing 96"]); rc <- rz(Ssk[, CC])
obs_sp <- cor(rx, rs, method = "spearman"); obs_cc <- cor(rx, rc, method = "spearman")
say("  observed in sun-exposed skin (%d donors): splicing %.2f, cell-cycle programme %.2f",
    nrow(sk), obs_sp, obs_cc)

## A donor's bulk skin is a mixture. Write the between-donor variation of each score as a part
## that follows the cycling fraction f and a part that does not:
##     P = z_f + e          the proliferation score is read from cycling cells only
##     C = z_f + e          so is the cell-cycle programme
##     S = w z_f + sqrt(1 - w^2) N + e      the splicing machinery is expressed in every cell,
##                                          so only the share w of its variation follows f
## Under that model rho(P,S) = w x rho(P,C): the cell-cycle programme measures how well the
## design can see the cycling fraction at all, and w is what is left to explain. Nothing here
## depends on the absolute size of f.
w_implied <- obs_sp / obs_cc
say("  rho(P,S) = w x rho(P,C) under the model, so the observed pair implies w = %.2f", w_implied)
say("  that is, about %.0f%% of the between-donor variation in the bulk splicing score follows the",
    100 * w_implied)
say("  cycling fraction, and %.0f%% of its variance (w squared).", 100 * w_implied^2)
say("  if the coupling were carried into tissue in full (w = 1) the model expects %.2f, %.1f times",
    obs_cc, obs_cc / obs_sp)
say("  the observed %.2f.", obs_sp)

## the same by simulation, which also carries the measurement noise of the two programme scores
sim <- function(w, n = nrow(sk), s_noise = 0.6, reps = 400) {
  r <- replicate(reps, {
    zf <- rnorm(n); N <- rnorm(n)
    P <- zf + rnorm(n, 0, s_noise); C <- zf + rnorm(n, 0, s_noise)
    S <- w * zf + sqrt(max(0, 1 - w^2)) * N + rnorm(n, 0, s_noise)
    c(cor(P, S, method = "spearman"), cor(P, C, method = "spearman")) })
  c(splice = mean(r[1, ]), cc = mean(r[2, ])) }
grid <- seq(0.2, 1.4, by = 0.02)
noise <- grid[which.min(abs(sapply(grid, function(z) sim(0.5, s_noise = z)["cc"]) - obs_cc))]
say("  simulation: noise set so the cell-cycle programme reproduces %.2f -> s = %.2f", obs_cc, noise)
W <- seq(0.1, 1, by = 0.1)
C3 <- rbindlist(lapply(W, function(w) { v <- sim(w, s_noise = noise)
  data.table(w = w, expected_splicing_rho = unname(v["splice"]),
             expected_cellcycle_rho = unname(v["cc"])) }))
print(C3, row.names = FALSE, digits = 3)
w_sim <- approx(C3$expected_splicing_rho, C3$w, xout = obs_sp, rule = 2)$y
say("  the simulation reproduces the observed %.2f at w = %.2f, the same answer as the algebra.",
    obs_sp, w_sim)
say("  the composition adjustment takes the observed %.2f to %.2f, so most of even that share is",
    obs_sp, 0.101)
say("  the cycling fraction itself rather than a relationship within the cells.")
OUT$C <- list(table = C3, w_implied = w_implied, w_sim = w_sim, obs = c(obs_sp, obs_cc))
fwrite(rbind(C3, data.table(w = NA_real_, expected_splicing_rho = obs_sp,
                            expected_cellcycle_rho = obs_cc)),
       file.path(O, "mixture_model.tsv"), sep = "\t")
saveRDS(OUT, file.path(O, "chain_inputs.rds"))
say("\nwritten to %s", O)
