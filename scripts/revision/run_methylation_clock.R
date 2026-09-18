## Methylation clocks in the same cultures that carry a measured division rate.
## Analysis plan fixed in preregistration_methylation_clock_ko.md (SHA-256 recorded
## before the data were downloaded); this script checks the digest before it runs.
##   M1  sanity: DNAm age rises with cumulative population doublings
##   M2  clock acceleration is associated with the measured division rate
##   M3  in the same cultures, the transcriptomic proliferation score tracks DNAm age
##   M4  the splicing score's association with DNAm age shrinks when division rate is
##       held constant
## Outputs: public_data_tierA/derived/methylation_clock/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table)})
D <- "public_data_tierA/derived"; MET <- "public_data_tierA/methylation"
O <- file.path(D, "methylation_clock"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

PREREG <- "preregistration_methylation_clock_ko.md"
LOCKED <- "be51a7550c982e2071ebdb80f36ab7098965351c6e14a65434e8df275505d929"
dig <- system(sprintf("sha256sum %s | cut -d' ' -f1", PREREG), intern = TRUE)
if (!identical(dig, LOCKED)) stop("preregistration digest changed: ", dig)
say("preregistration digest verified: %s", substr(dig, 1, 16))

## ---- betas, with the detection filter the preregistration fixed ---------------
B <- fread(cmd = paste("zcat", file.path(MET, "clock_cpg_betas.tsv.gz")))
cpg <- B[[1]]
nm <- names(B)[-1]
is_beta <- grepl(" beta$", nm); is_p <- grepl(" Detection Pval$", nm)
sent <- sub(" beta$", "", nm[is_beta])
BETA <- as.matrix(B[, nm[is_beta], with = FALSE]); rownames(BETA) <- cpg; colnames(BETA) <- sent
PVAL <- as.matrix(B[, nm[is_p], with = FALSE]); colnames(PVAL) <- sub(" Detection Pval$", "", nm[is_p])
PVAL <- PVAL[, sent, drop = FALSE]
say("clock CpGs: %d | samples: %d", nrow(BETA), ncol(BETA))
BETA[PVAL > 0.01] <- NA
say("measurements set to missing by the detection filter: %.2f%%", 100 * mean(is.na(BETA)))
frac_missing <- colMeans(is.na(BETA))
drop <- frac_missing > 0.30
say("samples dropped for >30%% missing clock CpGs: %d", sum(drop))
BETA <- BETA[, !drop, drop = FALSE]
## remaining gaps filled with the cohort mean of that CpG, as preregistered
nfill <- sum(is.na(BETA))
rmean <- rowMeans(BETA, na.rm = TRUE)
for (i in which(rowSums(is.na(BETA)) > 0)) BETA[i, is.na(BETA[i, ])] <- rmean[i]
say("single measurements filled with the CpG's cohort mean: %d (%.3f%%)", nfill, 100 * nfill / length(BETA))

## ---- clocks ------------------------------------------------------------------
coefdir <- file.path(MET, "clock_coefficients")
readcoef <- function(f) { x <- fread(file.path(coefdir, f)); setnames(x, 1:2, c("cpg", "coef")); x[] }
anti_trafo <- function(x, adult = 20) ifelse(x < 0, (1 + adult) * exp(x) - 1, (1 + adult) * x + adult)
clock <- function(tab, transform = FALSE) {
  ic <- tab[grepl("intercept", cpg, ignore.case = TRUE), coef]
  ic <- if (length(ic)) as.numeric(ic[1]) else 0
  tb <- tab[grepl("^cg", cpg)]
  g <- intersect(tb$cpg, rownames(BETA))
  lp <- ic + as.vector(t(BETA[g, , drop = FALSE]) %*% tb$coef[match(g, tb$cpg)])
  list(age = if (transform) anti_trafo(lp) else lp, n = length(g), n_total = nrow(tb))
}
H1 <- clock(readcoef("Horvath2013_AdditionalFile3.csv"), transform = TRUE)
H2 <- clock(readcoef("Horvath2.csv"), transform = TRUE)
PA <- clock(readcoef("PhenoAge.csv"))
HA <- clock(readcoef("Hannum.csv"))
for (x in list(c("Horvath 2013", H1$n, H1$n_total), c("skin & blood", H2$n, H2$n_total),
               c("PhenoAge", PA$n, PA$n_total), c("Hannum", HA$n, HA$n_total)))
  say("  %-14s CpGs available: %s of %s", x[1], x[2], x[3])

M <- data.table(sentrix = colnames(BETA), horvath = H1$age, skinblood = H2$age,
                phenoage = PA$age, hannum = HA$age)

## ---- sample identity and culture metadata ------------------------------------
gs <- fread(file.path(MET, "gsm_sentrix.tsv"), header = FALSE, col.names = c("gsm", "sentrix"))
ti <- fread(file.path(MET, "gsm_titles.tsv"), header = FALSE, col.names = c("gsm", "title"))
key <- merge(gs, ti, by = "gsm")
key[, uvn := sub("^Methylation_", "", title)]
M <- merge(M, key[, .(sentrix, gsm, uvn)], by = "sentrix")

brief <- readLines(file.path(MET, "GSE179847_samples_brief.txt"))
idx <- grep("^\\^SAMPLE", brief)
fields <- c("cell_line", "clinical_condition", "treatments", "percent_oxygen", "passage",
            "days_grown_udays", "population_doublings_utotal_divisions",
            "population_doubling_time_uhours_per_division", "unique_variable_name")
rec <- rbindlist(lapply(seq_along(idx), function(i) {
  blk <- brief[idx[i]:(if (i < length(idx)) idx[i + 1] - 1 else length(brief))]
  g <- function(f) { h <- grep(paste0("= ", f, ": "), blk, fixed = TRUE, value = TRUE)
    if (length(h)) sub(paste0(".*= ", f, ": "), "", h[1]) else NA_character_ }
  as.list(setNames(c(sub("^\\^SAMPLE = ", "", blk[1]), vapply(fields, g, character(1))), c("gsm", fields)))
}))
num <- c("percent_oxygen", "days_grown_udays", "population_doublings_utotal_divisions",
         "population_doubling_time_uhours_per_division")
rec[, (num) := lapply(.SD, as.numeric), .SDcols = num]
M <- merge(M, rec, by = "gsm")
M[, rate := 24 / population_doubling_time_uhours_per_division]
M <- M[is.finite(rate) & rate > 0]
say("\nmethylation samples with a measured division rate: %d", nrow(M))

## ---- RNA scores for the same cultures ----------------------------------------
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z; meta <- IN$meta
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
rna <- data.table(sample = colnames(z),
                  prolif = colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE]),
                  splice = colMeans(z[intersect(rownames(z), CORE), , drop = FALSE]),
                  rate_rna = meta$rate)
gm <- fread(file.path(D, "audit/all_geo_samples.tsv"))[dataset == "GSE179848"]
gm[, sample := paste0("Sample_", rnaseq_sampleid)]
rna <- merge(rna, gm[, .(sample, uvn = unique_variable_name)], by = "sample")
P <- merge(M, rna, by = "uvn")
say("cultures with both methylation and RNA: %d", nrow(P))

## ---- the four preregistered tests ---------------------------------------------
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(x) & is.finite(y))) }
res <- list()
HC <- M[grepl("^HC", cell_line) & treatments == "Control" & percent_oxygen == 21]
res$M1 <- sp(HC$horvath, HC$population_doublings_utotal_divisions)
## acceleration: residual of DNAm age on days in culture, within cell line
M[, accel := residuals(lm(horvath ~ days_grown_udays + cell_line, data = .SD, na.action = na.exclude))]
res$M2 <- sp(M$accel, M$rate)
res$M3 <- sp(P$prolif, P$horvath)
r_sp <- sp(P$splice, P$horvath)
pres <- residuals(lm(splice ~ rate, data = P)); hres <- residuals(lm(horvath ~ rate, data = P))
r_partial <- sp(pres, hres)
res$M4 <- c(rho = unname(r_partial["rho"]), p = unname(r_partial["p"]), n = unname(r_partial["n"]))

TAB <- data.table(
  test = c("M1 DNAm age vs cumulative doublings (untreated healthy)",
           "M2 clock acceleration vs measured division rate",
           "M3 proliferation score vs DNAm age (paired cultures)",
           "M4 splicing vs DNAm age, division rate held constant"),
  rho = sapply(res, `[`, "rho"), p = sapply(res, `[`, "p"), n = sapply(res, `[`, "n"))
TAB[, FDR := p.adjust(p, "BH")]
## the criteria are the ones fixed in the preregistration: M1-M3 need |rho| >= 0.3
## with P < 0.05; M4 needs the association to shrink by at least half when the
## division rate is held constant, which is a different question from significance
attenuation <- 1 - abs(r_partial["rho"]) / abs(r_sp["rho"])
TAB[, supported := c(abs(rho[1]) >= 0.3 & FDR[1] < 0.05,
                     abs(rho[2]) >= 0.3 & FDR[2] < 0.05,
                     abs(rho[3]) >= 0.3 & FDR[3] < 0.05,
                     attenuation >= 0.5)]
say("M4 attenuation when division rate is held constant: %.0f%% (criterion: at least 50%%)",
    100 * attenuation)
print(TAB, digits = 3)
fwrite(TAB, file.path(O, "preregistered_tests.tsv"), sep = "\t")

say("\nfor M4, the unadjusted association: rho = %+.2f (P = %.2g); holding division rate constant: rho = %+.2f (P = %.2g)",
    r_sp["rho"], r_sp["p"], r_partial["rho"], r_partial["p"])
say("other clocks against the measured division rate (acceleration):")
for (cl in c("skinblood", "phenoage", "hannum")) {
  M[, acc2 := residuals(lm(get(cl) ~ days_grown_udays + cell_line, data = .SD, na.action = na.exclude))]
  s <- sp(M$acc2, M$rate); say("   %-10s rho = %+.2f (P = %.2g)", cl, s["rho"], s["p"]) }
say("paired cultures: splicing vs proliferation rho = %+.2f | DNAm age vs measured rate rho = %+.2f",
    sp(P$splice, P$prolif)["rho"], sp(P$horvath, P$rate)["rho"])
fwrite(M, file.path(O, "methylation_samples.tsv"), sep = "\t")
fwrite(P, file.path(O, "paired_rna_methylation.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
