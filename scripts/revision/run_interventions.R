## Interventions (Supplementary Fig. S6; was main-text Fig. 8, before that Fig. 9)
## recomputed on the revised contrast table.
## Changes from the published version:
##   - GSE149694 "Fibroblast-D7 vs Fibroblast-D3" leaves the reprogramming class:
##     both arms are Sendai-OSKM-transduced cells still in fibroblast medium
##     (the series switched to reprogramming media on day 8), so the contrast is
##     four days of culture, not a reprogramming condition. Reprogramming: 8 -> 7.
##   - the donor-age contrast is the one recomputed in the redefined cohort.
##   - 2026-09-17: the "predicted" change of an intervention was the fitted value
##     of lm(spl ~ cc) over all 63 contrasts, which include the 21 interventions
##     themselves (in-sample). Each intervention class is now also predicted from
##     a fit that excludes that class (leave-class-out), from a fit that excludes
##     all 21 interventions at once, and from a fit that also excludes every contrast
##     of the same GEO series. Mean class residuals are given per contrast and per
##     series with the shift each could have detected. The 63-contrast coupling is
##     reported for the public contrasts only (the unreplicated local contrast was
##     and the named-factor means with and without it.
## Outputs: public_data_tierA/derived/interventions_revised/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; O <- file.path(D, "interventions_revised")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
fisher_ci <- function(r, n) tanh(atanh(r) + c(-1.96, 1.96) / sqrt(n - 3))
SFFILE <- "scripts/revision/named_splicing_factors.txt"
SF <- if (file.exists(SFFILE)) readLines(SFFILE) else
  c("SRSF1","SRSF2","SRSF3","SRSF5","SRSF6","SRSF7","HNRNPA1","HNRNPA2B1","HNRNPD",
    "HNRNPK","HNRNPM","SF1","SF3B1","SF3A3","SF3B5","U2AF1","U2AF2","PRPF8","PRPF3",
    "SNRPA","SNRPB","SNRPD2","SNRPF","RBM39","TRA2B","PTBP1")
say("named splicing factors: %d (%s)", length(SF), if (file.exists(SFFILE)) "sourced list" else "provisional list")

R <- read.delim(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
R$class[is.na(R$class)] <- "metabolic / culture"
R$name <- N$contrast[match(R$id, N$id)]
R$name[is.na(R$name)] <- CP$contrast[match(R$id[is.na(R$name)], CP$id)]
## series of origin, to count independent studies within a class
R$series <- N$dataset[match(R$id, N$id)]
R$series[is.na(R$series)] <- sub("_.*$", "", R$id[is.na(R$series)])

## (1) in-sample: fitted values of the 63-contrast regression, which contains the
##     interventions it is then said to predict
f <- lm(spl ~ cc, R); R$resid <- resid(f); R$pred <- predict(f)

## (2) leave-class-out: each intervention class is predicted from a fit to the
##     contrasts outside that class. (3) all 21 interventions left out at once.
INT <- c("reprogramming", "secretome")
R$pred_lco <- NA_real_
for (k in INT) { fk <- lm(spl ~ cc, R[R$class != k, ])
  R$pred_lco[R$class == k] <- predict(fk, R[R$class == k, ]) }
f_lao <- lm(spl ~ cc, R[!R$class %in% INT, ])
R$pred_lao <- ifelse(R$class %in% INT, predict(f_lao, R), NA_real_)
## (4) strictest: the class AND every contrast of the same GEO series are left out
##     (GSE149694 Fibroblast-D7, a culture contrast, shares its reference arm with
##     the five GSE149694 reprogramming media)
R$pred_lcso <- NA_real_
for (i in which(R$class %in% INT))
  R$pred_lcso[i] <- predict(lm(spl ~ cc, R[R$class != R$class[i] & R$series != R$series[i], ]), R[i, ])
R$resid_lco <- R$spl - R$pred_lco; R$resid_lao <- R$spl - R$pred_lao; R$resid_lcso <- R$spl - R$pred_lcso

S <- R[R$class %in% c("reprogramming", "secretome"), ]
S <- S[order(-S$spl), ]
print(data.frame(class = substr(S$class, 1, 13), cellcycle = round(S$cc, 3),
                 splicing = round(S$spl, 3), resid_insample = round(S$resid, 3),
                 resid_leave_class_out = round(S$resid_lco, 3), resid_all_22_out = round(S$resid_lao, 3),
                 series = S$series, name = substr(S$name, 1, 40)), row.names = FALSE)
write.table(S, file.path(O, "reprog_secretome_splicing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

CLS <- do.call(rbind, lapply(c("reprogramming", "secretome"), function(k) {
  s <- R[R$class == k, ]; ct <- suppressWarnings(cor.test(s$cc, s$spl, method = "spearman", exact = FALSE))
  ci <- fisher_ci(unname(ct$estimate), nrow(s))
  m <- lm(spl ~ cc, s); mo <- lm(spl ~ cc, R[R$class != k, ])
  data.frame(class = k, n = nrow(s), rho = unname(ct$estimate), lo = ci[1], hi = ci[2], p = ct$p.value,
             median_splicing = median(s$spl), median_residual = median(s$resid),
             n_series = length(unique(s$series)),
             slope_within_class = unname(coef(m)[2]), slope_lo = confint(m)["cc", 1], slope_hi = confint(m)["cc", 2],
             slope_outside_class = unname(coef(mo)[2]), slope_outside_lo = confint(mo)["cc", 1],
             slope_outside_hi = confint(mo)["cc", 2]) }))
print(CLS, digits = 3, row.names = FALSE)
write.table(CLS, file.path(O, "class_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- the fits behind each prediction scheme --------------------------------
fit_row <- function(lbl, d) { m <- lm(spl ~ cc, d); ci <- confint(m)["cc", ]
  ct <- suppressWarnings(cor.test(d$cc, d$spl, method = "spearman", exact = FALSE))
  rc <- fisher_ci(unname(ct$estimate), nrow(d))
  data.frame(fit = lbl, n_contrasts = nrow(d), n_series = length(unique(d$series)),
             intercept = unname(coef(m)[1]), slope = unname(coef(m)[2]), slope_lo = ci[1], slope_hi = ci[2],
             r2 = summary(m)$r.squared, rho = unname(ct$estimate), rho_lo = rc[1], rho_hi = rc[2],
             rho_p = ct$p.value, row.names = NULL) }
FITS <- rbind(fit_row("all 63 contrasts (in-sample)", R),
              fit_row("without reprogramming class", R[R$class != "reprogramming", ]),
              fit_row("without secretome class", R[R$class != "secretome", ]),
              fit_row("without all 21 interventions", R[!R$class %in% INT, ]))
say("\n=== the regression of splicing-set change on cell-cycle change, by what is left out ===")
print(FITS, digits = 4, row.names = FALSE)
write.table(FITS, file.path(O, "prediction_fits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## the same fit with the 23 genes shared with the cell-cycle set removed
Rn <- transform(R, spl = spl_noshare)
SHARED <- rbind(cbind(splicing_set = "177 genes", fit_row("all 63 contrasts (in-sample)", R)),
                cbind(splicing_set = "genes shared with the cell-cycle set removed",
                      fit_row("all 63 contrasts (in-sample)", Rn)))
say("\n63 contrasts: rho = %.4f [%.3f, %.3f], P = %.2e, R2 = %.4f, slope = %.4f",
    SHARED$rho[1], SHARED$rho_lo[1], SHARED$rho_hi[1], SHARED$rho_p[1], SHARED$r2[1], SHARED$slope[1])
say("   shared genes removed:     rho = %.4f [%.3f, %.3f], P = %.2e, R2 = %.4f, slope = %.4f",
    SHARED$rho[2], SHARED$rho_lo[2], SHARED$rho_hi[2], SHARED$rho_p[2], SHARED$r2[2], SHARED$slope[2])
write.table(SHARED, file.path(O, "coupling_shared_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- how much of an intervention effect the cell-cycle response predicts ----
## median absolute residual against median absolute change, for each scheme
SCHEME <- c(`in-sample (63-contrast fit)` = "resid", `leave-class-out` = "resid_lco",
            `all 21 interventions left out` = "resid_lao",
            `class and same-series contrasts left out` = "resid_lcso")
GROUPS <- list(`all interventions` = INT, reprogramming = "reprogramming", secretome = "secretome")
ACC <- do.call(rbind, lapply(names(SCHEME), function(sc) do.call(rbind, lapply(names(GROUPS), function(g) {
  s <- R[R$class %in% GROUPS[[g]], ]; e <- s[[SCHEME[[sc]]]]
  data.frame(scheme = sc, group = g, n_contrasts = nrow(s), n_series = length(unique(s$series)),
             median_abs_change = median(abs(s$spl)), median_abs_residual = median(abs(e)),
             ratio_residual_to_change = median(abs(e)) / median(abs(s$spl)),
             share_predicted = 1 - median(abs(e)) / median(abs(s$spl)),
             n_within_0.05 = sum(abs(e) <= 0.05), n_beyond_0.10 = sum(abs(e) > 0.10),
             max_abs_shift_vs_insample = max(abs(e - s$resid))) }))))
say("\n=== median absolute residual against median absolute change ===")
print(ACC, digits = 3, row.names = FALSE)
write.table(ACC, file.path(O, "prediction_accuracy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## the same question as a rank correlation: across the intervention contrasts alone,
## does the change in the splicing set follow the change in cell-cycle genes?
WITHIN <- do.call(rbind, lapply(names(GROUPS), function(g) { s <- R[R$class %in% GROUPS[[g]], ]
  ct <- suppressWarnings(cor.test(s$cc, s$spl, method = "spearman", exact = FALSE)); ci <- fisher_ci(unname(ct$estimate), nrow(s))
  pr <- cor.test(s$pred_lco, s$spl)
  data.frame(group = g, n_contrasts = nrow(s), rho_cellcycle_vs_splicing = unname(ct$estimate), lo = ci[1], hi = ci[2],
             p = ct$p.value, pearson_observed_vs_leave_class_out_prediction = unname(pr$estimate),
             pearson_lo = pr$conf.int[1], pearson_hi = pr$conf.int[2], pearson_p = pr$p.value) }))
say("\n=== within the intervention contrasts only ===")
print(WITHIN, digits = 3, row.names = FALSE)
write.table(WITHIN, file.path(O, "coupling_within_interventions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- mean class residual, its 95% CI, and the shift it could have detected ----
## unit = contrast treats the contrasts of a class as independent; unit = series
## first averages the contrasts of one GEO series (5 of the 7 reprogramming contrasts
## share one reference arm; the 14 secretome contrasts come from 8 series)
mde80 <- function(n, s) if (n < 2 || !is.finite(s) || s == 0) NA_real_ else
  tryCatch(power.t.test(n = n, sd = s, sig.level = 0.05, power = 0.80, type = "one.sample")$delta,
           error = function(e) NA_real_)
CR <- do.call(rbind, lapply(names(SCHEME), function(sc) do.call(rbind, lapply(c("reprogramming", "secretome"), function(k) {
  s <- R[R$class == k, ]; s$e <- s[[SCHEME[[sc]]]]
  do.call(rbind, lapply(c("contrast", "series"), function(u) {
    v <- if (u == "contrast") s$e else tapply(s$e, s$series, mean)
    tt <- t.test(v)
    data.frame(scheme = sc, class = k, unit = u, n = length(v), mean_residual = mean(v),
               lo = tt$conf.int[1], hi = tt$conf.int[2], half_width = diff(tt$conf.int) / 2,
               sd = sd(v), p = tt$p.value, detectable_shift_80pct_power = mde80(length(v), sd(v)),
               class_median_abs_change = median(abs(s$spl))) })) }))))
say("\n=== mean class residual with 95%% CI (one-sample t) ===")
print(CR, digits = 3, row.names = FALSE)
write.table(CR, file.path(O, "class_residuals_in_and_out_of_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (k in c("senescence", "donor age", "reprogramming", "secretome"))
  say("median splicing-factor change, %-14s %+.3f (n = %d)", k, median(R$spl[R$class == k]), sum(R$class == k))

## interventions that lie furthest from what their proliferative effect predicts,
## in both directions (reviewer 1 asked for the dissociations that go the other way)
say("\nlargest residuals in each direction, all 63 contrasts:")
print(data.frame(id = R$id, class = substr(R$class, 1, 13), cc = round(R$cc, 3), spl = round(R$spl, 3),
                 resid = round(R$resid, 3), name = substr(R$name, 1, 40))[order(-R$resid), ][c(1:6, (nrow(R) - 5):nrow(R)), ],
      row.names = FALSE)

## ---- named splicing factors across the reprogramming-related contrasts -----
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (fl in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
  d <- read.delim(fl); vec[[sub("\\.tsv$", "", basename(fl))]] <- setNames(d$logFC, toupper(d$gene)) }
rid <- S$id[S$class == "reprogramming"]
M <- sapply(rid, function(i) { v <- vec[[i]]; if (is.null(v)) return(rep(NA, length(SF))); v[match(SF, names(v))] })
rownames(M) <- SF
say("\n=== named splicing factors across %d reprogramming-related contrasts (log2 FC) ===", ncol(M))
keep <- rowSums(is.finite(M)) >= ncol(M) * 0.6
print(round(M[keep, , drop = FALSE], 2))
mean_resp <- sort(rowMeans(M, na.rm = TRUE))
say("consistently falling: %s", paste(sprintf("%s %+.2f", names(head(mean_resp, 3)), head(mean_resp, 3)), collapse = ", "))
say("consistently rising:  %s", paste(sprintf("%s %+.2f", names(tail(mean_resp, 5)), tail(mean_resp, 5)), collapse = ", "))
write.table(data.frame(gene = rownames(M), M, check.names = FALSE),
            file.path(O, "named_SF_reprogramming.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
MEANS <- data.frame(gene = rownames(M), n_measured = rowSums(is.finite(M)), mean_response = rowMeans(M, na.rm = TRUE),
                    n_rising = rowSums(M > 0, na.rm = TRUE))
MEANS <- MEANS[order(-MEANS$mean_response), ]
say("\n=== mean response of each named factor across the %d reprogramming contrasts ===", ncol(M))
print(MEANS, digits = 2, row.names = FALSE)
write.table(MEANS, file.path(O, "named_SF_mean_response.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

