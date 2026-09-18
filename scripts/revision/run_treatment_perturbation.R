## Do treatments that slow division lower the splicing score with it?
##
## The study's strongest limitation is that division rate was measured, not
## manipulated. The Cellular Lifespan Study does manipulate it: nine treatments
## slow or alter division in the same cell lines that provide the untreated
## trajectory. This is a post hoc analysis and is labelled as such; the plan below
## was fixed before any output was read.
##
##  Samples      all GSE179848 libraries with a finite division rate above zero
##  Treatments   every `cond` label with at least five libraries in at least two
##               cell lines; the others are reported but not modelled
##  Primary      splice ~ rate + days + line + oxygen + treatment
##               The paper's claim predicts treatment coefficients near zero: once
##               the measured rate is in the model, the treatment should add nothing.
##  Secondary    the same model for the 20-gene proliferation score, and
##               splice ~ proliferation + ... , which separates two failures:
##               splicing leaving the proliferation programme, or the programme
##               leaving the counted divisions
##  Descriptive  each treated culture against control cultures of the same line and
##               oxygen at the nearest time in culture: delta rate, delta splice,
##               delta proliferation
##  Reported     every treatment, in both directions, whatever the outcome
##
##  Oxygen (added 2026-09-17, post hoc). The treatment analysis holds oxygen fixed
##  by design, so the effect of 3% against 21% oxygen was never reported by it. It is
##  reported here in the same two ways as a treatment:
##    descriptive  each control culture at 3% oxygen against the control culture of
##                 the same cell line at 21% oxygen nearest in days grown
##    model        among control cultures, splice ~ rate + days + line + oxygen and
##                 splice ~ proliferation + days + line + oxygen
##  The older manuscript sentence (-0.02 divisions per day, splicing score 0.40
##  lower) came from run_division_rate.R, which took the median of every 3% library
##  minus the median of every 21% library of the same lines, treatments pooled and
##  time in culture unmatched. That contrast is recomputed here so the table shows
##  what it was made of.
##  oxygen_contrast.tsv: n = cultures at 3% oxygen; n_21 = distinct 21% cultures they
##  were matched to (unmatched rows: 21% libraries in the comparison); d_* = median of
##  3% minus 21% (unmatched rows: difference of medians); d_days = median absolute
##  difference in days grown; p_* = Wilcoxon signed-rank on the matched differences
##  (unmatched rows: rank-sum), per culture, cultures being serial passages of 3 lines.
##
##  Splicing score (corrected 2026-09-17). The score is the 177-gene pre-mRNA
##  processing score of Fig. 2, read from division_rate/sample_division_rate.tsv.
##  The run of 2026-09-16 took the `splice` column stored inside
##  hallmark_audit/audit_inputs.rds, which was written on 2026-09-13, before the
##  177-gene set existed, and is the earlier 96-gene score. Only the preregistered
##  GTEx figure uses the 96-gene set. The whole analysis is repeated on the 96-gene
##  score as a sensitivity (sensitivity_96gene_score/), and the outputs of the
##  2026-09-16 run are kept in previous_run_20260916_96gene_score/.
## Output: public_data_tierA/derived/treatment_perturbation/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table)})
D <- "public_data_tierA/derived"
O <- file.path(D, "treatment_perturbation"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z; M0 <- as.data.table(IN$meta)
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
M0[, prolif := colMeans(z[intersect(rownames(z), PROLIF), M0$sample, drop = FALSE])]

## ---- the splicing score: the one plotted in Fig. 2 ---------------------------
FIG2 <- fread(file.path(D, "division_rate/sample_division_rate.tsv"))
stopifnot(all(M0$sample %in% FIG2$sample))
M0[, splice96 := splice]
M0[, splice177 := FIG2$splice[match(sample, FIG2$sample)]]
stopifnot(max(abs(M0$rate - FIG2$rate[match(M0$sample, FIG2$sample)])) < 1e-12)
## 2026-09-18: the per-library scores are written out so Fig. 2 can show the 20-marker
## transcriptional proliferation score against the counted division rate directly.  The
## proxy is the only proliferation measure available outside this resource, so how well
## it tracks counting, and where it does not, belongs in the main figure.
fwrite(M0[, .(sample, line, cond, oxy, days, rate, prolif,
              splice = splice177, splice96)][order(cond, rate)],
       file.path(O, "library_scores_vs_counted_rate.tsv"), sep = "\t")

r177 <- cor(M0$rate, M0$splice177, method = "spearman")
r96  <- cor(M0$rate, M0$splice96,  method = "spearman")
say("splicing score against measured division rate, n = %d libraries: 177-gene rho = %.4f | stored 96-gene rho = %.4f | r between the two scores = %.4f",
    nrow(M0), r177, r96, cor(M0$splice177, M0$splice96))
stopifnot(round(r177, 2) == 0.71)          # the Fig. 2a value; stops if the score file is not the 177-gene one

sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = sum(is.finite(x) & is.finite(y))) }
spf <- function(s) sprintf("%+.2f (P = %.2g, n = %d)", s$rho, s$p, s$n)
wil <- function(v) if (length(v) < 2) NA_real_ else suppressWarnings(wilcox.test(v)$p.value)

analyse <- function(score, O, tag) {
  dir.create(O, showWarnings = FALSE, recursive = TRUE)
  say("\n\n######## %s ########", tag)
  M <- copy(M0); M[, splice := get(score)]
  M <- M[is.finite(rate) & rate > 0 & is.finite(splice)]
  say("libraries: %d | cell lines: %d | treatments: %d", nrow(M), uniqueN(M$line), uniqueN(M$cond))

  ## ---- which treatments can be modelled -------------------------------------
  cnt <- M[, .(n = .N, lines = uniqueN(line)), by = cond][order(-n)]
  cnt[, modelled := n >= 5 & lines >= 2]
  cnt[cond == "Control", modelled := TRUE]
  print(cnt)
  fwrite(cnt, file.path(O, "treatment_counts.tsv"), sep = "\t")
  say("treatments modelled: %d of %d (dropped: %s)", sum(cnt$modelled), nrow(cnt),
      paste(cnt[modelled == FALSE, cond], collapse = ", "))
  MM <- M[cond %in% cnt[modelled == TRUE, cond]]
  MM[, cond := relevel(factor(cond), ref = "Control")]
  MM[, line := factor(line)]; MM[, oxy := factor(oxy, levels = c(21, 3))]   # oxygen term reads 3% minus 21%

  ## ---- primary: residual treatment effect once the measured rate is in -------
  coefs <- function(fit, what) {
    s <- summary(fit)$coefficients
    k <- grep("^cond", rownames(s))
    data.table(response = what, term = sub("^cond", "", rownames(s)[k]),
               beta = s[k, 1], se = s[k, 2], p = s[k, 4])
  }
  f_splice_rate <- lm(splice ~ rate + days + line + oxy + cond, MM)
  f_splice_raw  <- lm(splice ~ days + line + oxy + cond, MM)
  f_prolif_rate <- lm(prolif ~ rate + days + line + oxy + cond, MM)
  f_splice_prol <- lm(splice ~ prolif + days + line + oxy + cond, MM)

  TAB <- rbindlist(list(coefs(f_splice_raw,  "splice, rate not in model"),
                        coefs(f_splice_rate, "splice, measured rate held constant"),
                        coefs(f_prolif_rate, "proliferation score, measured rate held constant"),
                        coefs(f_splice_prol, "splice, proliferation score held constant")))
  TAB[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  TAB[, FDR := p.adjust(p, "BH"), by = response]
  TAB[, explained_by_rate := lo <= 0 & hi >= 0]
  say("\n-- treatment coefficients (reference: Control) --")
  print(TAB[, .(response, term, beta = round(beta, 2), lo = round(lo, 2), hi = round(hi, 2),
                FDR = signif(FDR, 2))], nrows = 60)
  fwrite(TAB, file.path(O, "treatment_coefficients.tsv"), sep = "\t")
  say("treatments whose 95%% CI includes zero: %s",
      paste(TAB[, sprintf("%s %d of %d", response, sum(explained_by_rate), .N), by = response]$V1, collapse = " | "))

  say("\nvariance explained by the measured rate, with line/time/oxygen already in: %.3f",
      summary(f_splice_rate)$r.squared - summary(lm(splice ~ days + line + oxy, MM))$r.squared)

  ## ---- descriptive: each treated culture against its nearest control ---------
  ctl <- M[cond == "Control"]
  tre <- M[cond != "Control"]
  near <- rbindlist(lapply(seq_len(nrow(tre)), function(i) {
    x <- tre[i]
    c2 <- ctl[line == x$line & oxy == x$oxy]
    if (!nrow(c2)) return(NULL)
    j <- which.min(abs(c2$days - x$days))
    data.table(sample = x$sample, line = x$line, cond = x$cond, days = x$days,
               days_ctrl = c2$days[j], d_days = x$days - c2$days[j],
               d_rate = x$rate - c2$rate[j], d_splice = x$splice - c2$splice[j],
               d_prolif = x$prolif - c2$prolif[j])
  }))
  say("\nmatched treated cultures: %d of %d", nrow(near), nrow(tre))
  SUM <- near[, .(n = .N, d_rate = median(d_rate), d_splice = median(d_splice),
                  d_prolif = median(d_prolif), d_days = median(abs(d_days))), by = cond][order(d_rate)]
  print(SUM, digits = 2)
  fwrite(near, file.path(O, "matched_deltas.tsv"), sep = "\t")
  fwrite(SUM, file.path(O, "matched_summary.tsv"), sep = "\t")

  COR <- rbindlist(list(
    c(list(pair = "delta splice vs delta measured rate"),        sp(near$d_rate,   near$d_splice)),
    c(list(pair = "delta splice vs delta proliferation score"),  sp(near$d_prolif, near$d_splice)),
    c(list(pair = "delta proliferation score vs delta measured rate"), sp(near$d_rate, near$d_prolif))))
  COR[, unit := "matched treated cultures"]
  say("\nacross matched treated cultures (Spearman):")
  print(COR)
  fwrite(COR, file.path(O, "matched_correlations.tsv"), sep = "\t")

  ## ---- oxygen: 3% against 21%, the contrast the matching above holds fixed ----
  say("\n==== oxygen: 3%% against 21%% ====")
  print(M[, .N, by = .(oxy, cond)][order(oxy, cond)])
  match_oxygen <- function(arm) {
    lo <- M[cond == arm & oxy == 3]; hi <- M[cond == arm & oxy == 21]
    rbindlist(lapply(seq_len(nrow(lo)), function(i) {
      x <- lo[i]; c2 <- hi[line == x$line]
      if (!nrow(c2)) return(NULL)
      j <- which.min(abs(c2$days - x$days))
      data.table(arm = arm, sample = x$sample, line = x$line, days = x$days,
                 sample_21 = c2$sample[j], days_21 = c2$days[j], d_days = x$days - c2$days[j],
                 rate_3 = x$rate, rate_21 = c2$rate[j], d_rate = x$rate - c2$rate[j],
                 splice_3 = x$splice, splice_21 = c2$splice[j], d_splice = x$splice - c2$splice[j],
                 prolif_3 = x$prolif, prolif_21 = c2$prolif[j], d_prolif = x$prolif - c2$prolif[j])
    }))
  }
  srow <- function(d, contrast, line = "all") data.table(
    contrast = contrast, matched = TRUE, line = line, n_lines = uniqueN(d$line),
    n = nrow(d), n_21 = uniqueN(d$sample_21),
    d_rate = median(d$d_rate), d_splice = median(d$d_splice), d_prolif = median(d$d_prolif),
    d_days = median(abs(d$d_days)), max_abs_d_days = max(abs(d$d_days)),
    median_days_3 = median(d$days), median_days_21 = median(d$days_21),
    n_rate_lower = sum(d$d_rate < 0), n_splice_lower = sum(d$d_splice < 0), n_prolif_lower = sum(d$d_prolif < 0),
    p_rate = wil(d$d_rate), p_splice = wil(d$d_splice), p_prolif = wil(d$d_prolif))
  urow <- function(a, b, contrast) data.table(
    contrast = contrast, matched = FALSE, line = "all", n_lines = uniqueN(a$line),
    n = nrow(a), n_21 = nrow(b),
    d_rate = median(a$rate) - median(b$rate), d_splice = median(a$splice) - median(b$splice),
    d_prolif = median(a$prolif) - median(b$prolif),
    d_days = NA_real_, max_abs_d_days = NA_real_,
    median_days_3 = median(a$days), median_days_21 = median(b$days),
    n_rate_lower = NA_integer_, n_splice_lower = NA_integer_, n_prolif_lower = NA_integer_,
    p_rate = suppressWarnings(wilcox.test(a$rate, b$rate)$p.value),
    p_splice = suppressWarnings(wilcox.test(a$splice, b$splice)$p.value),
    p_prolif = suppressWarnings(wilcox.test(a$prolif, b$prolif)$p.value))

  oc <- match_oxygen("Control"); oi <- match_oxygen("Contact_Inhibition")
  say("control cultures at 3%% oxygen matched to a 21%% control of the same line: %d of %d",
      nrow(oc), nrow(M[cond == "Control" & oxy == 3]))
  print(oc[, .(sample, line, days, days_21, d_rate = round(d_rate, 3), d_splice = round(d_splice, 2),
               d_prolif = round(d_prolif, 2))])
  lines3 <- unique(M[oxy == 3, line])
  a_all <- M[oxy == 3];                    b_all <- M[oxy == 21 & line %in% lines3]
  a_ctl <- M[oxy == 3 & cond == "Control"]; b_ctl <- M[oxy == 21 & cond == "Control" & line %in% lines3]
  b_win <- b_ctl[days >= min(a_ctl$days) - 7 & days <= max(a_ctl$days) + 7]
  OX <- rbindlist(c(
    list(srow(oc, "control, 3% minus 21% oxygen, same line, nearest days grown")),
    lapply(split(oc, by = "line"), function(d) srow(d, "control, 3% minus 21% oxygen, same line, nearest days grown", d$line[1])),
    list(srow(oc[abs(d_days) <= 7], "control, 3% minus 21% oxygen, matches within 7 days grown only"),
         srow(oi, "contact inhibition, 3% minus 21% oxygen, same line, nearest days grown"),
         urow(a_ctl, b_ctl, "control, unmatched: median of 3% minus median of 21%, same lines, all days"),
         urow(a_ctl, b_win, "control, unmatched, 21% cultures within 7 days of the 3% time span"),
         urow(a_all, b_all, "older sentence (run_division_rate.R): every 3% library minus every 21% library of the same lines, treatments pooled"))))
  print(OX[, .(contrast = substr(contrast, 1, 60), line, n, n_21, d_rate = round(d_rate, 4),
               d_splice = round(d_splice, 4), d_prolif = round(d_prolif, 4), d_days,
               p_rate = signif(p_rate, 2), p_splice = signif(p_splice, 2), p_prolif = signif(p_prolif, 2))])
  say("composition of the pooled contrast: 3%% libraries = %s | 21%% libraries = %s",
      paste(a_all[, .N, by = cond][, sprintf("%s %d", cond, N)], collapse = ", "),
      paste(b_all[, .N, by = cond][order(-N)][, sprintf("%s %d", cond, N)], collapse = ", "))
  say("days grown, pooled contrast: 3%% median %.1f (range %.1f-%.1f) | 21%% median %.1f (range %.1f-%.1f)",
      median(a_all$days), min(a_all$days), max(a_all$days), median(b_all$days), min(b_all$days), max(b_all$days))
  fwrite(OX, file.path(O, "oxygen_contrast.tsv"), sep = "\t")
  fwrite(rbind(oc, oi), file.path(O, "oxygen_matched_deltas.tsv"), sep = "\t")

  ## the model of the primary analysis, oxygen as the term of interest
  ocoef <- function(fit, response, adjusted, set) {
    s <- summary(fit)$coefficients; k <- "oxy3"
    data.table(sample_set = set, response = response, held_constant = adjusted,
               term = "3% oxygen (reference 21%)", n = nobs(fit), n_3 = sum(model.frame(fit)$oxy == "3"),
               beta = s[k, 1], se = s[k, 2], p = s[k, 4],
               lo = s[k, 1] - 1.96 * s[k, 2], hi = s[k, 1] + 1.96 * s[k, 2])
  }
  fit_set <- function(d, set, extra = "") {
    d <- copy(d); d[, line := factor(line)]; d[, oxy := factor(oxy, levels = c(21, 3))]
    f <- function(lhs, rhs) lm(as.formula(paste(lhs, "~", rhs, "+ days + line + oxy", extra)), d)
    rbindlist(list(
      ocoef(lm(as.formula(paste("rate ~ days + line + oxy", extra)), d),   "measured division rate", "none", set),
      ocoef(lm(as.formula(paste("prolif ~ days + line + oxy", extra)), d), "proliferation score", "none", set),
      ocoef(lm(as.formula(paste("splice ~ days + line + oxy", extra)), d), "splicing score", "none", set),
      ocoef(f("splice", "rate"),   "splicing score", "measured division rate", set),
      ocoef(f("splice", "prolif"), "splicing score", "proliferation score", set),
      ocoef(f("prolif", "rate"),   "proliferation score", "measured division rate", set)))
  }
  OXM <- rbindlist(list(
    fit_set(M[cond == "Control"], "control cultures, all seven lines"),
    fit_set(M[cond == "Control" & line %in% lines3], "control cultures, the three lines grown at both oxygen tensions"),
    fit_set(MM, "all modelled libraries (the primary treatment model)", "+ cond")))
  say("\n-- oxygen coefficient, 3%% minus 21%% --")
  print(OXM[, .(sample_set = substr(sample_set, 1, 34), response, held_constant, n, n_3,
                beta = round(beta, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 2))])
  fwrite(OXM, file.path(O, "oxygen_models.tsv"), sep = "\t")

  list(TAB = TAB, SUM = SUM, COR = COR, OX = OX, OXM = OXM, n_matched = nrow(near))
}

A <- analyse("splice177", O, "PRIMARY: 177-gene pre-mRNA processing score (the score of Fig. 2)")
B <- analyse("splice96", file.path(O, "sensitivity_96gene_score"),
             "SENSITIVITY: the 96-gene score stored in audit_inputs.rds (what the 2026-09-16 run used)")

## ---- every quantity the text cites, under both scores ------------------------
pick <- function(R) {
  tc <- function(resp, trt) R$TAB[response == resp & term == trt, beta]
  sm <- function(trt, col) R$SUM[cond == trt][[col]]
  ox <- R$OX[contrast == "control, 3% minus 21% oxygen, same line, nearest days grown" & line == "all"]
  old <- R$OX[grepl("^older sentence", contrast)]
  om <- function(set, resp, adj) R$OXM[grepl(set, sample_set) & response == resp & held_constant == adj, beta]
  HR <- "splice, measured rate held constant"; HP <- "splice, proliferation score held constant"
  c(`matched treated cultures (n)` = R$n_matched,
    `rho, delta splice vs delta measured rate` = R$COR$rho[1],
    `rho, delta splice vs delta proliferation score` = R$COR$rho[2],
    `rho, delta proliferation score vs delta measured rate` = R$COR$rho[3],
    `contact inhibition, median delta rate` = sm("Contact_Inhibition", "d_rate"),
    `contact inhibition, median delta splice` = sm("Contact_Inhibition", "d_splice"),
    `rate held constant: DEX` = tc(HR, "DEX"), `rate held constant: oligomycin` = tc(HR, "Oligomycin"),
    `rate held constant: oligomycin + DEX` = tc(HR, "Oligomycin+DEX"),
    `rate held constant: mitoNUITs + DEX` = tc(HR, "mitoNUITs+DEX"),
    `rate held constant: treatments with CI including zero` = R$TAB[response == HR, sum(explained_by_rate)],
    `proliferation held constant: oligomycin` = tc(HP, "Oligomycin"),
    `proliferation held constant: 2-deoxyglucose` = tc(HP, "2-Deoxyglucose"),
    `proliferation held constant: treatments with CI including zero` = R$TAB[response == HP, sum(explained_by_rate)],
    `proliferation score, rate held constant: DEX` = tc("proliferation score, measured rate held constant", "DEX"),
    `oxygen, matched controls (n)` = ox$n,
    `oxygen, matched controls, median delta rate` = ox$d_rate,
    `oxygen, matched controls, median delta splice` = ox$d_splice,
    `oxygen, matched controls, median delta proliferation score` = ox$d_prolif,
    `oxygen, control model, splice with rate held constant` = om("all seven", "splicing score", "measured division rate"),
    `oxygen, control model, splice with proliferation held constant` = om("all seven", "splicing score", "proliferation score"),
    `oxygen, older pooled contrast, delta rate` = old$d_rate,
    `oxygen, older pooled contrast, delta splice` = old$d_splice)
}
CMP <- data.table(quantity = names(pick(A)), score_177_gene = pick(A), score_96_gene = pick(B))
say("\n\n==== cited quantities under the two scores ====")
print(CMP, digits = 3)
fwrite(CMP, file.path(O, "score_sensitivity.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
