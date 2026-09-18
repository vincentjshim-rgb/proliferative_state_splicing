#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Held-out generalization, executed exactly as fixed in
#   repro_cm_preregistration_phase2_ko.md  (SHA256 41e4cf4b...)
# 6 primary tests, BH correction across those 6, all results reported.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db); library(edgeR); library(readxl)
})
root <- "public_data_tierA"; ho <- file.path(root, "heldout")
out  <- file.path(root, "derived", "heldout_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")

## ---- fixed local anchor (never re-fitted) ---------------------------------
anchor <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                     check.names = FALSE)
anchor <- anchor[anchor$comparator_consistent & anchor$Repro_specific_score != 0 &
                 is.finite(anchor$Repro_specific_score), ]
prolif <- read.delim(file.path(root, "derived/direction_probe/P3_proliferation_loading_GSE179848.tsv"))
anchor$proliferation_loading <- prolif$proliferation_loading[match(anchor$gene, prolif$gene)]
cat("anchor genes:", nrow(anchor), "\n")

## ---- project-standard expression-decile-matched permutation ---------------
matched_permutation <- function(x, y, expr, nperm = 10000L, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average"); ry <- rank(y, ties.method = "average")
  br <- unique(quantile(expr, probs = seq(0, 1, 0.1), na.rm = TRUE))
  idxs <- split(seq_along(ry), cut(expr, breaks = br, include.lowest = TRUE, labels = FALSE))
  observed <- cor(rx, ry)
  null <- numeric(nperm)
  for (b in seq_len(nperm)) {
    p <- ry; for (i in idxs) p[i] <- sample(p[i]); null[b] <- cor(rx, p)
  }
  c(observed = observed,
    p_two_sided = (1 + sum(abs(null) >= abs(observed))) / (nperm + 1),
    q025 = unname(quantile(null, 0.025)), q975 = unname(quantile(null, 0.975)))
}

run_test <- function(id, dataset, contrast, prediction, eff) {
  # eff: data.frame(gene, logFC)
  eff <- eff[!is.na(eff$gene) & nzchar(eff$gene) & is.finite(eff$logFC), ]
  eff <- eff[!duplicated(eff$gene), ]
  j <- merge(anchor, eff, by = "gene")
  if (nrow(j) < 1000) {
    return(data.frame(id = id, dataset = dataset, contrast = contrast, prediction = prediction,
                      n_shared = nrow(j), rho = NA, p_matched = NA, null_q025 = NA, null_q975 = NA,
                      rho_proliferation_adjusted = NA, verdict = "NOT EVALUABLE (n_shared < 1000)"))
  }
  pm <- matched_permutation(j$Repro_specific_score, j$logFC, j$mean_log2_expression)
  k <- j[is.finite(j$proliferation_loading), ]
  rho_adj <- if (nrow(k) > 500) {
    cor(residuals(lm(k$Repro_specific_score ~ k$proliferation_loading)),
        residuals(lm(k$logFC ~ k$proliferation_loading)), method = "spearman")
  } else NA
  data.frame(id = id, dataset = dataset, contrast = contrast, prediction = prediction,
             n_shared = nrow(j), rho = unname(pm["observed"]), p_matched = unname(pm["p_two_sided"]),
             null_q025 = unname(pm["q025"]), null_q975 = unname(pm["q975"]),
             rho_proliferation_adjusted = rho_adj, verdict = NA)
}

de_logfc <- function(counts, group, ref) {           # simple two-group edgeR
  y <- DGEList(counts, group = relevel(factor(group), ref = ref))
  y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  y <- estimateDisp(y, model.matrix(~ y$samples$group))
  tt <- topTags(glmQLFTest(glmQLFit(y, model.matrix(~ y$samples$group)), coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), logFC = tt$logFC)
}

RES <- list(); AUX <- list()

## =========================== HO1  GSE297233 ================================
m <- read.csv(gzfile(file.path(ho, "GSE297233_raw_counts_matrix.csv.gz")), row.names = 1, check.names = FALSE)
sy <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", rownames(m)),
                              keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
ok <- !is.na(sy) & nzchar(sy)
.k <- filterByExpr(DGEList(as.matrix(m)), group = ifelse(grepl("_D4_", colnames(m)), "a", "b"))
MAPRATE <- rbind(data.frame(id = "HO1", dataset = "GSE297233",
  rate_all_rows = round(100 * mean(ok), 1), n_all = nrow(m),
  rate_expressed_rows = round(100 * mean(ok[.k]), 1), n_expressed = sum(.k)))
cat("HO1 mapping rate: all rows", round(100 * mean(ok), 1),
    "% / expressed rows", round(100 * mean(ok[.k]), 1), "%\n")
mm <- rowsum(as.matrix(m[ok, ]), sy[ok])
osk <- mm[, grep("^OSK_", colnames(mm))]
RES$HO1 <- run_test("HO1", "GSE297233", "OSK +dox (D4) vs -dox (D0), day 4",
                    "no positive alignment (extrinsic != intrinsic)",
                    de_logfc(osk, ifelse(grepl("_D4_", colnames(osk)), "dox", "nodox"), "nodox"))
mut <- mm[, grep("^O4YRSK_", colnames(mm))]
AUX$HO1_mutant <- run_test("HO1s", "GSE297233", "O4YRSK mutant +dox vs -dox (sensitivity)",
                           "no positive alignment",
                           de_logfc(mut, ifelse(grepl("_D4_", colnames(mut)), "dox", "nodox"), "nodox"))

## =========================== HO2  GSE307377 ================================
g7 <- read.delim(gzfile(file.path(ho, "GSE307377_raw.txt.gz")), check.names = FALSE)
sym7 <- toupper(sub("\\|.*$", "", g7[["Annotation/Divergence"]]))
cnt7 <- as.matrix(g7[, grep("^tags/", names(g7))]); storage.mode(cnt7) <- "numeric"
colnames(cnt7) <- sub("^tags/", "", sub(" .*$", "", colnames(cnt7)))
keep7 <- !is.na(sym7) & nzchar(sym7) & sym7 != "NA"
MAPRATE <- rbind(MAPRATE, data.frame(id = "HO2", dataset = "GSE307377",
  rate_all_rows = round(100 * mean(keep7), 1), n_all = length(keep7),
  rate_expressed_rows = round(100 * mean(keep7), 1), n_expressed = sum(keep7)))
cat("HO2 gene-symbol mapping rate:", round(100 * mean(keep7), 1), "%\n")
cnt7 <- round(rowsum(cnt7[keep7, ], sym7[keep7]))
grp7 <- ifelse(grepl("_O_", colnames(cnt7)), "old", "young")
## sex is not mappable from GEO metadata to these columns -> infer from expression
cpm7 <- cpm(DGEList(cnt7), log = TRUE, prior.count = 1)
ygenes <- intersect(c("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY"), rownames(cpm7))
yscore <- colMeans(cpm7[ygenes, , drop = FALSE])
xist <- if ("XIST" %in% rownames(cpm7)) cpm7["XIST", ] else rep(NA, ncol(cpm7))
sex7 <- ifelse(yscore > mean(range(yscore)), "male", "female")
AUX$HO2_sex <- data.frame(sample = colnames(cnt7), group = grp7,
                          Y_gene_mean_logCPM = round(yscore, 2),
                          XIST_logCPM = round(xist, 2), inferred_sex = sex7)
y7 <- DGEList(cnt7); y7 <- y7[filterByExpr(y7, group = grp7), , keep.lib.sizes = FALSE]
y7 <- calcNormFactors(y7)
d7 <- model.matrix(~ factor(grp7, levels = c("young", "old")) + factor(sex7))
y7 <- estimateDisp(y7, d7)
tt7 <- topTags(glmQLFTest(glmQLFit(y7, d7), coef = 2), n = Inf, sort.by = "none")$table
RES$HO2 <- run_test("HO2", "GSE307377", "old vs young dermal fibroblast (sex-adjusted)",
                    "rho < 0 (anchor opposes ageing)",
                    data.frame(gene = rownames(tt7), logFC = tt7$logFC))

## =========================== HO3  GSE116968 ================================
signed_fc_to_log2 <- function(v) { v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v), NA, ifelse(v >= 0, log2(pmax(v, 1e-6)), -log2(pmax(-v, 1e-6)))) }
ho3 <- function(file, label, id) {
  x <- as.data.frame(read_excel(file, sheet = 1))
  nm <- names(x)
  rescue_col <- grep("^Pre-Red.*/UV.*\\.fc$", nm, value = TRUE)[1]   # 1 h and 4 h sheets differ
  injury_col <- grep("^UV.*/Control\\.fc$",   nm, value = TRUE)[1]
  cat(sprintf("  %s columns: rescue='%s'  injury='%s'\n", label, rescue_col, injury_col))
  d <- data.frame(gene = toupper(x$Gene_Symbol),
                  rescue = signed_fc_to_log2(x[[rescue_col]]),
                  injury = signed_fc_to_log2(x[[injury_col]]))
  d <- d[!is.na(d$gene) & nzchar(d$gene) & is.finite(d$rescue) & is.finite(d$injury), ]
  d <- d[!duplicated(d$gene), ]
  d$resid <- residuals(lm(rescue ~ injury, data = d))   # remove the UV injury component
  list(primary = run_test(id, "GSE116968", paste0("Pre-Red vs UV, ", label,
         ", residualised for UV injury"), "rho > 0 (repair-aligned)",
         data.frame(gene = d$gene, logFC = d$resid)),
       raw = run_test(paste0(id, "r"), "GSE116968", paste0("Pre-Red vs UV, ", label,
         ", not residualised"), "reported alongside",
         data.frame(gene = d$gene, logFC = d$rescue)),
       injury = run_test(paste0(id, "i"), "GSE116968", paste0("UV vs Control, ", label,
         " (injury axis)"), "reported alongside",
         data.frame(gene = d$gene, logFC = d$injury)))
}
a <- ho3(file.path(ho, "GSE116968_Processed_data_NHDF_RNAseq_1_post_1h_.xlsx"), "1 h", "HO3a")
b <- ho3(file.path(ho, "GSE116968_Processed_data_NHDF_RNAseq_2_Post_4h_.xlsx"), "4 h", "HO3b")
RES$HO3a <- a$primary; RES$HO3b <- b$primary
AUX$HO3_extra <- rbind(a$raw, a$injury, b$raw, b$injury)

## =========================== HO4  GSE149694 ================================
td <- file.path(out, "GSE149694_files"); dir.create(td, showWarnings = FALSE)
utils::untar(file.path(ho, "GSE149694_RAW.tar"), exdir = td)
map <- do.call(rbind, lapply(strsplit(system(
  paste("zcat", file.path(root, "heldout/metadata/GSE149694_family.soft.gz"),
        "| awk '/^\\^SAMPLE/{g=$3} /cell subtype\\/time point:/{sub(/.*point: /,\"\"); print g\"\\t\"$0}'"),
  intern = TRUE), "\t"), function(z) data.frame(gsm = z[1], cond = z[2])))
fl <- list.files(td, pattern = "\\.txt\\.gz$", full.names = TRUE)
names(fl) <- sub("_.*$", "", basename(fl))
mats <- lapply(fl, function(f) { z <- read.delim(gzfile(f)); setNames(z[[3]], z[[1]]) })
gg <- Reduce(intersect, lapply(mats, names))
cnt4 <- do.call(cbind, lapply(mats, function(v) v[gg]))
sy4 <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", gg),
                               keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
ok4 <- !is.na(sy4) & nzchar(sy4)
.k4 <- filterByExpr(DGEList(cnt4), group = map$cond[match(colnames(cnt4), map$gsm)])
MAPRATE <- rbind(MAPRATE, data.frame(id = "HO4", dataset = "GSE149694",
  rate_all_rows = round(100 * mean(ok4), 1), n_all = length(ok4),
  rate_expressed_rows = round(100 * mean(ok4[.k4]), 1), n_expressed = sum(.k4)))
cat("HO4 mapping rate: all rows", round(100 * mean(ok4), 1),
    "% / expressed rows", round(100 * mean(ok4[.k4]), 1), "%\n")
cnt4 <- rowsum(cnt4[ok4, ], sy4[ok4])
cond <- map$cond[match(colnames(cnt4), map$gsm)]
sub47 <- cond %in% c("Fibroblast-D3", "Fibroblast-D7")
RES$HO4a <- run_test("HO4a", "GSE149694", "Fibroblast day 7 vs day 3",
                     "no positive alignment",
                     de_logfc(cnt4[, sub47], cond[sub47], "Fibroblast-D3"))
sub413 <- cond %in% c("Fibroblast-D3", "Primed-D13")
RES$HO4b <- run_test("HO4b", "GSE149694", "Primed day 13 vs fibroblast day 3",
                     "no positive alignment",
                     de_logfc(cnt4[, sub413], cond[sub413], "Fibroblast-D3"))
for (alt in c("NHSM-D13", "5iLAF-D13", "RSeT-D13", "t2iLGoY-D13")) {
  s <- cond %in% c("Fibroblast-D3", alt)
  if (sum(cond == alt) >= 2) AUX[[paste0("HO4_", alt)]] <-
    run_test(paste0("HO4b:", alt), "GSE149694", paste(alt, "vs fibroblast day 3 (sensitivity)"),
             "no positive alignment", de_logfc(cnt4[, s], cond[s], "Fibroblast-D3"))
}

## =========================== verdicts ======================================
P <- do.call(rbind, RES[c("HO1", "HO2", "HO3a", "HO3b", "HO4a", "HO4b")])
P$BH_FDR <- p.adjust(P$p_matched, "BH")
P$verdict <- mapply(function(id, r, q) {
  if (is.na(r)) return("NOT EVALUABLE")
  if (id == "HO2")  return(if (r < 0 && q < 0.05) "SUCCESS" else if (r > 0 && q < 0.05) "REFUTES MODEL" else "NOT SUPPORTED")
  if (grepl("^HO3", id)) return(if (r > 0 && q < 0.05) "SUCCESS" else if (r < 0 && q < 0.05) "REFUTES MODEL" else "NOT SUPPORTED")
  if (r > 0.15 && q < 0.05) return("REFUTES BOUNDARY") else return("BOUNDARY HELD")
}, P$id, P$rho, P$BH_FDR)
w(P, "heldout_primary_tests.tsv")
w(do.call(rbind, AUX[grep("HO2_sex", names(AUX), invert = TRUE)]), "heldout_sensitivity_tests.tsv")
w(AUX$HO2_sex, "GSE307377_inferred_sex.tsv")
w(MAPRATE, "heldout_gene_mapping_rates.tsv")

overall <- if (any(P$verdict == "REFUTES MODEL")) "REFUTED" else {
  s <- sum(P$verdict[P$id %in% c("HO2", "HO3a", "HO3b")] == "SUCCESS")
  boundary_ok <- !any(P$verdict == "REFUTES BOUNDARY")
  if (P$verdict[P$id == "HO2"] == "SUCCESS" &&
      any(P$verdict[P$id %in% c("HO3a", "HO3b")] == "SUCCESS") && boundary_ok) "CONFIRMED" else
  if (s >= 1) "PARTIAL" else "NOT SUPPORTED" }
w(data.frame(overall_verdict = overall,
             n_primary = nrow(P), n_success = sum(P$verdict == "SUCCESS"),
             preregistration_sha256 = "41e4cf4b548f74904db7ee8d6fefcf9e18cd801972bcdc1cec95de59a8788bd7"),
  "heldout_overall_verdict.tsv")
cat("\n================ PRIMARY TESTS ================\n")
print(P[, c("id", "dataset", "rho", "rho_proliferation_adjusted", "p_matched", "BH_FDR", "n_shared", "verdict")],
      row.names = FALSE, digits = 3)
cat("\nOVERALL:", overall, "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
