#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Held-out senescence test (HO5), executed exactly as fixed in
#   repro_cm_preregistration_heldout2_senescence_ko.md
# BH family = {HO5a, HO5b} only. The first held-out family (6 tests) is closed
# and is NOT re-corrected.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db); library(edgeR)
})
root <- "public_data_tierA"; h2 <- file.path(root, "heldout2")
out  <- file.path(root, "derived", "heldout2_senescence")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")

## ---- fixed local anchor + fixed proliferation loading ---------------------
anchor <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                     check.names = FALSE)
anchor <- anchor[anchor$comparator_consistent & anchor$Repro_specific_score != 0 &
                 is.finite(anchor$Repro_specific_score), ]
pl <- read.delim(file.path(root, "derived/direction_probe/P3_proliferation_loading_GSE179848.tsv"))
anchor$proliferation_loading <- pl$proliferation_loading[match(anchor$gene, pl$gene)]
cat("anchor genes:", nrow(anchor), "\n")

matched_permutation <- function(x, y, expr, nperm = 10000L, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average"); ry <- rank(y, ties.method = "average")
  br <- unique(quantile(expr, probs = seq(0, 1, 0.1), na.rm = TRUE))
  idxs <- split(seq_along(ry), cut(expr, breaks = br, include.lowest = TRUE, labels = FALSE))
  observed <- cor(rx, ry); null <- numeric(nperm)
  for (b in seq_len(nperm)) { p <- ry; for (i in idxs) p[i] <- sample(p[i]); null[b] <- cor(rx, p) }
  c(observed = observed,
    p_two_sided = (1 + sum(abs(null) >= abs(observed))) / (nperm + 1),
    q025 = unname(quantile(null, 0.025)), q975 = unname(quantile(null, 0.975)))
}

## ---- build the MRC5 count matrix ------------------------------------------
meta <- local({
  L <- readLines(gzfile(file.path(h2, "GSE306957_family.soft.gz")))
  gsm <- sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", L, value = TRUE))
  blocks <- split(L, cumsum(grepl("^\\^SAMPLE", L)))[-1]
  do.call(rbind, lapply(blocks, function(b) {
    g1 <- function(pat) { z <- grep(pat, b, value = TRUE)[1]
                          if (is.na(z)) NA else sub(".*= ", "", z) }
    data.frame(gsm = sub("^\\^SAMPLE = ", "", b[1]),
               title = g1("!Sample_title"),
               cell_line = g1("cell line:"),
               genotype = g1("genotype:"),
               treatment = g1("treatment:"),
               file = sub(".*/", "", g1("!Sample_supplementary_file_1")))
  }))
})
meta <- meta[!is.na(meta$cell_line) & grepl("MRC5", meta$cell_line), ]      # mouse liver excluded
meta$state <- ifelse(grepl("^Proliferating", meta$title), "proliferating", "senescent")
meta$arm <- ifelse(grepl("siScramble", meta$title), "siScramble",
            ifelse(grepl("siMAVS", meta$title), "siMAVS",
            ifelse(grepl("pLenti", meta$title), "pLenti",
            ifelse(grepl("IFIH1", meta$title), "IFIH1_KO", "DDX58_KO"))))
meta$batch <- ifelse(meta$arm %in% c("siScramble", "siMAVS"), "siRNA", "CRISPR")
cat("MRC5 samples:", nrow(meta), "\n"); print(table(meta$arm, meta$state))

read_counts <- function(f) { x <- read.delim(gzfile(file.path(h2, "files", f)), check.names = FALSE)
                             setNames(x[[ncol(x)]], x$Geneid) }
mats <- lapply(meta$file, read_counts)
gg <- Reduce(intersect, lapply(mats, names))
M <- do.call(cbind, lapply(mats, function(v) v[gg])); colnames(M) <- meta$gsm
sy <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", gg),
                              keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
ok <- !is.na(sy) & nzchar(sy)
kexp <- filterByExpr(DGEList(M), group = meta$state)
maprate <- data.frame(dataset = "GSE306957", n_all = length(ok),
                      rate_all_rows = round(100 * mean(ok), 1),
                      n_expressed = sum(kexp),
                      rate_expressed_rows = round(100 * mean(ok[kexp]), 1))
w(maprate, "HO5_gene_mapping_rates.tsv"); print(maprate, row.names = FALSE)
stopifnot(maprate$rate_expressed_rows >= 70)     # prespecified exclusion rule
Msym <- rowsum(M[ok, ], sy[ok])

de <- function(keep, design_batch) {
  m <- meta[keep, ]; y <- DGEList(Msym[, keep])
  y <- y[filterByExpr(y, group = m$state), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  st <- factor(m$state, levels = c("proliferating", "senescent"))
  d <- if (design_batch && length(unique(m$batch)) > 1) model.matrix(~ st + factor(m$batch)) else model.matrix(~ st)
  y <- estimateDisp(y, d)
  tt <- topTags(glmQLFTest(glmQLFit(y, d), coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), logFC = tt$logFC)
}

test_one <- function(id, label, eff, adjust_prolif) {
  j <- merge(anchor, eff, by = "gene")
  j <- j[is.finite(j$logFC), ]
  if (adjust_prolif) {
    j <- j[is.finite(j$proliferation_loading), ]
    x <- residuals(lm(j$Repro_specific_score ~ j$proliferation_loading))
    y <- residuals(lm(j$logFC ~ j$proliferation_loading))
  } else { x <- j$Repro_specific_score; y <- j$logFC }
  if (nrow(j) < 1000) return(data.frame(id = id, contrast = label, n_shared = nrow(j),
                                        rho = NA, p_matched = NA, null_q025 = NA, null_q975 = NA))
  pm <- matched_permutation(x, y, j$mean_log2_expression)
  data.frame(id = id, contrast = label, n_shared = nrow(j), rho = unname(pm["observed"]),
             p_matched = unname(pm["p_two_sided"]),
             null_q025 = unname(pm["q025"]), null_q975 = unname(pm["q975"]))
}

## ---- PRIMARY: pooled non-targeting controls, batch-adjusted ---------------
ctrl <- meta$arm %in% c("pLenti", "siScramble")
cat("\nprimary control set: proliferating", sum(ctrl & meta$state == "proliferating"),
    "| senescent", sum(ctrl & meta$state == "senescent"), "\n")
eff_primary <- de(ctrl, TRUE)
P <- rbind(
  test_one("HO5a", "senescent vs proliferating control (pooled, batch-adjusted)", eff_primary, FALSE),
  test_one("HO5b", "same contrast, both axes residualised on proliferation loading", eff_primary, TRUE))
P$BH_FDR <- p.adjust(P$p_matched, "BH")
P$verdict <- ifelse(is.na(P$rho), "NOT EVALUABLE",
             ifelse(P$rho < 0 & P$BH_FDR < 0.05, "SUCCESS",
             ifelse(P$rho > 0 & P$BH_FDR < 0.05, "REFUTES MODEL", "NOT SUPPORTED")))
w(P, "HO5_primary_tests.tsv")

## ---- SENSITIVITY ----------------------------------------------------------
S <- list()
for (a in c("pLenti", "siScramble", "IFIH1_KO", "DDX58_KO")) {
  k <- meta$arm == a
  if (sum(k & meta$state == "senescent") >= 2 && sum(k & meta$state == "proliferating") >= 2) {
    e <- de(k, FALSE)
    S[[a]] <- rbind(test_one(paste0("HO5s:", a), paste(a, "senescent vs proliferating"), e, FALSE),
                    test_one(paste0("HO5s:", a, ":prolifAdj"), paste(a, "senescent vs proliferating, proliferation-adjusted"), e, TRUE))
  }
}
SEN <- do.call(rbind, S)
w(SEN, "HO5_sensitivity_tests.tsv")

cat("\n================ HO5 PRIMARY (BH family of 2) ================\n")
print(P[, c("id", "n_shared", "rho", "p_matched", "BH_FDR", "verdict")], row.names = FALSE, digits = 3)
overall <- if (P$verdict[1] == "REFUTES MODEL") "REFUTED" else
           if (all(P$verdict == "SUCCESS")) "CONFIRMED" else
           if (P$verdict[1] == "SUCCESS") "PARTIAL" else "NOT SUPPORTED"
cat("\nOVERALL:", overall, "\n")
cat("\n================ SENSITIVITY ================\n")
print(SEN[, c("id", "n_shared", "rho", "p_matched")], row.names = FALSE, digits = 3)
w(data.frame(overall_verdict = overall,
             preregistration_sha256 = "6e5a7703cc7a5ac6c695b4bd10ef53e872c709371e18ef9135cfde33f7be5ff4",
             bh_family = "HO5a + HO5b only; the first held-out family of 6 is closed and not re-corrected"),
  "HO5_overall_verdict.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
