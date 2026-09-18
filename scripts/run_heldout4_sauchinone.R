#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Fourth held-out family - a single prespecified test.
#   repro_cm_preregistration_heldout4_sauchinone_ko.md (SHA256 c542d32f...)
# HO8 : GSE326951, HFF dermal fibroblasts, UVB + sauchinone vs UVB.
#       No untreated arm exists, so the contrast CANNOT be residualised on
#       injury; the non-residualised value is used, as fixed in advance.
# One test in this family, so no multiplicity correction is needed.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(edgeR))
root <- "public_data_tierA"; h3 <- file.path(root, "heldout3")
out  <- file.path(root, "derived", "heldout4_sauchinone")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)

AX <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
anchor <- AX[AX$comparator_consistent & AX$Repro_specific_score != 0 &
             is.finite(AX$Repro_specific_score), ]
AXf <- AX[is.finite(AX$UVA_meta) & is.finite(AX$senescence_meta), ]

matched_perm <- function(x, y, expr, nperm = 10000L, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average"); ry <- rank(y, ties.method = "average")
  br <- unique(quantile(expr, seq(0, 1, .1), na.rm = TRUE))
  idx <- split(seq_along(ry), cut(expr, br, include.lowest = TRUE, labels = FALSE))
  obs <- cor(rx, ry); null <- numeric(nperm)
  for (b in seq_len(nperm)) { p <- ry; for (i in idx) p[i] <- sample(p[i]); null[b] <- cor(rx, p) }
  c(rho = obs, p = (1 + sum(abs(null) >= abs(obs))) / (nperm + 1),
    q025 = unname(quantile(null, .025)), q975 = unname(quantile(null, .975)))
}

## ---- load ------------------------------------------------------------------
x <- read.csv(gzfile(file.path(h3, "GSE326951_expression_matrix.csv.gz")), check.names = FALSE)
names(x)[1] <- sub("^\\ufeff", "", names(x)[1])
cnt_cols <- grep("_count$", names(x), value = TRUE)
UVB  <- grep("^2-", cnt_cols, value = TRUE)      # soft: expression_matrix.csv column 2-x = UVB
SAU  <- grep("^3-", cnt_cols, value = TRUE)      # 3-x = UVB + sauchinone
cat("UVB columns      :", paste(UVB, collapse = ", "), "\n")
cat("sauchinone columns:", paste(SAU, collapse = ", "), "\n")
stopifnot(length(UVB) == 3, length(SAU) == 3)
M <- as.matrix(x[, c(UVB, SAU)]); storage.mode(M) <- "numeric"
sym <- toupper(trimws(x$GeneSymbol))
ok <- !is.na(sym) & nzchar(sym) & sym != "NA" & sym != "-"
M <- round(rowsum(M[ok, ], sym[ok]))
grp <- factor(c(rep("UVB", 3), rep("sauchinone", 3)), levels = c("UVB", "sauchinone"))
kexp <- filterByExpr(DGEList(M), group = grp)
maprate <- data.frame(dataset = "GSE326951", n_all = length(ok),
                      rate_all_rows = round(100 * mean(ok), 1),
                      n_expressed = sum(kexp),
                      rate_expressed_rows = round(100 * mean(ok[match(rownames(M), sym)][kexp]), 1))
maprate$rate_expressed_rows <- 100   # symbols are supplied directly; every retained row has one
w(maprate, "HO8_gene_mapping_rates.tsv"); print(maprate, row.names = FALSE)
cat("library sizes (M):", paste(round(colSums(M) / 1e6, 2), collapse = " "), "\n")

y <- DGEList(M, group = grp); y <- y[kexp, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
d <- model.matrix(~ grp); y <- estimateDisp(y, d)
tt <- topTags(glmQLFTest(glmQLFit(y, d), coef = 2), n = Inf, sort.by = "none")$table
eff <- data.frame(gene = rownames(tt), logFC = tt$logFC)
cat("genes tested:", nrow(eff), " | FDR<0.05:", sum(tt$FDR < 0.05), "\n")
w(data.frame(dataset = "GSE326951", contrast = "UVB + sauchinone vs UVB (rescue)",
             gene = eff$gene, logFC = eff$logFC),
  "HO8_signature.tsv")

## ---- PRIMARY ---------------------------------------------------------------
j <- merge(anchor, eff, by = "gene"); j <- j[is.finite(j$logFC), ]
stopifnot(nrow(j) >= 1000)
pm <- matched_perm(j$Repro_specific_score, j$logFC, j$mean_log2_expression)
P <- data.frame(id = "HO8", dataset = "GSE326951",
                contrast = "UVB + sauchinone vs UVB (non-residualised, no control arm exists)",
                n_shared = nrow(j), rho = unname(pm["rho"]), p_matched = unname(pm["p"]),
                null_q025 = unname(pm["q025"]), null_q975 = unname(pm["q975"]))
P$verdict <- if (P$rho > 0 && P$p_matched < 0.05) "SUCCESS" else
             if (P$rho < 0 && P$p_matched < 0.05) "REFUTES" else "NOT SUPPORTED"
w(P, "HO8_primary_test.tsv")
cat("\n================ HO8 (family of one, no correction needed) ================\n")
print(P[, c("id", "n_shared", "rho", "p_matched", "null_q025", "null_q975", "verdict")],
      row.names = FALSE, digits = 3)

## ---- SENSITIVITY -----------------------------------------------------------
k <- merge(AXf[, c("gene", "UVA_meta", "senescence_meta")], eff, by = "gene")
S <- data.frame(n_genes = nrow(k),
                rho_UVA = cor(k$logFC, k$UVA_meta, method = "spearman"),
                rho_senescence = cor(k$logFC, k$senescence_meta, method = "spearman"))
S$D <- S$rho_UVA - S$rho_senescence
w(S, "HO8_axis_loading.tsv")
cat("\naxis loading: rho_UVA =", round(S$rho_UVA, 3),
    "| rho_senescence =", round(S$rho_senescence, 3), "| D =", round(S$D, 3), "\n")

## agreement with the authors' own differential results (mislabelled .xlsx, actually TSV)
pub <- tryCatch(read.delim(file.path(h3, "GSE326951_DEG_results.xlsx"), check.names = FALSE),
                error = function(e) NULL)
if (!is.null(pub)) {
  lc <- grep("log2|logFC", names(pub), ignore.case = TRUE, value = TRUE)[1]
  gc_ <- grep("GeneSymbol|gene_name|symbol", names(pub), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(lc) && !is.na(gc_)) {
    pp <- data.frame(gene = toupper(pub[[gc_]]), pub_logFC = suppressWarnings(as.numeric(pub[[lc]])))
    pp <- pp[!duplicated(pp$gene) & is.finite(pp$pub_logFC), ]
    m <- merge(eff, pp, by = "gene")
    cat("agreement with the deposited DEG table:", nrow(m), "genes, Spearman rho =",
        round(cor(m$logFC, m$pub_logFC, method = "spearman"), 3),
        "(sign flips if the authors used the opposite reference)\n")
    w(data.frame(n = nrow(m), rho = cor(m$logFC, m$pub_logFC, method = "spearman"),
                 column_used = lc), "HO8_agreement_with_deposited_DEG.tsv")
  }
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
