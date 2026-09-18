#!/usr/bin/env Rscript
# Add the third held-out family's contrasts to the perturbation compendium.
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(edgeR); library(readxl) })
root <- "public_data_tierA"; h3 <- file.path(root, "heldout3")
sig  <- file.path(root, "derived", "compendium", "signatures")
dir.create(sig, recursive = TRUE, showWarnings = FALSE)
save_sig <- function(id, dataset, contrast, gene, logFC) {
  d <- data.frame(dataset = dataset, contrast = contrast, gene = toupper(gene), logFC = logFC)
  d <- d[!is.na(d$gene) & nzchar(d$gene) & is.finite(d$logFC), ]
  d <- d[!duplicated(d$gene), ]
  write.table(d, file.path(sig, paste0(id, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  %-52s %d genes\n", contrast, nrow(d)))
}
de <- function(cnt, grp, ref) {
  y <- DGEList(cnt); g <- factor(grp, levels = c(ref, setdiff(unique(grp), ref)))
  y <- y[filterByExpr(y, group = g), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  d <- model.matrix(~ g); y <- estimateDisp(y, d)
  tt <- topTags(glmQLFTest(glmQLFit(y, d), coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), logFC = tt$logFC)
}
## GSE240486
g1 <- read.delim(gzfile(file.path(h3, "GSE240486_all_sample_count.txt.gz")), check.names = FALSE)
m1 <- as.matrix(g1[, -(1:2)]); storage.mode(m1) <- "numeric"
s1 <- toupper(g1$SYMBOL); k1 <- !is.na(s1) & nzchar(s1); m1 <- rowsum(m1[k1, ], s1[k1])
gp1 <- sub("[0-9]+$", "", colnames(m1))
pick1 <- function(a, b) { k <- gp1 %in% c(a, b); de(m1[, k], gp1[k], a) }
save_sig("GSE240486_UVinjury", "GSE240486", "UV vs control (injury)",
         pick1("Control", "UV_Model")$gene, pick1("Control", "UV_Model")$logFC)
for (z in list(c("UV_Osmoter_L", "UV + 0.02% Osmoter vs UV (rescue)"),
               c("UV_Osmoter_H", "UV + 0.06% Osmoter vs UV (rescue)"))) {
  e <- pick1("UV_Model", z[1]); save_sig(paste0("GSE240486_", z[1]), "GSE240486", z[2], e$gene, e$logFC) }
e <- pick1("Control", "noUV_Osmoter_H")
save_sig("GSE240486_OsmoterAlone", "GSE240486", "0.06% Osmoter alone, no UV", e$gene, e$logFC)

## GSE222414
g2 <- read.delim(gzfile(file.path(h3, "GSE222414_counts.txt.gz")), check.names = FALSE)
m2 <- as.matrix(g2[, -1]); storage.mode(m2) <- "numeric"
s2 <- toupper(g2[[1]]); k2 <- !is.na(s2) & nzchar(s2); m2 <- rowsum(m2[k2, ], s2[k2])
lb2 <- c(A = "control", B = "UVB", C = "UVB_UIFSP5", D = "UVB_UIFSP6")[sub("-[0-9]+$", "", colnames(m2))]
pick2 <- function(a, b) { k <- lb2 %in% c(a, b); de(m2[, k], lb2[k], a) }
e <- pick2("control", "UVB"); save_sig("GSE222414_UVBinjury", "GSE222414", "UVB vs control (injury)", e$gene, e$logFC)
for (z in list(c("UVB_UIFSP5", "UVB + TAT-UIFSP5 vs UVB (rescue)"),
               c("UVB_UIFSP6", "UVB + TAT-UIFSP6 vs UVB (rescue)"))) {
  e <- pick2("UVB", z[1]); save_sig(paste0("GSE222414_", z[1]), "GSE222414", z[2], e$gene, e$logFC) }

## GSE329475
g3 <- as.data.frame(read_excel(file.path(h3, "GSE329475_all_gene_count_matrix.xlsx")))
m3 <- as.matrix(g3[, -1]); storage.mode(m3) <- "numeric"
s3 <- toupper(g3$gene_name); k3 <- !is.na(s3) & nzchar(s3); m3 <- rowsum(m3[k3, ], s3[k3])
gp3 <- ifelse(grepl("UVA", colnames(m3)), "UVA", "NC")
e <- de(m3, gp3, "NC")
save_sig("GSE329475_UVA3d", "GSE329475", "UVA, 3 consecutive days vs NC (injury)", e$gene, e$logFC)
cat("\nsignatures now in compendium:", length(list.files(sig)), "\n")
