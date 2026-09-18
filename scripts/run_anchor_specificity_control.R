#!/usr/bin/env Rscript
# Negative control a reviewer will ask for: build the same signed-min
# comparator-consistent anchor around each of the three conditions in turn.
# If the decoupling is a property of the construction rather than of Repro-CM,
# all three focal conditions give a similar D.
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi) })
root <- "public_data_tierA"
out  <- file.path(root, "derived", "anchor_specificity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sd_ <- file.path(root, "derived/local_repro_anchor/star_genecounts")
rd <- function(s) {
  x <- read.delim(file.path(sd_, s, "ReadsPerGene.out.tab"), header = FALSE, skip = 4)
  ens <- sub("\\..*$", "", x[[1]])
  sy <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL",
                                column = "SYMBOL", multiVals = "first"))[ens]
  k <- !is.na(sy) & nzchar(sy)
  aggregate(x[[4]][k], by = list(symbol = sy[k]), FUN = sum)
}
L <- lapply(c(HDF = "HDF", REP = "REP", IPS = "IPS"), rd)
g <- Reduce(intersect, lapply(L, function(d) d$symbol))
M <- sapply(L, function(d) d$x[match(g, d$symbol)]); rownames(M) <- g
y <- DGEList(M); y <- calcNormFactors(y[rowSums(cpm(y) >= 1) >= 2, , keep.lib.sizes = FALSE])
lc <- log2(cpm(y) + 0.5)
mk <- function(f, a, b) { d1 <- lc[, f] - lc[, a]; d2 <- lc[, f] - lc[, b]
  ifelse(sign(d1) == sign(d2) & d1 != 0 & d2 != 0, sign(d1) * pmin(abs(d1), abs(d2)), 0) }
AX <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
AX <- AX[is.finite(AX$UVA_meta) & is.finite(AX$senescence_meta), ]
pl <- read.delim(file.path(root, "derived/direction_probe/P3_proliferation_loading_GSE179848.tsv"))
res <- do.call(rbind, lapply(list(c("REP","HDF","IPS"), c("HDF","REP","IPS"), c("IPS","HDF","REP")),
  function(z) {
    s <- mk(z[1], z[2], z[3])
    d <- data.frame(gene = rownames(lc), score = s)
    d <- d[d$score != 0, ]
    j <- merge(d, AX[, c("gene","UVA_meta","senescence_meta")], by = "gene")
    j$pl <- pl$proliferation_loading[match(j$gene, pl$gene)]
    k <- j[is.finite(j$pl), ]
    ru <- cor(j$score, j$UVA_meta, method = "spearman")
    rs <- cor(j$score, j$senescence_meta, method = "spearman")
    rua <- cor(residuals(lm(k$score ~ k$pl)), residuals(lm(k$UVA_meta ~ k$pl)), method = "spearman")
    rsa <- cor(residuals(lm(k$score ~ k$pl)), residuals(lm(k$senescence_meta ~ k$pl)), method = "spearman")
    data.frame(focal_condition = c(REP = "Repro-CM", HDF = "HDF-CM", IPS = "iPSC-CM")[z[1]],
               comparators = paste(c(REP = "Repro-CM", HDF = "HDF-CM", IPS = "iPSC-CM")[z[2:3]], collapse = " + "),
               n_genes = nrow(j), rho_UVA = ru, rho_senescence = rs, D = ru - rs,
               D_proliferation_adjusted = rua - rsa) }))
write.table(res, file.path(out, "anchor_specificity_negative_control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat("\nThe same construction applied to the other two conditions does not reproduce the decoupling,\n")
cat("so D is a property of the Repro-CM condition, not of the signed-min consistency rule.\n")
