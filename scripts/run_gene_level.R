## Gene level: within the 96-gene module, is the loss of the age association
## under adjustment predicted by each gene's own coupling to proliferation?
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
D <- "public_data_tierA/derived"; O <- file.path(D, "outcome_vs_expression")
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
M <- fread(file.path(O, "GSE113957_sample_metrics.tsv"))

## rebuild the per-gene expression matrix for the same donors
suppressPackageStartupMessages(library(edgeR))
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
GM <- GM[, M$srr, drop = FALSE]
y <- DGEList(GM); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
g <- intersect(rownames(lc), CORE); say("core genes measured: %d of %d", length(g), length(CORE))

res <- do.call(rbind, lapply(g, function(gn) {
  v  <- lc[gn, ]
  f0 <- lm(v ~ I(M$age/10) + M$log_depth)
  f1 <- lm(v ~ I(M$age/10) + M$prolif + M$log_depth)
  b0 <- summary(f0)$coef[2, ]; b1 <- summary(f1)$coef[2, ]
  data.frame(gene = gn, r_prolif = cor(v, M$prolif),
             beta = b0[1], p = b0[4], beta_adj = b1[1], p_adj = b1[4],
             retained = b1[1]/b0[1]) }))
res$FDR <- p.adjust(res$p, "BH"); res$FDR_adj <- p.adjust(res$p_adj, "BH")
res$retained[!is.finite(res$retained)] <- NA
say("\ngenes with an age association at FDR<0.05:  unadjusted %d, adjusted %d",
    sum(res$FDR < 0.05), sum(res$FDR_adj < 0.05))
ct <- cor.test(res$r_prolif, res$retained, method = "spearman", exact = FALSE)
say("rho(gene's coupling to proliferation, fraction of age effect retained) = %+.3f, P = %.2e",
    ct$estimate, ct$p.value)
say("median coupling to proliferation across the 96 genes: %+.3f (range %+.2f to %+.2f)",
    median(res$r_prolif), min(res$r_prolif), max(res$r_prolif))
write.table(res, file.path(O, "gene_level_age_adjustment.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
## expression matrix of the core genes ordered by age, for the heatmap
ordm <- order(M$age)
Z <- t(scale(t(lc[g, ordm])))
write.table(data.frame(gene = rownames(Z), Z, check.names = FALSE),
            file.path(O, "core_gene_heatmap_matrix.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(M[ordm, .(srr, age, prolif, machinery, unannot_reads)],
            file.path(O, "heatmap_sample_order.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("heatmap matrix: %d genes x %d donors", nrow(Z), ncol(Z))
