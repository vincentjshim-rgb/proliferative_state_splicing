#!/usr/bin/env Rscript
# 로컬 4개 시료의 intron retention burden.
# n=1이므로 이것은 검정이 아니라 기술적 기술(descriptive)이며 p-value를 붙이지 않는다.
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"
ir   <- file.path(root, "derived", "intron_retention")
w <- function(x, f) write.table(x, file.path(ir, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
SAMPLES <- c(HDF = "HDF-CM", REP = "Repro-CM", IPS = "iPSC-CM", UVA0 = "UVA 15J, 0 h")

load_counts <- function(s, kind) {
  f <- file.path(ir, sprintf("%s.%s_unique.counts", s, kind))
  x <- read.delim(f, header = FALSE)
  names(x) <- c("chr", "start", "end", "gene", "count")
  data.frame(gene = x$gene, count = x$count, bp = x$end - x$start)
}
agg <- function(s, kind) {
  d <- load_counts(s, kind)
  a <- aggregate(cbind(count, bp) ~ gene, data = d, FUN = sum)
  names(a)[2:3] <- paste0(kind, c("_count", "_bp")); a
}
per_sample <- lapply(names(SAMPLES), function(s) {
  m <- merge(agg(s, "exonic"), agg(s, "intronic"), by = "gene")
  m$sample <- s
  m$exonic_rpkb   <- 1000 * m$exonic_count   / m$exonic_bp
  m$intronic_rpkb <- 1000 * m$intronic_count / m$intronic_bp
  m$IR_ratio <- m$intronic_rpkb / m$exonic_rpkb
  m
})
names(per_sample) <- names(SAMPLES)

## sample-level burden over genes that are adequately covered in ALL samples
common <- Reduce(intersect, lapply(per_sample, function(d) d$gene[d$exonic_count >= 50]))
cat("genes with exonic count >= 50 in all four samples:", length(common), "\n")
burden <- do.call(rbind, lapply(names(SAMPLES), function(s) {
  d <- per_sample[[s]]; d <- d[match(common, d$gene), ]
  tot_ex <- sum(d$exonic_count); tot_in <- sum(d$intronic_count)
  data.frame(
    sample = s, label = SAMPLES[[s]], n_genes = nrow(d),
    exonic_reads = tot_ex, intronic_reads = tot_in,
    intronic_read_fraction = tot_in / (tot_in + tot_ex),
    length_normalised_IR = (1000 * tot_in / sum(d$intronic_bp)) /
                           (1000 * tot_ex / sum(d$exonic_bp)),
    median_gene_IR_ratio = median(d$IR_ratio, na.rm = TRUE),
    q25_gene_IR_ratio = quantile(d$IR_ratio, 0.25, na.rm = TRUE),
    q75_gene_IR_ratio = quantile(d$IR_ratio, 0.75, na.rm = TRUE))
}))
w(burden, "local_IR_burden.tsv")
cat("\n================ SAMPLE-LEVEL IR BURDEN (descriptive, n = 1 per condition) ================\n")
print(burden[, c("label", "n_genes", "intronic_read_fraction",
                 "length_normalised_IR", "median_gene_IR_ratio")], row.names = FALSE, digits = 4)

## per-gene matrix + the conservative Repro-specific contrast used elsewhere
mat <- Reduce(function(a, b) merge(a, b, by = "gene"),
  lapply(names(SAMPLES), function(s) {
    d <- per_sample[[s]][, c("gene", "IR_ratio", "exonic_count")]
    names(d)[2:3] <- paste0(c("IR_", "exon_"), s); d }))
mat <- mat[mat$gene %in% common, ]
lg <- function(v) log2(v + 0.01)
mat$dH <- lg(mat$IR_REP) - lg(mat$IR_HDF)
mat$dI <- lg(mat$IR_REP) - lg(mat$IR_IPS)
mat$conservative_log2FC_IR <- ifelse(sign(mat$dH) == sign(mat$dI),
                                     sign(mat$dH) * pmin(abs(mat$dH), abs(mat$dI)), 0)
w(mat[order(-abs(mat$conservative_log2FC_IR)), ], "local_IR_per_gene.tsv")
cat("\ncomparator-consistent direction:",
    sprintf("%.1f%% of %d genes\n", 100 * mean(sign(mat$dH) == sign(mat$dI) & mat$dH != 0), nrow(mat)))
cat("genes with |conservative log2FC IR| >= log2(1.5):",
    sum(abs(mat$conservative_log2FC_IR) >= log2(1.5)), "\n")

## does the IR shift track the gene-expression anchor? (interpretation guard)
ax <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
j <- merge(mat, ax[ax$comparator_consistent, c("gene", "Repro_specific_score")], by = "gene")
cat("\ncor(IR conservative log2FC, expression anchor score) over", nrow(j), "genes: rho =",
    round(cor(j$conservative_log2FC_IR, j$Repro_specific_score, method = "spearman"), 4), "\n")
cat("\nNOTE: biological replicate = 1 per condition. These are descriptive placements,\n")
cat("      not differential intron-retention calls. No p-values are reported.\n")
writeLines(capture.output(sessionInfo()), file.path(ir, "sessionInfo.txt"))
