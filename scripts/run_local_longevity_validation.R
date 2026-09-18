#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(BiocParallel)
  library(fgsea)
  library(readxl)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "local_repro_anchor")

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

signed_z <- function(effect, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(effect) * qnorm(pvalue / 2, lower.tail = FALSE)
}

anchor <- read.delim(file.path(out, "local_Repro_specific_gene_rank.tsv"), check.names = FALSE)
anchor <- anchor[anchor$comparator_consistent & anchor$Repro_specific_score != 0, ]

reference_file <- file.path(root, "longevity", "Tyshkovskiy_2026_Table_S2_signatures.xlsx")
reference_sheets <- c(
  rodent_chronological_age = "(H) Rodents Chronological age",
  rodent_mortality = "(K) Rodents Mortality rate",
  rodent_mortality_age_adjusted = "(L) Rodents Mortality rate, adj",
  rodent_max_lifespan = "(M) Rodents Max lifespan",
  rodent_max_lifespan_age_adjusted = "(N) Rodents Max lifespan, adj",
  human_multitissue_age = "(R) Human aging multi-tissue"
)

matched_permutation <- function(x, y, expression_value, nperm = 5000, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  breaks <- unique(quantile(expression_value, probs = seq(0, 1, 0.1), na.rm = TRUE))
  bins <- cut(expression_value, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  observed <- cor(rx, ry)
  indices <- split(seq_along(ry), bins)
  null <- numeric(nperm)
  for (b in seq_len(nperm)) {
    permuted <- ry
    for (idx in indices) permuted[idx] <- sample(permuted[idx])
    null[b] <- cor(rx, permuted)
  }
  c(rho = observed,
    empirical_p = (1 + sum(abs(null) >= abs(observed))) / (nperm + 1),
    null_q025 = unname(quantile(null, 0.025)), null_q975 = unname(quantile(null, 0.975)))
}

rows <- lapply(names(reference_sheets), function(name) {
  x <- as.data.frame(read_excel(reference_file, sheet = reference_sheets[[name]]))
  reference <- data.frame(
    gene = toupper(x$Gene.symbol), slope = as.numeric(x$Slope),
    z = signed_z(as.numeric(x$Slope), as.numeric(x$P.Value))
  )
  joined <- merge(anchor[, c("gene", "Repro_specific_score", "mean_log2_expression")], reference, by = "gene")
  test <- matched_permutation(joined$Repro_specific_score, joined$z, joined$mean_log2_expression)
  expected <- if (grepl("max_lifespan", name)) "same" else "opposite"
  data.frame(
    reference = name, expected_relation = expected, n_shared = nrow(joined),
    rho_signed_z = test["rho"], rho_slope = cor(joined$Repro_specific_score, joined$slope, method = "spearman"),
    empirical_p_expression_matched = test["empirical_p"],
    null_q025 = test["null_q025"], null_q975 = test["null_q975"],
    direction_matches_hypothesis = ifelse(expected == "same", test["rho"] > 0, test["rho"] < 0)
  )
})
correlations <- do.call(rbind, rows)
correlations$FDR_expression_matched <- p.adjust(correlations$empirical_p_expression_matched, "BH")
write_tsv(correlations, "local_anchor_vs_Tyshkovskiy2026.tsv")

gene_file <- file.path(root, "longevity", "Tyshkovskiy_2019_Table_S6_genes.xlsx")
common <- as.data.frame(read_excel(gene_file, sheet = "A (Significant Common Genes)"))
associated <- as.data.frame(read_excel(gene_file, sheet = "B (Common and Lifespan Genes)"))
sets <- list(
  longevity_intervention_up = unique(toupper(common$Gene[common$`Sign of change` == "Upregulated"])),
  longevity_intervention_down = unique(toupper(common$Gene[common$`Sign of change` == "Downregulated"])),
  lifespan_positive = unique(toupper(associated$Gene[associated$`Sign of association` == "Positive"])),
  lifespan_negative = unique(toupper(associated$Gene[associated$`Sign of association` == "Negative"]))
)
ranks <- anchor$Repro_specific_score
names(ranks) <- anchor$gene
ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
set.seed(260911)
gsea <- fgseaMultilevel(sets, ranks, minSize = 5, maxSize = 500, eps = 0, BPPARAM = SerialParam())
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write_tsv(as.data.frame(gsea), "local_anchor_vs_Tyshkovskiy2019_gene_sets.tsv")

cat("Local longevity validation completed\n")
