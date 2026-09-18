#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(BiocParallel)
  library(fgsea)
  library(readxl)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "figure5_longevity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

read_effect <- function(folder, filename) {
  x <- read.delim(file.path(root, "derived", folder, filename), check.names = FALSE)
  x$gene <- toupper(x$gene)
  x[, c("gene", "logFC", "PValue", "FDR")]
}

targets <- list(
  GSE125429_chronic_UVA = read_effect("figure1_public_uva", "GSE125429_UVA_vs_control_DE.tsv"),
  GSE240226_acute_UVA = read_effect("figure1_public_uva", "GSE240226_UVA_vs_control_DE.tsv"),
  GSE302943_cumulative_UVA = read_effect("figure1_public_uva", "GSE302943_UVA_vs_control_DE.tsv"),
  GSE113957_adult_age = read_effect("figure2_public_aging", "GSE113957_age22_89_per_decade_DE.tsv"),
  GSE226189_age = read_effect("figure2_public_aging", "GSE226189_age_per_decade_DE.tsv"),
  GSE179848_late_replicative = read_effect("figure2_public_aging", "GSE179848_late_vs_early_donor_level_DE.tsv"),
  GSE109700_early_senescence = read_effect("figure2_public_aging", "GSE109700_early_vs_proliferating_DE.tsv"),
  GSE109700_deep_senescence = read_effect("figure2_public_aging", "GSE109700_deep_vs_proliferating_DE.tsv"),
  GSE191055_P27 = read_effect("figure2_public_aging", "GSE191055_P27_vs_P4_DE.tsv"),
  GSE93535_SIPS = read_effect("figure2_public_aging", "GSE93535_SIPS_vs_Q_DE.tsv"),
  GSE165177_MPTR_all_days = read_effect("figure2_mptr", "GSE165177_MPTR_paired_DE.tsv")
)

reference_file <- file.path(root, "longevity", "Tyshkovskiy_2026_Table_S2_signatures.xlsx")
reference_sheets <- c(
  rodent_chronological_age = "(H) Rodents Chronological age",
  rodent_mortality = "(K) Rodents Mortality rate",
  rodent_mortality_age_adjusted = "(L) Rodents Mortality rate, adj",
  rodent_max_lifespan = "(M) Rodents Max lifespan",
  rodent_max_lifespan_age_adjusted = "(N) Rodents Max lifespan, adj",
  human_multitissue_age = "(R) Human aging multi-tissue"
)
references <- lapply(reference_sheets, function(sheet) {
  x <- as.data.frame(read_excel(reference_file, sheet = sheet))
  data.frame(
    gene = toupper(x$Gene.symbol), slope = as.numeric(x$Slope),
    PValue = as.numeric(x$P.Value), FDR = as.numeric(x$P.Adjusted)
  )
})

compare_effects <- function(target, reference, target_name, reference_name) {
  joined <- merge(target[, c("gene", "logFC", "PValue")],
                  reference[, c("gene", "slope", "PValue")], by = "gene",
                  suffixes = c("_target", "_reference"))
  x <- signed_z(joined$logFC, joined$PValue_target)
  y <- signed_z(joined$slope, joined$PValue_reference)
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  data.frame(
    target = target_name, reference = reference_name, n_shared = nrow(joined),
    rho_signed_z = unname(test$estimate), p_value = test$p.value,
    same_direction_fraction = mean(sign(joined$logFC) == sign(joined$slope))
  )
}

correlations <- do.call(rbind, lapply(names(targets), function(target_name) {
  do.call(rbind, lapply(names(references), function(reference_name) {
    compare_effects(targets[[target_name]], references[[reference_name]], target_name, reference_name)
  }))
}))
write_tsv(correlations, "public_effects_vs_Tyshkovskiy2026_signatures.tsv")

gene_file <- file.path(root, "longevity", "Tyshkovskiy_2019_Table_S6_genes.xlsx")
common <- as.data.frame(read_excel(gene_file, sheet = "A (Significant Common Genes)"))
associated <- as.data.frame(read_excel(gene_file, sheet = "B (Common and Lifespan Genes)"))
gene_sets <- list(
  longevity_intervention_up = unique(toupper(common$Gene[common$`Sign of change` == "Upregulated"])),
  longevity_intervention_down = unique(toupper(common$Gene[common$`Sign of change` == "Downregulated"])),
  lifespan_positive = unique(toupper(associated$Gene[associated$`Sign of association` == "Positive"])),
  lifespan_negative = unique(toupper(associated$Gene[associated$`Sign of association` == "Negative"]))
)

run_gene_set_test <- function(target, target_name) {
  ranks <- signed_z(target$logFC, target$PValue)
  names(ranks) <- target$gene
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  set.seed(260911)
  result <- fgseaMultilevel(gene_sets, ranks, minSize = 5, maxSize = 500, eps = 0,
                           BPPARAM = SerialParam())
  result$target <- target_name
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  as.data.frame(result)[, c("target", "pathway", "size", "NES", "pval", "padj", "leadingEdge")]
}

gene_set_results <- do.call(rbind, lapply(names(targets), function(name) run_gene_set_test(targets[[name]], name)))
write_tsv(gene_set_results, "public_effects_vs_Tyshkovskiy2019_gene_sets.tsv")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Longevity convergence audit completed:", out, "\n")
