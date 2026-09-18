#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(BiocParallel)
  library(fgsea)
  library(limma)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "figure2_mptr")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
meta <- read.delim(file.path(root, "derived", "audit", "all_geo_samples.tsv"), check.names = FALSE)

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_mptr <- function(path) {
  x <- read.delim(gzfile(path), check.names = FALSE)
  list(
    symbol = toupper(trimws(x[[1]])),
    key = paste(x[[1]], x[[7]], x[[2]], x[[3]], x[[4]], sep = "|"),
    expression = as.matrix(x[, 13:ncol(x), drop = FALSE])
  )
}

part1 <- read_mptr(file.path(root, "reprogramming", "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz"))
part2 <- read_mptr(file.path(root, "reprogramming", "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz"))
part2_order <- match(part1$key, part2$key)
stopifnot(!anyNA(part2_order), !anyDuplicated(part1$key), !anyDuplicated(part2$key))
expression <- cbind(part1$expression, part2$expression[part2_order, , drop = FALSE])
storage.mode(expression) <- "numeric"
keep_gene <- !is.na(part1$symbol) & nzchar(part1$symbol) & part1$symbol != "--"
expression <- limma::avereps(expression[keep_gene, , drop = FALSE], ID = part1$symbol[keep_gene])

mptr_meta <- meta[meta$dataset == "GSE165177", ]
normalize_title <- function(x) tolower(gsub("[ _]+", "_", trimws(x)))
mptr_meta <- mptr_meta[match(normalize_title(colnames(expression)), normalize_title(mptr_meta$title)), ]
stopifnot(!anyNA(mptr_meta$accession))

transient_idx <- which(grepl("_transiently_reprogrammed_", mptr_meta$title) &
                         !grepl("intermediate", mptr_meta$title, ignore.case = TRUE))
negative_idx <- which(grepl("_negative_control_", mptr_meta$title) &
                        !grepl("intermediate", mptr_meta$title, ignore.case = TRUE))

pair_key <- function(title) {
  donor <- sub("_.*", "", title)
  day <- sub(".*_(\\d+)days_.*", "\\1", title)
  experiment <- sub(".*_(exp[12])$", "\\1", title)
  paste(donor, day, experiment, sep = "_")
}
transient_keys <- vapply(mptr_meta$title[transient_idx], pair_key, character(1))
negative_keys <- vapply(mptr_meta$title[negative_idx], pair_key, character(1))
matched_negative <- negative_idx[match(transient_keys, negative_keys)]
if (anyNA(matched_negative)) stop("A transiently reprogrammed sample lacks a time-matched negative control")

selected <- as.vector(rbind(matched_negative, transient_idx))
sample_data <- data.frame(
  pair = factor(rep(transient_keys, each = 2), levels = transient_keys),
  condition = factor(rep(c("negative", "transient"), length(transient_keys)),
                     levels = c("negative", "transient")),
  row.names = colnames(expression)[selected]
)
design <- model.matrix(~ pair + condition, data = sample_data)
fit <- eBayes(lmFit(expression[, selected, drop = FALSE], design), trend = TRUE, robust = TRUE)
tab <- topTable(fit, coef = "conditiontransient", number = Inf, sort.by = "none")
mptr <- data.frame(
  dataset = "GSE165177", contrast = "successful_MPTR_vs_time_matched_negative_control_13_pairs",
  gene = rownames(tab), logFC = tab$logFC, AveExpr = tab$AveExpr,
  statistic = tab$t, PValue = tab$P.Value, FDR = tab$adj.P.Val
)
write_tsv(mptr[order(mptr$PValue), ], "GSE165177_MPTR_paired_DE.tsv")
write_tsv(
  data.frame(
    pair = transient_keys,
    negative_sample = mptr_meta$title[matched_negative],
    transient_sample = mptr_meta$title[transient_idx],
    donor = sub("_.*", "", transient_keys),
    day = as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", transient_keys)),
    experiment = sub(".*_", "", transient_keys)
  ),
  "selected_13_pairs.tsv"
)

# The published MPTR paper treats day 13 as the principal rejuvenation window.
# Pooling days 10/13/15/17 estimates a generic reprogramming response instead,
# so retain the paired all-day model above but also report each day separately.
pair_delta <- expression[, transient_idx, drop = FALSE] - expression[, matched_negative, drop = FALSE]
colnames(pair_delta) <- transient_keys
pair_day <- as.integer(sub("^[^_]+_([^_]+)_.*", "\\1", transient_keys))
pair_donor <- sub("_.*", "", transient_keys)

day_effects <- list()
day_inventory <- list()
for (day_value in sort(unique(pair_day))) {
  day_columns <- which(pair_day == day_value)
  donors <- unique(pair_donor[day_columns])
  donor_delta <- sapply(donors, function(donor) {
    donor_columns <- day_columns[pair_donor[day_columns] == donor]
    rowMeans(pair_delta[, donor_columns, drop = FALSE])
  })
  if (is.null(dim(donor_delta))) donor_delta <- matrix(donor_delta, ncol = 1,
                                                       dimnames = list(rownames(pair_delta), donors))
  effect <- data.frame(
    dataset = "GSE165177",
    contrast = paste0("MPTR_day", day_value, "_vs_time_matched_negative"),
    gene = rownames(donor_delta),
    logFC = rowMeans(donor_delta),
    donor_SD = if (ncol(donor_delta) > 1) apply(donor_delta, 1, sd) else NA_real_,
    n_pairs = length(day_columns),
    n_unique_donors = length(donors)
  )
  day_effects[[as.character(day_value)]] <- effect
  write_tsv(effect[order(-abs(effect$logFC)), ],
            paste0("GSE165177_MPTR_day", day_value, "_donor_aggregated_descriptive.tsv"))
  day_inventory[[as.character(day_value)]] <- data.frame(
    day = day_value,
    n_pairs = length(day_columns),
    n_unique_donors = length(donors),
    donors = paste(donors, collapse = ";"),
    inference_status = if (length(donors) >= 3) "descriptive_small_n" else "descriptive_only_n_donors_lt3"
  )
}
write_tsv(do.call(rbind, day_inventory), "MPTR_day_specific_sample_inventory.tsv")

signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

read_effect <- function(filename) {
  x <- read.delim(file.path(root, "derived", "figure2_public_aging", filename), check.names = FALSE)
  x[, c("gene", "logFC", "PValue", "FDR")]
}

reference <- list(
  GSE113957_adult_age = read_effect("GSE113957_age22_89_per_decade_DE.tsv"),
  GSE226189_age = read_effect("GSE226189_age_per_decade_DE.tsv"),
  GSE179848_late_replicative = read_effect("GSE179848_late_vs_early_donor_level_DE.tsv"),
  GSE109700_early_senescence = read_effect("GSE109700_early_vs_proliferating_DE.tsv"),
  GSE109700_deep_senescence = read_effect("GSE109700_deep_vs_proliferating_DE.tsv"),
  GSE191055_P27 = read_effect("GSE191055_P27_vs_P4_DE.tsv"),
  GSE93535_SIPS = read_effect("GSE93535_SIPS_vs_Q_DE.tsv")
)

correlation <- function(target, label) {
  joined <- merge(mptr[, c("gene", "logFC", "PValue")], target[, c("gene", "logFC", "PValue")],
                  by = "gene", suffixes = c("_MPTR", "_reference"))
  x <- signed_z(joined$logFC_MPTR, joined$PValue_MPTR)
  y <- signed_z(joined$logFC_reference, joined$PValue_reference)
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  data.frame(
    reference = label, n_shared = nrow(joined), rho_MPTR_vs_reference = unname(test$estimate),
    anti_concordance_score = -unname(test$estimate), p_value = test$p.value,
    opposite_direction_fraction = mean(sign(joined$logFC_MPTR) != sign(joined$logFC_reference))
  )
}
gene_correlations <- do.call(rbind, lapply(names(reference), function(name) correlation(reference[[name]], name)))
write_tsv(gene_correlations, "MPTR_vs_aging_senescence_gene_correlations.tsv")

day_correlation <- function(day_effect, target, day_value, label) {
  joined <- merge(day_effect[, c("gene", "logFC")], target[, c("gene", "logFC")],
                  by = "gene", suffixes = c("_MPTR", "_reference"))
  test <- cor.test(joined$logFC_MPTR, joined$logFC_reference, method = "spearman", exact = FALSE)
  data.frame(
    day = day_value, reference = label, n_shared = nrow(joined),
    rho_logFC_MPTR_vs_reference = unname(test$estimate), p_value = test$p.value,
    opposite_direction_fraction = mean(sign(joined$logFC_MPTR) != sign(joined$logFC_reference))
  )
}
day_gene_correlations <- do.call(rbind, lapply(names(day_effects), function(day_value) {
  do.call(rbind, lapply(names(reference), function(name) {
    day_correlation(day_effects[[day_value]], reference[[name]], as.integer(day_value), name)
  }))
}))
write_tsv(day_gene_correlations, "MPTR_day_specific_vs_aging_senescence_logFC_correlations.tsv")

gmt_connection <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt_lines <- readLines(gmt_connection, warn = FALSE)
close(gmt_connection)
reactome <- lapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) unique(x[-c(1, 2)]))
names(reactome) <- vapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) paste(x[1], x[2], sep = "__"), character(1))

run_fgsea <- function(tab) {
  ranks <- signed_z(tab$logFC, tab$PValue)
  names(ranks) <- tab$gene
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  set.seed(260911)
  ans <- fgseaMultilevel(
    reactome, ranks, minSize = 15, maxSize = 500, eps = 1e-10,
    BPPARAM = SerialParam()
  )
  ans$leadingEdge <- vapply(ans$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  as.data.frame(ans)
}

fg_mptr <- run_fgsea(mptr)
write_tsv(fg_mptr[order(fg_mptr$padj, -abs(fg_mptr$NES)), ], "GSE165177_MPTR_Reactome_fgsea.tsv")
pathway_rows <- lapply(names(reference), function(name) {
  fg_ref <- run_fgsea(reference[[name]])
  joined <- merge(fg_mptr[, c("pathway", "NES")], fg_ref[, c("pathway", "NES")],
                  by = "pathway", suffixes = c("_MPTR", "_reference"))
  test <- cor.test(joined$NES_MPTR, joined$NES_reference, method = "spearman", exact = FALSE)
  data.frame(
    reference = name, n_shared_pathways = nrow(joined),
    rho_MPTR_vs_reference = unname(test$estimate), anti_concordance_score = -unname(test$estimate),
    p_value = test$p.value
  )
})
pathway_correlations <- do.call(rbind, pathway_rows)
write_tsv(pathway_correlations, "MPTR_vs_aging_senescence_Reactome_correlations.tsv")

summary <- data.frame(
  metric = c("matched_pairs_all_days", "unique_donors_all_days", "tested_genes", "FDR05_genes_all_days"),
  value = c(length(transient_keys), length(unique(pair_donor)), nrow(mptr), sum(mptr$FDR < 0.05, na.rm = TRUE))
)
write_tsv(summary, "analysis_summary.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("MPTR matched-pair analysis completed:", out, "\n")
