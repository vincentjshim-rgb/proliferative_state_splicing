#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(BiocParallel)
  library(edgeR)
  library(fgsea)
  library(org.Hs.eg.db)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "local_repro_anchor")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(out, "star_genecounts")

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_star_counts <- function(sample, count_column = 4) {
  x <- read.delim(file.path(source_dir, sample, "ReadsPerGene.out.tab"), header = FALSE,
                  skip = 4, check.names = FALSE)
  x$ensembl <- sub("\\..*$", "", x[[1]])
  x$symbol <- mapIds(org.Hs.eg.db, keys = unique(x$ensembl), keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first")[x$ensembl]
  x <- x[!is.na(x$symbol) & nzchar(x$symbol), ]
  aggregate(x[[count_column]], by = list(symbol = x$symbol), FUN = sum) |>
    setNames(c("symbol", "count"))
}

make_expression <- function(count_column) {
  samples <- lapply(c(HDF = "HDF", IPS = "IPS", REP = "REP"), read_star_counts,
                    count_column = count_column)
  expression <- Reduce(function(x, y) merge(x, y, by = "symbol", all = TRUE),
                       lapply(names(samples), function(name) {
                         x <- samples[[name]]
                         names(x)[2] <- paste0("count_", name)
                         x
                       }))
  expression[is.na(expression)] <- 0
  expression
}

expression_reverse <- make_expression(4)
expression_unstranded <- make_expression(2)

score_anchor <- function(pseudocount = 0.5, expression_filter = 1,
                         count_mode = c("reverse_stranded", "unstranded")) {
  count_mode <- match.arg(count_mode)
  expression <- if (count_mode == "reverse_stranded") expression_reverse else expression_unstranded
  count_matrix <- as.matrix(expression[, c("count_HDF", "count_IPS", "count_REP")])
  rownames(count_matrix) <- toupper(expression$symbol)
  dge <- calcNormFactors(DGEList(count_matrix))
  values <- cpm(dge, log = FALSE)
  colnames(values) <- c("HDF", "IPS", "REP")
  keep <- rowSums(values >= expression_filter) >= 2
  x <- expression[keep, ]
  v <- values[keep, , drop = FALSE]
  log_hdf <- log2(v[, "HDF"] + pseudocount)
  log_ips <- log2(v[, "IPS"] + pseudocount)
  log_rep <- log2(v[, "REP"] + pseudocount)
  d_hdf <- log_rep - log_hdf
  d_ips <- log_rep - log_ips
  consistent <- sign(d_hdf) == sign(d_ips) & sign(d_hdf) != 0
  score <- ifelse(consistent, sign(d_hdf) * pmin(abs(d_hdf), abs(d_ips)), 0)
  data.frame(
    gene = toupper(x$symbol), count_HDF = x$count_HDF, count_iPSC = x$count_IPS, count_Repro = x$count_REP,
    normalized_HDF = v[, "HDF"], normalized_iPSC = v[, "IPS"], normalized_Repro = v[, "REP"],
    logFC_Repro_vs_HDF = d_hdf, logFC_Repro_vs_iPSC = d_ips,
    comparator_consistent = consistent, Repro_specific_score = score,
    mean_log2_expression = rowMeans(cbind(log_hdf, log_ips, log_rep)),
    normalization = paste0("STAR_", count_mode, "_TMM_CPM")
  )
}

anchor <- score_anchor()
write_tsv(anchor[order(-abs(anchor$Repro_specific_score)), ], "local_Repro_specific_gene_rank.tsv")
anchor_unstranded <- score_anchor(count_mode = "unstranded")
write_tsv(anchor_unstranded[order(-abs(anchor_unstranded$Repro_specific_score)), ],
          "local_Repro_specific_gene_rank_unstranded_sensitivity.tsv")
normalization_joined <- merge(
  anchor[anchor$comparator_consistent, c("gene", "Repro_specific_score")],
  anchor_unstranded[anchor_unstranded$comparator_consistent, c("gene", "Repro_specific_score")],
  by = "gene", suffixes = c("_reverse", "_unstranded")
)
write_tsv(data.frame(
  n_shared = nrow(normalization_joined),
  spearman_rho = cor(normalization_joined$Repro_specific_score_reverse,
                     normalization_joined$Repro_specific_score_unstranded, method = "spearman"),
  top250_overlap_fraction = length(intersect(
    head(normalization_joined$gene[order(-abs(normalization_joined$Repro_specific_score_reverse))], 250),
    head(normalization_joined$gene[order(-abs(normalization_joined$Repro_specific_score_unstranded))], 250)
  )) / 250
), "reverse_stranded_vs_unstranded_stability.tsv")

settings <- expand.grid(pseudocount = c(0.1, 0.5, 1), expression_filter = c(0.5, 1, 2))
stability_rows <- list()
base_rank <- anchor[anchor$comparator_consistent, c("gene", "Repro_specific_score")]
for (i in seq_len(nrow(settings))) {
  alt <- score_anchor(settings$pseudocount[i], settings$expression_filter[i])
  alt <- alt[alt$comparator_consistent, c("gene", "Repro_specific_score")]
  joined <- merge(base_rank, alt, by = "gene", suffixes = c("_base", "_alt"))
  stability_rows[[i]] <- data.frame(
    pseudocount = settings$pseudocount[i], expression_filter = settings$expression_filter[i],
    n_genes = nrow(alt), n_shared = nrow(joined),
    spearman_rho = cor(joined$Repro_specific_score_base, joined$Repro_specific_score_alt, method = "spearman"),
    top250_overlap_fraction = length(intersect(
      head(alt$gene[order(-abs(alt$Repro_specific_score))], 250),
      head(base_rank$gene[order(-abs(base_rank$Repro_specific_score))], 250)
    )) / 250
  )
}
write_tsv(do.call(rbind, stability_rows), "anchor_parameter_stability.tsv")

signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

read_effect <- function(folder, filename, label, axis, expected_relation, has_p = TRUE) {
  x <- read.delim(file.path(root, "derived", folder, filename), check.names = FALSE)
  x$gene <- toupper(x$gene)
  stat <- if (has_p) signed_z(x$logFC, x$PValue) else x$logFC
  data.frame(gene = x$gene, external_stat = stat, external_logFC = x$logFC,
             dataset_contrast = label, axis = axis, expected_relation = expected_relation)
}

effects <- list(
  read_effect("figure1_public_uva", "GSE125429_UVA_vs_control_DE.tsv", "GSE125429_chronic_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE240226_UVA_vs_control_DE.tsv", "GSE240226_acute_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE302943_UVA_vs_control_DE.tsv", "GSE302943_cumulative_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE89005_single_6h_UVA_vs_sham_DE.tsv", "GSE89005_single6h_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE89005_single_24h_UVA_vs_sham_DE.tsv", "GSE89005_single24h_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE89005_repeated_24h_UVA_vs_sham_DE.tsv", "GSE89005_repeated24h_UVA", "UVA_injury", "opposite"),
  read_effect("figure1_public_uva", "GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv", "GSE240226_Maifuyin_rescue", "UVA_rescue", "same"),
  read_effect("figure1_public_uva", "GSE240226_succinate_rescue_vs_UVA_DE.tsv", "GSE240226_succinate_rescue", "UVA_rescue", "same"),
  read_effect("figure2_public_aging", "GSE113957_age22_89_per_decade_DE.tsv", "GSE113957_adult_age", "chronological_age", "opposite"),
  read_effect("figure2_public_aging", "GSE226189_age_per_decade_DE.tsv", "GSE226189_age", "chronological_age", "opposite"),
  read_effect("figure2_public_aging", "GSE179848_late_vs_early_donor_level_DE.tsv", "GSE179848_late_replicative", "senescence", "opposite"),
  read_effect("figure2_public_aging", "GSE109700_early_vs_proliferating_DE.tsv", "GSE109700_early_senescence", "senescence", "opposite"),
  read_effect("figure2_public_aging", "GSE109700_deep_vs_proliferating_DE.tsv", "GSE109700_deep_senescence", "senescence", "opposite"),
  read_effect("figure2_public_aging", "GSE191055_P27_vs_P4_DE.tsv", "GSE191055_P27", "senescence", "opposite"),
  read_effect("figure2_public_aging", "GSE93535_SIPS_vs_Q_DE.tsv", "GSE93535_SIPS", "senescence", "opposite"),
  read_effect("figure2_public_aging", "GSE93535_SIPS1201_vs_SIPS_DE.tsv", "GSE93535_rescue", "senescence_rescue", "same"),
  read_effect("figure2_mptr", "GSE165177_MPTR_paired_DE.tsv", "GSE165177_MPTR_all_days", "partial_reprogramming", "same"),
  read_effect("figure2_mptr", "GSE165177_MPTR_day13_donor_aggregated_descriptive.tsv", "GSE165177_MPTR_day13", "partial_reprogramming", "same", FALSE)
)

matched_permutation <- function(x, y, expression_value, nperm = 5000, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  breaks <- unique(quantile(expression_value, probs = seq(0, 1, 0.1), na.rm = TRUE))
  bins <- cut(expression_value, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  observed <- cor(rx, ry)
  null <- numeric(nperm)
  indices <- split(seq_along(ry), bins)
  for (b in seq_len(nperm)) {
    permuted <- ry
    for (idx in indices) permuted[idx] <- sample(permuted[idx])
    null[b] <- cor(rx, permuted)
  }
  c(observed = observed,
    empirical_p_two_sided = (1 + sum(abs(null) >= abs(observed))) / (nperm + 1),
    null_q025 = unname(quantile(null, 0.025)), null_q975 = unname(quantile(null, 0.975)))
}

test_effect <- function(effect) {
  local <- anchor[anchor$comparator_consistent & anchor$Repro_specific_score != 0, ]
  joined <- merge(local, effect, by = "gene")
  perm <- matched_permutation(joined$Repro_specific_score, joined$external_stat, joined$mean_log2_expression)
  rho_logfc <- cor(joined$Repro_specific_score, joined$external_logFC, method = "spearman")
  data.frame(
    dataset_contrast = effect$dataset_contrast[1], axis = effect$axis[1],
    expected_relation = effect$expected_relation[1], n_shared = nrow(joined),
    rho_rank_stat = perm["observed"], rho_rank_logFC = rho_logfc,
    empirical_p_expression_matched = perm["empirical_p_two_sided"],
    null_q025 = perm["null_q025"], null_q975 = perm["null_q975"],
    direction_matches_hypothesis = ifelse(effect$expected_relation[1] == "same", perm["observed"] > 0, perm["observed"] < 0),
    same_direction_fraction = mean(sign(joined$Repro_specific_score) == sign(joined$external_logFC))
  )
}

association <- do.call(rbind, lapply(effects, test_effect))
association$FDR_expression_matched <- p.adjust(association$empirical_p_expression_matched, "BH")
write_tsv(association, "local_anchor_vs_external_association.tsv")

gmt_connection <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt_lines <- readLines(gmt_connection, warn = FALSE)
close(gmt_connection)
reactome <- lapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) unique(toupper(x[-c(1, 2)])))
names(reactome) <- vapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) paste(x[1], x[2], sep = "__"), character(1))

local_ranks <- anchor$Repro_specific_score
names(local_ranks) <- anchor$gene
local_ranks <- sort(local_ranks[local_ranks != 0 & is.finite(local_ranks)], decreasing = TRUE)
set.seed(260911)
local_fgsea <- fgseaMultilevel(reactome, local_ranks, minSize = 15, maxSize = 500, eps = 1e-10,
                               BPPARAM = SerialParam())
local_fgsea$leadingEdge <- vapply(local_fgsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write_tsv(as.data.frame(local_fgsea)[order(local_fgsea$padj, -abs(local_fgsea$NES)), ],
          "local_Repro_specific_Reactome_fgsea.tsv")

gene_set_rows <- list()
for (effect in effects) {
  ranks <- effect$external_stat
  names(ranks) <- effect$gene
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  for (n_top in c(100, 250, 500, 1000)) {
    up <- head(anchor$gene[order(-anchor$Repro_specific_score)], n_top)
    down <- head(anchor$gene[order(anchor$Repro_specific_score)], n_top)
    sets <- list(local_Repro_up = up, local_Repro_down = down)
    set.seed(260911 + n_top)
    ans <- fgseaMultilevel(sets, ranks, minSize = 10, maxSize = 2000, eps = 0,
                           BPPARAM = SerialParam())
    nes_up <- ans$NES[match("local_Repro_up", ans$pathway)]
    nes_down <- ans$NES[match("local_Repro_down", ans$pathway)]
    gene_set_rows[[length(gene_set_rows) + 1]] <- data.frame(
      dataset_contrast = effect$dataset_contrast[1], axis = effect$axis[1],
      expected_relation = effect$expected_relation[1], top_n = n_top,
      NES_local_up = nes_up, NES_local_down = nes_down,
      directional_score = (nes_up - nes_down) / 2,
      padj_local_up = ans$padj[match("local_Repro_up", ans$pathway)],
      padj_local_down = ans$padj[match("local_Repro_down", ans$pathway)]
    )
  }
}
write_tsv(do.call(rbind, gene_set_rows), "local_top_gene_sets_in_external_effects.tsv")

summary <- data.frame(
  metric = c("mapped_symbols", "filtered_genes", "comparator_consistent_genes",
             "consistent_fraction", "score_abs_ge_log2_1.25", "score_abs_ge_log2_1.5", "score_abs_ge_1"),
  value = c(nrow(expression_reverse), nrow(anchor), sum(anchor$comparator_consistent),
            mean(anchor$comparator_consistent),
            sum(abs(anchor$Repro_specific_score) >= log2(1.25)),
            sum(abs(anchor$Repro_specific_score) >= log2(1.5)),
            sum(abs(anchor$Repro_specific_score) >= 1))
)
write_tsv(summary, "analysis_summary.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Local Repro-CM anchor analysis completed:", out, "\n")
