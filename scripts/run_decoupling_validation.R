#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

root <- "public_data_tierA"
anchor_dir <- file.path(root, "derived", "local_repro_anchor")
out <- file.path(root, "derived", "decoupling_validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
}

signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

rank_z <- function(x) {
  ok <- is.finite(x)
  ans <- rep(NA_real_, length(x))
  ans[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok))
  ans
}

read_effect <- function(folder, filename, label, axis) {
  x <- read.delim(file.path(root, "derived", folder, filename), check.names = FALSE)
  x$gene <- toupper(x$gene)
  x <- x[!duplicated(x$gene), ]
  data.frame(gene = x$gene,
             stat = signed_z(x$logFC, x$PValue),
             logFC = x$logFC,
             dataset = label, axis = axis)
}

effects <- list(
  GSE240226 = read_effect("figure1_public_uva", "GSE240226_UVA_vs_control_DE.tsv",
                          "GSE240226 acute UVA", "UVA"),
  GSE302943 = read_effect("figure1_public_uva", "GSE302943_UVA_vs_control_DE.tsv",
                          "GSE302943 cumulative UVA", "UVA"),
  GSE179848 = read_effect("figure2_public_aging", "GSE179848_late_vs_early_donor_level_DE.tsv",
                          "GSE179848 longitudinal RS", "senescence"),
  GSE109700 = read_effect("figure2_public_aging", "GSE109700_deep_vs_proliferating_DE.tsv",
                          "GSE109700 deep RS", "senescence"),
  GSE191055 = read_effect("figure2_public_aging", "GSE191055_P27_vs_P4_DE.tsv",
                          "GSE191055 passage RS", "senescence"),
  GSE93535 = read_effect("figure2_public_aging", "GSE93535_SIPS_vs_Q_DE.tsv",
                         "GSE93535 SIPS", "senescence")
)

anchor <- read.delim(file.path(anchor_dir, "local_Repro_specific_gene_rank.tsv"),
                     check.names = FALSE)
anchor$gene <- toupper(anchor$gene)
anchor <- anchor[!duplicated(anchor$gene), ]

make_matrix <- function(effect_names = names(effects)) {
  result <- anchor[, c("gene", "Repro_specific_score", "logFC_Repro_vs_HDF",
                       "logFC_Repro_vs_iPSC", "mean_log2_expression",
                       "comparator_consistent", "count_HDF", "count_iPSC", "count_Repro")]
  for (nm in effect_names) {
    z <- effects[[nm]][, c("gene", "stat")]
    names(z)[2] <- nm
    result <- merge(result, z, by = "gene")
  }
  result
}

primary <- make_matrix()
primary$UVA_meta <- rowMeans(cbind(rank_z(primary$GSE240226), rank_z(primary$GSE302943)))
primary$senescence_meta <- rowMeans(cbind(rank_z(primary$GSE179848), rank_z(primary$GSE109700),
                                         rank_z(primary$GSE191055), rank_z(primary$GSE93535)))
write_tsv(primary, "primary_gene_axis_matrix.tsv")

score_metrics <- function(data, score, label = "primary") {
  ok <- is.finite(score) & score != 0 & data$comparator_consistent
  u <- cor(score[ok], data$UVA_meta[ok], method = "spearman")
  s <- cor(score[ok], data$senescence_meta[ok], method = "spearman")
  residual_s <- resid(lm(senescence_meta ~ UVA_meta, data = data[ok, ]))
  sr <- cor(score[ok], residual_s, method = "spearman")
  data.frame(analysis = label, n_genes = sum(ok), rho_UVA = u,
             rho_senescence = s, decoupling_index = u - s,
             rho_senescence_residualized_for_UVA = sr)
}

primary_metrics <- score_metrics(primary, primary$Repro_specific_score)

# Expression-decile-matched empirical null for the central contrast.
set.seed(260911)
ok <- primary$comparator_consistent & primary$Repro_specific_score != 0 &
      is.finite(primary$Repro_specific_score)
pdat <- primary[ok, ]
breaks <- unique(quantile(pdat$mean_log2_expression, probs = seq(0, 1, 0.1), na.rm = TRUE))
pdat$expression_bin <- cut(pdat$mean_log2_expression, breaks = breaks,
                           include.lowest = TRUE, labels = FALSE)
bin_index <- split(seq_len(nrow(pdat)), pdat$expression_bin)
nperm <- 10000L
null <- matrix(NA_real_, nrow = nperm, ncol = 3,
               dimnames = list(NULL, c("rho_UVA", "rho_senescence", "decoupling_index")))
for (i in seq_len(nperm)) {
  perm <- pdat$Repro_specific_score
  for (idx in bin_index) perm[idx] <- sample(perm[idx])
  null[i, 1] <- cor(perm, pdat$UVA_meta, method = "spearman")
  null[i, 2] <- cor(perm, pdat$senescence_meta, method = "spearman")
  null[i, 3] <- null[i, 1] - null[i, 2]
}
null_df <- as.data.frame(null)
write_tsv(null_df, "decoupling_expression_matched_null_10000.tsv")
primary_metrics$empirical_p_decoupling <-
  (1 + sum(null_df$decoupling_index >= primary_metrics$decoupling_index)) / (nperm + 1)
primary_metrics$null_q025 <- quantile(null_df$decoupling_index, 0.025)
primary_metrics$null_q975 <- quantile(null_df$decoupling_index, 0.975)
write_tsv(primary_metrics, "primary_decoupling_metrics.tsv")

# Per-dataset effect sizes, with the expression-matched interval inherited from
# the earlier 5,000-permutation analysis.
assoc <- read.delim(file.path(anchor_dir, "local_anchor_vs_external_association.tsv"),
                    check.names = FALSE)
name_map <- c(
  GSE240226_acute_UVA = "GSE240226 acute UVA",
  GSE302943_cumulative_UVA = "GSE302943 cumulative UVA",
  GSE179848_late_replicative = "GSE179848 longitudinal RS",
  GSE109700_deep_senescence = "GSE109700 deep RS",
  GSE191055_P27 = "GSE191055 passage RS",
  GSE93535_SIPS = "GSE93535 SIPS"
)
dataset_metrics <- assoc[match(names(name_map), assoc$dataset_contrast), ]
dataset_metrics$display_label <- unname(name_map[dataset_metrics$dataset_contrast])
dataset_metrics$primary_axis <- ifelse(dataset_metrics$axis == "UVA_injury", "UVA response", "Senescence")
dataset_metrics$expected_in_decoupling_model <- ifelse(dataset_metrics$axis == "UVA_injury", "same", "opposite")
dataset_metrics$matches_decoupling_model <- ifelse(dataset_metrics$axis == "UVA_injury",
                                                    dataset_metrics$rho_rank_stat > 0,
                                                    dataset_metrics$rho_rank_stat < 0)
write_tsv(dataset_metrics, "primary_dataset_metrics.tsv")

# Leave one dataset out from either meta-axis.
loo_rows <- list()
for (drop in names(effects)) {
  keep <- setdiff(names(effects), drop)
  dat <- make_matrix(keep)
  uva_names <- intersect(c("GSE240226", "GSE302943"), keep)
  sen_names <- intersect(c("GSE179848", "GSE109700", "GSE191055", "GSE93535"), keep)
  dat$UVA_meta <- rowMeans(sapply(uva_names, function(nm) rank_z(dat[[nm]])))
  dat$senescence_meta <- rowMeans(sapply(sen_names, function(nm) rank_z(dat[[nm]])))
  loo_rows[[drop]] <- cbind(dropped_dataset = drop,
                            score_metrics(dat, dat$Repro_specific_score, paste0("drop_", drop)))
}
loo <- do.call(rbind, loo_rows)
rownames(loo) <- NULL
write_tsv(loo, "leave_one_dataset_out.tsv")

# Comparator sensitivity: the central model must be visible against HDF-CM and
# iPSC-CM separately, not only after taking their conservative minimum.
comparator_rows <- rbind(
  score_metrics(transform(primary, comparator_consistent = TRUE),
                primary$logFC_Repro_vs_HDF, "Repro_vs_HDF"),
  score_metrics(transform(primary, comparator_consistent = TRUE),
                primary$logFC_Repro_vs_iPSC, "Repro_vs_iPSC"),
  primary_metrics[, names(score_metrics(primary, primary$Repro_specific_score))]
)
comparator_rows$analysis[3] <- "conservative_AND_anchor"
write_tsv(comparator_rows, "comparator_sensitivity.tsv")

# Remove gene families/programs most likely to dominate the association.
gmt_con <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt_lines <- readLines(gmt_con, warn = FALSE)
close(gmt_con)
gmt_split <- strsplit(gmt_lines, "\t", fixed = TRUE)
pathway_names <- vapply(gmt_split, `[`, character(1), 1)
program_idx <- grepl("RIBOSOM|TRANSLAT|RRNA|NONSENSE.MEDIATED|INTERFERON|EIF2|RNA PROCESS",
                     pathway_names, ignore.case = TRUE)
program_genes <- unique(toupper(unlist(lapply(gmt_split[program_idx], function(x) x[-c(1, 2)]))))
cell_cycle_idx <- grepl("CELL CYCLE|MITOTIC|DNA REPLICATION|G1.S|G2.M",
                        pathway_names, ignore.case = TRUE)
cell_cycle_genes <- unique(toupper(unlist(lapply(gmt_split[cell_cycle_idx], function(x) x[-c(1, 2)]))))
filters <- list(
  none = rep(TRUE, nrow(primary)),
  no_RPL_RPS_MT = !grepl("^(RPL|RPS|MRPL|MRPS|MT-)", primary$gene),
  no_translation_rRNA_IFN = !(primary$gene %in% program_genes),
  no_cell_cycle_DNA_replication = !(primary$gene %in% cell_cycle_genes),
  no_translation_IFN_or_cell_cycle = !(primary$gene %in% union(program_genes, cell_cycle_genes)),
  no_top250_local = !(primary$gene %in% head(anchor$gene[order(-abs(anchor$Repro_specific_score))], 250)),
  no_top500_local = !(primary$gene %in% head(anchor$gene[order(-abs(anchor$Repro_specific_score))], 500))
)
filter_rows <- lapply(names(filters), function(nm) {
  dat <- primary[filters[[nm]], ]
  score_metrics(dat, dat$Repro_specific_score, nm)
})
filter_rows <- do.call(rbind, filter_rows)
write_tsv(filter_rows, "gene_program_removal_sensitivity.tsv")

# Technical read-thinning only: this assesses computational stability, not
# biological reproducibility. Counts are independently thinned to 50%.
base_score <- setNames(anchor$Repro_specific_score, anchor$gene)
thin_rows <- vector("list", 100)
set.seed(260912)
for (i in seq_len(100)) {
  counts <- cbind(
    HDF = rbinom(nrow(anchor), anchor$count_HDF, 0.5),
    iPSC = rbinom(nrow(anchor), anchor$count_iPSC, 0.5),
    Repro = rbinom(nrow(anchor), anchor$count_Repro, 0.5)
  )
  rownames(counts) <- anchor$gene
  cpm_mat <- cpm(calcNormFactors(DGEList(counts)), log = FALSE)
  keep <- rowSums(cpm_mat >= 1) >= 2
  dh <- log2(cpm_mat[, "Repro"] + 0.5) - log2(cpm_mat[, "HDF"] + 0.5)
  di <- log2(cpm_mat[, "Repro"] + 0.5) - log2(cpm_mat[, "iPSC"] + 0.5)
  consistent <- sign(dh) == sign(di) & sign(dh) != 0 & keep
  score <- ifelse(consistent, sign(dh) * pmin(abs(dh), abs(di)), 0)
  names(score) <- rownames(counts)
  common <- intersect(primary$gene, names(score)[score != 0])
  d <- primary[match(common, primary$gene), ]
  d$comparator_consistent <- TRUE
  met <- score_metrics(d, score[common], paste0("thin_", i))
  rank_common <- intersect(names(base_score)[base_score != 0], names(score)[score != 0])
  met$rho_with_full_anchor <- cor(base_score[rank_common], score[rank_common], method = "spearman")
  met$top250_overlap <- length(intersect(
    head(names(sort(abs(base_score), decreasing = TRUE)), 250),
    head(names(sort(abs(score), decreasing = TRUE)), 250))) / 250
  thin_rows[[i]] <- met
}
thin <- do.call(rbind, thin_rows)
write_tsv(thin, "read_thinning_50pct_100runs.tsv")

# Pathway-level replication using existing independent Reactome GSEA results.
read_nes <- function(folder, filename, label, axis) {
  x <- read.delim(file.path(root, "derived", folder, filename), check.names = FALSE)
  x <- x[, c("pathway", "NES")]
  names(x)[2] <- label
  x
}
local_nes <- read.delim(file.path(anchor_dir, "local_Repro_specific_Reactome_fgsea.tsv"),
                        check.names = FALSE)[, c("pathway", "NES", "padj")]
names(local_nes)[2:3] <- c("local_NES", "local_FDR")
nes_list <- list(
  GSE240226 = read_nes("figure1_public_uva", "GSE240226_Reactome_fgsea.tsv", "GSE240226", "UVA"),
  GSE302943 = read_nes("figure1_public_uva", "GSE302943_Reactome_fgsea.tsv", "GSE302943", "UVA"),
  GSE179848 = read_nes("figure2_public_aging", "GSE179848_late_Reactome_fgsea.tsv", "GSE179848", "senescence"),
  GSE109700 = read_nes("figure2_public_aging", "GSE109700_deep_Reactome_fgsea.tsv", "GSE109700", "senescence"),
  GSE191055 = read_nes("figure2_public_aging", "GSE191055_P27_Reactome_fgsea.tsv", "GSE191055", "senescence"),
  GSE93535 = read_nes("figure2_public_aging", "GSE93535_SIPS_Reactome_fgsea.tsv", "GSE93535", "senescence")
)
pathway <- local_nes
for (nm in names(nes_list)) pathway <- merge(pathway, nes_list[[nm]], by = "pathway")
pathway$UVA_meta_NES <- rowMeans(pathway[, c("GSE240226", "GSE302943")])
pathway$senescence_meta_NES <- rowMeans(pathway[, c("GSE179848", "GSE109700", "GSE191055", "GSE93535")])
write_tsv(pathway[order(-abs(pathway$local_NES)), ], "pathway_axis_matrix.tsv")
pathway_metrics <- data.frame(
  level = "Reactome_pathway",
  n_pathways = nrow(pathway),
  rho_local_vs_UVA = cor(pathway$local_NES, pathway$UVA_meta_NES, method = "spearman"),
  rho_local_vs_senescence = cor(pathway$local_NES, pathway$senescence_meta_NES, method = "spearman")
)
pathway_metrics$decoupling_index <- pathway_metrics$rho_local_vs_UVA - pathway_metrics$rho_local_vs_senescence
write_tsv(pathway_metrics, "pathway_decoupling_metrics.tsv")

# A rescue effect is expected to oppose the original UVA contrast. Therefore,
# test whether Repro-CM still aligns with two rescue contrasts after removing
# the component linearly predictable from the UVA-injury signature.
injury <- effects$GSE240226[, c("gene", "stat")]
names(injury)[2] <- "injury"
rescue_effects <- list(
  Maifuyin = read_effect("figure1_public_uva", "GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv",
                         "Maifuyin rescue", "rescue"),
  succinate = read_effect("figure1_public_uva", "GSE240226_succinate_rescue_vs_UVA_DE.tsv",
                          "succinate rescue", "rescue")
)
rescue_rows <- list()
for (nm in names(rescue_effects)) {
  rr <- rescue_effects[[nm]][, c("gene", "stat")]
  names(rr)[2] <- "rescue"
  dat <- merge(anchor[anchor$comparator_consistent & anchor$Repro_specific_score != 0,
                      c("gene", "Repro_specific_score", "mean_log2_expression")], injury, by = "gene")
  dat <- merge(dat, rr, by = "gene")
  residual_rescue <- resid(lm(rank_z(rescue) ~ rank_z(injury), data = dat))
  observed <- cor(dat$Repro_specific_score, residual_rescue, method = "spearman")
  set.seed(260913)
  breaks <- unique(quantile(dat$mean_log2_expression, probs = seq(0, 1, 0.1), na.rm = TRUE))
  bins <- cut(dat$mean_log2_expression, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  idxs <- split(seq_len(nrow(dat)), bins)
  null_rescue <- numeric(5000)
  for (i in seq_along(null_rescue)) {
    perm <- dat$Repro_specific_score
    for (idx in idxs) perm[idx] <- sample(perm[idx])
    null_rescue[i] <- cor(perm, residual_rescue, method = "spearman")
  }
  rescue_rows[[nm]] <- data.frame(
    rescue = nm, n_genes = nrow(dat),
    rho_injury_vs_rescue = cor(dat$injury, dat$rescue, method = "spearman"),
    rho_local_vs_rescue = cor(dat$Repro_specific_score, dat$rescue, method = "spearman"),
    rho_local_vs_rescue_residualized_for_injury = observed,
    empirical_p_expression_matched = (1 + sum(abs(null_rescue) >= abs(observed))) / (length(null_rescue) + 1),
    null_q025 = quantile(null_rescue, 0.025), null_q975 = quantile(null_rescue, 0.975)
  )
}
rescue_rows <- do.call(rbind, rescue_rows)
rescue_rows$FDR_expression_matched <- p.adjust(rescue_rows$empirical_p_expression_matched, "BH")
write_tsv(rescue_rows, "rescue_residualized_validation.tsv")

# Central evidence figures. These are diagnostic manuscript drafts; panel labels
# explicitly show which statement each plot tests.
p1 <- ggplot(dataset_metrics,
             aes(x = rho_rank_stat, y = reorder(display_label, rho_rank_stat),
                 colour = primary_axis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = null_q025, xmax = null_q975), height = 0.18, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("UVA response" = "#D55E00", "Senescence" = "#0072B2")) +
  labs(x = "Spearman rho with Repro-CM anchor", y = NULL, colour = NULL,
       title = "A. Independent datasets show opposite axis relationships",
       subtitle = "Bars: 95% expression-matched permutation interval") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

observed_d <- primary_metrics$decoupling_index
p2 <- ggplot(null_df, aes(x = decoupling_index)) +
  geom_histogram(bins = 60, fill = "grey75", colour = "white") +
  geom_vline(xintercept = observed_d, colour = "#CC0000", linewidth = 1.1) +
  annotate("text", x = observed_d, y = Inf, vjust = 1.6, hjust = 1.05,
           label = sprintf("observed = %.3f\nempirical p = %.5f",
                           observed_d, primary_metrics$empirical_p_decoupling), colour = "#990000") +
  labs(x = "Decoupling index: rho(UVA) - rho(senescence)", y = "Permutations",
       title = "B. Central contrast exceeds the matched null") +
  theme_bw(base_size = 11)

p3dat <- rbind(
  data.frame(test = loo$dropped_dataset, value = loo$decoupling_index, class = "Leave-one-study-out"),
  data.frame(test = filter_rows$analysis, value = filter_rows$decoupling_index, class = "Gene-removal sensitivity"),
  data.frame(test = comparator_rows$analysis, value = comparator_rows$decoupling_index, class = "Comparator sensitivity")
)
p3 <- ggplot(p3dat, aes(x = value, y = reorder(test, value), colour = class)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") + geom_point(size = 2.6) +
  labs(x = "Decoupling index", y = NULL, colour = NULL,
       title = "C. The contrast survives prespecified sensitivity tests") +
  theme_bw(base_size = 10) + theme(legend.position = "top")

ggsave(file.path(out, "Figure_decoupling_A_dataset_forest.png"), p1, width = 8.2, height = 4.8, dpi = 300)
ggsave(file.path(out, "Figure_decoupling_B_null.png"), p2, width = 6.3, height = 4.6, dpi = 300)
ggsave(file.path(out, "Figure_decoupling_C_robustness.png"), p3, width = 8.2, height = 6.2, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Decoupling validation completed:", out, "\n")
