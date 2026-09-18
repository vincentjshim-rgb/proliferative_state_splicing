#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(edgeR)
  library(fgsea)
  library(ggplot2)
  library(limma)
  library(org.Hs.eg.db)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "figure1_public_uva")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
meta <- read.delim(file.path(root, "derived", "audit", "all_geo_samples.tsv"), check.names = FALSE)

write_tsv <- function(x, filename) {
  write.table(x, file.path(out, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

collapse_matrix <- function(x, symbols, method = c("sum", "average")) {
  method <- match.arg(method)
  symbols <- toupper(trimws(symbols))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "--"
  x <- x[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (method == "sum") {
    rowsum(x, group = symbols, reorder = FALSE)
  } else {
    limma::avereps(x, ID = symbols)
  }
}

map_ensembl <- function(ids) {
  clean <- sub("\\..*$", "", ids)
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(clean),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  unname(mapped[clean])
}

run_edger <- function(counts, groups, contrasts, dataset) {
  groups <- factor(groups)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  y <- DGEList(counts = round(counts), group = groups)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)

  results <- list()
  for (contrast_name in names(contrasts)) {
    contrast <- makeContrasts(contrasts = contrasts[[contrast_name]], levels = design)
    test <- glmQLFTest(fit, contrast = contrast)
    tab <- topTags(test, n = Inf, sort.by = "none")$table
    tab$gene <- rownames(tab)
    tab$dataset <- dataset
    tab$contrast <- contrast_name
    tab <- tab[, c("dataset", "contrast", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
    results[[contrast_name]] <- tab
    write_tsv(tab[order(tab$PValue), ], paste0(dataset, "_", contrast_name, "_DE.tsv"))
  }

  logcpm <- cpm(y, log = TRUE, prior.count = 0.5)
  list(results = results, logcpm = logcpm, keep = keep, design = design)
}

# GSE125429: use the submitted paired-donor log2 block (columns 2-9).
gse125_path <- file.path(root, "uva", "GSE125429_log2_and_linear_values.txt.gz")
header_lines <- readLines(gzfile(gse125_path), n = 9, warn = FALSE)
sample_names <- strsplit(header_lines[8], "\t", fixed = TRUE)[[1]][2:9]
gse125_raw <- read.delim(
  gzfile(gse125_path), skip = 9, header = FALSE, fill = TRUE,
  quote = "\"", check.names = FALSE
)
gse125_raw <- gse125_raw[nzchar(trimws(gse125_raw[[1]])), , drop = FALSE]
gse125_expr <- as.matrix(gse125_raw[, 2:9])
storage.mode(gse125_expr) <- "numeric"
colnames(gse125_expr) <- sample_names
gse125_expr <- collapse_matrix(gse125_expr, gse125_raw[[1]], method = "average")

gse125_design_data <- data.frame(
  donor = factor(sub(" .*", "", sample_names)),
  condition = factor(ifelse(grepl("UVA", sample_names, ignore.case = TRUE), "UVA", "Control"),
                     levels = c("Control", "UVA"))
)
gse125_design <- model.matrix(~ donor + condition, data = gse125_design_data)
gse125_fit <- eBayes(lmFit(gse125_expr, gse125_design), trend = TRUE, robust = TRUE)
gse125_tab <- topTable(gse125_fit, coef = "conditionUVA", number = Inf, sort.by = "none")
gse125_tab$gene <- rownames(gse125_tab)
gse125_tab$dataset <- "GSE125429"
gse125_tab$contrast <- "UVA_vs_control"
gse125_tab <- gse125_tab[, c("dataset", "contrast", "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
names(gse125_tab)[names(gse125_tab) == "P.Value"] <- "PValue"
names(gse125_tab)[names(gse125_tab) == "adj.P.Val"] <- "FDR"
write_tsv(gse125_tab[order(gse125_tab$PValue), ], "GSE125429_UVA_vs_control_DE.tsv")

# GSE240226: reconstruct a raw-count matrix; labels come from GEO metadata.
gse240_files <- sort(list.files(file.path(root, "extracted", "GSE240226"), pattern = "\\.txt\\.gz$", full.names = TRUE))
gse240_tables <- lapply(gse240_files, function(path) read.delim(gzfile(path), check.names = FALSE))
gse240_ids <- gse240_tables[[1]][[1]]
stopifnot(all(vapply(gse240_tables, function(x) identical(x[[1]], gse240_ids), logical(1))))
gse240_symbols <- gse240_tables[[1]][[9]]
gse240_counts <- do.call(cbind, lapply(gse240_tables, function(x) x[[2]]))
gse240_accessions <- sub("_.*", "", basename(gse240_files))
colnames(gse240_counts) <- gse240_accessions
gse240_counts <- collapse_matrix(gse240_counts, gse240_symbols, method = "sum")
gse240_meta <- meta[match(gse240_accessions, meta$accession), ]
stopifnot(!anyNA(gse240_meta$accession))
gse240_group <- gsub(" ", "_", gse240_meta$treatment)
gse240_group <- gsub("\\+", "_plus_", gse240_group)
gse240_fit <- run_edger(
  gse240_counts,
  gse240_group,
  contrasts = c(
    UVA_vs_control = "UVA-control",
    Maifuyin_rescue_vs_UVA = "Maifuyin_plus_UVA-UVA",
    succinate_rescue_vs_UVA = "succinic_acid_plus_UVA-UVA"
  ),
  dataset = "GSE240226"
)

# GSE302943: GEO sample titles/treatments are reversed relative to submitted filenames.
# Publication Figure S6 labels the first five columns untreated and the next five
# UVA, matching filenames; reported MMP1-up/COL1A1-down directions also match only
# this assignment. Use filename labels and retain the conflict as an explicit audit.
gse302_files <- sort(list.files(file.path(root, "extracted", "GSE302943"), pattern = "\\.txt\\.gz$", full.names = TRUE))
gse302_tables <- lapply(gse302_files, function(path) read.delim(gzfile(path), check.names = FALSE))
gse302_ids <- gse302_tables[[1]][[1]]
stopifnot(all(vapply(gse302_tables, function(x) identical(x[[1]], gse302_ids), logical(1))))
gse302_counts <- do.call(cbind, lapply(gse302_tables, function(x) x[[2]]))
gse302_accessions <- sub("_.*", "", basename(gse302_files))
colnames(gse302_counts) <- gse302_accessions
gse302_symbols <- map_ensembl(gse302_ids)
gse302_counts <- collapse_matrix(gse302_counts, gse302_symbols, method = "sum")
gse302_meta <- meta[match(gse302_accessions, meta$accession), ]
stopifnot(!anyNA(gse302_meta$accession))
gse302_group <- ifelse(grepl("Untreated", basename(gse302_files), ignore.case = TRUE), "control", "UVA")
gse302_label_audit <- data.frame(
  accession = gse302_accessions,
  submitted_filename = basename(gse302_files),
  GEO_treatment = gse302_meta$treatment,
  resolved_group = gse302_group,
  metadata_conflict = tolower(gse302_meta$treatment) != tolower(gse302_group),
  adjudication_evidence = "publication_Figure_S6_column_labels_and_reported_marker_directions"
)
write_tsv(gse302_label_audit, "GSE302943_label_conflict_audit.tsv")
gse302_fit <- run_edger(
  gse302_counts,
  gse302_group,
  contrasts = c(UVA_vs_control = "UVA-control"),
  dataset = "GSE302943"
)

# GSE89005 sensitivity analysis: only conservatively sequence-mapped probes.
gse890_lines <- readLines(gzfile(file.path(root, "uva", "GSE89005_series_matrix.txt.gz")), warn = FALSE)
gse890_begin <- which(gse890_lines == "!series_matrix_table_begin") + 1
gse890_end <- which(gse890_lines == "!series_matrix_table_end") - 1
gse890_raw <- read.delim(
  textConnection(paste(gse890_lines[gse890_begin:gse890_end], collapse = "\n")),
  check.names = FALSE, quote = "\""
)
probe_map <- read.delim(
  file.path(root, "derived", "GPL16956_probe_mapping", "GPL16956_probe_to_GENCODE_v50.tsv"),
  check.names = FALSE
)
probe_map <- probe_map[probe_map$status == "PASS" & nzchar(probe_map$gene_symbol), c("probe", "gene_symbol")]
probe_index <- match(gse890_raw[[1]], probe_map$probe)
keep_probe <- !is.na(probe_index)
gse890_expr <- as.matrix(gse890_raw[keep_probe, -1])
storage.mode(gse890_expr) <- "numeric"
gse890_expr <- collapse_matrix(gse890_expr, probe_map$gene_symbol[probe_index[keep_probe]], method = "average")
gse890_accessions <- colnames(gse890_expr)
gse890_meta <- meta[match(gse890_accessions, meta$accession), ]
stopifnot(!anyNA(gse890_meta$accession))
gse890_group <- rep(NA_character_, length(gse890_accessions))
gse890_group[grepl("sham, parallel to single", gse890_meta$title, ignore.case = TRUE)] <- "single_sham"
gse890_group[grepl("6 h.*single UVA", gse890_meta$title, ignore.case = TRUE)] <- "single_UVA_6h"
gse890_group[grepl("24 h.*single UVA", gse890_meta$title, ignore.case = TRUE)] <- "single_UVA_24h"
gse890_group[grepl("sham, parallel to repeated", gse890_meta$title, ignore.case = TRUE)] <- "repeated_sham"
gse890_group[grepl("24 h.*repeated UVA", gse890_meta$title, ignore.case = TRUE)] <- "repeated_UVA_24h"
keep_sample <- !is.na(gse890_group)
gse890_design <- model.matrix(~ 0 + factor(gse890_group[keep_sample]))
colnames(gse890_design) <- levels(factor(gse890_group[keep_sample]))
gse890_fit0 <- lmFit(gse890_expr[, keep_sample, drop = FALSE], gse890_design)
gse890_contrasts <- c(
  single_6h_UVA_vs_sham = "single_UVA_6h-single_sham",
  single_24h_UVA_vs_sham = "single_UVA_24h-single_sham",
  repeated_24h_UVA_vs_sham = "repeated_UVA_24h-repeated_sham"
)
gse890_results <- list()
for (contrast_name in names(gse890_contrasts)) {
  contrast <- makeContrasts(contrasts = gse890_contrasts[[contrast_name]], levels = gse890_design)
  fit <- eBayes(contrasts.fit(gse890_fit0, contrast), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, number = Inf, sort.by = "none")
  tab$gene <- rownames(tab)
  tab$dataset <- "GSE89005"
  tab$contrast <- contrast_name
  tab <- tab[, c("dataset", "contrast", "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  names(tab)[names(tab) == "P.Value"] <- "PValue"
  names(tab)[names(tab) == "adj.P.Val"] <- "FDR"
  gse890_results[[contrast_name]] <- tab
  write_tsv(tab[order(tab$PValue), ], paste0("GSE89005_", contrast_name, "_DE.tsv"))
}

# Equal-study-weight signed Stouffer discovery meta-score.
signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

discovery <- merge(
  gse125_tab[, c("gene", "logFC", "PValue")],
  gse240_fit$results$UVA_vs_control[, c("gene", "logFC", "PValue")],
  by = "gene", suffixes = c("_GSE125429", "_GSE240226")
)
discovery$z_GSE125429 <- signed_z(discovery$logFC_GSE125429, discovery$PValue_GSE125429)
discovery$z_GSE240226 <- signed_z(discovery$logFC_GSE240226, discovery$PValue_GSE240226)
discovery$z_meta <- (discovery$z_GSE125429 + discovery$z_GSE240226) / sqrt(2)
discovery$P_meta <- 2 * pnorm(abs(discovery$z_meta), lower.tail = FALSE)
discovery$FDR_meta <- p.adjust(discovery$P_meta, method = "BH")
discovery$direction_consistent <- sign(discovery$logFC_GSE125429) == sign(discovery$logFC_GSE240226)
discovery$injury_rank <- ifelse(discovery$direction_consistent, discovery$z_meta, 0)
discovery <- discovery[order(-discovery$injury_rank), ]
write_tsv(discovery, "discovery_UVA_meta_signature.tsv")

# Held-out validation against GSE302943.
heldout <- merge(
  discovery,
  gse302_fit$results$UVA_vs_control[, c("gene", "logFC", "PValue", "FDR")],
  by = "gene", suffixes = c("", "_GSE302943")
)
heldout$z_GSE302943 <- signed_z(heldout$logFC, heldout$PValue)
names(heldout)[names(heldout) == "logFC"] <- "logFC_GSE302943"
names(heldout)[names(heldout) == "PValue"] <- "PValue_GSE302943"
names(heldout)[names(heldout) == "FDR"] <- "FDR_GSE302943"
heldout$heldout_same_direction <- sign(heldout$z_meta) == sign(heldout$z_GSE302943)
write_tsv(heldout[order(-abs(heldout$z_meta)), ], "heldout_GSE302943_validation.tsv")

# GSE240226 rescue relationships (shared-UVA comparator; interpreted descriptively).
injury240 <- gse240_fit$results$UVA_vs_control[, c("gene", "logFC")]
names(injury240)[2] <- "injury_logFC"
rescue <- merge(
  injury240,
  gse240_fit$results$Maifuyin_rescue_vs_UVA[, c("gene", "logFC")],
  by = "gene"
)
names(rescue)[3] <- "Maifuyin_rescue_logFC"
rescue <- merge(
  rescue,
  gse240_fit$results$succinate_rescue_vs_UVA[, c("gene", "logFC")],
  by = "gene"
)
names(rescue)[4] <- "succinate_rescue_logFC"
write_tsv(rescue, "GSE240226_injury_rescue_effects.tsv")

# Reactome enrichment of the discovery injury rank.
gmt_connection <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt_lines <- readLines(gmt_connection, warn = FALSE)
close(gmt_connection)
reactome <- lapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) unique(x[-c(1, 2)]))
names(reactome) <- vapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) paste(x[1], x[2], sep = "__"), character(1))
ranks <- discovery$injury_rank
names(ranks) <- discovery$gene
ranks <- sort(ranks[ranks != 0 & !duplicated(names(ranks))], decreasing = TRUE)
set.seed(260910)
fg <- fgseaMultilevel(
  reactome,
  ranks,
  minSize = 15,
  maxSize = 500,
  eps = 1e-10,
  BPPARAM = BiocParallel::SerialParam()
)
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
write_tsv(as.data.frame(fg), "discovery_UVA_meta_Reactome_fgsea.tsv")

run_reactome <- function(tab, dataset_name) {
  values <- signed_z(tab$logFC, tab$PValue)
  names(values) <- tab$gene
  values <- sort(values[is.finite(values) & !duplicated(names(values))], decreasing = TRUE)
  set.seed(260910)
  ans <- fgseaMultilevel(
    reactome, values, minSize = 15, maxSize = 500, eps = 1e-10,
    BPPARAM = BiocParallel::SerialParam()
  )
  ans$leadingEdge <- vapply(ans$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  ans$dataset <- dataset_name
  ans <- ans[, c("dataset", "pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")]
  write_tsv(as.data.frame(ans[order(ans$padj, -abs(ans$NES)), ]), paste0(dataset_name, "_Reactome_fgsea.tsv"))
  ans
}

fg125 <- run_reactome(gse125_tab, "GSE125429")
fg240 <- run_reactome(gse240_fit$results$UVA_vs_control, "GSE240226")
fg302 <- run_reactome(gse302_fit$results$UVA_vs_control, "GSE302943")
fg890 <- lapply(names(gse890_results), function(contrast_name) {
  run_reactome(gse890_results[[contrast_name]], paste0("GSE89005_", contrast_name))
})
names(fg890) <- names(gse890_results)
pathway_compare <- Reduce(
  function(x, y) merge(x, y, by = "pathway"),
  list(
    as.data.frame(fg125)[, c("pathway", "NES", "padj")],
    as.data.frame(fg240)[, c("pathway", "NES", "padj")],
    as.data.frame(fg302)[, c("pathway", "NES", "padj")]
  )
)
names(pathway_compare) <- c(
  "pathway", "NES_GSE125429", "FDR_GSE125429",
  "NES_GSE240226", "FDR_GSE240226", "NES_GSE302943", "FDR_GSE302943"
)
write_tsv(pathway_compare, "Reactome_NES_cross_study.tsv")

cross_enrichment <- function(source, target, source_name, target_name, top_n = 250) {
  source_z <- signed_z(source$logFC, source$PValue)
  names(source_z) <- source$gene
  source_z <- sort(source_z[is.finite(source_z) & !duplicated(names(source_z))], decreasing = TRUE)
  target_z <- signed_z(target$logFC, target$PValue)
  names(target_z) <- target$gene
  target_z <- sort(target_z[is.finite(target_z) & !duplicated(names(target_z))], decreasing = TRUE)
  sets <- list(source_up = names(head(source_z, top_n)), source_down = names(tail(source_z, top_n)))
  set.seed(260910)
  ans <- fgseaMultilevel(
    sets, target_z, minSize = 15, maxSize = 500, eps = 1e-10,
    BPPARAM = BiocParallel::SerialParam()
  )
  ans$leadingEdge <- vapply(ans$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  ans$source <- source_name
  ans$target <- target_name
  ans$top_n <- top_n
  as.data.frame(ans[, c("source", "target", "top_n", "pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")])
}

cross_gsea <- do.call(rbind, list(
  cross_enrichment(gse125_tab, gse240_fit$results$UVA_vs_control, "GSE125429", "GSE240226"),
  cross_enrichment(gse240_fit$results$UVA_vs_control, gse125_tab, "GSE240226", "GSE125429"),
  cross_enrichment(gse125_tab, gse302_fit$results$UVA_vs_control, "GSE125429", "GSE302943"),
  cross_enrichment(gse240_fit$results$UVA_vs_control, gse302_fit$results$UVA_vs_control, "GSE240226", "GSE302943")
))
write_tsv(cross_gsea, "top250_cross_study_fgsea.tsv")

effect_correlation <- function(x, y, x_name, y_name) {
  joined <- merge(
    x[, c("gene", "logFC", "PValue")],
    y[, c("gene", "logFC", "PValue")],
    by = "gene", suffixes = c("_x", "_y")
  )
  zx <- signed_z(joined$logFC_x, joined$PValue_x)
  zy <- signed_z(joined$logFC_y, joined$PValue_y)
  test <- cor.test(zx, zy, method = "spearman", exact = FALSE)
  data.frame(
    analysis_level = "gene_signed_z", comparison_x = x_name, comparison_y = y_name,
    n_shared = nrow(joined), rho = unname(test$estimate), p_value = test$p.value
  )
}

pathway_correlation <- function(x, y, x_name, y_name) {
  joined <- merge(as.data.frame(x)[, c("pathway", "NES")], as.data.frame(y)[, c("pathway", "NES")],
                  by = "pathway", suffixes = c("_x", "_y"))
  test <- cor.test(joined$NES_x, joined$NES_y, method = "spearman", exact = FALSE)
  data.frame(
    analysis_level = "Reactome_NES", comparison_x = x_name, comparison_y = y_name,
    n_shared = nrow(joined), rho = unname(test$estimate), p_value = test$p.value
  )
}

reference_tabs <- list(
  GSE125429_chronic_5week = gse125_tab,
  GSE240226_acute_3h = gse240_fit$results$UVA_vs_control,
  GSE302943_cumulative_5day = gse302_fit$results$UVA_vs_control
)
reference_fg <- list(
  GSE125429_chronic_5week = fg125,
  GSE240226_acute_3h = fg240,
  GSE302943_cumulative_5day = fg302
)
context_rows <- list()
for (contrast_name in names(gse890_results)) {
  for (reference_name in names(reference_tabs)) {
    context_rows[[length(context_rows) + 1]] <- effect_correlation(
      gse890_results[[contrast_name]], reference_tabs[[reference_name]],
      paste0("GSE89005_", contrast_name), reference_name
    )
    context_rows[[length(context_rows) + 1]] <- pathway_correlation(
      fg890[[contrast_name]], reference_fg[[reference_name]],
      paste0("GSE89005_", contrast_name), reference_name
    )
  }
}
context_correlations <- do.call(rbind, context_rows)
write_tsv(context_correlations, "GSE89005_context_correlations.tsv")

# Prespecified summary statistics.
cor_discovery <- cor.test(discovery$z_GSE125429, discovery$z_GSE240226, method = "spearman", exact = FALSE)
cor_heldout <- cor.test(heldout$z_meta, heldout$z_GSE302943, method = "spearman", exact = FALSE)
cor_heldout_125 <- cor.test(heldout$z_GSE125429, heldout$z_GSE302943, method = "spearman", exact = FALSE)
cor_heldout_240 <- cor.test(heldout$z_GSE240226, heldout$z_GSE302943, method = "spearman", exact = FALSE)
cor_path_125_240 <- cor.test(pathway_compare$NES_GSE125429, pathway_compare$NES_GSE240226, method = "spearman", exact = FALSE)
cor_path_125_302 <- cor.test(pathway_compare$NES_GSE125429, pathway_compare$NES_GSE302943, method = "spearman", exact = FALSE)
cor_path_240_302 <- cor.test(pathway_compare$NES_GSE240226, pathway_compare$NES_GSE302943, method = "spearman", exact = FALSE)
cor_maifuyin <- cor.test(rescue$injury_logFC, rescue$Maifuyin_rescue_logFC, method = "spearman", exact = FALSE)
cor_succinate <- cor.test(rescue$injury_logFC, rescue$succinate_rescue_logFC, method = "spearman", exact = FALSE)
summary <- data.frame(
  metric = c(
    "discovery_genes_shared", "discovery_direction_consistent_fraction",
    "discovery_spearman_rho", "discovery_spearman_p",
    "heldout_genes_shared", "heldout_same_direction_fraction",
    "heldout_spearman_rho", "heldout_spearman_p",
    "heldout_GSE125429_spearman_rho", "heldout_GSE125429_spearman_p",
    "heldout_GSE240226_spearman_rho", "heldout_GSE240226_spearman_p",
    "Reactome_NES_GSE125429_vs_GSE240226_rho", "Reactome_NES_GSE125429_vs_GSE240226_p",
    "Reactome_NES_GSE125429_vs_GSE302943_rho", "Reactome_NES_GSE125429_vs_GSE302943_p",
    "Reactome_NES_GSE240226_vs_GSE302943_rho", "Reactome_NES_GSE240226_vs_GSE302943_p",
    "GSE125429_FDR05_genes", "GSE240226_injury_FDR05_genes", "GSE302943_FDR05_genes",
    "GSE240226_injury_vs_Maifuyin_rescue_rho", "GSE240226_injury_vs_Maifuyin_rescue_p",
    "GSE240226_injury_vs_succinate_rescue_rho", "GSE240226_injury_vs_succinate_rescue_p"
  ),
  value = c(
    nrow(discovery), mean(discovery$direction_consistent),
    unname(cor_discovery$estimate), cor_discovery$p.value,
    nrow(heldout), mean(heldout$heldout_same_direction),
    unname(cor_heldout$estimate), cor_heldout$p.value,
    unname(cor_heldout_125$estimate), cor_heldout_125$p.value,
    unname(cor_heldout_240$estimate), cor_heldout_240$p.value,
    unname(cor_path_125_240$estimate), cor_path_125_240$p.value,
    unname(cor_path_125_302$estimate), cor_path_125_302$p.value,
    unname(cor_path_240_302$estimate), cor_path_240_302$p.value,
    sum(gse125_tab$FDR < 0.05, na.rm = TRUE),
    sum(gse240_fit$results$UVA_vs_control$FDR < 0.05, na.rm = TRUE),
    sum(gse302_fit$results$UVA_vs_control$FDR < 0.05, na.rm = TRUE),
    unname(cor_maifuyin$estimate), cor_maifuyin$p.value,
    unname(cor_succinate$estimate), cor_succinate$p.value
  )
)
write_tsv(summary, "analysis_summary.tsv")

# Diagnostic plots; these are intermediate and not yet manuscript figures.
p1 <- ggplot(discovery, aes(z_GSE125429, z_GSE240226, color = direction_consistent)) +
  geom_point(alpha = 0.22, size = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey75")) +
  labs(x = "GSE125429 signed z", y = "GSE240226 signed z", color = "Same direction") +
  theme_classic(base_size = 11)
ggsave(file.path(out, "discovery_concordance.png"), p1, width = 6, height = 5, dpi = 180)

p2 <- ggplot(heldout, aes(z_meta, z_GSE302943, color = heldout_same_direction)) +
  geom_point(alpha = 0.22, size = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) +
  scale_color_manual(values = c(`TRUE` = "#1B7837", `FALSE` = "grey75")) +
  labs(x = "Discovery UVA meta z", y = "Held-out GSE302943 signed z", color = "Same direction") +
  theme_classic(base_size = 11)
ggsave(file.path(out, "heldout_validation.png"), p2, width = 6, height = 5, dpi = 180)

rescue_long <- rbind(
  data.frame(injury = rescue$injury_logFC, rescue = rescue$Maifuyin_rescue_logFC, treatment = "Maifuyin"),
  data.frame(injury = rescue$injury_logFC, rescue = rescue$succinate_rescue_logFC, treatment = "Succinate")
)
p3 <- ggplot(rescue_long, aes(injury, rescue)) +
  geom_point(alpha = 0.18, size = 0.6, color = "#762A83") +
  geom_hline(yintercept = 0, linewidth = 0.25) + geom_vline(xintercept = 0, linewidth = 0.25) +
  facet_wrap(~ treatment) +
  labs(x = "UVA injury logFC", y = "Rescue vs UVA logFC") +
  theme_classic(base_size = 11)
ggsave(file.path(out, "GSE240226_injury_rescue.png"), p3, width = 9, height = 4.5, dpi = 180)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Figure 1 public UVA analysis completed:", out, "\n")
