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
  library(readxl)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "figure2_public_aging")
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
  if (method == "sum") rowsum(x, symbols, reorder = FALSE) else limma::avereps(x, ID = symbols)
}

map_ids <- function(ids, keytype) {
  clean <- if (keytype == "ENSEMBL") sub("\\..*$", "", ids) else as.character(ids)
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = unique(clean), keytype = keytype,
    column = "SYMBOL", multiVals = "first"
  )
  unname(mapped[clean])
}

run_edger_groups <- function(counts, groups, contrasts, dataset) {
  groups <- factor(groups)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  y <- DGEList(round(counts))
  keep <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  ans <- list()
  for (contrast_name in names(contrasts)) {
    contrast <- makeContrasts(contrasts = contrasts[[contrast_name]], levels = design)
    tab <- topTags(glmQLFTest(fit, contrast = contrast), n = Inf, sort.by = "none")$table
    tab$gene <- rownames(tab)
    tab$dataset <- dataset
    tab$contrast <- contrast_name
    tab <- tab[, c("dataset", "contrast", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
    ans[[contrast_name]] <- tab
    write_tsv(tab[order(tab$PValue), ], paste0(dataset, "_", contrast_name, "_DE.tsv"))
  }
  ans
}

run_edger_continuous <- function(counts, sample_data, dataset) {
  sample_data$age_decade <- as.numeric(sample_data$age) / 10
  sample_data$sex <- relevel(factor(tolower(sample_data$sex)), ref = "male")
  design <- model.matrix(~ age_decade + sex, data = sample_data)
  y <- DGEList(round(counts))
  keep <- filterByExpr(y, design = design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  tab <- topTags(glmQLFTest(fit, coef = which(colnames(design) == "age_decade")), n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$dataset <- dataset
  tab$contrast <- "per_10_years_adjusted_for_sex"
  tab <- tab[, c("dataset", "contrast", "gene", "logFC", "logCPM", "F", "PValue", "FDR")]
  write_tsv(tab[order(tab$PValue), ], paste0(dataset, "_age_per_decade_DE.tsv"))
  tab
}

signed_z <- function(logfc, pvalue) {
  pvalue <- pmax(pmin(pvalue, 1), .Machine$double.xmin)
  sign(logfc) * qnorm(pvalue / 2, lower.tail = FALSE)
}

# GSE109700: proliferating, early senescent and deep/late senescent, n=3 each.
gse109_raw <- read.delim(
  gzfile(file.path(root, "senescence", "GSE109700_LF1_Counts_HiSat2_FeatureCounts.txt.gz")),
  skip = 1, check.names = FALSE
)
gse109_counts <- as.matrix(gse109_raw[, -1])
storage.mode(gse109_counts) <- "numeric"
gse109_symbols <- map_ids(gse109_raw[[1]], "ENTREZID")
gse109_counts <- collapse_matrix(gse109_counts, gse109_symbols, "sum")
gse109_groups <- rep(c("proliferating", "early", "deep"), each = 3)
gse109 <- run_edger_groups(
  gse109_counts, gse109_groups,
  c(early_vs_proliferating = "early-proliferating", deep_vs_proliferating = "deep-proliferating"),
  "GSE109700"
)

# GSE191055: submitted DESeq2 table for HDF P27 (n=4) versus P4 (n=4).
gse191_raw <- as.data.frame(read_excel(file.path(root, "senescence", "GSE191055_gene.P27_vs_P4.xlsx")))
gse191 <- data.frame(
  dataset = "GSE191055", contrast = "P27_vs_P4", gene = toupper(gse191_raw$external_gene_name),
  logFC = gse191_raw$log2FoldChange, mean_expression = gse191_raw$baseMean,
  statistic = gse191_raw$log2FoldChange / gse191_raw$lfcSE,
  PValue = gse191_raw$pvalue, FDR = gse191_raw$padj
)
gse191 <- gse191[!is.na(gse191$gene) & nzchar(gse191$gene), ]
gse191 <- gse191[order(gse191$PValue), ]
gse191 <- gse191[!duplicated(gse191$gene), ]
write_tsv(gse191, "GSE191055_P27_vs_P4_DE.tsv")

# GSE93535: extract the six side-by-side Cuffdiff blocks from publication Table S1.
gse935_sheet <- as.data.frame(read_excel(
  file.path(root, "extracted", "PMC5895844", "41514_2018_23_MOESM4_ESM.xlsx"),
  sheet = "comparing sign. diff. genes", col_names = FALSE
))
extract_gse935_block <- function(start_col, contrast_name) {
  x <- data.frame(
    gene = toupper(as.character(gse935_sheet[-c(1, 2), start_col])),
    sample_1 = as.character(gse935_sheet[-c(1, 2), start_col + 2]),
    sample_2 = as.character(gse935_sheet[-c(1, 2), start_col + 3]),
    logFC = suppressWarnings(as.numeric(gse935_sheet[-c(1, 2), start_col + 6])),
    statistic = suppressWarnings(as.numeric(gse935_sheet[-c(1, 2), start_col + 7])),
    PValue = suppressWarnings(as.numeric(gse935_sheet[-c(1, 2), start_col + 8])),
    FDR = suppressWarnings(as.numeric(gse935_sheet[-c(1, 2), start_col + 9]))
  )
  x <- x[!is.na(x$gene) & nzchar(x$gene) & !grepl(",", x$gene) & !is.na(x$logFC), ]
  x$dataset <- "GSE93535"
  x$contrast <- contrast_name
  x
}
combine_gse935 <- function(starts, contrast_name) {
  x <- do.call(rbind, lapply(starts, extract_gse935_block, contrast_name = contrast_name))
  x <- x[order(x$PValue), ]
  x[!duplicated(x$gene), c("dataset", "contrast", "gene", "logFC", "statistic", "PValue", "FDR", "sample_1", "sample_2")]
}
gse935 <- list(
  SIPS_vs_Q = combine_gse935(c(1, 13), "SIPS_vs_Q"),
  Q1201_vs_Q = combine_gse935(c(25, 37), "Q1201_vs_Q"),
  SIPS1201_vs_SIPS = combine_gse935(c(49, 61), "SIPS1201_vs_SIPS")
)
for (name in names(gse935)) write_tsv(gse935[[name]], paste0("GSE93535_", name, "_DE.tsv"))

# GSE226189: 82 healthy donors, age 22-89, raw gene counts.
gse226_files <- sort(list.files(
  file.path(root, "extracted", "GSE226189"), pattern = "geneCOUNT\\.txt\\.gz$", full.names = TRUE
))
gse226_tables <- lapply(gse226_files, function(path) read.delim(gzfile(path), check.names = FALSE))
gse226_ids <- gse226_tables[[1]][[1]]
stopifnot(all(vapply(gse226_tables, function(x) identical(x[[1]], gse226_ids), logical(1))))
gse226_counts <- do.call(cbind, lapply(gse226_tables, function(x) x[[2]]))
gse226_accessions <- sub("_.*", "", basename(gse226_files))
colnames(gse226_counts) <- gse226_accessions
gse226_counts <- collapse_matrix(gse226_counts, map_ids(gse226_ids, "ENSEMBL"), "sum")
gse226_meta <- meta[match(gse226_accessions, meta$accession), ]
stopifnot(!anyNA(gse226_meta$accession))
gse226_sample_data <- data.frame(age = as.numeric(gse226_meta$age_years), sex = gse226_meta$sex)
gse226 <- run_edger_continuous(gse226_counts, gse226_sample_data, "GSE226189")

# GSE113957: healthy donors only; HGPS is retained outside the primary age model.
gse113_raw <- read.delim(gzfile(file.path(root, "aging", "GSE113957_fpkm.txt.gz")), check.names = FALSE)
gse113_sample_names <- colnames(gse113_raw)[9:ncol(gse113_raw)]
gse113_meta <- meta[match(gse113_sample_names, meta$title), ]
stopifnot(!anyNA(gse113_meta$accession))
gse113_symbol <- sub("\\|.*", "", gse113_raw[["Annotation/Divergence"]])
gse113_fpkm <- as.matrix(gse113_raw[, 9:ncol(gse113_raw)])
storage.mode(gse113_fpkm) <- "numeric"
colnames(gse113_fpkm) <- gse113_sample_names
gse113_fpkm <- collapse_matrix(gse113_fpkm, gse113_symbol, "sum")
healthy <- tolower(gse113_meta$disease) == "normal"
gse113_ages <- as.numeric(gse113_meta$age)
run_gse113_age <- function(sample_keep, contrast_label, output_name) {
  expression <- log2(gse113_fpkm[, sample_keep, drop = FALSE] + 0.1)
  sample_data <- data.frame(
    age_decade = gse113_ages[sample_keep] / 10,
    sex = relevel(factor(tolower(gse113_meta$sex[sample_keep])), ref = "male"),
    platform = factor(gse113_meta$platform_id[sample_keep])
  )
  design <- model.matrix(~ age_decade + sex + platform, data = sample_data)
  fit <- eBayes(lmFit(expression, design), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = "age_decade", number = Inf, sort.by = "none")
  ans <- data.frame(
    dataset = "GSE113957", contrast = contrast_label,
    gene = rownames(tab), logFC = tab$logFC, AveExpr = tab$AveExpr,
    statistic = tab$t, PValue = tab$P.Value, FDR = tab$adj.P.Val
  )
  write_tsv(ans[order(ans$PValue), ], output_name)
  ans
}
adult_overlap <- healthy & gse113_ages >= 22 & gse113_ages <= 89
gse113 <- run_gse113_age(
  adult_overlap,
  "per_10_years_age22_89_adjusted_for_sex_and_platform",
  "GSE113957_age22_89_per_decade_DE.tsv"
)
gse113_all_age <- run_gse113_age(
  healthy,
  "per_10_years_age1_96_developmental_sensitivity",
  "GSE113957_age1_96_sensitivity_DE.tsv"
)

# GSE179848: donor is the biological unit. For each of four healthy cell lines
# under untreated 21% oxygen culture, average the last two minus first two
# longitudinal time points, then test the four independent donor differences.
gse179_raw <- read.csv(
  gzfile(file.path(root, "aging", "GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz")),
  check.names = FALSE
)
gse179_ids <- gse179_raw[[1]]
gse179_counts <- as.matrix(gse179_raw[, -1])
storage.mode(gse179_counts) <- "numeric"
gse179_symbols <- map_ids(gse179_ids, "REFSEQ")
gse179_counts <- collapse_matrix(gse179_counts, gse179_symbols, "sum")
gse179_sample_names <- colnames(gse179_counts)
gse179_meta <- meta[meta$dataset == "GSE179848", ]
gse179_meta$matrix_name <- paste0("Sample_", gse179_meta$rnaseq_sampleid)
gse179_meta <- gse179_meta[match(gse179_sample_names, gse179_meta$matrix_name), ]
stopifnot(!anyNA(gse179_meta$accession))
gse179_keep <- gse179_meta$clinical_condition == "Normal" &
  gse179_meta$treatments == "Control" &
  as.numeric(gse179_meta$percent_oxygen) == 21 &
  gse179_meta$cell_line %in% paste0("HC", 1:4)
gse179_y <- DGEList(gse179_counts[, gse179_keep, drop = FALSE])
gse179_filter <- filterByExpr(gse179_y, group = factor(gse179_meta$cell_line[gse179_keep]))
gse179_y <- calcNormFactors(gse179_y[gse179_filter, , keep.lib.sizes = FALSE])
gse179_logcpm <- cpm(gse179_y, log = TRUE, prior.count = 0.5)
gse179_submeta <- gse179_meta[gse179_keep, ]
donors <- paste0("HC", 1:4)
gse179_delta <- sapply(donors, function(donor) {
  idx <- which(gse179_submeta$cell_line == donor)
  idx <- idx[order(as.numeric(gse179_submeta$days_grown_udays[idx]))]
  if (length(idx) < 4) stop("Fewer than four longitudinal samples for ", donor)
  rowMeans(gse179_logcpm[, tail(idx, 2), drop = FALSE]) -
    rowMeans(gse179_logcpm[, head(idx, 2), drop = FALSE])
})
gse179_design <- matrix(1, nrow = ncol(gse179_delta), ncol = 1,
                        dimnames = list(colnames(gse179_delta), "late_minus_early"))
gse179_fit <- eBayes(lmFit(gse179_delta, gse179_design), robust = TRUE)
gse179_tab <- topTable(gse179_fit, coef = "late_minus_early", number = Inf, sort.by = "none")
gse179 <- data.frame(
  dataset = "GSE179848", contrast = "late2_minus_early2_within_donor_n4",
  gene = rownames(gse179_tab), logFC = gse179_tab$logFC,
  donor_SD = apply(gse179_delta, 1, sd), statistic = gse179_tab$t,
  PValue = gse179_tab$P.Value, FDR = gse179_tab$adj.P.Val,
  direction_donors = apply(gse179_delta, 1, function(x) paste(sign(x), collapse = ",")),
  same_direction_donors = apply(gse179_delta, 1, function(x) max(sum(x > 0), sum(x < 0)))
)
write_tsv(gse179[order(gse179$PValue), ], "GSE179848_late_vs_early_donor_level_DE.tsv")
write_tsv(
  data.frame(
    donor = donors,
    n_timepoints = vapply(donors, function(d) sum(gse179_submeta$cell_line == d), integer(1)),
    first_day = vapply(donors, function(d) min(as.numeric(gse179_submeta$days_grown_udays[gse179_submeta$cell_line == d])), numeric(1)),
    last_day = vapply(donors, function(d) max(as.numeric(gse179_submeta$days_grown_udays[gse179_submeta$cell_line == d])), numeric(1))
  ),
  "GSE179848_selected_trajectory_summary.tsv"
)

# Correlations test whether age and senescence can legitimately be called one axis.
effect_correlation <- function(x, y, x_name, y_name) {
  joined <- merge(x[, c("gene", "logFC", "PValue")], y[, c("gene", "logFC", "PValue")],
                  by = "gene", suffixes = c("_x", "_y"))
  zx <- signed_z(joined$logFC_x, joined$PValue_x)
  zy <- signed_z(joined$logFC_y, joined$PValue_y)
  test <- cor.test(zx, zy, method = "spearman", exact = FALSE)
  data.frame(
    comparison_x = x_name, comparison_y = y_name, n_shared = nrow(joined),
    rho = unname(test$estimate), p_value = test$p.value,
    same_direction_fraction = mean(sign(joined$logFC_x) == sign(joined$logFC_y))
  )
}

comparisons <- do.call(rbind, list(
  effect_correlation(gse113, gse226, "GSE113957_age", "GSE226189_age"),
  effect_correlation(gse109$early_vs_proliferating, gse109$deep_vs_proliferating, "GSE109700_early", "GSE109700_deep"),
  effect_correlation(gse109$deep_vs_proliferating, gse191, "GSE109700_deep", "GSE191055_P27"),
  effect_correlation(gse109$deep_vs_proliferating, gse935$SIPS_vs_Q, "GSE109700_deep", "GSE93535_SIPS"),
  effect_correlation(gse191, gse935$SIPS_vs_Q, "GSE191055_P27", "GSE93535_SIPS"),
  effect_correlation(gse113, gse109$deep_vs_proliferating, "GSE113957_age", "GSE109700_deep"),
  effect_correlation(gse226, gse109$deep_vs_proliferating, "GSE226189_age", "GSE109700_deep"),
  effect_correlation(gse113, gse191, "GSE113957_age", "GSE191055_P27"),
  effect_correlation(gse226, gse191, "GSE226189_age", "GSE191055_P27")
  ,effect_correlation(gse179, gse113, "GSE179848_late", "GSE113957_age")
  ,effect_correlation(gse179, gse226, "GSE179848_late", "GSE226189_age")
  ,effect_correlation(gse179, gse109$deep_vs_proliferating, "GSE179848_late", "GSE109700_deep")
  ,effect_correlation(gse179, gse191, "GSE179848_late", "GSE191055_P27")
))
write_tsv(comparisons, "cross_dataset_effect_correlations.tsv")

# Rescue should oppose the SIPS direction; shared comparator makes this descriptive.
rescue_compare <- effect_correlation(
  gse935$SIPS_vs_Q, gse935$SIPS1201_vs_SIPS,
  "GSE93535_SIPS_vs_Q", "GSE93535_SIPS1201_vs_SIPS"
)
write_tsv(rescue_compare, "GSE93535_injury_rescue_correlation.tsv")

# Reactome pathway profiles for each independent axis.
gmt_connection <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt_lines <- readLines(gmt_connection, warn = FALSE)
close(gmt_connection)
reactome <- lapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) unique(x[-c(1, 2)]))
names(reactome) <- vapply(strsplit(gmt_lines, "\t", fixed = TRUE), function(x) paste(x[1], x[2], sep = "__"), character(1))

run_fgsea <- function(tab, label) {
  ranks <- signed_z(tab$logFC, tab$PValue)
  names(ranks) <- tab$gene
  ranks <- sort(ranks[is.finite(ranks) & !duplicated(names(ranks))], decreasing = TRUE)
  set.seed(260910)
  ans <- fgseaMultilevel(
    reactome, ranks, minSize = 15, maxSize = 500, eps = 1e-10,
    BPPARAM = BiocParallel::SerialParam()
  )
  ans$leadingEdge <- vapply(ans$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  ans$dataset_contrast <- label
  ans <- as.data.frame(ans[, c("dataset_contrast", "pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge")])
  write_tsv(ans[order(ans$padj, -abs(ans$NES)), ], paste0(label, "_Reactome_fgsea.tsv"))
  ans
}

axis_tables <- list(
  GSE113957_age = gse113,
  GSE226189_age = gse226,
  GSE109700_early = gse109$early_vs_proliferating,
  GSE109700_deep = gse109$deep_vs_proliferating,
  GSE191055_P27 = gse191,
  GSE93535_SIPS = gse935$SIPS_vs_Q,
  GSE93535_rescue = gse935$SIPS1201_vs_SIPS,
  GSE179848_late = gse179
)
fg_results <- lapply(names(axis_tables), function(name) run_fgsea(axis_tables[[name]], name))
names(fg_results) <- names(axis_tables)

pathway_correlation <- function(x, y, x_name, y_name) {
  joined <- merge(x[, c("pathway", "NES")], y[, c("pathway", "NES")], by = "pathway", suffixes = c("_x", "_y"))
  test <- cor.test(joined$NES_x, joined$NES_y, method = "spearman", exact = FALSE)
  data.frame(comparison_x = x_name, comparison_y = y_name, n_shared = nrow(joined),
             rho = unname(test$estimate), p_value = test$p.value)
}
pathway_pairs <- list(
  c("GSE113957_age", "GSE226189_age"),
  c("GSE109700_deep", "GSE191055_P27"),
  c("GSE109700_deep", "GSE93535_SIPS"),
  c("GSE191055_P27", "GSE93535_SIPS"),
  c("GSE113957_age", "GSE109700_deep"),
  c("GSE226189_age", "GSE109700_deep"),
  c("GSE93535_SIPS", "GSE93535_rescue")
  ,c("GSE179848_late", "GSE113957_age")
  ,c("GSE179848_late", "GSE226189_age")
  ,c("GSE179848_late", "GSE109700_deep")
  ,c("GSE179848_late", "GSE191055_P27")
)
pathway_comparisons <- do.call(rbind, lapply(pathway_pairs, function(pair) {
  pathway_correlation(fg_results[[pair[1]]], fg_results[[pair[2]]], pair[1], pair[2])
}))
write_tsv(pathway_comparisons, "cross_dataset_Reactome_NES_correlations.tsv")

summary <- data.frame(
  dataset_contrast = names(axis_tables),
  tested_genes = vapply(axis_tables, nrow, integer(1)),
  FDR05_genes = vapply(axis_tables, function(x) sum(x$FDR < 0.05, na.rm = TRUE), integer(1))
)
write_tsv(summary, "analysis_summary.tsv")

p <- ggplot(comparisons, aes(reorder(paste(comparison_x, comparison_y, sep = " vs "), rho), rho)) +
  geom_col(fill = "#4D4D4D") + geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() + labs(x = NULL, y = "Spearman rho of signed gene statistics") +
  theme_classic(base_size = 10)
ggsave(file.path(out, "cross_dataset_gene_concordance.png"), p, width = 9, height = 6, dpi = 180)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("Figure 2 public aging/senescence analysis completed:", out, "\n")
