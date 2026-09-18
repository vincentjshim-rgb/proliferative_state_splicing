## Scratch only: rebuild the ageing core with an explicit list of series, using exactly the
## rule of run_conserved_core.R + run_splicing_secretome_link.R. Writes nothing into the project.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(fgsea)})
D <- "public_data_tierA/derived"; OUT <- commandArgs(TRUE)[1]
Z <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); M <- Z$meta
N  <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D,"compendium/signatures"), full.names=TRUE)) {
  d <- read.delim(f); vec[[sub("\\.tsv$","",basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
lg <- list(GSE179848_late="figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
 GSE113957_age="figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
 GSE226189_age="figure2_public_aging/GSE226189_age_per_decade_DE.tsv")
for (id in names(lg)) { d <- read.delim(file.path(D, lg[[id]])); vec[[id]] <- setNames(d$logFC, toupper(d$gene)) }
## same per-study averaging and MAD scaling as study_vec() in run_conserved_core.R
study_vec_ids <- function(ids) { out <- list()
  for (ds in unique(M$dataset[match(ids, M$id)])) {
    ii <- ids[M$dataset[match(ids, M$id)] == ds]; vs <- vec[ii]; vs <- vs[!sapply(vs, is.null)]
    gg <- Reduce(union, lapply(vs, names)); mat <- sapply(vs, function(v) v[gg]); if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
    rownames(mat) <- gg; v <- rowMeans(mat, na.rm = TRUE)
    sc <- mad(v, na.rm = TRUE); if (!is.finite(sc) || sc == 0) sc <- sd(v, na.rm = TRUE); out[[ds]] <- v / sc }
  out }
td <- tempfile(); dir.create(td); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
spl <- unique(unlist(P[grep("mRNA Splicing|Processing of Capped Intron|3'-end processing", names(P))]))
capped <- P[[grep("^Processing of Capped Intron-Containing Pre-mRNA", names(P))[1]]]
CORE96 <- readLines(file.path(D, "conserved_core/age_down_splicing_core.txt"))
build <- function(ids, label) {
  sv <- study_vec_ids(ids); ns <- length(sv); min_s <- min(4, ns)
  gg <- Reduce(union, lapply(sv, names)); mat <- sapply(sv, function(v) v[gg]); rownames(mat) <- gg
  nm <- rowSums(is.finite(mat)); mat <- mat[nm >= max(min_s, ceiling(0.75 * ns)), , drop = FALSE]
  n_pos <- rowSums(mat > 0, na.rm = TRUE); n_neg <- rowSums(mat < 0, na.rm = TRUE); nm <- rowSums(is.finite(mat))
  dn <- rownames(mat)[pmax(n_pos, n_neg) == nm & n_pos < n_neg]
  core <- intersect(dn, spl); uni <- rownames(mat)
  ## enrichment of the capped-intron pathway among age-down genes, and its rank among all pathways
  hyp <- function(set) { s <- intersect(set, uni); k <- length(intersect(dn, s))
    c(k = k, size = length(s), p = phyper(k - 1, length(s), length(uni) - length(s), length(dn), lower.tail = FALSE)) }
  PP <- P[sapply(P, function(s) length(intersect(s, uni)) >= 20 & length(intersect(s, uni)) <= 500)]
  pv <- sapply(PP, function(s) hyp(s)["p"]); fdr <- p.adjust(pv, "BH"); names(fdr) <- sub("\\.p$", "", names(fdr))
  cap <- hyp(capped); cap_fdr <- fdr[grep("^Processing of Capped Intron", names(fdr))[1]]
  top <- head(sort(fdr), 3)
  cat(sprintf("\n== %s ==\nseries: %s\nstudies: %d | genes tested: %d | age-down unanimous: %d | expected by chance: %.0f\n",
      label, paste(ids, collapse = ", "), ns, length(uni), length(dn), sum(0.5^nm)))
  cat(sprintf("splicing core: %d genes | shared with the published 96: %d | new: %d | lost: %d\n",
      length(core), length(intersect(core, CORE96)), length(setdiff(core, CORE96)), length(setdiff(CORE96, core))))
  cat(sprintf("capped-intron pre-mRNA processing: %d of %d, BH-FDR %.1e, rank %d of %d pathways\n",
      cap["k"], cap["size"], cap_fdr, match(names(cap_fdr), names(sort(fdr))), length(fdr)))
  cat("top 3 pathways:", paste(sprintf("%s (%.0e)", substr(sub("__.*", "", names(top)), 1, 45), top), collapse = " | "), "\n")
  writeLines(core, file.path(OUT, paste0("core_", label, ".txt")))
  invisible(core) }
AGE4 <- c("GSE307377_age", "GSE179848_late", "GSE113957_age", "GSE226189_age")
c4 <- build(AGE4, "published_four_series")
cat(sprintf("\nREPRODUCES THE PUBLISHED 96: %s\n", setequal(c4, CORE96)))
build(setdiff(AGE4, "GSE226189_age"), "without_GSE226189")
build(setdiff(AGE4, "GSE179848_late"), "donor_cohorts_only")
