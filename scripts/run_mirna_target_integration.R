#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(BiocParallel)
  library(fgsea)
})

root <- "public_data_tierA"
out <- file.path(root, "derived", "local_mirna_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                        row.names = FALSE, na = "")

# Stream the large OmniPath archive and retain only experimentally curated
# miRTarBase edges for the two reprogramming-associated miRNA families.
network_file <- file.path(root, "network", "omnipath_webservice_interactions__latest.tsv.gz")
con <- gzfile(network_file, "rt")
header <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
idx <- setNames(match(c("source_genesymbol", "target_genesymbol", "sources", "references", "mirnatarget"),
                      header), c("source", "target", "sources", "references", "mirnatarget"))
edges <- list()
repeat {
  lines <- readLines(con, n = 50000, warn = FALSE)
  if (!length(lines)) break
  fields <- strsplit(lines, "\t", fixed = TRUE)
  keep <- vapply(fields, function(x) length(x) >= max(idx) && x[idx["mirnatarget"]] == "True" &&
                   grepl("miRTarBase", x[idx["sources"]], fixed = TRUE) &&
                   nzchar(x[idx["references"]]) &&
                   grepl("^hsa-(miR-(302[a-z]?|367|371[ab]?|372|373))$",
                         x[idx["source"]], ignore.case = TRUE), logical(1))
  if (any(keep)) {
    rows <- fields[keep]
    edges[[length(edges) + 1]] <- data.frame(
      miRNA = vapply(rows, `[`, character(1), idx["source"]),
      target = toupper(vapply(rows, `[`, character(1), idx["target"])),
      sources = vapply(rows, `[`, character(1), idx["sources"]),
      references = vapply(rows, `[`, character(1), idx["references"])
    )
  }
}
close(con)
edges <- unique(do.call(rbind, edges))
edges$family <- ifelse(grepl("hsa-miR-(302|367)", edges$miRNA, ignore.case = TRUE),
                       "miR-302/367", "miR-371/372/373")
write_tsv(edges, "OmniPath_miRTarBase_reprogramming_miRNA_edges.tsv")

sets <- split(edges$target, edges$family)
sets <- lapply(sets, unique)
set_summary <- data.frame(family = names(sets), validated_targets = lengths(sets))
write_tsv(set_summary, "validated_target_set_sizes.tsv")

anchor <- read.delim(file.path(root, "derived", "local_repro_anchor",
                               "local_Repro_specific_gene_rank.tsv"), check.names = FALSE)
anchor$gene <- toupper(anchor$gene)
axes <- read.delim(file.path(root, "derived", "decoupling_validation",
                             "primary_gene_axis_matrix.tsv"), check.names = FALSE)
axes$gene <- toupper(axes$gene)

ranks <- list(
  Repro_CM_anchor = setNames(anchor$Repro_specific_score, anchor$gene),
  UVA_meta_axis = setNames(axes$UVA_meta, axes$gene),
  senescence_meta_axis = setNames(axes$senescence_meta, axes$gene)
)
rows <- list()
for (nm in names(ranks)) {
  stat <- ranks[[nm]]
  stat <- sort(stat[is.finite(stat) & stat != 0 & !duplicated(names(stat))], decreasing = TRUE)
  set.seed(260914)
  ans <- fgseaMultilevel(sets, stat, minSize = 10, maxSize = 2000, eps = 0,
                         BPPARAM = SerialParam())
  ans$axis <- nm
  ans$leadingEdge <- vapply(ans$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  rows[[nm]] <- as.data.frame(ans)
}
result <- do.call(rbind, rows)
rownames(result) <- NULL
result$FDR_within_axis <- ave(result$pval, result$axis, FUN = function(x) p.adjust(x, "BH"))
write_tsv(result[, c("axis", "pathway", "size", "ES", "NES", "pval", "FDR_within_axis", "leadingEdge")],
          "validated_miRNA_target_GSEA.tsv")

cat("miRNA target integration completed:", out, "\n")
