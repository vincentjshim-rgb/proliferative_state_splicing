#!/usr/bin/env Rscript

project_lib <- normalizePath("analysis_r_lib", mustWork = FALSE)
dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_lib, .libPaths()))

packages <- c(
  "limma",
  "edgeR",
  "fgsea",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "metafor",
  "readxl",
  "data.table"
)

missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  BiocManager::install(missing, lib = project_lib, ask = FALSE, update = FALSE)
}

status <- vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
write.table(
  data.frame(package = packages, installed = unname(status)),
  file = "public_data_tierA/logs/r_package_status.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(status)) {
  stop("One or more required packages were not installed; see r_package_status.tsv")
}

writeLines(capture.output(sessionInfo()), "public_data_tierA/logs/R_sessionInfo.txt")
cat("All required R packages are available in", project_lib, "\n")
