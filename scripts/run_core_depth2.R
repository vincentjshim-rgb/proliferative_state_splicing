.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(edgeR))
D <- "public_data_tierA/derived"; O <- file.path(D, "core_depth")
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core.txt"))
modscore <- function(mat, set) { mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop=FALSE]
  mat <- mat[apply(mat, 1, sd) > 0, , drop=FALSE]; z <- t(scale(t(mat)))
  g <- intersect(rownames(z), set); list(score = colMeans(z[g,,drop=FALSE]), n = length(g)) }

## ---- is GSE226189's batch confounded with age? -----------------------------
A  <- read.delim(file.path(O, "donor_age_core_scores.tsv")); d <- A[A$study == "GSE226189", ]
ff <- list.files(file.path(O, "GSE226189"), pattern = "geneCOUNT")
bn <- ff[match(d$sample, sub("^(GSM[0-9]+)_.*$", "\\1", ff))]
d$grp <- sub("^.*_[MF]_(GR[0-9]+)_.*$", "\\1", bn)
say("GSE226189 batch vs age: Kruskal-Wallis P = %.3f", kruskal.test(age ~ factor(grp), d)$p.value)
print(aggregate(age ~ grp, d, function(x) c(n = length(x), median = median(x), range = diff(range(x)))))

## ---- GSE179848: 345 samples, longitudinal, control condition only ----------
say("\n=== GSE179848 longitudinal ===")
cnt <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                row.names = 1, check.names = FALSE)
ln <- readLines(gzfile("public_data_tierA/metadata/GSE179848_family.soft.gz"))
gsm <- ttl <- character(0); ch <- list(); cur <- NULL; acc <- character(0)
for (l in ln) {
  if (grepl("^\\^SAMPLE", l)) { if (!is.null(cur)) ch[[cur]] <- acc
    cur <- sub(".*= ", "", l); acc <- character(0); gsm <- c(gsm, cur) }
  if (grepl("^!Sample_title", l)) ttl <- c(ttl, sub(".*= ", "", l))
  if (grepl("^!Sample_characteristics_ch1", l)) acc <- c(acc, sub(".*= ", "", l)) }
if (!is.null(cur)) ch[[cur]] <- acc
getf <- function(v, k) { i <- grep(paste0("^", k, ":"), v); if (!length(i)) NA else sub(".*: ", "", v[i[1]]) }
meta <- data.frame(gsm = gsm, title = ttl,
  line   = sapply(ch[gsm], getf, "cell_line"),
  cond   = sapply(ch[gsm], getf, "treatments"),
  days   = suppressWarnings(as.numeric(sapply(ch[gsm], getf, "days_grown"))),
  pass   = suppressWarnings(as.numeric(sapply(ch[gsm], getf, "passage"))),
  oxy    = sapply(ch[gsm], getf, "percent_oxygen"),
  clin   = sapply(ch[gsm], getf, "clinical_condition"),
  sid    = sapply(ch[gsm], getf, "study_sample_id"), stringsAsFactors = FALSE)
say("metadata rows %d | fields with values: days %d, passage %d, condition %d",
    nrow(meta), sum(is.finite(meta$days)), sum(is.finite(meta$pass)), sum(!is.na(meta$cond)))
print(head(sort(table(meta$cond), decreasing = TRUE), 6))
## map count columns (Sample_N) to GSM via study_sample_id order in the title
say("count matrix: %d genes x %d samples; first cols: %s", nrow(cnt), ncol(cnt),
    paste(head(colnames(cnt), 3), collapse = ", "))
saveRDS(list(meta = meta), file.path(O, "GSE179848_meta.rds"))
