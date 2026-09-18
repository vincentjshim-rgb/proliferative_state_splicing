#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Fibroblast perturbation compendium.
# Every public fibroblast contrast we can build from data already on disk is
# scored on the same two reference axes (UVA meta, senescence meta) and on the
# decoupling index D.  The Repro-CM anchor becomes one point on that map.
# Signatures that contributed to a reference axis are flagged as circular.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi); library(readxl)
})
root <- "public_data_tierA"
der  <- file.path(root, "derived")
out  <- file.path(der, "compendium")
dir.create(file.path(out, "signatures"), recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
sig_path <- function(id) file.path(out, "signatures", paste0(id, ".tsv"))
save_sig <- function(id, dataset, contrast, gene, logFC) {
  d <- data.frame(dataset = dataset, contrast = contrast, gene = toupper(gene), logFC = logFC)
  d <- d[!is.na(d$gene) & nzchar(d$gene) & is.finite(d$logFC), ]
  d <- d[!duplicated(d$gene), ]
  write.table(d, sig_path(id), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(nrow(d))
}
map_sym <- function(ids, kt) {
  s <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", ids),
                               keytype = kt, column = "SYMBOL", multiVals = "first"))
  unname(s)
}

## ===========================================================================
## PART A - GSE179848: nine treatments, hypoxia, and a mitochondrial mutation
## ===========================================================================
cat("== PART A: GSE179848 perturbations ==\n")
meta <- read.delim(file.path(der, "audit", "all_geo_samples.tsv"), check.names = FALSE)
gm <- meta[meta$dataset == "GSE179848", ]
raw <- read.csv(gzfile(file.path(root, "aging",
        "GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz")), check.names = FALSE)
sy <- map_sym(raw[[1]], "REFSEQ")
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy)
cnt <- rowsum(cnt[ok, ], sy[ok])
gm$matrix_name <- paste0("Sample_", gm$rnaseq_sampleid)
gm <- gm[match(colnames(cnt), gm$matrix_name), ]
gm$days <- as.numeric(gm$days_grown_udays)
gm$ox   <- as.character(gm$percent_oxygen)
stopifnot(!anyNA(gm$accession))

de_179 <- function(keep, grp, ref, label, id, use_cellline = TRUE) {
  if (file.exists(sig_path(id))) { cat(sprintf("  %-34s cached\n", label))
    return(data.frame(id = id, dataset = "GSE179848", contrast = label,
                      n_samples = sum(keep), class = NA)) }
  m <- gm[keep, ]; y <- DGEList(cnt[, keep])
  g <- factor(grp[keep], levels = c(ref, setdiff(unique(grp[keep]), ref)))
  y <- y[filterByExpr(y, group = g), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  cl <- factor(m$cell_line)
  ## cell line is fully confounded with genotype in the SURF1 comparison
  dsn <- if (use_cellline && nlevels(cl) > 1) model.matrix(~ g + cl + scale(m$days))
         else                                model.matrix(~ g + scale(m$days))
  if (qr(dsn)$rank < ncol(dsn)) dsn <- model.matrix(~ g + scale(m$days))
  y <- estimateDisp(y, dsn)
  tt <- topTags(glmQLFTest(glmQLFit(y, dsn), coef = 2), n = Inf, sort.by = "none")$table
  n <- save_sig(id, "GSE179848", label, rownames(tt), tt$logFC)
  cat(sprintf("  %-34s n_samples=%3d  genes=%d\n", label, sum(keep), n))
  data.frame(id = id, dataset = "GSE179848", contrast = label,
             n_samples = sum(keep), class = NA)
}
normal21 <- gm$clinical_condition == "Normal" & gm$ox == "21"
A <- list()
TRT <- c("DEX", "Contact_Inhibition", "2-Deoxyglucose", "betahydroxybutyrate",
         "Galactose", "mitoNUITs", "mitoNUITs+DEX", "Oligomycin", "Oligomycin+DEX")
for (t in TRT) {
  k <- normal21 & gm$treatments %in% c("Control", t)
  if (sum(gm$treatments[k] == t) >= 5)
    A[[t]] <- de_179(k, gm$treatments, "Control",
                     paste(t, "vs Control (21% O2)"), paste0("GSE179848_", gsub("[^A-Za-z0-9]", "", t)))
}
k <- gm$clinical_condition == "Normal" & gm$treatments == "Control"
A[["hypoxia"]] <- de_179(k, ifelse(gm$ox == "3", "O2_3", "O2_21"), "O2_21",
                         "3% vs 21% oxygen (Control)", "GSE179848_hypoxia")
k <- gm$treatments == "Control" & gm$ox == "21"
A[["SURF1"]] <- de_179(k, gm$clinical_condition, "Normal",
                       "SURF1 mutation vs normal (Control)", "GSE179848_SURF1",
                       use_cellline = FALSE)

## ===========================================================================
## PART B - held-out datasets, effect sizes only (no permutation here)
## ===========================================================================
cat("== PART B: held-out datasets ==\n")
ho <- file.path(root, "heldout"); h2 <- file.path(root, "heldout2")
simple_de <- function(counts, group, ref, covar = NULL) {
  y <- DGEList(counts); g <- factor(group, levels = c(ref, setdiff(unique(group), ref)))
  y <- y[filterByExpr(y, group = g), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  d <- if (is.null(covar)) model.matrix(~ g) else model.matrix(~ g + factor(covar))
  y <- estimateDisp(y, d)
  tt <- topTags(glmQLFTest(glmQLFit(y, d), coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(tt), logFC = tt$logFC)
}
B <- list()
add <- function(id, dataset, label, eff, n_samples) {
  if (!file.exists(sig_path(id))) save_sig(id, dataset, label, eff$gene, eff$logFC)
  cat(sprintf("  %-40s genes=%d\n", label, nrow(eff)))
  B[[id]] <<- data.frame(id = id, dataset = dataset, contrast = label,
                         n_samples = n_samples, class = NA)
}
## GSE297233 OSK
m297 <- read.csv(gzfile(file.path(ho, "GSE297233_raw_counts_matrix.csv.gz")), row.names = 1, check.names = FALSE)
s297 <- map_sym(rownames(m297), "ENSEMBL"); o <- !is.na(s297) & nzchar(s297)
mm297 <- rowsum(as.matrix(m297[o, ]), s297[o])
for (arm in c("OSK", "O4YRSK")) {
  cc <- mm297[, grep(paste0("^", arm, "_"), colnames(mm297))]
  add(paste0("GSE297233_", arm), "GSE297233", paste(arm, "+dox vs -dox, day 4"),
      simple_de(cc, ifelse(grepl("_D4_", colnames(cc)), "dox", "nodox"), "nodox"), ncol(cc))
}
## GSE307377 old vs young
g7 <- read.delim(gzfile(file.path(ho, "GSE307377_raw.txt.gz")), check.names = FALSE)
s7 <- toupper(sub("\\|.*$", "", g7[["Annotation/Divergence"]]))
c7 <- as.matrix(g7[, grep("^tags/", names(g7))]); storage.mode(c7) <- "numeric"
colnames(c7) <- sub("^tags/", "", sub(" .*$", "", colnames(c7)))
k7 <- !is.na(s7) & nzchar(s7) & s7 != "NA"; c7 <- round(rowsum(c7[k7, ], s7[k7]))
add("GSE307377_age", "GSE307377", "old vs young dermal fibroblast",
    simple_de(c7, ifelse(grepl("_O_", colnames(c7)), "old", "young"), "young"), ncol(c7))
## GSE116968 red light / UV
sfc <- function(v) { v <- suppressWarnings(as.numeric(v))
  ifelse(is.na(v), NA, ifelse(v >= 0, log2(pmax(v, 1e-6)), -log2(pmax(-v, 1e-6)))) }
for (fi in 1:2) {
  f <- file.path(ho, c("GSE116968_Processed_data_NHDF_RNAseq_1_post_1h_.xlsx",
                       "GSE116968_Processed_data_NHDF_RNAseq_2_Post_4h_.xlsx")[fi])
  tp <- c("1h", "4h")[fi]
  x <- as.data.frame(read_excel(f, sheet = 1)); nm <- names(x)
  gset <- toupper(x$Gene_Symbol)
  cols <- list(red_vs_control = grep("^Red/Control\\.fc$", nm, value = TRUE)[1],
               uv_vs_control  = grep("^UV.*/Control\\.fc$", nm, value = TRUE)[1],
               prered_vs_uv   = grep("^Pre-Red.*/UV.*\\.fc$", nm, value = TRUE)[1])
  inj <- sfc(x[[cols$uv_vs_control]])
  for (cn in names(cols)) {
    v <- sfc(x[[cols[[cn]]]])
    add(paste0("GSE116968_", cn, "_", tp), "GSE116968", paste(cn, tp),
        data.frame(gene = gset, logFC = v)[is.finite(v), ], 6)
  }
  d <- data.frame(gene = gset, r = sfc(x[[cols$prered_vs_uv]]), i = inj)
  d <- d[is.finite(d$r) & is.finite(d$i), ]
  d$res <- residuals(lm(r ~ i, data = d))
  add(paste0("GSE116968_prered_vs_uv_injuryresid_", tp), "GSE116968",
      paste("Pre-Red vs UV, injury-residualised,", tp),
      data.frame(gene = d$gene, logFC = d$res), 6)
}
## GSE149694 reprogramming media
td <- file.path(der, "heldout_validation", "GSE149694_files")
if (dir.exists(td)) {
  map <- do.call(rbind, lapply(strsplit(system(paste(
    "zcat", file.path(ho, "metadata/GSE149694_family.soft.gz"),
    "| awk '/^\\^SAMPLE/{g=$3} /cell subtype\\/time point:/{sub(/.*point: /,\"\"); print g\"\\t\"$0}'"),
    intern = TRUE), "\t"), function(z) data.frame(gsm = z[1], cond = z[2])))
  fl <- list.files(td, pattern = "\\.txt\\.gz$", full.names = TRUE); names(fl) <- sub("_.*$", "", basename(fl))
  mats <- lapply(fl, function(f) { z <- read.delim(gzfile(f)); setNames(z[[3]], z[[1]]) })
  gg <- Reduce(intersect, lapply(mats, names))
  c4 <- do.call(cbind, lapply(mats, function(v) v[gg]))
  s4 <- map_sym(gg, "ENSEMBL"); o4 <- !is.na(s4) & nzchar(s4)
  c4 <- rowsum(c4[o4, ], s4[o4]); cond <- map$cond[match(colnames(c4), map$gsm)]
  for (tgt in c("Fibroblast-D7", "Primed-D13", "NHSM-D13", "5iLAF-D13", "RSeT-D13", "t2iLGoY-D13")) {
    k <- cond %in% c("Fibroblast-D3", tgt)
    if (sum(cond == tgt) >= 2)
      add(paste0("GSE149694_", gsub("[^A-Za-z0-9]", "", tgt)), "GSE149694",
          paste(tgt, "vs Fibroblast-D3"),
          simple_de(c4[, k], cond[k], "Fibroblast-D3"), sum(k))
  }
}
## GSE306957 senescence arms
f2 <- file.path(h2, "files")
if (dir.exists(f2)) {
  L <- readLines(gzfile(file.path(h2, "GSE306957_family.soft.gz")))
  bl <- split(L, cumsum(grepl("^\\^SAMPLE", L)))[-1]
  md <- do.call(rbind, lapply(bl, function(b) { g1 <- function(p) { z <- grep(p, b, value = TRUE)[1]
        if (is.na(z)) NA else sub(".*= ", "", z) }
    data.frame(title = g1("!Sample_title"), cell = g1("cell line:"),
               file = sub(".*/", "", g1("!Sample_supplementary_file_1"))) }))
  md <- md[!is.na(md$cell) & grepl("MRC5", md$cell), ]
  md$state <- ifelse(grepl("^Proliferating", md$title), "proliferating", "senescent")
  md$arm <- ifelse(grepl("siScramble", md$title), "siScramble",
             ifelse(grepl("siMAVS", md$title), "siMAVS",
             ifelse(grepl("pLenti", md$title), "pLenti",
             ifelse(grepl("IFIH1", md$title), "IFIH1_KO", "DDX58_KO"))))
  mats <- lapply(md$file, function(f) { x <- read.delim(gzfile(file.path(f2, f)), check.names = FALSE)
                                        setNames(x[[ncol(x)]], x$Geneid) })
  gg <- Reduce(intersect, lapply(mats, names))
  M6 <- do.call(cbind, lapply(mats, function(v) v[gg])); colnames(M6) <- md$title
  s6 <- map_sym(gg, "ENSEMBL"); o6 <- !is.na(s6) & nzchar(s6)
  M6 <- rowsum(M6[o6, ], s6[o6])
  for (a in c("pLenti", "siScramble", "IFIH1_KO", "DDX58_KO")) {
    k <- md$arm == a
    if (sum(k & md$state == "senescent") >= 2 && sum(k & md$state == "proliferating") >= 2)
      add(paste0("GSE306957_", a), "GSE306957", paste(a, "senescent vs proliferating"),
          simple_de(M6[, k], md$state[k], "proliferating"), sum(k))
  }
}
w(do.call(rbind, c(A, B)), "compendium_inventory_raw.tsv")
cat("PART A+B done\n")
