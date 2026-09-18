#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Phase 2, H1b : is splicing OUTCOME an axis of donor age that survives
#                adjustment for proliferative state?
#   GSE113957 / SRP144355, 143 human dermal fibroblast donors, age 1-96.
#   recount3 junction counts + recount3 gene sums (for the proliferation score).
# Added to the Phase 2 preregistration as H1b when recount3 coverage was
# confirmed (see repro_cm_preregistration_phase2_ko.md, section 1).
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- "public_data_tierA"; rc <- file.path(root, "recount3")
out  <- file.path(root, "derived", "phase2_splicing")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
SRP <- "SRP144355"

## ---- sample map -----------------------------------------------------------
soft <- readLines(gzfile(file.path(root, "metadata", "GSE113957_family.soft.gz")))
blk <- split(soft, cumsum(grepl("^\\^SAMPLE", soft)))[-1]
gmap <- do.call(rbind, lapply(blk, function(b) data.frame(
  gsm = sub("^\\^SAMPLE = ", "", b[1]),
  title = sub(".*= ", "", grep("!Sample_title", b, value = TRUE)[1]),
  srx = sub(".*term=", "", grep("!Sample_relation = SRA", b, value = TRUE)[1]))))
md <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", SRP))),
            select = c("rail_id", "external_id", "experiment_acc"), data.table = FALSE)
gmap$rail_id <- md$rail_id[match(gmap$srx, md$experiment_acc)]
gmap$srr    <- md$external_id[match(gmap$srx, md$experiment_acc)]   # gene_sums columns are SRR
cl <- tolower(gmap$title)
gmap$age <- as.numeric(sub("^[^_]+_([0-9]+)(yr|ys)?([0-9]*mos)?_.*$", "\\1", cl))
gmap$sex <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", cl)
gmap$sex[gmap$sex == "f"] <- "female"; gmap$sex[gmap$sex == "m"] <- "male"
gmap$sex[!gmap$sex %in% c("male", "female")] <- NA
gmap$hgps <- grepl("hgps|progeria", cl)
use <- !gmap$hgps & !is.na(gmap$age) & !is.na(gmap$sex) & gmap$age >= 22 & gmap$age <= 96 &
       !is.na(gmap$rail_id)
cat("total GEO samples:", nrow(gmap), " | healthy adults 22-96 with rail_id:", sum(use), "\n")
G <- gmap[use, ]

## ---- junction matrix ------------------------------------------------------
ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))),
            what = character(), quiet = TRUE); ids <- ids[ids != "rail_id"]
RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))),
            data.table = FALSE)
trip <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)),
                          "| tail -n +4"), header = FALSE, data.table = FALSE)
cat("junctions:", nrow(RR), " samples:", length(ids), " nonzero:", nrow(trip), "\n")
J <- matrix(0L, nrow = nrow(RR), ncol = length(ids), dimnames = list(NULL, ids))
J[cbind(trip[[1]], trip[[2]])] <- as.integer(trip[[3]])
rm(trip); gc()
J <- J[, match(as.character(G$rail_id), colnames(J)), drop = FALSE]
colnames(J) <- G$gsm
depth <- colSums(J)
cat("junction reads per donor: median", round(median(depth)/1e6, 2), "M, range",
    paste(round(range(depth)/1e6, 2), collapse = "-"), "M\n")

## ---- proliferation score from recount3 gene sums --------------------------
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)),
                        "| grep -v '^##'"), data.table = FALSE)
rownames(gs) <- gs[[1]]; gs <- as.matrix(gs[, -1])
stopifnot(all(G$srr %in% colnames(gs)))
gs <- gs[, match(G$srr, colnames(gs)), drop = FALSE]
cat("gene_sums matched for", ncol(gs), "donors\n")
con <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt <- readLines(con, warn = FALSE); close(con); sp <- strsplit(gmt, "\t", fixed = TRUE)
sets <- lapply(sp, function(x) unique(toupper(x[-c(1, 2)])))
names(sets) <- vapply(sp, function(x) x[1], character(1))
prolif_genes <- unique(unlist(sets[c("Cell Cycle, Mitotic", "DNA Replication")]))
suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
sy <- suppressMessages(mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", rownames(gs)),
                              keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
ok <- !is.na(sy) & nzchar(sy)
gsx <- rowsum(gs[ok, ], sy[ok])
lcpm <- cpm(calcNormFactors(DGEList(gsx)), log = TRUE, prior.count = 1)
Z <- t(scale(t(lcpm[rowMeans(lcpm) > 0, ])))
P <- colMeans(Z[rownames(Z) %in% prolif_genes, , drop = FALSE])
cat("proliferation meta-gene from", sum(rownames(Z) %in% prolif_genes), "genes; sd =",
    round(sd(P), 3), "; cor(prolif, age) =", round(cor(P, G$age), 3), "\n")

## ---- sample-level splicing metrics ---------------------------------------
ann <- RR$annotated == 1
metrics <- data.frame(
  gsm = G$gsm, age = G$age, sex = G$sex, prolif = P,
  log_depth = log10(depth),
  unannotated_read_fraction = colSums(J[!ann, ]) / depth,
  unannotated_junction_fraction = colSums(J[!ann, ] > 0) / colSums(J > 0))
w(metrics, "H1b_sample_metrics.tsv")
mod <- function(y, adj) {
  f <- if (adj) lm(y ~ I(age/10) + prolif + sex + log_depth, data = metrics)
       else      lm(y ~ I(age/10) + sex + log_depth, data = metrics)
  c(summary(f)$coef["I(age/10)", c(1, 4)])
}
sm <- do.call(rbind, lapply(c("unannotated_read_fraction", "unannotated_junction_fraction"),
  function(v) { a <- mod(metrics[[v]], FALSE); b <- mod(metrics[[v]], TRUE)
    data.frame(metric = v, r_with_prolif = cor(metrics[[v]], metrics$prolif),
               beta_per_decade = a[1], p = a[2],
               beta_prolif_adjusted = b[1], p_adjusted = b[2]) }))
cat("\n=== sample-level splicing-noise metrics vs donor age (n =", nrow(metrics), ") ===\n")
print(sm, row.names = FALSE, digits = 3)
w(sm, "H1b_sample_metric_age_models.tsv")

## ---- per-junction PSI vs age ---------------------------------------------
psi_age <- function(anchor_cols, label, mincov = 20, min_frac = 0.9) {
  key <- do.call(paste, c(RR[anchor_cols], sep = "_"))
  nalt <- ave(rep(1, nrow(RR)), key, FUN = sum)
  sel <- which(nalt >= 2)
  Js <- J[sel, , drop = FALSE]; ks <- key[sel]
  tot <- apply(Js, 2, function(v) ave(v, ks, FUN = sum))
  ok <- rowMeans(tot >= mincov) >= min_frac
  Js <- Js[ok, , drop = FALSE]; tot <- tot[ok, , drop = FALSE]
  PSI <- Js / pmax(tot, 1)
  keep <- apply(PSI, 1, function(r) !(min(r) > 0.98)) & apply(PSI, 1, sd) > 0
  PSI <- PSI[keep, , drop = FALSE]
  lg <- log((PSI + 0.01) / (1 - PSI + 0.01))
  fit <- function(adj) {
    d <- if (adj) model.matrix(~ I(age/10) + prolif + sex + log_depth, data = metrics)
         else      model.matrix(~ I(age/10) + sex + log_depth, data = metrics)
    e <- eBayes(lmFit(lg, d))
    topTable(e, coef = "I(age/10)", number = Inf, sort.by = "none")
  }
  A <- fit(FALSE); B <- fit(TRUE)
  list(summary = data.frame(
      analysis = label, n_junctions = nrow(lg),
      n_FDR05_unadjusted = sum(A$adj.P.Val < 0.05),
      n_FDR05_prolif_adjusted = sum(B$adj.P.Val < 0.05),
      median_abs_logFC_per_decade = median(abs(A$logFC)),
      cor_t_unadj_vs_adj = cor(A$t, B$t, method = "spearman")),
    A = A, B = B)
}
R1 <- psi_age(c("chromosome", "start", "strand"), "5'-anchored groups")
R2 <- psi_age(c("chromosome", "end", "strand"),   "3'-anchored groups")
SS <- rbind(R1$summary, R2$summary)
cat("\n=== per-junction PSI vs donor age (limma, n =", nrow(metrics), "donors) ===\n")
print(SS, row.names = FALSE, digits = 3)
w(SS, "H1b_junction_psi_age_summary.tsv")
w(data.frame(analysis = "5prime", logFC = R1$A$logFC, t = R1$A$t, FDR = R1$A$adj.P.Val,
             t_prolif_adj = R1$B$t, FDR_prolif_adj = R1$B$adj.P.Val),
  "H1b_junction_psi_age_5prime.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_H1b.txt"))
