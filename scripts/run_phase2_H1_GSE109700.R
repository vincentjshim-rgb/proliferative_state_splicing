#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Phase 2, H1 re-run on a properly powered dataset.
#   GSE109700 / SRP131506, LF1 lung fibroblasts, ribodepleted total RNA,
#   proliferating / early senescent / deep senescent, n = 3 each,
#   19.5-36.6 M junction reads per sample (5-10x GSE93535).
# Moderated t (limma) on logit-PSI, plus a label-permutation empirical null.
# NOTE: this dataset has NO quiescence arm, so it addresses
#   "does senescence change splicing outcome at all", NOT
#   "is that separable from proliferative arrest".
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma) })
root <- "public_data_tierA"; rc <- file.path(root, "recount3")
out  <- file.path(root, "derived", "phase2_splicing")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
SRP <- "SRP131506"

## ---- sample map ------------------------------------------------------------
soft <- readLines(gzfile(file.path(root, "metadata", "GSE109700_family.soft.gz")))
bl <- split(soft, cumsum(grepl("^\\^SAMPLE", soft)))[-1]
gm <- do.call(rbind, lapply(bl, function(b) { g1 <- function(p) { z <- grep(p, b, value = TRUE)[1]
    if (is.na(z)) NA else sub(".*= ", "", z) }
  data.frame(gsm = sub("^\\^SAMPLE = ", "", b[1]), title = g1("!Sample_title"),
             srx = sub(".*term=", "", g1("!Sample_relation = SRA"))) }))
md <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", SRP))),
            select = c("rail_id", "external_id", "experiment_acc"), data.table = FALSE)
gm$rail_id <- md$rail_id[match(gm$srx, md$experiment_acc)]
gm$group <- ifelse(grepl("^Proliferating", gm$title), "proliferating",
            ifelse(grepl("^Early", gm$title), "early", "deep"))
print(table(gm$group))

## ---- junctions -------------------------------------------------------------
ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))),
            what = character(), quiet = TRUE); ids <- ids[ids != "rail_id"]
RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))),
            data.table = FALSE)
tr <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)),
                        "| tail -n +4"), header = FALSE, data.table = FALSE)
S <- matrix(0L, nrow = nrow(RR), ncol = length(ids), dimnames = list(NULL, ids))
S[cbind(tr[[1]], tr[[2]])] <- as.integer(tr[[3]]); rm(tr); gc()
S <- S[, match(as.character(gm$rail_id), colnames(S)), drop = FALSE]
colnames(S) <- gm$title
cat("junction reads per sample (M):\n"); print(round(colSums(S) / 1e6, 2))
cat("annotated junction fraction:", round(mean(RR$annotated == 1), 3), "\n")

## ---- PSI + moderated t -----------------------------------------------------
build <- function(smp, mincov, min_frac, cols) {
  key <- do.call(paste, c(RR[cols], sep = "_"))
  nalt <- ave(rep(1, nrow(RR)), key, FUN = sum); sel <- which(nalt >= 2)
  Js <- S[sel, smp, drop = FALSE]; ks <- key[sel]
  tot <- apply(Js, 2, function(v) ave(v, ks, FUN = sum))
  ok <- rowMeans(tot >= mincov) >= min_frac
  Js <- Js[ok, , drop = FALSE]; tot <- tot[ok, , drop = FALSE]
  P <- Js / pmax(tot, 1)
  keep <- apply(P, 1, function(r) !(min(r) > 0.98)) & apply(P, 1, sd) > 0
  list(psi = P[keep, , drop = FALSE], idx = sel[ok][keep])
}
lg <- function(x) log((x + 0.01) / (1 - x + 0.01))
test <- function(b, ga, gb) {
  L <- lg(b$psi)
  grp <- factor(ifelse(colnames(L) %in% gb, "b", "a"), levels = c("a", "b"))
  tt <- topTable(eBayes(lmFit(L, model.matrix(~ grp))), coef = 2, number = Inf, sort.by = "none")
  d <- rowMeans(b$psi[, gb, drop = FALSE]) - rowMeans(b$psi[, ga, drop = FALSE])
  c(n = nrow(L), FDR05 = sum(tt$adj.P.Val < 0.05),
    FDR05_dPSI10 = sum(tt$adj.P.Val < 0.05 & abs(d) >= 0.10),
    FDR01 = sum(tt$adj.P.Val < 0.01), med_absdPSI = median(abs(d)))
}
PRO <- gm$title[gm$group == "proliferating"]
EAR <- gm$title[gm$group == "early"]
DEE <- gm$title[gm$group == "deep"]

cat("\n############ filter sweep, deep senescence vs proliferating ############\n")
grid <- expand.grid(mincov = c(10, 20, 40, 80), min_frac = c(0.834, 1.0))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  b <- build(c(PRO, DEE), grid$mincov[i], grid$min_frac[i], c("chromosome", "start", "strand"))
  r <- test(b, PRO, DEE)
  data.frame(mincov = grid$mincov[i], min_frac = grid$min_frac[i],
             n_junctions = unname(r["n"]), FDR05 = unname(r["FDR05"]),
             FDR05_dPSI10 = unname(r["FDR05_dPSI10"]), FDR01 = unname(r["FDR01"]),
             med_absdPSI = round(unname(r["med_absdPSI"]), 4))
}))
print(res, row.names = FALSE)
w(res, "H1_GSE109700_filter_sweep.tsv")

MC <- 20; MF <- 1.0
cat("\n############ contrasts at mincov =", MC, ", all samples ############\n")
ctr <- list(`deep vs proliferating` = list(PRO, DEE),
            `early vs proliferating` = list(PRO, EAR),
            `deep vs early (both senescent)` = list(EAR, DEE))
CT <- do.call(rbind, lapply(names(ctr), function(n) {
  b <- build(unlist(ctr[[n]]), MC, MF, c("chromosome", "start", "strand"))
  r <- test(b, ctr[[n]][[1]], ctr[[n]][[2]])
  data.frame(contrast = n, n_junctions = unname(r["n"]), FDR05 = unname(r["FDR05"]),
             FDR05_dPSI10 = unname(r["FDR05_dPSI10"]), FDR01 = unname(r["FDR01"]),
             med_absdPSI = round(unname(r["med_absdPSI"]), 4))
}))
print(CT, row.names = FALSE); w(CT, "H1_GSE109700_contrasts.tsv")

## ---- empirical null: every other 3-vs-3 split of the same six samples -------
cat("\n############ label-permutation null (deep vs proliferating samples) ############\n")
six <- c(PRO, DEE)
b6 <- build(six, MC, MF, c("chromosome", "start", "strand"))
combs <- combn(6, 3, simplify = FALSE)
seen <- c(); NUL <- c()
for (k in combs) {
  ga <- six[k]; gb <- six[-k]
  tag <- paste(sort(ga), collapse = "|")
  if (tag %in% seen) next
  seen <- c(seen, tag, paste(sort(gb), collapse = "|"))
  is_true <- setequal(ga, PRO) || setequal(ga, DEE)
  r <- test(b6, ga, gb)
  NUL <- rbind(NUL, data.frame(split = paste(sub(" rep", "", ga), collapse = ","),
                               true_split = is_true, FDR05 = unname(r["FDR05"]),
                               FDR05_dPSI10 = unname(r["FDR05_dPSI10"])))
}
print(NUL, row.names = FALSE); w(NUL, "H1_GSE109700_label_permutation_null.tsv")
t_obs <- NUL$FDR05[NUL$true_split][1]
nul <- NUL$FDR05[!NUL$true_split]
cat(sprintf("\nTRUE split: %d events at FDR<0.05\nOther 3v3 splits: median %d, max %d (n = %d splits)\n",
            t_obs, median(nul), max(nul), length(nul)))
cat(sprintf("empirical p = %.3f\n", (1 + sum(nul >= t_obs)) / (length(nul) + 1)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_GSE109700.txt"))
