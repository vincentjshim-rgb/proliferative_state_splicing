#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Phase 2, H1 : is splicing OUTCOME a senescence axis that is separable from
#               proliferative arrest?
#   Quiescent (Q) vs stress-induced senescent (SIPS), BOTH at population
#   doubling 15, n = 3 vs 3.  GSE93535 / SRP096629, recount3 junction counts.
# Prespecified in repro_cm_preregistration_phase2_ko.md (SHA256 41e4cf4b...).
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"; rc <- file.path(root, "recount3")
out  <- file.path(root, "derived", "phase2_splicing")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
SRP <- "SRP096629"

## ---- sample map: GSM -> title -> SRX -> SRR -> rail_id --------------------
soft <- readLines(gzfile(file.path(root, "metadata", "GSE93535_family.soft.gz")))
gsm <- sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", soft, value = TRUE))
blk <- split(soft, cumsum(grepl("^\\^SAMPLE", soft)))[-1]
gmap <- do.call(rbind, lapply(blk, function(b) data.frame(
  gsm = sub("^\\^SAMPLE = ", "", b[1]),
  title = sub(".*= ", "", grep("!Sample_title", b, value = TRUE)[1]),
  srx = sub(".*term=", "", grep("!Sample_relation = SRA", b, value = TRUE)[1]))))
ena <- read.delim(file.path(root, "metadata", "tierB_runinfo", "PRJNA361062.tsv"))
rid <- read.delim(gzfile(file.path(rc, sprintf("sra.recount_project.%s.MD.gz", SRP))))
srr2rail <- setNames(rid$rail_id, rid$external_id)
gmap$runs <- sapply(gmap$srx, function(x) paste(ena$run_accession[ena$experiment_accession == x], collapse = ","))
gmap$arm <- sub("_[0-9]+$", "", gmap$title)
gmap$group <- ifelse(gmap$arm %in% c("Q", "SIPS"), gmap$arm, paste0(gmap$arm, "_compound"))
cat("samples:", nrow(gmap), "\n"); print(table(gmap$arm))

## ---- load recount3 junction matrix ----------------------------------------
ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))),
            what = character(), quiet = TRUE)
ids <- ids[ids != "rail_id"]
RRf <- read.delim(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))))
mm  <- gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)), "rt")
hdr <- readLines(mm, n = 3); dims <- as.integer(strsplit(trimws(hdr[3]), "\\s+")[[1]])
trip <- scan(mm, what = integer(), quiet = TRUE); close(mm)
trip <- matrix(trip, ncol = 3, byrow = TRUE)
cat("junction matrix:", dims[1], "junctions x", dims[2], "samples;", nrow(trip), "nonzero\n")
stopifnot(dims[1] == nrow(RRf), dims[2] == length(ids))
J <- matrix(0L, nrow = dims[1], ncol = dims[2], dimnames = list(NULL, ids))
J[cbind(trip[, 1], trip[, 2])] <- trip[, 3]

## collapse technical runs to biological samples
colrail <- as.character(ids)
S <- sapply(seq_len(nrow(gmap)), function(i) {
  r <- as.character(srr2rail[strsplit(gmap$runs[i], ",")[[1]]])
  k <- which(colrail %in% r); stopifnot(length(k) >= 1)
  if (length(k) == 1) J[, k] else rowSums(J[, k, drop = FALSE])
})
colnames(S) <- gmap$title
cat("collapsed to", ncol(S), "biological samples; total junction reads per sample:\n")
print(round(colSums(S) / 1e6, 2))

## ---- PSI within splice-site groups ----------------------------------------
psi_table <- function(anchor_cols, mincov, grp_a, grp_b) {
  key <- do.call(paste, c(RRf[anchor_cols], sep = "_"))
  nalt <- ave(rep(1, nrow(RRf)), key, FUN = sum)
  sel <- which(nalt >= 2)
  Ssel <- S[sel, , drop = FALSE]; ksel <- key[sel]
  tot <- apply(Ssel, 2, function(v) ave(v, ksel, FUN = sum))
  smp <- c(grp_a, grp_b)
  ok <- rowSums(tot[, smp, drop = FALSE] >= mincov) == length(smp)
  P <- Ssel[ok, smp, drop = FALSE] / tot[ok, smp, drop = FALSE]
  keep <- !(apply(P, 1, min) > 0.98)
  list(psi = P[keep, , drop = FALSE], idx = sel[ok][keep],
       cov = tot[ok, smp, drop = FALSE][keep, , drop = FALSE])
}
diff_psi <- function(pt, grp_a, grp_b, label) {
  P <- pt$psi
  a <- P[, grp_a, drop = FALSE]; b <- P[, grp_b, drop = FALSE]
  d <- rowMeans(b) - rowMeans(a)
  lg <- function(x) log((x + 0.01) / (1 - x + 0.01))
  la <- lg(a); lb <- lg(b)
  pv <- suppressWarnings(vapply(seq_len(nrow(P)), function(i) {
    x <- la[i, ]; y <- lb[i, ]
    if (stats::sd(c(x, y)) == 0) return(NA_real_)
    tryCatch(stats::t.test(y, x)$p.value, error = function(e) NA_real_)
  }, numeric(1)))
  fdr <- p.adjust(pv, "BH")
  list(summary = data.frame(
         label = label, n_testable = nrow(P),
         median_abs_dPSI = median(abs(d), na.rm = TRUE),
         q90_abs_dPSI = unname(quantile(abs(d), 0.9, na.rm = TRUE)),
         n_dPSI_ge_0.10 = sum(abs(d) >= 0.10, na.rm = TRUE),
         n_sig_FDR05 = sum(fdr < 0.05, na.rm = TRUE),
         n_sig_and_dPSI_ge_0.10 = sum(fdr < 0.05 & abs(d) >= 0.10, na.rm = TRUE)),
       detail = data.frame(junction_row = pt$idx, dPSI = d, p = pv, FDR = fdr))
}

run_contrast <- function(grp_a, grp_b, label, mincov = 20) {
  res <- lapply(list(c("chromosome", "start", "strand"), c("chromosome", "end", "strand")),
                function(cols) {
                  pt <- psi_table(cols, mincov, grp_a, grp_b)
                  diff_psi(pt, grp_a, grp_b, paste(label, cols[2], "anchored"))
                })
  list(summary = do.call(rbind, lapply(res, `[[`, "summary")),
       detail = do.call(rbind, lapply(res, `[[`, "detail")))
}

## depth is strongly unequal between arms; keep the raw matrix and a
## depth-matched copy obtained by binomial downsampling to the smallest library.
S_raw <- S
set.seed(93535)
target <- min(colSums(S_raw))
S_ds <- S_raw
for (j in seq_len(ncol(S_raw))) {
  p <- target / sum(S_raw[, j])
  if (p < 1) S_ds[, j] <- rbinom(nrow(S_raw), S_raw[, j], p)
}
cat("\ndepth-matched target per sample:", round(target / 1e6, 2), "M junction reads\n")

Q    <- gmap$title[gmap$arm == "Q"]
SIPS <- gmap$title[gmap$arm == "SIPS"]
Q12  <- gmap$title[gmap$arm == "Q1201"]
S12  <- gmap$title[gmap$arm == "SIPS1201"]
cat("\nPRIMARY H1 contrast: Q (", paste(Q, collapse = ","), ") vs SIPS (",
    paste(SIPS, collapse = ","), ") -- both at PD15\n")

H1 <- run_contrast(Q, SIPS, "SIPS vs Q (primary, raw depth)")
S <- S_ds
H1d <- run_contrast(Q, SIPS, "SIPS vs Q (depth-matched)")
S <- S_raw
w(rbind(H1$summary, H1d$summary), "H1_SIPS_vs_Q_summary.tsv")
cat("\n================ H1 PRIMARY: SIPS vs quiescent, matched PD ================\n")
print(rbind(H1$summary, H1d$summary), row.names = FALSE, digits = 3)

## negative control: the same test between two arms that differ only by compound
NC <- run_contrast(Q, Q12, "Q1201 vs Q (compound-only negative control)")
S <- S_ds; NCd <- run_contrast(Q, Q12, "Q1201 vs Q (negative control, depth-matched)"); S <- S_raw
cat("\n================ negative control: compound within quiescent ================\n")
print(rbind(NC$summary, NCd$summary), row.names = FALSE, digits = 3)
w(rbind(NC$summary, NCd$summary), "H1_negative_control_summary.tsv")

## exploratory (not prespecified): unannotated junction fraction per sample
ann <- RRf$annotated == 1
frac <- data.frame(sample = colnames(S), arm = gmap$arm[match(colnames(S), gmap$title)],
                   total_junction_reads = colSums(S),
                   unannotated_read_fraction = colSums(S[!ann, ]) / colSums(S))
w(frac, "H1_unannotated_junction_fraction_exploratory.tsv")
cat("\n=== exploratory: unannotated junction read fraction (not prespecified) ===\n")
print(frac[order(frac$arm), c("sample", "arm", "unannotated_read_fraction")],
      row.names = FALSE, digits = 4)
a <- frac$unannotated_read_fraction[frac$arm == "Q"]
b <- frac$unannotated_read_fraction[frac$arm == "SIPS"]
cat(sprintf("Q mean = %.5f ; SIPS mean = %.5f ; Welch p = %.3g\n",
            mean(a), mean(b), t.test(b, a)$p.value))

saveRDS(list(gmap = gmap, S = S, RR = RRf), file.path(out, "H1_inputs.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
