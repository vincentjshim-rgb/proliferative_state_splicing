#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# The decisive question: is the senescence splicing change separable from
# proliferative arrest?  GSE93535 is the only dataset on disk with a quiescence
# arm at MATCHED population doubling (PD15).  Use the full 2x2 design
#   state (quiescent | SIPS)  x  compound (none | 1201)
# to get n = 6 vs 6 for the state main effect instead of 3 vs 3.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(limma) })
root <- "public_data_tierA"
out  <- file.path(root, "derived", "phase2_splicing")
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
d <- readRDS(file.path(out, "H1_inputs.rds"))
S0 <- d$S; RR <- d$RR; gm <- d$gmap
gm$state <- ifelse(grepl("^SIPS", gm$arm), "senescent", "quiescent")
gm$compound <- ifelse(grepl("1201$", gm$arm), "c1201", "none")
info <- data.frame(sample = colnames(S0),
                   state = gm$state[match(colnames(S0), gm$title)],
                   compound = gm$compound[match(colnames(S0), gm$title)],
                   junction_reads_M = round(colSums(S0) / 1e6, 2))
print(info, row.names = FALSE)
cat("\nstate x compound:\n"); print(table(info$state, info$compound))

build <- function(S, smp, mincov, min_frac) {
  key <- do.call(paste, c(RR[c("chromosome", "start", "strand")], sep = "_"))
  nalt <- ave(rep(1, nrow(RR)), key, FUN = sum); sel <- which(nalt >= 2)
  Js <- S[sel, smp, drop = FALSE]; ks <- key[sel]
  tot <- apply(Js, 2, function(v) ave(v, ks, FUN = sum))
  ok <- rowMeans(tot >= mincov) >= min_frac
  Js <- Js[ok, , drop = FALSE]; tot <- tot[ok, , drop = FALSE]
  P <- Js / pmax(tot, 1)
  k <- apply(P, 1, function(r) !(min(r) > 0.98)) & apply(P, 1, sd) > 0
  P[k, , drop = FALSE]
}
fit2 <- function(P, st, cp) {
  L <- log((P + 0.01) / (1 - P + 0.01))
  dsn <- model.matrix(~ factor(st, levels = c("quiescent", "senescent")) + factor(cp))
  tt <- topTable(eBayes(lmFit(L, dsn)), coef = 2, number = Inf, sort.by = "none")
  dp <- rowMeans(P[, st == "senescent", drop = FALSE]) - rowMeans(P[, st == "quiescent", drop = FALSE])
  c(n = nrow(L), FDR05 = sum(tt$adj.P.Val < 0.05),
    FDR05_dPSI10 = sum(tt$adj.P.Val < 0.05 & abs(dp) >= 0.10),
    minP = min(tt$P.Value))
}
smp <- info$sample; st <- info$state; cp <- info$compound
cat("\n############ senescence main effect, n = 6 vs 6, filter sweep ############\n")
res <- do.call(rbind, lapply(c(5, 10, 20, 40), function(mc) {
  do.call(rbind, lapply(c(0.834, 1.0), function(mf) {
    P <- build(S0, smp, mc, mf); if (nrow(P) < 100) return(NULL)
    r <- fit2(P, st, cp)
    data.frame(mincov = mc, min_frac = mf, n_junctions = unname(r["n"]),
               FDR05 = unname(r["FDR05"]), FDR05_dPSI10 = unname(r["FDR05_dPSI10"]),
               minP = signif(unname(r["minP"]), 3)) })) }))
print(res, row.names = FALSE); w(res, "H2_GSE93535_2x2_filter_sweep.tsv")

## depth-matched copy (the two shallow runs dominate, so also try dropping them)
set.seed(93535); tgt <- min(colSums(S0)); Sd <- S0
for (j in seq_len(ncol(S0))) { p <- tgt / sum(S0[, j]); if (p < 1) Sd[, j] <- rbinom(nrow(S0), S0[, j], p) }
cat("\n############ depth-matched to", round(tgt / 1e6, 2), "M ############\n")
P <- build(Sd, smp, 5, 1.0); print(round(fit2(P, st, cp), 5))

drop <- info$junction_reads_M >= 2          # drop the two shallow runs
cat("\n############ dropping runs below 2 M junction reads (n =", sum(drop), ") ############\n")
print(table(info$state[drop], info$compound[drop]))
for (mc in c(10, 20, 40)) {
  P <- build(S0, smp[drop], mc, 1.0)
  cat(sprintf("mincov=%2d  ", mc)); print(round(fit2(P, st[drop], cp[drop]), 5))
}

## empirical null: permute the state label within compound strata
cat("\n############ label-permutation null (state permuted within compound) ############\n")
P <- build(S0, smp[drop], 20, 1.0)
st2 <- st[drop]; cp2 <- cp[drop]
obs <- fit2(P, st2, cp2)[["FDR05"]]
set.seed(11); nul <- replicate(50, {
  s <- st2; for (lv in unique(cp2)) { i <- which(cp2 == lv); s[i] <- sample(s[i]) }
  if (all(s == st2) || all(s != st2)) return(NA); fit2(P, s, cp2)[["FDR05"]] })
nul <- nul[!is.na(nul)]
cat(sprintf("observed: %d | permuted: median %d, max %d (n=%d)\nempirical p = %.3f\n",
            obs, median(nul), max(nul), length(nul), (1 + sum(nul >= obs)) / (length(nul) + 1)))
