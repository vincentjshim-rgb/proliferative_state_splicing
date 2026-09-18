## Does a transcriptomic age predictor built in cultured fibroblasts read proliferation?
##
## Plan, fixed before running (post hoc with respect to the study as a whole, and
## labelled as such in the manuscript):
##   cohort    the 107 normal adult donors of GSE113957, the same matrix and the same
##             covariates (log depth, repository, instrument, sex) as everywhere else
##   features  genes passing filterByExpr, TMM-normalised log2 CPM, covariates
##             regressed out within each training fold, the 3,000 most variable genes
##   model     ridge regression solved in the dual (n << p), lambda chosen by inner
##             5-fold CV on the training fold only
##   outcome   chronological age in years; 10-fold CV repeated 5 times
##   readouts  (1) accuracy of the full predictor, r and MAE
##             (2) correlation of the cross-validated prediction with the 20-marker
##                 proliferation score
##             (3) accuracy of a predictor that uses the proliferation score alone
##             (4) accuracy when proliferation is regressed out of every gene inside
##                 each training fold
##             (5) a permuted-age control, to show the pipeline does not overfit
##   rule      whatever comes out is reported; if the predictor is unrelated to
##             proliferation that is the result
## Outputs: public_data_tierA/derived/age_predictor/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "age_predictor")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
set.seed(20260916)

## ---- data -------------------------------------------------------------------
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))
DON <- as.data.table(merge(DON, MET[, .(srr, log_depth)], by = "srr"))
P <- DON[gse113957_keep(DON, "primary")]
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
GM <- GM[, P$srr, drop = FALSE]
y <- DGEList(GM); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); LC <- cpm(y, log = TRUE)
say("matrix: %d genes x %d donors | age %d-%d", nrow(LC), ncol(LC), min(P$age), max(P$age))

PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
z <- t(scale(t(LC[apply(LC, 1, sd) > 0, ])))
P[, prolif := colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE])]
age <- P$age
COV <- model.matrix(~ log_depth + repo + instr + sex, data = P)

## ---- ridge in the dual, lambda by inner CV -----------------------------------
ridge_fit <- function(X, yy, lambdas = 10^seq(-1, 7, length.out = 33), inner = 5) {
  n <- nrow(X); folds <- sample(rep_len(1:inner, n))
  err <- sapply(lambdas, function(lam) {
    e <- 0
    for (k in 1:inner) {
      tr <- folds != k; te <- !tr
      mx <- colMeans(X[tr, , drop = FALSE]); my <- mean(yy[tr])
      Xt <- sweep(X[tr, , drop = FALSE], 2, mx); Xe <- sweep(X[te, , drop = FALSE], 2, mx)
      K <- Xt %*% t(Xt); a <- solve(K + lam * diag(nrow(K)), yy[tr] - my)
      e <- e + sum(((Xe %*% t(Xt)) %*% a + my - yy[te])^2)
    }
    e })
  lam <- lambdas[which.min(err)]
  mx <- colMeans(X); my <- mean(yy); Xt <- sweep(X, 2, mx)
  a <- solve(Xt %*% t(Xt) + lam * diag(nrow(X)), yy - my)
  list(a = a, Xt = Xt, mx = mx, my = my, lambda = lam)
}
ridge_pred <- function(f, Xnew) as.vector((sweep(Xnew, 2, f$mx) %*% t(f$Xt)) %*% f$a + f$my)

## features are built inside the training fold only
cv_predict <- function(yy, remove_prolif = FALSE, remove_cov = FALSE, ngene = 5000, folds = 10, reps = 5) {
  n <- length(yy); preds <- matrix(NA_real_, n, reps)
  for (r in 1:reps) {
    fold <- sample(rep_len(1:folds, n))
    for (k in 1:folds) {
      tr <- fold != k; te <- !tr
      ## what is removed from the features, fitted on the training donors only
      if (remove_cov || remove_prolif) {
        des    <- if (remove_prolif) cbind(COV[tr, ], prolif = P$prolif[tr]) else COV[tr, ]
        des_te <- if (remove_prolif) cbind(COV[te, , drop = FALSE], prolif = P$prolif[te]) else COV[te, , drop = FALSE]
        B <- qr.coef(qr(des), t(LC[, tr, drop = FALSE])); B[is.na(B)] <- 0
        Rtr <- t(LC[, tr, drop = FALSE]) - des %*% B
        Rte <- t(LC[, te, drop = FALSE]) - des_te %*% B
      } else {
        Rtr <- t(LC[, tr, drop = FALSE]); Rte <- t(LC[, te, drop = FALSE])
      }
      v <- apply(Rtr, 2, sd); g <- order(v, decreasing = TRUE)[1:min(ngene, sum(v > 0))]
      ## ridge is scale sensitive: standardise genes using the training donors
      mu <- colMeans(Rtr[, g, drop = FALSE]); sdv <- apply(Rtr[, g, drop = FALSE], 2, sd)
      Xtr <- sweep(sweep(Rtr[, g, drop = FALSE], 2, mu), 2, sdv, "/")
      Xte <- sweep(sweep(Rte[, g, drop = FALSE], 2, mu), 2, sdv, "/")
      f <- ridge_fit(Xtr, yy[tr])
      preds[te, r] <- ridge_pred(f, Xte)
    }
  }
  rowMeans(preds)
}

res <- list()
res$full <- cv_predict(age)
res$cov <- cv_predict(age, remove_cov = TRUE)
res$noprolif <- cv_predict(age, remove_prolif = TRUE, remove_cov = TRUE)
## proliferation score alone, same CV structure
prolif_only <- local({
  n <- length(age); pr <- matrix(NA_real_, n, 5)
  for (r in 1:5) { fold <- sample(rep_len(1:10, n))
    for (k in 1:10) { tr <- fold != k; te <- !tr
      d <- data.frame(age = age, prolif = P$prolif)
      m <- lm(age ~ prolif, d[tr, ]); pr[te, r] <- predict(m, d[te, ]) } }
  rowMeans(pr) })
res$prolif_only <- prolif_only
## permuted-age control: one permutation is not a control, because cross-validated
## predictions are not independent across donors; the null is run 20 times
PERM <- replicate(20, { p <- cv_predict(sample(age), reps = 1); cor(p, age) })
res$permuted <- cv_predict(sample(age), reps = 1)

report <- function(nm, p) {
  ct <- cor.test(p, age)
  data.frame(model = nm, r = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
             p = ct$p.value, MAE = mean(abs(p - age)),
             r_with_proliferation = cor(p, P$prolif), row.names = NULL) }
TAB <- do.call(rbind, list(
  report("all genes", res$full),
  report("all genes, cohort covariates removed", res$cov),
  report("proliferation score alone", res$prolif_only),
  report("all genes, proliferation removed", res$noprolif),
  report("permuted age (one draw)", res$permuted)))
print(TAB, digits = 3, row.names = FALSE)
write.table(TAB, file.path(O, "age_predictor_cv.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

say("\npermuted-age null over 20 draws: mean r = %+.3f, 95%% range %+.2f to %+.2f, none above %.2f",
    mean(PERM), quantile(PERM, .025), quantile(PERM, .975), max(PERM))
write.table(data.frame(permutation = seq_along(PERM), r = PERM),
            file.path(O, "age_predictor_permutation_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("correlation between the full predictor's output and the proliferation score: %+.2f", cor(res$full, P$prolif))
say("correlation between proliferation score and chronological age: %+.2f", cor(P$prolif, age))
## how much of the predictor's output is explained by proliferation
m <- summary(lm(res$full ~ P$prolif))
say("proliferation explains %.0f%% of the variance in the predicted age", 100 * m$r.squared)

fwrite(data.table(srr = P$srr, age = age, prolif = P$prolif,
                  pred_full = res$full, pred_prolif_only = res$prolif_only,
                  pred_no_prolif = res$noprolif, pred_permuted = res$permuted),
       file.path(O, "age_predictor_per_donor.tsv"), sep = "\t")

## which genes carry the predictor: fit once on all donors and look at the top weights
R <- t(LC)
v <- apply(R, 2, sd); g <- order(v, decreasing = TRUE)[1:5000]
Xall <- scale(R[, g, drop = FALSE])
f <- ridge_fit(Xall, age)
w <- as.vector(t(f$Xt) %*% f$a); names(w) <- colnames(R)[g]
top <- sort(w)[c(1:25, (length(w) - 24):length(w))]
suppressPackageStartupMessages(library(fgsea))
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
RE <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
CC <- unique(unlist(RE[grep("^Cell Cycle|^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(RE))]))
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
big <- names(sort(abs(w), decreasing = TRUE))[1:200]
say("\namong the 200 genes with the largest weights: %d cell-cycle genes (expected %.0f), %d splicing-set genes (expected %.0f)",
    sum(big %in% CC), 200 * mean(names(w) %in% CC), sum(big %in% CORE), 200 * mean(names(w) %in% CORE))
write.table(data.frame(gene = names(w), weight = w, in_cell_cycle = names(w) %in% CC,
                       in_splicing_set = names(w) %in% CORE),
            file.path(O, "age_predictor_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
