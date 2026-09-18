## What is conserved across laboratories when the global signature is not?
## Per class, per gene: does the direction of change agree across independent studies?
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; O <- file.path(D, "conserved_core")
dir.create(O, showWarnings = FALSE); say <- function(...) cat(sprintf(...), "\n"); set.seed(42)
Z <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); M <- Z$meta
N  <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))

## ---- rebuild all vectors -----------------------------------------------------
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D,"compendium/signatures"), full.names=TRUE)) {
  d <- read.delim(f); vec[[sub("\\.tsv$","",basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
lg <- list(GSE240226_UVA="figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv",
 GSE240226_Maifuyin="figure1_public_uva/GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv",
 GSE240226_succinate="figure1_public_uva/GSE240226_succinate_rescue_vs_UVA_DE.tsv",
 GSE302943_UVA="figure1_public_uva/GSE302943_UVA_vs_control_DE.tsv",
 GSE125429_UVA="figure1_public_uva/GSE125429_UVA_vs_control_DE.tsv",
 GSE89005_single6h="figure1_public_uva/GSE89005_single_6h_UVA_vs_sham_DE.tsv",
 GSE89005_single24h="figure1_public_uva/GSE89005_single_24h_UVA_vs_sham_DE.tsv",
 GSE89005_repeat24h="figure1_public_uva/GSE89005_repeated_24h_UVA_vs_sham_DE.tsv",
 GSE109700_deep="figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
 GSE109700_early="figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv",
 GSE191055_P27="figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
 GSE93535_SIPS="figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
 GSE93535_rescueSIPS="figure2_public_aging/GSE93535_SIPS1201_vs_SIPS_DE.tsv",
 GSE93535_rescueQ="figure2_public_aging/GSE93535_Q1201_vs_Q_DE.tsv",
 GSE179848_late="figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
 GSE113957_age="figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
 GSE226189_age="figure2_public_aging/GSE226189_age_per_decade_DE.tsv",
 GSE165177_MPTR="figure2_mptr/GSE165177_MPTR_paired_DE.tsv")
for (id in names(lg)) { f <- file.path(D, lg[[id]]); if (!file.exists(f)) next
  d <- read.delim(f); vec[[id]] <- setNames(d$logFC, toupper(d$gene)) }
loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))

## ---- one vector per STUDY per class (average within study) -------------------
study_vec <- function(k) {
  idx <- which(M$class == k); out <- list()
  for (ds in unique(M$dataset[idx])) {
    ii <- idx[M$dataset[idx] == ds]
    vs <- vec[M$id[ii]]; vs <- vs[!sapply(vs, is.null)]
    gg <- Reduce(union, lapply(vs, names))
    mat <- sapply(vs, function(v) v[gg]); if (is.null(dim(mat))) mat <- matrix(mat, ncol=1)
    rownames(mat) <- gg
    v <- rowMeans(mat, na.rm = TRUE)
    ## robust scaling: studies differ up to 30-fold in fold-change scale, so each
    ## study is put on a common footing before the cross-study model is fitted
    sc <- mad(v, na.rm = TRUE); if (!is.finite(sc) || sc == 0) sc <- sd(v, na.rm = TRUE)
    out[[ds]] <- v / sc }
  out }

## per gene, treat independent studies as replicates and fit an intercept-only
## model with empirical-Bayes moderation (limma); this uses effect size, not just sign
suppressPackageStartupMessages(library(limma))
core_of <- function(k, min_studies = 4) {
  sv <- study_vec(k); ns <- length(sv)
  if (ns < min_studies) return(NULL)
  gg  <- Reduce(union, lapply(sv, names))
  mat <- sapply(sv, function(v) v[gg]); rownames(mat) <- gg
  n_meas <- rowSums(is.finite(mat))
  mat <- mat[n_meas >= max(min_studies, ceiling(0.75*ns)), , drop = FALSE]
  fit <- eBayes(lmFit(mat, design = matrix(1, ncol(mat), 1)))
  tt  <- topTable(fit, number = Inf, sort.by = "none")
  n_pos <- rowSums(mat > 0, na.rm = TRUE); n_neg <- rowSums(mat < 0, na.rm = TRUE)
  nm    <- rowSums(is.finite(mat))
  data.frame(class = k, gene = rownames(mat), n_studies = nm,
             k_agree = pmax(n_pos, n_neg), dir = ifelse(n_pos >= n_neg, 1, -1),
             mean_lfc = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
             row.names = NULL) }

## how many genes are unanimous, against the binomial expectation
unanimity <- function(k, min_studies = 4) {
  sv <- study_vec(k); ns <- length(sv); if (ns < min_studies) return(NULL)
  gg <- Reduce(union, lapply(sv, names))
  mat <- sapply(sv, function(v) v[gg]); rownames(mat) <- gg
  nm <- rowSums(is.finite(mat)); mat <- mat[nm >= max(min_studies, ceiling(0.75*ns)), , drop=FALSE]
  nm <- nm[nm >= max(min_studies, ceiling(0.75*ns))]
  u  <- (rowSums(mat > 0, na.rm=TRUE) == nm) | (rowSums(mat < 0, na.rm=TRUE) == nm)
  exp_u <- sum(2 * 0.5^nm)
  data.frame(class = k, n_studies = ns, genes = length(nm), unanimous = sum(u),
             expected = exp_u, fold = sum(u)/exp_u,
             p = poisson.test(sum(u), r = exp_u)$p.value) }

CLS <- c("senescence","UV injury","photoprotection / rescue","secretome",
         "reprogramming","metabolic / culture","donor age")
res <- do.call(rbind, lapply(CLS, core_of))
write.table(res, file.path(O, "gene_sign_consistency.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

say("=== conserved core per class (genes with cross-study sign consistency, BH-FDR < 0.05) ===")
UN <- do.call(rbind, lapply(CLS, unanimity))
tab <- do.call(rbind, lapply(CLS, function(k) {
  s <- res[res$class == k, ]; if (!nrow(s)) return(NULL)
  core <- s[s$FDR < 0.05, ]; u <- UN[UN$class == k, ]
  data.frame(class = k, n_studies = u$n_studies, genes_tested = nrow(s),
             core = nrow(core), pct = 100*nrow(core)/nrow(s),
             up = sum(core$dir > 0), down = sum(core$dir < 0),
             unanimous = u$unanimous, expected = round(u$expected, 1),
             fold = round(u$fold, 2)) }))
print(tab, digits = 3, row.names = FALSE)
write.table(UN, file.path(O, "unanimity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(tab, file.path(O, "core_sizes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- what are the UV core genes? --------------------------------------------
for (k in c("UV injury","senescence","secretome")) {
  s <- res[res$class == k & res$FDR < 0.05, ]
  s <- s[order(-abs(s$mean_lfc)), ]
  say("\n--- %s core, %d genes; strongest 18 ---", k, nrow(s))
  if (nrow(s)) print(head(data.frame(gene=s$gene, k=sprintf("%d/%d", s$k_agree, s$n_studies),
        lfc=round(s$mean_lfc,3), FDR=signif(s$FDR,2)), 18), row.names = FALSE) }
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo.txt"))
