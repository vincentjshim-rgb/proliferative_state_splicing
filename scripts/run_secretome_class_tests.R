## Preregistered tests HS1-HS3 (SHA-256 55bf016d...), executed after the data were retrieved.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; OUT <- file.path(D, "secretome_class")
say <- function(...) cat(sprintf(...), "\n"); set.seed(20260912)

S  <- readRDS(file.path(D, "compendium_structure/structure.rds"))
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
NEW <- read.delim(file.path(OUT, "secretome_signatures.tsv"))
vecs <- lapply(split(NEW, NEW$id), function(d) setNames(d$logFC, d$gene))
meta_new <- unique(NEW[, c("id","dataset","contrast","recipient")])

## ---- existing 60 signature vectors -------------------------------------------
old <- readRDS(file.path(D, "compendium_structure/structure.rds"))
## rebuild the vectors exactly as in run_compendium_structure.R
source_vecs <- local({
  sigs <- list()
  for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
    d <- read.delim(f); sigs[[sub("\\.tsv$","",basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
  legacy <- list(
    GSE240226_UVA="figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv",
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
  for (id in names(legacy)) { f <- file.path(D, legacy[[id]]); if (!file.exists(f)) next
    d <- read.delim(f); sigs[[id]] <- setNames(d$logFC, toupper(d$gene)) }
  loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
  sigs[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))
  sigs })
allv <- c(source_vecs, vecs)
allv <- lapply(allv, function(v) { v <- v[is.finite(v)]; v[!duplicated(names(v))] })
say("total signatures: %d  (was %d, added %d)", length(allv), length(source_vecs), length(vecs))

## ---- similarity matrix -------------------------------------------------------
n <- length(allv); ids <- names(allv)
R <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
for (i in 1:n) for (j in i:n) {
  g <- intersect(names(allv[[i]]), names(allv[[j]]))
  if (length(g) >= 1000) R[i,j] <- R[j,i] <- cor(allv[[i]][g], allv[[j]][g], method = "spearman") }
diag(R) <- 1
saveRDS(list(R = R), file.path(OUT, "similarity_with_secretome.rds"))
write.table(round(R,4), file.path(OUT, "similarity_with_secretome.tsv"), sep="\t", quote=FALSE)

## ---- class labels ------------------------------------------------------------
ml <- S$meta[, c("id","dataset","class")]
ml2 <- data.frame(id = meta_new$id, dataset = meta_new$dataset,
  class = ifelse(meta_new$id == "UVX", "UV injury", "secretome"))
M <- rbind(ml, ml2); M <- M[match(ids, M$id), ]
M$dataset[is.na(M$dataset)] <- "?"
print(table(M$class))

## between-study coherence per class (same-study pairs excluded)
between <- function(k, drop = character()) {
  idx <- which(M$class == k & !M$id %in% drop)
  sub <- R[idx, idx]; ds <- M$dataset[idx]
  same <- outer(ds, ds, "==")
  v <- sub[upper.tri(sub) & !same]; v[is.finite(v)] }
CL <- do.call(rbind, lapply(sort(unique(M$class)), function(k) {
  v <- between(k); if (length(v) < 2) return(NULL)
  ci <- quantile(replicate(3000, median(sample(v, replace = TRUE))), c(.025,.975))
  data.frame(class = k, n_sig = sum(M$class == k), n_pairs = length(v),
             median_r = median(v), lo = ci[1], hi = ci[2]) }))
CL <- CL[order(-CL$median_r), ]
say("\n=== between-study coherence, all classes ===")
print(CL, digits = 3, row.names = FALSE)
write.table(CL, file.path(OUT, "between_study_coherence.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ================== HS1  secretome between-study coherence ==================
sec_main <- between("secretome", drop = "SC11")     # gingival excluded per prereg
say("\n### HS1  secretome, between-study pairs (gingival excluded)")
say("    n pairs = %d,  median r = %+.4f  [%.3f, %.3f]",
    length(sec_main), median(sec_main),
    quantile(replicate(3000, median(sample(sec_main, replace=TRUE))), .025),
    quantile(replicate(3000, median(sample(sec_main, replace=TRUE))), .975))
verdict <- if (median(sec_main) > 0.35) "REFUTED (coheres like senescence)" else
           if (abs(median(sec_main)) < 0.15) "PREDICTION HELD (no coherent class)" else
           "PARTIAL (0.15-0.35, not interpreted either way)"
say("    prediction |median r| < 0.15, refuting > 0.35  ->  %s", verdict)
sec_all <- between("secretome")
say("    with the gingival boundary contrast included: median r = %+.4f (n = %d)",
    median(sec_all), length(sec_all))

## ================== HS2  positive control ===================================
axg <- setNames(AX$senescence_meta, toupper(AX$gene))
g <- intersect(names(allv[["SC1"]]), names(axg))
r_hs2 <- cor(allv[["SC1"]][g], axg[g], method = "spearman")
say("\n### HS2  paracrine-senescence CM (SC1) on the senescence reference axis")
say("    r = %+.4f over %s shared genes   (prediction r > +0.20)  ->  %s",
    r_hs2, format(length(g), big.mark=","), ifelse(r_hs2 > 0.20, "HELD", "FAILED"))
g2 <- intersect(names(allv[["SC2"]]), names(axg))
say("    SC2 (day 10) r = %+.4f", cor(allv[["SC2"]][g2], axg[g2], method = "spearman"))

## ================== HS3  Repro-CM position ==================================
REG <- c("SC3","SC3b","SC3c","SC3d","SC4","SC5","SC6","SC7","SC8","SC9","SC10")
SEN <- c("SC1","SC2")
rc  <- R["LOCAL_ReproCM", ]
say("\n### HS3  Repro-CM similarity")
say("    regenerative secretomes (n=%d): median r = %+.4f", length(REG), median(rc[REG], na.rm=TRUE))
say("    senescence-transmitting CM (n=%d): median r = %+.4f", length(SEN), median(rc[SEN], na.rm=TRUE))
say("    direction predicted: regenerative > senescence-transmitting  ->  %s",
    ifelse(median(rc[REG],na.rm=TRUE) > median(rc[SEN],na.rm=TRUE), "HELD", "FAILED"))
tab <- data.frame(id = c(REG, SEN, "SC11"), r = as.numeric(rc[c(REG, SEN, "SC11")]))
tab$contrast <- meta_new$contrast[match(tab$id, meta_new$id)]
tab$dataset  <- meta_new$dataset[match(tab$id, meta_new$id)]
print(tab[order(-tab$r), c("id","dataset","r","contrast")], digits = 3, row.names = FALSE)
write.table(tab, file.path(OUT, "HS3_reproCM_vs_secretome.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- the added UV contrast ---------------------------------------------------
uv <- which(M$class == "UV injury")
say("\n### UV injury class with the added GSE134533 contrast (n = %d)", length(uv))
sub <- R[uv, uv]; ds <- M$dataset[uv]; same <- outer(ds, ds, "==")
v <- sub[upper.tri(sub) & !same]; v <- v[is.finite(v)]
say("    between-study median r = %+.4f (n = %d pairs)", median(v), length(v))
rx <- R["UVX", uv]; say("    GSE134533 HDF UVB vs the other UV contrasts: median r = %+.4f",
                        median(rx[names(rx) != "UVX"], na.rm = TRUE))
saveRDS(list(R = R, meta = M), file.path(OUT, "structure_with_secretome.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_tests.txt"))
