.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; O <- file.path(D, "secretome_class")
OUT <- file.path(D, "benchmark_figs"); dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n"); set.seed(3)
Z <- readRDS(file.path(O, "structure_with_secretome.rds")); R <- Z$R; M <- Z$meta
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
N  <- read.delim(file.path(O, "secretome_signatures.tsv"))

## ---- rebuild every signature vector (for magnitude + axis loading) ----------
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
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

M$mag <- sapply(M$id, function(i) { v <- vec[[i]]; if (is.null(v)) NA else median(abs(v[is.finite(v)]), na.rm=TRUE) })
axs <- setNames(AX$senescence_meta, toupper(AX$gene))
M$rho_sen <- sapply(M$id, function(i) { v <- vec[[i]]; if (is.null(v)) return(NA)
  g <- intersect(names(v), names(axs)); if (length(g) < 1000) return(NA)
  cor(v[g], axs[g], method = "spearman") })
M$n_sen <- sapply(M$id, function(i) { v <- vec[[i]]; if (is.null(v)) return(NA)
  length(intersect(names(v), names(axs))) })

## ---- MDS on the 76 x 76 matrix ---------------------------------------------
Rf <- R; Rf[is.na(Rf)] <- 0
md <- cmdscale(as.dist(1 - Rf), k = 2, eig = TRUE)
ev <- md$eig[md$eig > 0]; ev <- 100*ev/sum(ev)
M$MDS1 <- md$points[,1]; M$MDS2 <- md$points[,2]
say("MDS eigen: %s", paste(sprintf("%.1f%%", head(ev,4)), collapse = "  "))
ok <- is.finite(M$rho_sen)
say("MDS1 vs rho_senescence: r = %+.3f", cor(M$MDS1[ok], M$rho_sen[ok]))
writeLines(sprintf("%.2f", head(ev, 5)), file.path(OUT, "mds_eigen.txt"))

## ---- within-study vs between-study coherence per class ---------------------
WB <- do.call(rbind, lapply(sort(unique(M$class)), function(k) {
  idx <- which(M$class == k); if (length(idx) < 2) return(NULL)
  s <- R[idx, idx]; ds <- M$dataset[idx]; sm <- outer(ds, ds, "==")
  w <- s[upper.tri(s) & sm];  w <- w[is.finite(w)]
  b <- s[upper.tri(s) & !sm]; b <- b[is.finite(b)]
  bc <- if (length(b) > 1) quantile(replicate(3000, median(sample(b, length(b), TRUE))), c(.025,.975)) else c(NA,NA)
  wc <- if (length(w) > 1) quantile(replicate(3000, median(sample(w, length(w), TRUE))), c(.025,.975)) else c(NA,NA)
  data.frame(class = k, n_sig = length(idx),
             n_within = length(w), within = ifelse(length(w) > 0, median(w), NA),
             w_lo = wc[1], w_hi = wc[2],
             n_between = length(b), between = ifelse(length(b) > 0, median(b), NA),
             b_lo = bc[1], b_hi = bc[2],
             mag = median(M$mag[idx], na.rm = TRUE)) }))
WB <- WB[order(-WB$between), ]
print(WB, digits = 3, row.names = FALSE)
write.table(WB, file.path(OUT, "within_between_by_class.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- per-dataset within-study coherence for the secretome class ------------
DS <- do.call(rbind, lapply(unique(M$dataset[M$class == "secretome"]), function(ds) {
  idx <- which(M$dataset == ds & M$class == "secretome"); if (length(idx) < 2) return(NULL)
  s <- R[idx, idx]; v <- s[upper.tri(s)]; v <- v[is.finite(v)]
  data.frame(dataset = ds, n = length(idx), n_pairs = length(v),
             median_r = median(v), min = min(v), max = max(v)) }))
DS <- DS[order(-DS$median_r), ]; print(DS, digits = 3, row.names = FALSE)
write.table(DS, file.path(OUT, "secretome_within_study.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- UV-injury pairwise detail (the single coherent pair) ------------------
uv <- which(M$class == "UV injury"); s <- R[uv, uv]
rownames(s) <- colnames(s) <- M$id[uv]
ds <- M$dataset[uv]; sm <- outer(ds, ds, "==")
pp <- data.frame(a = rownames(s)[row(s)[upper.tri(s) & !sm]],
                 b = colnames(s)[col(s)[upper.tri(s) & !sm]],
                 r = s[upper.tri(s) & !sm])
pp <- pp[is.finite(pp$r), ]; pp <- pp[order(-pp$r), ]
say("\nUV injury between-study pairs: %d; top 3:", nrow(pp)); print(head(pp, 3), digits = 3, row.names = FALSE)
write.table(pp, file.path(OUT, "uv_pairwise.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
sen <- which(M$class == "senescence"); s2 <- R[sen, sen]
ds2 <- M$dataset[sen]; sm2 <- outer(ds2, ds2, "==")
ps <- data.frame(r = s2[upper.tri(s2) & !sm2]); ps <- ps[is.finite(ps$r), , drop = FALSE]
write.table(ps, file.path(OUT, "sen_pairwise.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

saveRDS(list(R = R, meta = M), file.path(OUT, "figdata.rds"))
write.table(M, file.path(OUT, "signature_metadata.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("\nprepared -> %s", OUT)
