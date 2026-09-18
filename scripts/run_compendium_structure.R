## What actually structures the 60-perturbation space?
## Builds the full 60 x 60 signature-similarity matrix and asks whether the
## structure is captured by (a) perturbation class, (b) module content, or
## (c) the two hand-built reference axes.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D   <- "public_data_tierA/derived"
OUT <- file.path(D, "compendium_structure"); dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
set.seed(20260912)

## ---------------------------------------------------------- 1. all 60 sigs
AX  <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
sigs <- list()
for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
  d <- read.delim(f)
  sigs[[sub("\\.tsv$", "", basename(f))]] <-
    list(dataset = d$dataset[1], contrast = d$contrast[1],
         eff = data.frame(gene = toupper(d$gene), logFC = d$logFC)) }
legacy <- list(
  GSE240226_UVA=c("figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv","GSE240226","UVA vs control"),
  GSE240226_Maifuyin=c("figure1_public_uva/GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv","GSE240226","Maifuyin rescue"),
  GSE240226_succinate=c("figure1_public_uva/GSE240226_succinate_rescue_vs_UVA_DE.tsv","GSE240226","succinate rescue"),
  GSE302943_UVA=c("figure1_public_uva/GSE302943_UVA_vs_control_DE.tsv","GSE302943","cumulative UVA"),
  GSE125429_UVA=c("figure1_public_uva/GSE125429_UVA_vs_control_DE.tsv","GSE125429","chronic UVA"),
  GSE89005_single6h=c("figure1_public_uva/GSE89005_single_6h_UVA_vs_sham_DE.tsv","GSE89005","single UVA 6 h"),
  GSE89005_single24h=c("figure1_public_uva/GSE89005_single_24h_UVA_vs_sham_DE.tsv","GSE89005","single UVA 24 h"),
  GSE89005_repeat24h=c("figure1_public_uva/GSE89005_repeated_24h_UVA_vs_sham_DE.tsv","GSE89005","repeated UVA 24 h"),
  GSE109700_deep=c("figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv","GSE109700","deep replicative senescence"),
  GSE109700_early=c("figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv","GSE109700","early replicative senescence"),
  GSE191055_P27=c("figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv","GSE191055","P27 vs P4"),
  GSE93535_SIPS=c("figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv","GSE93535","SIPS vs quiescent"),
  GSE93535_rescueSIPS=c("figure2_public_aging/GSE93535_SIPS1201_vs_SIPS_DE.tsv","GSE93535","compound 1201 in SIPS"),
  GSE93535_rescueQ=c("figure2_public_aging/GSE93535_Q1201_vs_Q_DE.tsv","GSE93535","compound 1201 in quiescent"),
  GSE179848_late=c("figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv","GSE179848","late vs early passage"),
  GSE113957_age=c("figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv","GSE113957","donor age per decade"),
  GSE226189_age=c("figure2_public_aging/GSE226189_age_per_decade_DE.tsv","GSE226189","donor age per decade"),
  GSE165177_MPTR=c("figure2_mptr/GSE165177_MPTR_paired_DE.tsv","GSE165177","transient reprogramming (MPTR)"))
for (id in names(legacy)) {
  f <- file.path(D, legacy[[id]][1]); if (!file.exists(f)) next
  d <- read.delim(f)
  sigs[[id]] <- list(dataset = legacy[[id]][2], contrast = legacy[[id]][3],
                     eff = data.frame(gene = toupper(d$gene), logFC = d$logFC)) }
loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
sigs[["LOCAL_ReproCM"]] <- list(dataset = "this study", contrast = "Repro-CM anchor",
  eff = data.frame(gene = toupper(loc$gene), logFC = loc$Repro_specific_score))
say("signatures assembled: %d", length(sigs))

## ------------------------------------------------- 2. 60 x 60 similarity
ids <- names(sigs); n <- length(ids)
R <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
NSH <- matrix(NA_integer_, n, n, dimnames = list(ids, ids))
vec <- lapply(sigs, function(s) { e <- s$eff; e <- e[is.finite(e$logFC), ]
  e <- e[!duplicated(e$gene), ]; setNames(e$logFC, e$gene) })
for (i in seq_len(n)) for (j in seq_len(n)) {
  if (j < i) next
  g <- intersect(names(vec[[i]]), names(vec[[j]]))
  NSH[i,j] <- NSH[j,i] <- length(g)
  if (length(g) >= 1000) {
    r <- suppressWarnings(cor(vec[[i]][g], vec[[j]][g], method = "spearman"))
    R[i,j] <- R[j,i] <- r } }
diag(R) <- 1
say("pairs with >=1,000 shared genes: %d of %d (%.0f%%)",
    sum(!is.na(R[upper.tri(R)])), sum(upper.tri(R)),
    100*mean(!is.na(R[upper.tri(R)])))
say("median shared genes: %s", format(median(NSH[upper.tri(NSH)]), big.mark = ","))
write.table(round(R, 4), file.path(OUT, "signature_similarity_60x60.tsv"), sep = "\t", quote = FALSE)

## ------------------------------------------------------- 3. class labels
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
cls <- function(id, con) {
  if (id == "LOCAL_ReproCM") "Repro-CM"
  else if (grepl("senesc|SIPS|_deep|_early|P27|_KO|siScramble|pLenti", paste(id,con), ignore.case=TRUE)) "senescence"
  else if (grepl("age|old_vs_young", paste(id,con), ignore.case=TRUE)) "donor age"
  else if (grepl("rescue|red_vs|prered|Osmoter|sauchinone|UIFSP|Maifuyin|succinate|1201",
                 paste(id,con), ignore.case=TRUE)) "photoprotection / rescue"
  else if (grepl("UV", paste(id,con))) "UV injury"
  else if (grepl("D13|D7|dox|OSK|MPTR|reprogram", paste(id,con), ignore.case=TRUE)) "reprogramming"
  else "metabolic / culture" }
meta <- data.frame(id = ids,
  dataset  = sapply(sigs, `[[`, "dataset"),
  contrast = sapply(sigs, `[[`, "contrast"), stringsAsFactors = FALSE)
meta$class <- mapply(cls, meta$id, meta$contrast)
meta$D     <- CP$D[match(meta$id, CP$id)]
meta$circ  <- CP$circularity[match(meta$id, CP$id)]
print(table(meta$class))

## ----------------------------------------- 4. is the structure class-driven?
Rf <- R; Rf[is.na(Rf)] <- 0
dmat <- as.dist(1 - Rf)
hc   <- hclust(dmat, method = "average")
mds  <- cmdscale(dmat, k = 3, eig = TRUE)
ev   <- mds$eig[mds$eig > 0]; ev <- 100 * ev / sum(ev)
say("classical MDS variance: axis1 %.1f%%  axis2 %.1f%%  axis3 %.1f%%", ev[1], ev[2], ev[3])
meta <- cbind(meta, MDS1 = mds$points[,1], MDS2 = mds$points[,2], MDS3 = mds$points[,3])

## within- vs between-class similarity, permutation test on the class labels
ut <- upper.tri(R)
ci <- outer(meta$class, meta$class, "==")
within  <- R[ut & ci]; between <- R[ut & !ci]
within  <- within[is.finite(within)]; between <- between[is.finite(between)]
obs <- median(within) - median(between)
perm <- replicate(9999, {
  s  <- sample(meta$class); cs <- outer(s, s, "==")
  w  <- R[ut & cs]; b <- R[ut & !cs]
  median(w[is.finite(w)]) - median(b[is.finite(b)]) })
pv <- (1 + sum(perm >= obs)) / (1 + length(perm))
say("within-class median r = %+.3f (n=%d), between-class %+.3f (n=%d)",
    median(within), length(within), median(between), length(between))
say("difference = %+.3f, label-permutation P = %.4f", obs, pv)
writeLines(sprintf("within=%.4f between=%.4f diff=%.4f perm_p=%.4f",
                   median(within), median(between), obs, pv),
           file.path(OUT, "class_coupling_permutation.txt"))

## per-class coherence
coh <- do.call(rbind, lapply(sort(unique(meta$class)), function(k) {
  idx <- which(meta$class == k); if (length(idx) < 2) return(NULL)
  sub <- R[idx, idx]; v <- sub[upper.tri(sub)]; v <- v[is.finite(v)]
  out <- R[idx, -idx]; out <- out[is.finite(out)]
  data.frame(class = k, n = length(idx), within_median = median(v),
             out_median = median(out), delta = median(v) - median(out)) }))
coh <- coh[order(-coh$delta), ]; print(coh, digits = 3)
write.table(coh, file.path(OUT, "class_coherence.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## -------------------------- 5. does module content explain pairwise similarity?
MO <- read.delim(file.path(D, "what_drives_alignment/module_scores_vs_alignment.tsv"))
mods <- c("HSF1","premRNA","splicing","UPR_PERK","prolif","ECM")
MOm  <- MO[match(meta$id, MO$id), mods]
ok   <- which(rowSums(is.na(MOm)) == 0)
say("signatures with complete module scores: %d", length(ok))
Z    <- scale(as.matrix(MOm[ok, ]))
Dm   <- as.matrix(dist(Z))                      # module-content distance
Rs   <- R[ok, ok]
utk  <- upper.tri(Rs) & is.finite(Rs)
ct   <- suppressWarnings(cor.test(Dm[utk], Rs[utk], method = "spearman", exact = FALSE))
say("module-content distance vs signature similarity: rho = %+.3f, P = %.2e",
    ct$estimate, ct$p.value)
## Mantel permutation over signature labels
mp <- replicate(9999, { p <- sample(seq_len(length(ok)))
  suppressWarnings(cor(Dm[utk], Rs[p, p][utk], method = "spearman")) })
mantel_p <- (1 + sum(mp <= ct$estimate)) / (1 + length(mp))
say("Mantel permutation P = %.4f", mantel_p)

## the same for the two hand-built axes only
AXd <- as.matrix(dist(scale(cbind(CP$rho_UVA, CP$rho_senescence)[match(meta$id, CP$id), ])))
AXd <- AXd[ok, ok]
ct2 <- suppressWarnings(cor.test(AXd[utk], Rs[utk], method = "spearman", exact = FALSE))
say("two-axis distance vs signature similarity:      rho = %+.3f, P = %.2e",
    ct2$estimate, ct2$p.value)
writeLines(c(sprintf("module_distance_rho=%.4f p=%.3g mantel_p=%.4f", ct$estimate, ct$p.value, mantel_p),
             sprintf("two_axis_distance_rho=%.4f p=%.3g", ct2$estimate, ct2$p.value)),
           file.path(OUT, "mantel_module_vs_axis.txt"))

## ------------------------------- 6. what is Repro-CM's nearest neighbourhood?
rc <- R["LOCAL_ReproCM", ]; rc <- rc[names(rc) != "LOCAL_ReproCM"]
nb <- data.frame(id = names(rc), r = as.numeric(rc))
nb$class <- meta$class[match(nb$id, meta$id)]
nb$contrast <- meta$contrast[match(nb$id, meta$id)]
nb <- nb[order(-nb$r), ]
say("Repro-CM nearest neighbours:")
print(head(nb[, c("id","class","r")], 10), digits = 3)
say("Repro-CM furthest:")
print(tail(nb[, c("id","class","r")], 5), digits = 3)
write.table(nb, file.path(OUT, "reproCM_neighbours.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
## rank-sum: is Repro-CM closer to photoprotection than to UV injury?
w <- wilcox.test(nb$r[nb$class == "photoprotection / rescue"], nb$r[nb$class == "UV injury"])
say("Repro-CM vs photoprotection (median r = %+.3f) vs UV injury (%+.3f): Wilcoxon P = %.4f",
    median(nb$r[nb$class=="photoprotection / rescue"]), median(nb$r[nb$class=="UV injury"]), w$p.value)
w2 <- wilcox.test(nb$r[nb$class == "photoprotection / rescue"], nb$r[nb$class == "senescence"])
say("            vs senescence (%+.3f): Wilcoxon P = %.4f",
    median(nb$r[nb$class=="senescence"]), w2$p.value)

write.table(meta, file.path(OUT, "compendium_metadata_mds.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(R = R, NSH = NSH, meta = meta, hc = hc), file.path(OUT, "structure.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
say("done -> %s", OUT)
