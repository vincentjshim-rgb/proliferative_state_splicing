.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; OUT <- file.path(D, "compendium_structure")
S <- readRDS(file.path(OUT, "structure.rds")); R <- S$R; meta <- S$meta
say <- function(...) cat(sprintf(...), "\n"); set.seed(11)

## ---- 1. class coherence with bootstrap CIs -------------------------------
boot_med <- function(v, B = 2000) { v <- v[is.finite(v)]
  if (length(v) < 2) return(c(NA, NA))
  quantile(replicate(B, median(sample(v, replace = TRUE))), c(.025, .975)) }
coh <- do.call(rbind, lapply(sort(unique(meta$class)), function(k) {
  idx <- which(meta$class == k); if (length(idx) < 2) return(NULL)
  sub <- R[idx, idx]; w <- sub[upper.tri(sub)]; w <- w[is.finite(w)]
  o   <- R[idx, -idx]; o <- o[is.finite(o)]
  ci  <- boot_med(w)
  wt  <- wilcox.test(w, o)
  data.frame(class = k, n_sig = length(idx), n_pairs = length(w),
             within = median(w), lo = ci[1], hi = ci[2], between = median(o),
             delta = median(w) - median(o), p = wt$p.value) }))
coh$FDR <- p.adjust(coh$p, "BH"); coh <- coh[order(-coh$delta), ]
print(coh, digits = 3, row.names = FALSE)
write.table(coh, file.path(OUT, "class_coherence_CI.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- 2. Repro-CM against every class ------------------------------------
rc <- R["LOCAL_ReproCM", ]; rc <- rc[names(rc) != "LOCAL_ReproCM"]
cl <- meta$class[match(names(rc), meta$id)]
nb <- data.frame(id = names(rc), r = as.numeric(rc), class = cl)
nb <- nb[is.finite(nb$r), ]
say("\nRepro-CM similarity by class (n usable = %d):", nrow(nb))
byc <- do.call(rbind, lapply(sort(unique(nb$class)), function(k) {
  v <- nb$r[nb$class == k]; o <- nb$r[nb$class != k]
  data.frame(class = k, n = length(v), median_r = median(v),
             p_vs_rest = if (length(v) > 1) wilcox.test(v, o)$p.value else NA) }))
byc$FDR <- p.adjust(byc$p_vs_rest, "BH"); byc <- byc[order(-byc$median_r), ]
print(byc, digits = 3, row.names = FALSE)
write.table(byc, file.path(OUT, "reproCM_by_class.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## exclude circular contrasts (those that define an axis) and repeat
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
nb$circ <- CP$circularity[match(nb$id, CP$id)]
ind <- nb[nb$circ == "independent", ]
say("\nafter removing axis-defining and same-study contrasts (n = %d):", nrow(ind))
byi <- do.call(rbind, lapply(sort(unique(ind$class)), function(k) {
  v <- ind$r[ind$class == k]
  data.frame(class = k, n = length(v), median_r = median(v), max_r = max(v)) }))
print(byi[order(-byi$median_r), ], digits = 3, row.names = FALSE)

## ---- 3. is the neighbourhood bioenergetic rather than photoprotective? ----
BIO <- c("GSE179848_Oligomycin","GSE179848_OligomycinDEX","GSE179848_2Deoxyglucose",
         "GSE179848_SURF1","GSE179848_Galactose","GSE179848_mitoNUITs",
         "GSE179848_mitoNUITsDEX","GSE179848_hypoxia")
bio <- nb$r[nb$id %in% BIO]; oth <- nb$r[!nb$id %in% BIO & nb$circ == "independent"]
say("\nbioenergetic / mitochondrial arms (n = %d): median r = %+.3f  [%s]",
    length(bio), median(bio), paste(sprintf("%+.3f", sort(bio, TRUE)), collapse = ", "))
say("all other independent contrasts (n = %d): median r = %+.3f", length(oth), median(oth))
say("Wilcoxon P = %.4f", wilcox.test(bio, oth)$p.value)

## ---- 4. how much of the space do the two hand-built axes explain? --------
Rf <- R; Rf[is.na(Rf)] <- 0
e  <- eigen((Rf + t(Rf))/2, symmetric = TRUE)$values
e  <- e[e > 0]; e <- 100*e/sum(e)
say("\neigenvalue spectrum of the 60x60 matrix: %s", paste(sprintf("%.1f%%", head(e, 6)), collapse = "  "))
ax <- cbind(CP$rho_UVA, CP$rho_senescence)[match(meta$id, CP$id), ]
keep <- which(rowSums(is.na(ax)) == 0)
pc <- cmdscale(as.dist(1 - Rf), k = 2)
say("MDS1 vs rho_senescence: r = %+.3f | MDS1 vs rho_UVA: r = %+.3f",
    cor(pc[keep,1], ax[keep,2]), cor(pc[keep,1], ax[keep,1]))
say("MDS2 vs rho_senescence: r = %+.3f | MDS2 vs rho_UVA: r = %+.3f",
    cor(pc[keep,2], ax[keep,2]), cor(pc[keep,2], ax[keep,1]))

## ---- 5. does the senescence axis alone beat the two-axis model? ----------
utk <- upper.tri(Rf) & is.finite(R)
d1 <- as.matrix(dist(scale(ax[,2, drop = FALSE])))         # senescence only
d2 <- as.matrix(dist(scale(ax)))                            # both axes
say("\nsenescence-axis distance vs similarity: rho = %+.3f",
    cor(d1[utk], R[utk], method = "spearman", use = "complete.obs"))
say("two-axis distance vs similarity:        rho = %+.3f",
    cor(d2[utk], R[utk], method = "spearman", use = "complete.obs"))
