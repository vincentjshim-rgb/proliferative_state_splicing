## Is between-study coherence explained by effect magnitude rather than by class?
## If a class's contrasts simply have small effects, attenuation alone lowers r.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; OUT <- file.path(D, "secretome_class")
say <- function(...) cat(sprintf(...), "\n"); set.seed(7)
Z <- readRDS(file.path(OUT, "structure_with_secretome.rds")); R <- Z$R; M <- Z$meta

## rebuild vectors to get magnitudes
NEW <- read.delim(file.path(OUT, "secretome_signatures.tsv"))
AX  <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
vec <- lapply(split(NEW, NEW$id), function(d) setNames(d$logFC, d$gene))
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

M$mag <- sapply(M$id, function(i) { v <- vec[[i]]; if (is.null(v)) NA else
  median(abs(v[is.finite(v)]), na.rm = TRUE) })
say("=== effect magnitude by class (median |logFC| of the contrast) ===")
ag <- aggregate(mag ~ class, M, function(x) c(median = median(x), min = min(x), max = max(x)))
ag <- data.frame(class = ag$class, as.data.frame(ag$mag))
CL <- read.delim(file.path(OUT, "between_study_coherence.tsv"))
ag$coh <- CL$median_r[match(ag$class, CL$class)]
ag <- ag[order(-ag$coh), ]; print(ag, digits = 3, row.names = FALSE)
say("\nclass-level: rho(magnitude, coherence) = %+.3f",
    cor(ag$median, ag$coh, method = "spearman", use = "complete.obs"))

## pair-level: does |logFC| of the two members predict their similarity?
pr <- list()
for (k in unique(M$class)) {
  idx <- which(M$class == k); if (length(idx) < 2) next
  for (a in seq_along(idx)) for (b in seq_along(idx)) { if (b <= a) next
    i <- idx[a]; j <- idx[b]
    if (M$dataset[i] == M$dataset[j] || !is.finite(R[i,j])) next
    pr[[length(pr)+1]] <- data.frame(class = k, r = R[i,j],
      gmag = sqrt(M$mag[i] * M$mag[j]), minmag = min(M$mag[i], M$mag[j])) } }
P <- do.call(rbind, pr)
say("\n=== pair level (n = %d between-study pairs) ===", nrow(P))
ct <- cor.test(P$gmag, P$r, method = "spearman", exact = FALSE)
say("rho(geometric-mean magnitude, similarity) = %+.3f, P = %.2e", ct$estimate, ct$p.value)
m1 <- lm(r ~ gmag, P); m2 <- lm(r ~ gmag + factor(class), P)
say("R2  magnitude only            = %.3f", summary(m1)$r.squared)
say("R2  magnitude + class         = %.3f", summary(m2)$r.squared)
say("class effect after magnitude: F = %.1f, P = %.2e",
    anova(m1, m2)$F[2], anova(m1, m2)$`Pr(>F)`[2])

## magnitude-matched comparison: senescence vs the rest within overlapping magnitude
rng <- range(P$gmag[P$class == "secretome"])
sub <- P[P$gmag >= rng[1] & P$gmag <= rng[2], ]
say("\nrestricted to the secretome magnitude range (%.3f-%.3f), n = %d pairs:", rng[1], rng[2], nrow(sub))
ags <- aggregate(r ~ class, sub, function(x) c(n = length(x), med = median(x)))
print(data.frame(class = ags$class, as.data.frame(ags$r)), digits = 3, row.names = FALSE)
sen <- sub$r[sub$class == "senescence"]
if (length(sen) > 1) { for (k in setdiff(unique(sub$class), "senescence")) {
  o <- sub$r[sub$class == k]; if (length(o) < 2) next
  say("  senescence vs %-24s Wilcoxon P = %.4f", k, wilcox.test(sen, o)$p.value) } }
write.table(P, file.path(OUT, "pairwise_magnitude_vs_similarity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(ag, file.path(OUT, "class_magnitude_vs_coherence.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## within-study positive control: does the pipeline detect secretome effects at all?
say("\n=== within-study control (same experiment, different secretome fraction) ===")
for (ds in c("GSE251807","GSE266052","GSE306748","GSE139563")) {
  idx <- which(M$dataset == ds & M$class == "secretome")
  if (length(idx) < 2) next
  s <- R[idx, idx]; v <- s[upper.tri(s)]
  say("  %-10s n contrasts = %d, within-study median r = %+.3f  [%s]", ds, length(idx),
      median(v, na.rm = TRUE), paste(sprintf("%+.2f", sort(v, TRUE)), collapse = ", ")) }
