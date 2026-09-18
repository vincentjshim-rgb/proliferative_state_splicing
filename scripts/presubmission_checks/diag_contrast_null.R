## Is the effect of a secretome on the pre-mRNA processing core separable from
## its effect on proliferation?  Contact inhibition shows the core also falls
## with reversible arrest, so this control is required.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(fgsea))
D <- "public_data_tierA/derived"; O <- file.path(D, "conserved_core")
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(O, "age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
gmtz <- "public_data_tierA/network/ReactomePathways.gmt.zip"; td <- tempdir(); unzip(gmtz, exdir = td)
P <- gmtPathways(list.files(td, pattern="\\.gmt$", full.names=TRUE)[1])
CC <- unique(unlist(P[grep("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(P))]))
say("cell-cycle module: %d genes | splicing core: %d genes | overlap: %d",
    length(CC), length(CORE), length(intersect(CC, CORE)))

Z <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); M <- Z$meta
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D,"compendium/signatures"), full.names=TRUE)) {
  d <- read.delim(f); vec[[sub("\\.tsv$","",basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
lgf <- list(GSE109700_deep="figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
 GSE179848_late="figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
 GSE113957_age="figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
 ## GSE226189 excluded from the study (processing batches coincide with age strata)
 GSE93535_SIPS="figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
 GSE191055_P27="figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
 GSE109700_early="figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv")
for (id in names(lgf)) { f <- file.path(D, lgf[[id]]); if (!file.exists(f)) next
  d <- read.delim(f); vec[[id]] <- setNames(d$logFC, toupper(d$gene)) }
loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))

dlt <- function(v, set) { g <- intersect(names(v), set); if (length(g) < 15) return(NA)
  bg <- v[is.finite(v)]; mean(bg[names(bg) %in% set]) - mean(bg) }
R <- do.call(rbind, lapply(names(vec), function(i)
  data.frame(id = i, spl = dlt(vec[[i]], CORE), cc = dlt(vec[[i]], CC))))
R$class <- M$class[match(R$id, M$id)]
## class corrections (2026-09-13). The compendium metadata placed two UVB-vs-control
## contrasts in "metabolic / culture" (empty description fell to the default) and the
## late-vs-early passage series in "donor age" (the pattern "age" matched "passage").
## ---------------- diagnostic: is the splicing-set coupling special? ----------------
set.seed(20260915)
R <- R[is.finite(R$spl) & is.finite(R$cc), ]
obs_rho <- cor(R$cc, R$spl, method = "spearman"); obs_r2 <- summary(lm(spl ~ cc, R))$r.squared
cat(sprintf("\nobserved: rho = %.3f, R2 = %.3f, n = %d contrasts\n", obs_rho, obs_r2, nrow(R)))
## expression level in cultured dermal fibroblasts (GSE113957 FPKM, adults) for matching
fp <- read.delim(gzfile("public_data_tierA/aging/GSE113957_fpkm.txt.gz"), check.names = FALSE)
sym <- toupper(sub("\\|.*$", "", fp[["Annotation/Divergence"]]))
X <- as.matrix(fp[, grepl("^[0-9]+_[0-9]+", names(fp))])
lev <- tapply(rowMeans(log2(X + 1)), sym, max); lev <- lev[names(lev) != "" & names(lev) != "-"]
cat(sprintf("expression levels for %d gene symbols (splicing set covered: %d of %d)\n", length(lev), length(intersect(CORE, names(lev))), length(CORE)))
meas <- table(unlist(lapply(vec[R$id], function(v) names(v)[is.finite(v)])))
univ <- intersect(names(meas)[meas >= 0.8 * nrow(R)], names(lev))
cat(sprintf("genes measured in >= 80%% of contrasts and in GSE113957: %d (splicing set: %d, cell-cycle set: %d)\n",
    length(univ), length(intersect(univ, CORE)), length(intersect(univ, CC))))
lv <- lev[univ]
bins <- cut(lv, quantile(lv, seq(0, 1, 0.05)), include.lowest = TRUE, labels = FALSE); names(bins) <- univ
core_b <- table(bins[intersect(CORE, univ)])
pool <- setdiff(univ, c(CORE, CC))
cat("pool size:", length(pool), "| smallest pool bin used:", min(table(bins[pool])[names(core_b)]), "| largest core bin:", max(core_b), "\n")
draw <- function() unlist(lapply(names(core_b), function(b) sample(pool[bins[pool] == as.integer(b)], core_b[[b]])))
nullstat <- t(replicate(1000, { g <- draw(); s <- sapply(R$id, function(i) dlt(vec[[i]], g))
  ok <- is.finite(s); c(rho = cor(R$cc[ok], s[ok], method = "spearman"), r2 = summary(lm(s[ok] ~ R$cc[ok]))$r.squared) }))
cat(sprintf("expression-matched random sets (177 genes, excluding splicing and cell-cycle genes), 1000 draws:\n  rho median %.2f [95%%: %.2f to %.2f], max %.2f\n  R2  median %.2f [95%%: %.2f to %.2f], max %.2f\n  empirical P(rho >= observed) = %.3f\n",
    median(nullstat[, 1]), quantile(nullstat[, 1], .025), quantile(nullstat[, 1], .975), max(nullstat[, 1]),
    median(nullstat[, 2]), quantile(nullstat[, 2], .025), quantile(nullstat[, 2], .975), max(nullstat[, 2]),
    mean(nullstat[, 1] >= obs_rho)))
## other programmes for context (same score, same contrasts)
ctx <- c("Translation", "Metabolism of RNA", "Collagen formation", "Extracellular matrix organization", "Unfolded Protein Response (UPR)", "Cellular Senescence")
for (nm in ctx) { g <- setdiff(P[[nm]], CC); s <- sapply(R$id, function(i) dlt(vec[[i]], g)); ok <- is.finite(s)
  cat(sprintf("  %-36s %3d genes (cell-cycle genes removed): rho = %+.2f, R2 = %.2f\n", nm, length(g), cor(R$cc[ok], s[ok], method = "spearman"), summary(lm(s[ok] ~ R$cc[ok]))$r.squared)) }
## leverage: without the senescence class, and one contrast per study
nos <- R$class != "senescence"
cat(sprintf("without senescence contrasts: rho = %.2f, R2 = %.2f, n = %d\n", cor(R$cc[nos], R$spl[nos], method = "spearman"), summary(lm(spl ~ cc, R[nos, ]))$r.squared, sum(nos)))
R$study <- sub("_.*$", "", R$id)
St <- aggregate(cbind(spl, cc) ~ study, data = R, FUN = mean)
cat(sprintf("one averaged value per study: rho = %.2f, R2 = %.2f, n = %d\n", cor(St$cc, St$spl, method = "spearman"), summary(lm(spl ~ cc, St))$r.squared, nrow(St)))
