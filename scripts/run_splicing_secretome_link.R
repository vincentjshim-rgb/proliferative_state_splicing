## Do secretome interventions move the age-associated pre-mRNA processing core?
## SUPERSEDED for the manuscript: this script wrote the version-1 96-gene set, which is kept only because the
## GTEx preregistration fixed it. The current set is written by define_splicing_core_v2.R.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(fgsea)})
D <- "public_data_tierA/derived"; O <- file.path(D, "conserved_core")
say <- function(...) cat(sprintf(...), "\n"); set.seed(5)
res <- read.delim(file.path(O, "gene_sign_consistency.tsv"))
Z   <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); M <- Z$meta
N   <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
AX  <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D,"compendium/signatures"), full.names=TRUE)) {
  d <- read.delim(f); vec[[sub("\\.tsv$","",basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
lgf <- list(GSE109700_deep="figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
 GSE179848_late="figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
 GSE113957_age="figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
 GSE226189_age="figure2_public_aging/GSE226189_age_per_decade_DE.tsv",
 GSE165177_MPTR="figure2_mptr/GSE165177_MPTR_paired_DE.tsv",
 GSE240226_UVA="figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv",
 GSE93535_SIPS="figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
 GSE191055_P27="figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv")
for (id in names(lgf)) { f <- file.path(D, lgf[[id]]); if (!file.exists(f)) next
  d <- read.delim(f); vec[[id]] <- setNames(d$logFC, toupper(d$gene)) }
loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))

## ---- the two cores ----------------------------------------------------------
ag <- res[res$class == "donor age" & res$k_agree == res$n_studies, ]
AGE_DN <- ag$gene[ag$dir < 0]; AGE_UP <- ag$gene[ag$dir > 0]
sc <- res[res$class == "secretome" & res$k_agree == res$n_studies, ]
SEC_UP <- sc$gene[sc$dir > 0]; SEC_DN <- sc$gene[sc$dir < 0]
say("age core: %d down, %d up | secretome core: %d up, %d down",
    length(AGE_DN), length(AGE_UP), length(SEC_UP), length(SEC_DN))
say("overlap age-down x secretome-up: %d genes", length(intersect(AGE_DN, SEC_UP)))

## splicing subset of the age-down core
gmtz <- "public_data_tierA/network/ReactomePathways.gmt.zip"
td <- tempdir(); unzip(gmtz, exdir = td)
P <- gmtPathways(list.files(td, pattern="\\.gmt$", full.names=TRUE)[1])
spl <- unique(unlist(P[grep("mRNA Splicing|Processing of Capped Intron|3'-end processing", names(P))]))
AGE_SPL <- intersect(AGE_DN, spl)
say("age-down core genes in pre-mRNA processing: %d", length(AGE_SPL))
writeLines(AGE_SPL, file.path(O, "age_down_splicing_core.txt"))

## ---- score every contrast on the age cores ----------------------------------
score <- function(v, set) { g <- intersect(names(v), set); if (length(g) < 15) return(c(NA, NA))
  bg <- v[is.finite(v)]; s <- bg[names(bg) %in% set]
  c(mean(s) - mean(bg), wilcox.test(s, bg[!names(bg) %in% set])$p.value) }
ids <- names(vec)
R <- do.call(rbind, lapply(ids, function(i) {
  a <- score(vec[[i]], AGE_SPL); b <- score(vec[[i]], SEC_UP)
  data.frame(id = i, spl_delta = a[1], spl_p = a[2], sec_delta = b[1], sec_p = b[2]) }))
R$class <- M$class[match(R$id, M$id)]
R$contrast <- N$contrast[match(R$id, N$id)]
R$contrast[R$id=="LOCAL_ReproCM"] <- "reprogramming-phase secretome"
R$spl_FDR <- p.adjust(R$spl_p, "BH")
write.table(R, file.path(O, "contrasts_on_age_cores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

say("\n=== movement of the age-associated pre-mRNA processing core (%d genes) ===", length(AGE_SPL))
say("positive = the contrast raises these genes, i.e. opposite to what ageing does\n")
for (k in c("donor age","senescence","secretome","Repro-CM","photoprotection / rescue",
            "UV injury","reprogramming","metabolic / culture")) {
  s <- R[R$class == k & is.finite(R$spl_delta), ]
  if (!nrow(s)) next
  say("%-24s n=%2d  median delta = %+.4f  (range %+.3f to %+.3f)  FDR<0.05: %d",
      k, nrow(s), median(s$spl_delta), min(s$spl_delta), max(s$spl_delta), sum(s$spl_FDR < 0.05, na.rm=TRUE)) }
say("\n--- secretome contrasts, ranked by how much they raise the age-down splicing core ---")
ss <- R[R$class %in% c("secretome","Repro-CM") & is.finite(R$spl_delta), ]
ss <- ss[order(-ss$spl_delta), ]
print(data.frame(id=ss$id, delta=round(ss$spl_delta,4), FDR=signif(ss$spl_FDR,2),
                 contrast=substr(ss$contrast,1,44)), row.names=FALSE)
