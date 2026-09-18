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
## Corrected here, where the contrast table is written, so that the 96-gene core
## (defined upstream from the original four ageing series) is not redefined.
FIX <- c(GSE116968_uv_vs_control_1h = "UV injury", GSE116968_uv_vs_control_4h = "UV injury",
         GSE179848_late = "senescence")
hit <- R$id %in% names(FIX); R$class[hit] <- FIX[R$id[hit]]
R$contrast <- N$contrast[match(R$id, N$id)]
R$contrast[R$id == "LOCAL_ReproCM"] <- "reprogramming-phase secretome"
R <- R[is.finite(R$spl) & is.finite(R$cc), ]

ct <- cor.test(R$cc, R$spl, method = "spearman", exact = FALSE)
say("\nacross all %d contrasts: rho(cell-cycle delta, splicing delta) = %+.3f, P = %.2e",
    nrow(R), ct$estimate, ct$p.value)
f <- lm(spl ~ cc, R); R$resid <- resid(f)
say("linear fit: splicing = %+.3f + %.3f x cell-cycle,  R2 = %.3f",
    coef(f)[1], coef(f)[2], summary(f)$r.squared)

S <- R[R$class %in% c("secretome","Repro-CM"), ]
S <- S[order(-S$resid), ]
say("\n=== secretome contrasts: splicing effect beyond what proliferation predicts ===")
print(data.frame(id = S$id, cellcycle = round(S$cc, 3), splicing = round(S$spl, 3),
                 residual = round(S$resid, 3), contrast = substr(S$contrast, 1, 40)),
      row.names = FALSE)
sen <- R[R$class == "senescence", ]
say("\nsenescence contrasts: median cell-cycle %+.3f, splicing %+.3f, residual %+.3f",
    median(sen$cc), median(sen$spl), median(sen$resid))
say("Repro-CM: cell-cycle %+.3f, splicing %+.3f, residual %+.3f  (rank %d of %d secretomes)",
    R$cc[R$id=="LOCAL_ReproCM"], R$spl[R$id=="LOCAL_ReproCM"], R$resid[R$id=="LOCAL_ReproCM"],
    which(S$id == "LOCAL_ReproCM"), nrow(S))
write.table(R, file.path(O, "splicing_vs_proliferation.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
