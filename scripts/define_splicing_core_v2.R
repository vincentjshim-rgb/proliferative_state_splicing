## Splicing gene set, version 2 (2026-09-13).
## Genes that fall in every one of three cultured-fibroblast ageing series, intersected with
## Reactome pre-mRNA processing. GSE226189 is excluded from the study because its processing
## batches coincide with age strata (run_core_depth2.R).
## Same rule as run_conserved_core.R + run_splicing_secretome_link.R (a gene must be measured in
## every series and fall in every one), but the series are named here explicitly instead of being
## taken from compendium class labels, which mislabelled the passage series as donor age.
## The 96-gene version-1 set (age_down_splicing_core.txt) is kept unchanged: the GTEx
## preregistration fixed that file.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(fgsea))
D <- "public_data_tierA/derived"; O <- file.path(D, "conserved_core")
say <- function(...) cat(sprintf(...), "\n")
SERIES <- c(GSE113957_age  = "figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
            GSE179848_late = "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
            GSE307377_age  = "compendium/signatures/GSE307377_age.tsv")
vec <- lapply(SERIES, function(f) { d <- read.delim(file.path(D, f)); setNames(d$logFC, toupper(d$gene)) })
gg  <- Reduce(union, lapply(vec, names))
mat <- sapply(vec, function(v) v[gg]); rownames(mat) <- gg
mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]          # measured in every series
DOWN <- rownames(mat)[rowSums(mat < 0) == ncol(mat)]                         # falls in every series
say("series: %s", paste(names(SERIES), collapse = ", "))
say("genes measured in all %d series: %d | falling in all: %d | expected by chance: %.0f (%.2f-fold)",
    ncol(mat), nrow(mat), length(DOWN), nrow(mat) * 0.5^ncol(mat), length(DOWN) / (nrow(mat) * 0.5^ncol(mat)))

td <- tempfile(); dir.create(td); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
spl <- unique(unlist(P[grep("mRNA Splicing|Processing of Capped Intron|3'-end processing", names(P))]))
CORE <- intersect(DOWN, spl)
say("splicing set v2: %d genes (v1 had 96; shared %d)", length(CORE),
    length(intersect(CORE, readLines(file.path(O, "age_down_splicing_core.txt")))))

## Reactome over-representation among genes falling in every series (sets of 15-500 genes)
uni <- rownames(mat)
PP  <- P[sapply(P, function(s) { n <- length(intersect(s, uni)); n >= 15 && n <= 500 })]
EN  <- do.call(rbind, lapply(names(PP), function(nm) { s <- intersect(PP[[nm]], uni); k <- length(intersect(DOWN, s))
  data.frame(pathway = sub("__.*", "", nm), k = k, size = length(s),
             fold = (k / length(DOWN)) / (length(s) / length(uni)),
             p = phyper(k - 1, length(s), length(uni) - length(s), length(DOWN), lower.tail = FALSE)) }))
EN$FDR <- p.adjust(EN$p, "BH"); EN <- EN[order(EN$p), ]; EN$rank <- seq_len(nrow(EN))
say("\ntop Reactome pathways among genes falling in all three series:")
print(head(transform(EN, fold = round(fold, 2), FDR = signif(FDR, 2))[, c("rank","pathway","k","size","fold","FDR")], 14), row.names = FALSE)

writeLines(CORE, file.path(O, "age_down_splicing_core_v2.txt"))
writeLines(DOWN, file.path(O, "age_down_unanimous_v2.txt"))
write.table(EN, file.path(O, "age_core_v2_reactome.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
