## GO / Reactome enrichment of genes whose splicing changes with donor age,
## against the universe of genes with a testable event.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(fgsea))
D <- "public_data_tierA/derived"; O <- file.path(D, "psi_age")
say <- function(...) cat(sprintf(...), "\n")
sel <- readLines(file.path(O, "genes_age_spliced.txt"))
uni <- readLines(file.path(O, "genes_universe.txt"))
say("selected %d genes from a universe of %d", length(sel), length(uni))
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
P <- P[sapply(P, length) >= 15 & sapply(P, length) <= 600]
ora <- do.call(rbind, lapply(names(P), function(nm) {
  g <- intersect(P[[nm]], uni); if (length(g) < 12) return(NULL)
  a <- length(intersect(sel, g)); if (a < 4) return(NULL)
  data.frame(pathway = sub("__.*", "", nm), hits = a, set = length(g),
             expected = length(g)*length(sel)/length(uni),
             p = phyper(a-1, length(g), length(uni)-length(g), length(sel), lower.tail = FALSE)) }))
ora$FDR <- p.adjust(ora$p, "BH"); ora$fold <- ora$hits/ora$expected
ora <- ora[order(ora$p), ]
say("\npathways at BH-FDR < 0.05: %d", sum(ora$FDR < 0.05))
print(head(data.frame(pathway = substr(ora$pathway, 1, 56), hits = ora$hits, set = ora$set,
                      fold = round(ora$fold, 2), FDR = signif(ora$FDR, 2)), 20), row.names = FALSE)
write.table(ora, file.path(O, "reactome_ora_spliced.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
