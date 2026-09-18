.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(fgsea))
D <- "public_data_tierA/derived"; O <- file.path(D, "conserved_core")
say <- function(...) cat(sprintf(...), "\n")
res <- read.delim(file.path(O, "gene_sign_consistency.tsv"))
## Reactome gene sets
gmtz <- "public_data_tierA/network/ReactomePathways.gmt.zip"
td <- tempdir(); unzip(gmtz, exdir = td)
gmt <- list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1]
P <- gmtPathways(gmt)
P <- P[sapply(P, length) >= 15 & sapply(P, length) <= 500]
say("Reactome sets: %d", length(P))

for (k in c("senescence","secretome","UV injury","photoprotection / rescue","donor age")) {
  s <- res[res$class == k, ]; if (!nrow(s)) next
  u <- s[s$k_agree == s$n_studies, ]
  up <- u$gene[u$dir > 0]; dn <- u$gene[u$dir < 0]
  uni <- s$gene
  ora <- function(sel, lab) {
    if (length(sel) < 10) return(NULL)
    out <- do.call(rbind, lapply(names(P), function(nm) {
      g <- intersect(P[[nm]], uni); if (length(g) < 10) return(NULL)
      a <- length(intersect(sel, g)); if (a < 3) return(NULL)
      pv <- phyper(a-1, length(g), length(uni)-length(g), length(sel), lower.tail = FALSE)
      data.frame(pathway = sub("__.*","",nm), k = a, set = length(g), p = pv) }))
    if (is.null(out)) return(NULL)
    out$FDR <- p.adjust(out$p, "BH"); out$dir <- lab
    out[order(out$p), ][1:min(6, nrow(out)), ] }
  r <- rbind(ora(up, "up"), ora(dn, "down"))
  say("\n===== %s : %d unanimous of %d tested (%d up, %d down) =====",
      k, nrow(u), nrow(s), length(up), length(dn))
  if (!is.null(r)) print(data.frame(dir=r$dir, pathway=substr(r$pathway,1,52),
        hits=sprintf("%d/%d", r$k, r$set), FDR=signif(r$FDR,2)), row.names=FALSE)
  write.table(u, file.path(O, sprintf("unanimous_%s.tsv", gsub("[^a-z]","",tolower(k)))),
              sep="\t", quote=FALSE, row.names=FALSE)
}
