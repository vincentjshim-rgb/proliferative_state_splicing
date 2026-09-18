## Reviewer 1, major 7b: the senescence result should not rest on one curated set.
## SenMayo (Saul et al. 2022, MSigDB SAUL_SEN_MAYO) is the panel researchers use most
## often, and it is scored here exactly as the Reactome Cellular Senescence set was:
## against the measured division rate in the 328 counted libraries, and dissected by
## whether each gene also belongs to the Reactome cell-cycle union.
## Outputs: public_data_tierA/derived/reviewer_sensitivities/senmayo_*.tsv
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(fgsea)})
D <- "public_data_tierA/derived"; O <- file.path(D, "reviewer_sensitivities")
say <- function(...) cat(sprintf(...), "\n")

SM <- readLines("public_data_tierA/senescence/senmayo/SAUL_SEN_MAYO.v2026.1.Hs.grp")
SM <- toupper(SM[!grepl("^#|^SAUL_SEN_MAYO$", SM)]); SM <- SM[nzchar(SM)]
say("SenMayo: %d genes", length(SM))

IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z; rate <- IN$meta$rate; P <- IN$pathways
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
RE <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
CC <- unique(unlist(RE[grep("^Cell Cycle|^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(RE))]))
REACT <- unique(unlist(RE[grep("^Cellular Senescence", names(RE))]))

score <- function(g) { gg <- intersect(rownames(z), g); colMeans(z[gg, , drop = FALSE]) }
tab <- do.call(rbind, lapply(list(list("SenMayo (Saul 2022)", SM), list("Reactome Cellular Senescence", REACT)), function(x) {
  nm <- x[[1]]; g <- intersect(rownames(z), x[[2]])
  s <- score(g); ct <- cor.test(s, rate, method = "spearman", exact = FALSE)
  gnoCC <- setdiff(g, CC); s2 <- score(gnoCC); ct2 <- cor.test(s2, rate, method = "spearman", exact = FALSE)
  gw <- apply(z[g, , drop = FALSE], 1, function(v) cor(v, rate, method = "spearman"))
  inCC <- g %in% CC
  w <- suppressWarnings(wilcox.test(gw[inCC], gw[!inCC]))
  data.frame(set = nm, n_measured = length(g), n_cellcycle = sum(inCC),
             rho_set = unname(ct$estimate), p_set = ct$p.value,
             rho_set_noCC = unname(ct2$estimate), n_noCC = length(gnoCC),
             median_rho_cellcycle_genes = median(gw[inCC]), median_rho_other_genes = median(gw[!inCC]),
             wilcoxon_p = w$p.value, row.names = NULL) }))
print(tab, digits = 3, row.names = FALSE)
write.table(tab, file.path(O, "senmayo_division_rate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## the canonical effectors inside SenMayo, and the overlap between the two panels
say("\noverlap between the two panels: %d genes", length(intersect(SM, REACT)))
canon <- c("CDKN1A", "CDKN2A", "CDKN2B", "ATM", "LMNB1", "E2F1", "IL6", "IL1A", "IL1B", "CXCL8", "SERPINE1", "MMP3", "IGFBP3")
gw <- sapply(intersect(canon, rownames(z)), function(g) cor(z[g, ], rate, method = "spearman"))
say("canonical genes present in this matrix: %s",
    paste(sprintf("%s %+.2f", names(gw), gw), collapse = ", "))
say("of these, %d are SenMayo members: %s", sum(names(gw) %in% SM), paste(intersect(names(gw), SM), collapse = ", "))
write.table(data.frame(gene = names(gw), rho_rate = gw, in_senmayo = names(gw) %in% SM,
                       in_reactome = names(gw) %in% REACT),
            file.path(O, "senmayo_canonical_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
