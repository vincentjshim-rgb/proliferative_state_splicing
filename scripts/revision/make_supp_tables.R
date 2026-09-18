## Supplementary tables for the revised manuscript.
##   S1  datasets and contrasts, with full contrast descriptions (the published
##       version truncated them) and the source publication of each GEO series
##   S4  the GSE113957 cohort sensitivity table
## Outputs: public_data_tierA/derived/figures_sciadv/supp_table{1,4}.tsv and .md
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; O <- file.path(D, "figures_sciadv")
say <- function(...) cat(sprintf(...), "\n")

R <- read.delim(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
R$name <- N$contrast[match(R$id, N$id)]
R$name[is.na(R$name)] <- CP$contrast[match(R$id[is.na(R$name)], CP$id)]
R$series <- sub("^(GSE[0-9]+).*$", "\\1", R$id)
R$series[!grepl("^GSE", R$series)] <- NA
## the secretome contrasts carry their series in the signature table
R$series[is.na(R$series)] <- N$dataset[match(R$id[is.na(R$series)], N$id)]
say("contrasts: %d | series: %d", nrow(R), length(unique(na.omit(R$series))))

## source publications, verified against PubMed/Crossref (scripts/revision/references_v6.md)
SRC <- c(
 GSE93535  = "Lammermann, I. et al. npj Aging Mech. Dis. 4, 4 (2018)",
 GSE109700 = "De Cecco, M. et al. Nature 566, 73-78 (2019)",
 GSE113957 = "Fleischer, J. G. et al. Genome Biol. 19, 221 (2018)",
 GSE139563 = "Lopez-Antona, I. et al. Aging Cell 21, e13580 (2022)",
 GSE149694 = "Liu, X. et al. Nature 586, 101-107 (2020)",
 GSE179848 = "Sturm, G. et al. Sci. Data 9, 751 (2022); Commun. Biol. 6, 22 (2023)",
 GSE191055 = "Hasegawa, T. et al. Cell 186, 1417-1431.e20 (2023)",
 GSE222414 = "Yang, T. et al. MedComm 5, e625 (2024)",
 GSE251807 = "Lei, R., da Silva, T. B., Cui, Z. & Ye, H. Sci. Rep. 15, 19383 (2025)",
 GSE266052 = "Brizio, M. et al. Stem Cell Res. Ther. 15, 329 (2024)",
 GSE268248 = "Imani, A., Panahipour, L., Kuhtreiber, H., Mildner, M. & Gruber, R. Cells 13, 1308 (2024)",
 GSE279804 = "Yao, T. et al. J. Cell. Mol. Med. 29, e70877 (2025)",
 GSE282054 = "Abdelmohsen, K. et al. Aging Cell 25, e70368 (2026)",
 GSE293186 = "Yuan, H. et al. J. Invest. Dermatol. 146, 1369-1381.e11 (2026)",
 GSE297233 = "Lu, J. Y. et al. Cell 188, 5895-5911.e17 (2025)",
 GSE306957 = "Victorelli, S. et al. Nat. Commun. 16, 10992 (2025)",
 GSE307377 = "Tanaka, H. et al. bioRxiv https://doi.org/10.1101/2025.09.17.676776 (2025)",
 GSE326951 = "Zhang, X. et al. Pharm. Biol. 64, 764-782 (2026)",
 GSE329475 = "Zhong, J. et al. Cells 15, 1497 (2026)",
 GSE116968 = "no associated publication",
 GSE134533 = "no associated publication",
 GSE240486 = "no associated publication",
 GSE306748 = "no associated publication")
R$source <- SRC[R$series]
miss <- unique(R$series[is.na(R$source)])
if (length(miss)) say("WARNING: no source recorded for %s", paste(miss, collapse = ", "))

CLS <- c("donor age" = "donor age", "senescence" = "senescence", "UV injury" = "UV injury",
         "photoprotection / rescue" = "photoprotection", "reprogramming" = "reprogramming media",
         "metabolic / culture" = "culture and metabolic", "secretome" = "secretome",
         "Repro-CM" = "secretome")
R$class_label <- CLS[R$class]
R <- R[order(R$class_label, R$series, R$id), ]
S1 <- data.frame(accession = R$series, class = R$class_label, contrast = R$name,
                 cellcycle_change = round(R$cc, 3), splicing_change = round(R$spl, 3),
                 source_publication = R$source)
write.table(S1, file.path(O, "supp_table1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
say("S1: %d rows, %d distinct series", nrow(S1), length(unique(S1$accession)))

## ---- S4: cohort sensitivity -------------------------------------------------
M <- read.delim(file.path(D, "cohort_revised/sensitivity_metrics.tsv"))
G <- read.delim(file.path(D, "cohort_revised/sensitivity_genes.tsv"))
P <- read.delim(file.path(D, "psi_cohorts/psi_events_by_cohort.tsv"))
lab <- c(published = "deposited 142 samples", normal = "133 normal donors",
         primary = "107 normal adults 20+ (primary)", adult2082 = "76 normal adults 20-82")
mk <- function(metric, cols) { x <- M[M$metric == metric, ]; x$cohort <- lab[x$cohort]; x[, cols] }
S4 <- data.frame(
  cohort = lab[G$cohort], n = G$n,
  rho_age_proliferation = round(G$rho_age_prolif, 2),
  machinery_beta = round(M$beta[M$metric == "machinery"][match(G$cohort, M$cohort[M$metric == "machinery"])], 3),
  machinery_P = signif(M$p[M$metric == "machinery"][match(G$cohort, M$cohort[M$metric == "machinery"])], 2),
  machinery_beta_adj = round(M$beta_adj[M$metric == "machinery"][match(G$cohort, M$cohort[M$metric == "machinery"])], 3),
  machinery_P_adj = signif(M$p_adj[M$metric == "machinery"][match(G$cohort, M$cohort[M$metric == "machinery"])], 2),
  pct_of_age_effect_lost = round(M$pct_lost[M$metric == "machinery"][match(G$cohort, M$cohort[M$metric == "machinery"])]),
  unannotated_junction_beta = signif(M$beta[M$metric == "unannot_reads"][match(G$cohort, M$cohort[M$metric == "unannot_reads"])], 2),
  unannotated_junction_P = signif(M$p[M$metric == "unannot_reads"][match(G$cohort, M$cohort[M$metric == "unannot_reads"])], 2),
  genes_age_associated = G$age_assoc, genes_keeping_decline = G$decline_kept,
  psi_age_events = P$age_events[match(lab[G$cohort], c("deposited 142 samples", "133 normal donors",
     "107 normal adults 20+ (primary)", "76 normal adults 20-82")[match(P$cohort, P$cohort)])])
S4$psi_age_events <- P$age_events[c(1, 2, 3, 4)][match(G$cohort, c("published", "normal", "primary", "adult2082"))]
write.table(S4, file.path(O, "supp_table4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(S4[, c("cohort", "n", "machinery_P", "machinery_P_adj", "unannotated_junction_P",
             "genes_keeping_decline", "psi_age_events")], row.names = FALSE)

## markdown versions for the manuscript
md <- function(df) {
  h <- paste("|", paste(gsub("_", " ", names(df)), collapse = " | "), "|")
  s <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
  rows <- apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
  c(h, s, rows) }
writeLines(md(S1), file.path(O, "supp_table1.md"))
writeLines(md(S4), file.path(O, "supp_table4.md"))
say("wrote supp_table1.tsv/.md and supp_table4.tsv/.md")
