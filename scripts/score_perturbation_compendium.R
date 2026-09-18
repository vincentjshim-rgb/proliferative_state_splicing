#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Score every compendium signature on the two fixed reference axes and on D.
# Signatures that contributed to a reference axis are flagged; they are shown
# as positive/negative controls, never as independent evidence.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"; der <- file.path(root, "derived")
out  <- file.path(der, "compendium")
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")

AX <- read.delim(file.path(der, "decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
AX <- AX[is.finite(AX$UVA_meta) & is.finite(AX$senescence_meta) &
         is.finite(AX$mean_log2_expression), ]
cat("reference gene universe:", nrow(AX), "\n")

## ---- harvest ---------------------------------------------------------------
sigs <- list()
for (f in list.files(file.path(out, "signatures"), full.names = TRUE)) {
  d <- read.delim(f)
  sigs[[sub("\\.tsv$", "", basename(f))]] <-
    list(dataset = d$dataset[1], contrast = d$contrast[1],
         eff = data.frame(gene = d$gene, logFC = d$logFC))
}
legacy <- list(
  GSE240226_UVA       = c("figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv", "GSE240226", "UVA vs control"),
  GSE240226_Maifuyin  = c("figure1_public_uva/GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv", "GSE240226", "Maifuyin rescue vs UVA"),
  GSE240226_succinate = c("figure1_public_uva/GSE240226_succinate_rescue_vs_UVA_DE.tsv", "GSE240226", "succinate rescue vs UVA"),
  GSE302943_UVA       = c("figure1_public_uva/GSE302943_UVA_vs_control_DE.tsv", "GSE302943", "cumulative UVA vs control"),
  GSE125429_UVA       = c("figure1_public_uva/GSE125429_UVA_vs_control_DE.tsv", "GSE125429", "chronic UVA vs control"),
  GSE89005_single6h   = c("figure1_public_uva/GSE89005_single_6h_UVA_vs_sham_DE.tsv", "GSE89005", "single UVA 6 h"),
  GSE89005_single24h  = c("figure1_public_uva/GSE89005_single_24h_UVA_vs_sham_DE.tsv", "GSE89005", "single UVA 24 h"),
  GSE89005_repeat24h  = c("figure1_public_uva/GSE89005_repeated_24h_UVA_vs_sham_DE.tsv", "GSE89005", "repeated UVA 24 h"),
  GSE109700_deep      = c("figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv", "GSE109700", "deep replicative senescence"),
  GSE109700_early     = c("figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv", "GSE109700", "early replicative senescence"),
  GSE191055_P27       = c("figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv", "GSE191055", "P27 vs P4"),
  GSE93535_SIPS       = c("figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv", "GSE93535", "SIPS vs quiescent"),
  GSE93535_rescueSIPS = c("figure2_public_aging/GSE93535_SIPS1201_vs_SIPS_DE.tsv", "GSE93535", "compound 1201 in SIPS"),
  GSE93535_rescueQ    = c("figure2_public_aging/GSE93535_Q1201_vs_Q_DE.tsv", "GSE93535", "compound 1201 in quiescent"),
  GSE179848_late      = c("figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv", "GSE179848", "late vs early passage"),
  GSE113957_age       = c("figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv", "GSE113957", "donor age, per decade"),
  GSE226189_age       = c("figure2_public_aging/GSE226189_age_per_decade_DE.tsv", "GSE226189", "donor age, per decade"),
  GSE165177_MPTR      = c("figure2_mptr/GSE165177_MPTR_paired_DE.tsv", "GSE165177", "transient reprogramming (MPTR)"))
for (id in names(legacy)) {
  f <- file.path(der, legacy[[id]][1]); if (!file.exists(f)) next
  d <- read.delim(f)
  sigs[[id]] <- list(dataset = legacy[[id]][2], contrast = legacy[[id]][3],
                     eff = data.frame(gene = toupper(d$gene), logFC = d$logFC))
}
loc <- AX[AX$comparator_consistent & AX$Repro_specific_score != 0, ]
sigs[["LOCAL_ReproCM"]] <- list(dataset = "this study", contrast = "Repro-CM anchor (n = 1 per condition)",
                                eff = data.frame(gene = loc$gene, logFC = loc$Repro_specific_score))
cat("signatures harvested:", length(sigs), "\n")

## ---- scoring ---------------------------------------------------------------
IN_UVA <- c("GSE240226_UVA", "GSE302943_UVA")
IN_SEN <- c("GSE179848_late", "GSE109700_deep", "GSE191055_P27", "GSE93535_SIPS")
AXIS_STUDIES <- c("GSE240226", "GSE302943", "GSE179848", "GSE109700", "GSE191055", "GSE93535")
NPERM <- 2000L
score <- function(id) {
  s <- sigs[[id]]
  j <- merge(AX[, c("gene", "UVA_meta", "senescence_meta", "mean_log2_expression")], s$eff, by = "gene")
  j <- j[is.finite(j$logFC), ]
  if (nrow(j) < 500) return(NULL)
  ru <- cor(j$logFC, j$UVA_meta, method = "spearman")
  rs <- cor(j$logFC, j$senescence_meta, method = "spearman")
  set.seed(7)
  br <- unique(quantile(j$mean_log2_expression, seq(0, 1, 0.1)))
  idx <- split(seq_len(nrow(j)), cut(j$mean_log2_expression, br, include.lowest = TRUE, labels = FALSE))
  rl <- rank(j$logFC)
  nullD <- replicate(NPERM, { p <- rl; for (i in idx) p[i] <- sample(p[i])
    cor(p, j$UVA_meta, method = "spearman") - cor(p, j$senescence_meta, method = "spearman") })
  D <- ru - rs
  data.frame(id = id, dataset = s$dataset, contrast = s$contrast, n_genes = nrow(j),
             rho_UVA = ru, rho_senescence = rs, D = D,
             p_matched = (1 + sum(abs(nullD) >= abs(D))) / (NPERM + 1),
             null_q025 = quantile(nullD, 0.025), null_q975 = quantile(nullD, 0.975),
             circularity = ifelse(id %in% c(IN_UVA, IN_SEN), "IN reference axis",
                           ifelse(s$dataset %in% AXIS_STUDIES, "same study as an axis", "independent")))
}
R <- do.call(rbind, lapply(names(sigs), score))
R$FDR <- p.adjust(R$p_matched, "BH")
R <- R[order(-R$D), ]
w(R, "compendium_D_scores.tsv")
cat("\n################ PERTURBATION COMPENDIUM, ranked by D ################\n")
print(R[, c("dataset", "contrast", "n_genes", "rho_UVA", "rho_senescence", "D", "FDR", "circularity")],
      row.names = FALSE, digits = 3)
