#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Does the anchor simply align with whatever perturbs the cell hardest?
# GSE326951 (sauchinone) gave rho = +0.373 with a median |logFC| of 0.363,
# three times any other photoprotection contrast. If alignment tracks
# perturbation magnitude, the photoprotection claim is much weaker than it looks.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"; der <- file.path(root, "derived")
out  <- file.path(der, "magnitude_vs_alignment")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)

AX <- read.delim(file.path(der, "decoupling_validation/primary_gene_axis_matrix.tsv"), check.names = FALSE)
anchor <- AX[AX$comparator_consistent & AX$Repro_specific_score != 0 & is.finite(AX$Repro_specific_score), ]
D <- read.delim(file.path(der, "compendium", "compendium_D_scores.tsv"))

## every compendium signature, plus the legacy DE tables the scorer also used
sigdir <- file.path(der, "compendium", "signatures")
load_sig <- function(id) {
  f <- file.path(sigdir, paste0(id, ".tsv"))
  if (file.exists(f)) { d <- read.delim(f); return(data.frame(gene = d$gene, logFC = d$logFC)) }
  NULL
}
legacy <- c(GSE240226_UVA = "figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv",
  GSE240226_Maifuyin = "figure1_public_uva/GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv",
  GSE240226_succinate = "figure1_public_uva/GSE240226_succinate_rescue_vs_UVA_DE.tsv",
  GSE302943_UVA = "figure1_public_uva/GSE302943_UVA_vs_control_DE.tsv",
  GSE125429_UVA = "figure1_public_uva/GSE125429_UVA_vs_control_DE.tsv",
  GSE89005_single6h = "figure1_public_uva/GSE89005_single_6h_UVA_vs_sham_DE.tsv",
  GSE89005_single24h = "figure1_public_uva/GSE89005_single_24h_UVA_vs_sham_DE.tsv",
  GSE89005_repeat24h = "figure1_public_uva/GSE89005_repeated_24h_UVA_vs_sham_DE.tsv",
  GSE109700_deep = "figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
  GSE109700_early = "figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv",
  GSE191055_P27 = "figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
  GSE93535_SIPS = "figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
  GSE93535_rescueSIPS = "figure2_public_aging/GSE93535_SIPS1201_vs_SIPS_DE.tsv",
  GSE93535_rescueQ = "figure2_public_aging/GSE93535_Q1201_vs_Q_DE.tsv",
  GSE179848_late = "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
  GSE113957_age = "figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
  GSE226189_age = "figure2_public_aging/GSE226189_age_per_decade_DE.tsv",
  GSE165177_MPTR = "figure2_mptr/GSE165177_MPTR_paired_DE.tsv")
rows <- list()
for (i in seq_len(nrow(D))) {
  id <- D$id[i]
  s <- if (id == "LOCAL_ReproCM")
         data.frame(gene = anchor$gene, logFC = anchor$Repro_specific_score)
       else if (!is.null(load_sig(id))) load_sig(id)
       else if (id %in% names(legacy) && file.exists(file.path(der, legacy[[id]]))) {
         d <- read.delim(file.path(der, legacy[[id]])); data.frame(gene = toupper(d$gene), logFC = d$logFC)
       } else NULL
  if (is.null(s)) next
  s <- s[is.finite(s$logFC), ]; s <- s[!duplicated(s$gene), ]
  j <- merge(anchor[, c("gene", "Repro_specific_score")], s, by = "gene")
  if (nrow(j) < 500) next
  rows[[id]] <- data.frame(id = id, dataset = D$dataset[i], contrast = D$contrast[i],
    circularity = D$circularity[i], D = D$D[i], rho_UVA = D$rho_UVA[i],
    rho_senescence = D$rho_senescence[i],
    median_abs_logFC = median(abs(s$logFC)),
    frac_absFC_gt_05 = mean(abs(s$logFC) > 0.5),
    rho_anchor = cor(j$Repro_specific_score, j$logFC, method = "spearman"),
    n_shared = nrow(j))
}
R <- do.call(rbind, rows)
cls <- function(id, ds, ct) {
  if (id == "LOCAL_ReproCM") return("Repro-CM")
  if (grepl("senescen|P27|late vs early|SIPS vs", ct, ignore.case = TRUE)) return("senescence")
  if (grepl("donor age|old vs young", ct, ignore.case = TRUE)) return("donor age")
  if (grepl("rescue|Pre-Red|prered|red_vs|Osmoter|UIFSP|Maifuyin|succinate|sauchinone", ct, ignore.case = TRUE))
    return("photoprotection")
  if (grepl("UVA|UVB|UV |uv_vs|injury|solar", ct, ignore.case = TRUE)) return("UV injury")
  if (ds %in% c("GSE149694", "GSE297233", "GSE165177")) return("reprogramming")
  if (ds == "GSE179848") return("metabolic / culture")
  "other"
}
R$class <- mapply(cls, R$id, R$dataset, R$contrast)
R <- R[order(-R$median_abs_logFC), ]
w(R, "magnitude_vs_alignment.tsv")

cat("=========== is alignment explained by perturbation magnitude? ===========\n")
f <- function(sel, lab) {
  x <- R[sel, ]
  if (nrow(x) < 5) return(invisible())
  ct <- suppressWarnings(cor.test(x$median_abs_logFC, x$rho_anchor, method = "spearman"))
  cat(sprintf("%-34s n=%2d  rho(magnitude, anchor alignment) = %+.3f  p = %.3f\n",
              lab, nrow(x), unname(ct$estimate), ct$p.value))
}
f(rep(TRUE, nrow(R)), "all signatures")
f(R$class != "Repro-CM", "excluding Repro-CM itself")
f(R$class == "photoprotection", "photoprotection only")
f(R$class == "metabolic / culture", "metabolic / culture only")
f(R$class %in% c("UV injury", "senescence"), "injury + senescence")

cat("\n=========== the twelve largest perturbations on the map ===========\n")
top <- head(R[R$class != "Repro-CM", ], 12)
cat(sprintf("%-48s %10s %11s %8s  %s\n", "contrast", "med|logFC|", "rho_anchor", "D", "class"))
for (i in seq_len(nrow(top)))
  cat(sprintf("%-48s %10.3f %+11.3f %+8.3f  %s\n",
      substr(paste0(top$dataset[i], ": ", top$contrast[i]), 1, 48),
      top$median_abs_logFC[i], top$rho_anchor[i], top$D[i], top$class[i]))

cat("\n=========== photoprotection contrasts, magnitude vs alignment ===========\n")
pp <- R[R$class == "photoprotection", ]
pp <- pp[order(-pp$median_abs_logFC), ]
cat(sprintf("%-46s %10s %11s\n", "contrast", "med|logFC|", "rho_anchor"))
for (i in seq_len(nrow(pp)))
  cat(sprintf("%-46s %10.3f %+11.3f\n",
      substr(paste0(pp$dataset[i], ": ", pp$contrast[i]), 1, 46),
      pp$median_abs_logFC[i], pp$rho_anchor[i]))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
