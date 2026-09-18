#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# What makes a perturbation align with the Repro-CM anchor?
# Magnitude was ruled out. Among agent-alone arms the top alignments are
# bioenergetic stressors (oligomycin, 2-deoxyglucose, SURF1 deficiency) and the
# bottom is quiescence, which points at the adaptive stress-response module that
# Figure 4B identified as one of the two decoupling carriers.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"; der <- file.path(root, "derived")
out  <- file.path(der, "what_drives_alignment")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)

AX <- read.delim(file.path(der, "decoupling_validation/primary_gene_axis_matrix.tsv"), check.names = FALSE)
anchor <- AX[AX$comparator_consistent & AX$Repro_specific_score != 0 & is.finite(AX$Repro_specific_score), ]
R <- read.delim(file.path(der, "magnitude_vs_alignment", "magnitude_vs_alignment.tsv"))

con <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
L <- readLines(con, warn = FALSE); close(con); sp <- strsplit(L, "\t", fixed = TRUE)
sets <- lapply(sp, function(x) unique(toupper(x[-c(1, 2)])))
names(sets) <- vapply(sp, function(x) x[1], character(1))
pk <- function(nm) unique(unlist(sets[intersect(nm, names(sets))]))
MOD <- list(
  UPR_PERK = pk(c("PERK regulates gene expression", "Unfolded Protein Response (UPR)",
                  "ATF4 activates genes in response to endoplasmic reticulum stress")),
  HSF1     = pk(c("Regulation of HSF1-mediated heat shock response", "HSF1 activation",
                  "Cellular response to heat stress")),
  splicing = pk(c("mRNA Splicing", "mRNA Splicing - Major Pathway", "snRNP Assembly")),
  premRNA  = pk("Processing of Capped Intron-Containing Pre-mRNA"),
  prolif   = pk(c("Cell Cycle, Mitotic", "DNA Replication")),
  ECM      = pk("Extracellular matrix organization"))

sigdir <- file.path(der, "compendium", "signatures")
legacy <- c(GSE240226_UVA="figure1_public_uva/GSE240226_UVA_vs_control_DE.tsv",
  GSE240226_Maifuyin="figure1_public_uva/GSE240226_Maifuyin_rescue_vs_UVA_DE.tsv",
  GSE240226_succinate="figure1_public_uva/GSE240226_succinate_rescue_vs_UVA_DE.tsv",
  GSE302943_UVA="figure1_public_uva/GSE302943_UVA_vs_control_DE.tsv",
  GSE125429_UVA="figure1_public_uva/GSE125429_UVA_vs_control_DE.tsv",
  GSE89005_single6h="figure1_public_uva/GSE89005_single_6h_UVA_vs_sham_DE.tsv",
  GSE89005_single24h="figure1_public_uva/GSE89005_single_24h_UVA_vs_sham_DE.tsv",
  GSE89005_repeat24h="figure1_public_uva/GSE89005_repeated_24h_UVA_vs_sham_DE.tsv",
  GSE109700_deep="figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
  GSE109700_early="figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv",
  GSE191055_P27="figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
  GSE93535_SIPS="figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
  GSE93535_rescueSIPS="figure2_public_aging/GSE93535_SIPS1201_vs_SIPS_DE.tsv",
  GSE93535_rescueQ="figure2_public_aging/GSE93535_Q1201_vs_Q_DE.tsv",
  GSE179848_late="figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
  GSE113957_age="figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
  GSE226189_age="figure2_public_aging/GSE226189_age_per_decade_DE.tsv",
  GSE165177_MPTR="figure2_mptr/GSE165177_MPTR_paired_DE.tsv")
get_sig <- function(id) {
  f <- file.path(sigdir, paste0(id, ".tsv"))
  if (file.exists(f)) { d <- read.delim(f); return(data.frame(gene = d$gene, logFC = d$logFC)) }
  if (id %in% names(legacy) && file.exists(file.path(der, legacy[[id]]))) {
    d <- read.delim(file.path(der, legacy[[id]])); return(data.frame(gene = toupper(d$gene), logFC = d$logFC)) }
  NULL
}
## module score of a signature = expression-matched z of its module mean logFC
mz <- function(eff, members, nperm = 2000, seed = 5) {
  e <- eff[is.finite(eff$logFC), ]
  ins <- e$gene %in% members
  if (sum(ins) < 10) return(NA_real_)
  obs <- mean(e$logFC[ins]); set.seed(seed)
  null <- replicate(nperm, mean(sample(e$logFC, sum(ins))))
  (obs - mean(null)) / sd(null)
}
rows <- list()
for (i in seq_len(nrow(R))) {
  id <- R$id[i]; if (id == "LOCAL_ReproCM") next
  s <- get_sig(id); if (is.null(s)) next
  s <- s[is.finite(s$logFC), ]; s <- s[!duplicated(s$gene), ]
  z <- vapply(MOD, function(m) mz(s, m), numeric(1))
  rows[[id]] <- cbind(R[i, c("id","dataset","contrast","class","rho_anchor","D","median_abs_logFC")],
                      as.data.frame(as.list(z)))
}
M <- do.call(rbind, rows)
w(M, "module_scores_vs_alignment.tsv")

cat("========== what predicts alignment with the anchor? (n =", nrow(M), "signatures) ==========\n")
cat(sprintf("%-22s %10s %9s\n", "predictor", "rho", "p"))
for (v in c(names(MOD), "median_abs_logFC")) {
  x <- M[[v]]; k <- is.finite(x)
  ct <- suppressWarnings(cor.test(x[k], M$rho_anchor[k], method = "spearman"))
  cat(sprintf("%-22s %+10.3f %9.4f\n", v, unname(ct$estimate), ct$p.value))
}
cat("\n========== partial: adaptive stress modules, controlling for magnitude ==========\n")
for (v in c("UPR_PERK", "HSF1")) {
  k <- is.finite(M[[v]])
  r1 <- residuals(lm(M$rho_anchor[k] ~ M$median_abs_logFC[k]))
  r2 <- residuals(lm(M[[v]][k] ~ M$median_abs_logFC[k]))
  cat(sprintf("%-22s partial rho = %+.3f  (n = %d)\n", v, cor(r1, r2, method = "spearman"), sum(k)))
}
cat("\n========== agent-alone arms, ranked by UPR/PERK score ==========\n")
alone <- M[grepl("alone|red_vs_control", M$contrast) |
           (M$dataset == "GSE179848" & grepl("vs Control|oxygen|SURF1", M$contrast)), ]
alone <- alone[order(-alone$UPR_PERK), ]
cat(sprintf("%-44s %9s %9s %10s\n", "contrast", "UPR_PERK", "HSF1", "rho_anchor"))
for (i in seq_len(nrow(alone)))
  cat(sprintf("%-44s %+9.2f %+9.2f %+10.3f\n",
      substr(paste0(alone$dataset[i], ": ", alone$contrast[i]), 1, 44),
      alone$UPR_PERK[i], alone$HSF1[i], alone$rho_anchor[i]))
ct <- suppressWarnings(cor.test(alone$UPR_PERK, alone$rho_anchor, method = "spearman"))
cat(sprintf("\nwithin agent-alone arms: rho(UPR/PERK, alignment) = %+.3f, p = %.4f, n = %d\n",
            unname(ct$estimate), ct$p.value, nrow(alone)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
