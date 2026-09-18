#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# How far do the two reference axes generalise?
# For every compendium signature that did NOT contribute to an axis, ask how
# strongly it loads on that axis.  A reference axis that only reproduces inside
# the studies that built it is not a reference axis.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
root <- "public_data_tierA"
out  <- file.path(root, "derived", "axis_characterisation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)
R <- read.delim(file.path(root, "derived", "compendium", "compendium_D_scores.tsv"))

## ---- hand classification of every UV-family and senescence-family contrast --
uv <- data.frame(
  id = c("GSE240226_UVA","GSE302943_UVA","GSE89005_single6h","GSE89005_single24h",
         "GSE89005_repeat24h","GSE125429_UVA","GSE116968_uv_vs_control_1h",
         "GSE116968_uv_vs_control_4h","GSE116968_red_vs_control_1h","GSE116968_red_vs_control_4h",
         "GSE116968_prered_vs_uv_1h","GSE116968_prered_vs_uv_4h",
         "GSE116968_prered_vs_uv_injuryresid_1h","GSE116968_prered_vs_uv_injuryresid_4h",
         "GSE240226_Maifuyin","GSE240226_succinate"),
  kind = c("injury, acute (axis)","injury, cumulative (axis)","injury, acute 6 h","injury, 24 h",
           "injury, repeated","injury, chronic","injury, UVB 1 h","injury, UVB 4 h",
           "photoprotection, independent","photoprotection, independent",
           "photoprotection, independent","photoprotection, independent",
           "photoprotection, independent","photoprotection, independent",
           "rescue within the axis study","rescue within the axis study"),
  in_axis = c(TRUE, TRUE, rep(FALSE, 14)))
sen <- data.frame(
  id = c("GSE179848_late","GSE109700_deep","GSE191055_P27","GSE93535_SIPS",
         "GSE109700_early","GSE306957_pLenti","GSE306957_siScramble",
         "GSE306957_IFIH1_KO","GSE306957_DDX58_KO"),
  kind = c(rep("senescence (axis)", 4), "senescence, same study",
           rep("senescence, independent", 4)),
  in_axis = c(rep(TRUE, 4), rep(FALSE, 5)))

grab <- function(tab, col) {
  m <- merge(tab, R[, c("id", "dataset", "contrast", "rho_UVA", "rho_senescence", "D", "n_genes")], by = "id")
  m$rho_axis <- m[[col]]
  m[order(-m$rho_axis), ]
}
U <- grab(uv, "rho_UVA"); S <- grab(sen, "rho_senescence")
w(U, "UVA_axis_generalisation.tsv"); w(S, "senescence_axis_generalisation.tsv")

say <- function(x, label) {
  cat("\n########", label, "########\n")
  cat(sprintf("%-42s %-26s %9s\n", "signature", "kind", "rho_axis"))
  for (i in seq_len(nrow(x)))
    cat(sprintf("%-42s %-26s %+9.3f%s\n", substr(paste0(x$dataset[i], ": ", x$contrast[i]), 1, 42),
                x$kind[i], x$rho_axis[i], ifelse(x$in_axis[i], "  [axis]", "")))
}
say(U, "UVA axis"); say(S, "senescence axis")

cat("\n######## generalisation summary ########\n")
grp <- function(x, sel, lab) {
  v <- x$rho_axis[sel]
  data.frame(axis = lab, group = NA, n = length(v), median = round(median(v), 3),
             min = round(min(v), 3), max = round(max(v), 3))
}
sm <- rbind(
  transform(grp(U, U$in_axis, "UVA"), group = "defines the axis"),
  transform(grp(U, !U$in_axis & grepl("^injury", U$kind), "UVA"), group = "independent UV INJURY"),
  transform(grp(U, U$kind == "photoprotection, independent", "UVA"), group = "independent PHOTOPROTECTION"),
  transform(grp(U, U$kind == "rescue within the axis study", "UVA"), group = "rescue inside the axis study"),
  transform(grp(S, S$in_axis, "senescence"), group = "defines the axis"),
  transform(grp(S, S$kind == "senescence, independent", "senescence"), group = "independent SENESCENCE"))
print(sm[, c("axis", "group", "n", "median", "min", "max")], row.names = FALSE)
w(sm, "axis_generalisation_summary.tsv")

inj <- U$rho_axis[!U$in_axis & grepl("^injury", U$kind)]
pho <- U$rho_axis[U$kind == "photoprotection, independent"]
cat(sprintf("\nUVA axis, independent contrasts only: injury median %+.3f (n=%d) vs photoprotection median %+.3f (n=%d)\n",
            median(inj), length(inj), median(pho), length(pho)))
cat(sprintf("Wilcoxon rank-sum p = %.4f  (n = %d vs %d)\n",
            suppressWarnings(wilcox.test(pho, inj)$p.value), length(pho), length(inj)))
cat(sprintf("ranges overlap only at the earliest injury timepoint: injury max %+.3f (UVA 6 h), photoprotection min %+.3f\n",
            max(inj), min(pho)))
cat("\nNOTE: the two GSE240226 rescue contrasts are negative BY CONSTRUCTION -- that study defines\n")
cat("      the UVA axis, so undoing its own UVA response must anti-correlate with it.\n")
cat(sprintf("\nlocal Repro-CM anchor sits at rho_UVA = %+.3f -- with the photoprotection group, not the injury group.\n",
            R$rho_UVA[R$id == "LOCAL_ReproCM"]))
si <- S$rho_axis[S$kind == "senescence, independent"]
cat(sprintf("\nsenescence axis, 4 fully independent contrasts: %s (median %+.3f)\n",
            paste(sprintf("%+.3f", sort(si)), collapse = " "), median(si)))
