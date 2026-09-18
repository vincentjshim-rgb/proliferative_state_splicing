#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Third held-out family, executed exactly as fixed in
#   repro_cm_preregistration_heldout3_photoprotection_ko.md (SHA256 315e68ea...)
# HO6a-d : does the photoprotection alignment replicate in two independent studies?
# HO7    : does the UVA axis fail to load an independent ACCUMULATED UVA injury
#          contrast, as the post-hoc axis interpretation predicts?
# BH within this family of 5 only. Families 1 (6 tests) and 2 (2 tests) stay closed.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(edgeR); library(readxl) })
root <- "public_data_tierA"; h3 <- file.path(root, "heldout3")
out  <- file.path(root, "derived", "heldout3_photoprotection")
dir.create(file.path(out, "signatures"), recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE, row.names = FALSE)

AX <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
anchor <- AX[AX$comparator_consistent & AX$Repro_specific_score != 0 &
             is.finite(AX$Repro_specific_score), ]
AXf <- AX[is.finite(AX$UVA_meta) & is.finite(AX$senescence_meta), ]
cat("anchor genes:", nrow(anchor), " | axis universe:", nrow(AXf), "\n")

matched_perm <- function(x, y, expr, nperm = 10000L, seed = 260911) {
  set.seed(seed)
  rx <- rank(x, ties.method = "average"); ry <- rank(y, ties.method = "average")
  br <- unique(quantile(expr, seq(0, 1, 0.1), na.rm = TRUE))
  idx <- split(seq_along(ry), cut(expr, br, include.lowest = TRUE, labels = FALSE))
  obs <- cor(rx, ry); null <- numeric(nperm)
  for (b in seq_len(nperm)) { p <- ry; for (i in idx) p[i] <- sample(p[i]); null[b] <- cor(rx, p) }
  c(rho = obs, p = (1 + sum(abs(null) >= abs(obs))) / (nperm + 1),
    q025 = unname(quantile(null, .025)), q975 = unname(quantile(null, .975)))
}
de <- function(cnt, grp, ref) {
  y <- DGEList(cnt); g <- factor(grp, levels = c(ref, setdiff(unique(grp), ref)))
  y <- y[filterByExpr(y, group = g), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  d <- model.matrix(~ g); y <- estimateDisp(y, d)
  tt <- topTags(glmQLFTest(glmQLFit(y, d), coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = toupper(rownames(tt)), logFC = tt$logFC)
}
MAP <- list()

## ---------------- GSE240486 : NHDF, Osmoter ----------------
g1 <- read.delim(gzfile(file.path(h3, "GSE240486_all_sample_count.txt.gz")), check.names = FALSE)
m1 <- as.matrix(g1[, -(1:2)]); storage.mode(m1) <- "numeric"
sym1 <- toupper(g1$SYMBOL); k1 <- !is.na(sym1) & nzchar(sym1)
m1 <- rowsum(m1[k1, ], sym1[k1])
grp1 <- sub("[0-9]+$", "", colnames(m1))
cat("\nGSE240486 groups:\n"); print(table(grp1))
kexp1 <- filterByExpr(DGEList(m1), group = grp1)
MAP$GSE240486 <- data.frame(dataset = "GSE240486", n_all = nrow(m1),
  rate_all_rows = 100, n_expressed = sum(kexp1), rate_expressed_rows = 100)
sel <- function(g) colnames(m1)[grp1 %in% g]
inj1 <- de(m1[, sel(c("Control", "UV_Model"))], grp1[grp1 %in% c("Control", "UV_Model")], "Control")
res1 <- list()
for (nm in c("UV_Osmoter_L", "UV_Osmoter_H"))
  res1[[nm]] <- de(m1[, sel(c("UV_Model", nm))], grp1[grp1 %in% c("UV_Model", nm)], "UV_Model")
for (nm in c("noUV_Osmoter_L", "noUV_Osmoter_H"))
  res1[[nm]] <- de(m1[, sel(c("Control", nm))], grp1[grp1 %in% c("Control", nm)], "Control")

## ---------------- GSE222414 : WS1, TAT-UIFSP ----------------
g2 <- read.delim(gzfile(file.path(h3, "GSE222414_counts.txt.gz")), check.names = FALSE)
m2 <- as.matrix(g2[, -1]); storage.mode(m2) <- "numeric"
sym2 <- toupper(g2[[1]]); k2 <- !is.na(sym2) & nzchar(sym2)
m2 <- rowsum(m2[k2, ], sym2[k2])
grp2 <- sub("-[0-9]+$", "", colnames(m2))
lab2 <- c(A = "control", B = "UVB", C = "UVB_UIFSP5", D = "UVB_UIFSP6")[grp2]
cat("\nGSE222414 groups:\n"); print(table(lab2))
MAP$GSE222414 <- data.frame(dataset = "GSE222414", n_all = nrow(m2), rate_all_rows = 100,
  n_expressed = sum(filterByExpr(DGEList(m2), group = lab2)), rate_expressed_rows = 100)
s2 <- function(g) colnames(m2)[lab2 %in% g]
inj2 <- de(m2[, s2(c("control", "UVB"))], lab2[lab2 %in% c("control", "UVB")], "control")
res2 <- list()
for (nm in c("UVB_UIFSP5", "UVB_UIFSP6"))
  res2[[nm]] <- de(m2[, s2(c("UVB", nm))], lab2[lab2 %in% c("UVB", nm)], "UVB")

## ---------------- GSE329475 : HDF, 3-day UVA ----------------
g3 <- as.data.frame(read_excel(file.path(h3, "GSE329475_all_gene_count_matrix.xlsx")))
m3 <- as.matrix(g3[, -1]); storage.mode(m3) <- "numeric"
sym3 <- toupper(g3$gene_name); k3 <- !is.na(sym3) & nzchar(sym3)
m3 <- rowsum(m3[k3, ], sym3[k3])
grp3 <- ifelse(grepl("UVA", colnames(m3)), "UVA", "NC")
cat("\nGSE329475 groups:\n"); print(table(grp3))
MAP$GSE329475 <- data.frame(dataset = "GSE329475", n_all = nrow(m3), rate_all_rows = 100,
  n_expressed = sum(filterByExpr(DGEList(m3), group = grp3)), rate_expressed_rows = 100)
uva3 <- de(m3, grp3, "NC")
w(do.call(rbind, MAP), "HO67_gene_mapping_rates.tsv")

## ---------------- residualise rescue on injury ----------------
resid_on <- function(rescue, injury) {
  d <- merge(rescue, injury, by = "gene", suffixes = c("_r", "_i"))
  d <- d[is.finite(d$logFC_r) & is.finite(d$logFC_i), ]
  data.frame(gene = d$gene, logFC = residuals(lm(logFC_r ~ logFC_i, data = d)))
}
test_anchor <- function(id, ds, lab, eff) {
  j <- merge(anchor, eff, by = "gene"); j <- j[is.finite(j$logFC), ]
  if (nrow(j) < 1000) return(data.frame(id = id, dataset = ds, contrast = lab,
    n_shared = nrow(j), rho = NA, p_matched = NA, null_q025 = NA, null_q975 = NA))
  pm <- matched_perm(j$Repro_specific_score, j$logFC, j$mean_log2_expression)
  data.frame(id = id, dataset = ds, contrast = lab, n_shared = nrow(j),
             rho = unname(pm["rho"]), p_matched = unname(pm["p"]),
             null_q025 = unname(pm["q025"]), null_q975 = unname(pm["q975"]))
}
test_axis <- function(id, ds, lab, eff) {
  j <- merge(AXf[, c("gene","UVA_meta","senescence_meta","mean_log2_expression")], eff, by = "gene")
  j <- j[is.finite(j$logFC), ]
  pm <- matched_perm(j$logFC, j$UVA_meta, j$mean_log2_expression)
  data.frame(id = id, dataset = ds, contrast = lab, n_shared = nrow(j),
             rho_UVA = unname(pm["rho"]), p_matched = unname(pm["p"]),
             rho_senescence = cor(j$logFC, j$senescence_meta, method = "spearman"),
             D = unname(pm["rho"]) - cor(j$logFC, j$senescence_meta, method = "spearman"))
}

## ---------------- PRIMARY ----------------
P <- rbind(
  test_anchor("HO6a", "GSE240486", "UV+0.02% Osmoter vs UV, injury-residualised",
              resid_on(res1$UV_Osmoter_L, inj1)),
  test_anchor("HO6b", "GSE240486", "UV+0.06% Osmoter vs UV, injury-residualised",
              resid_on(res1$UV_Osmoter_H, inj1)),
  test_anchor("HO6c", "GSE222414", "UVB+TAT-UIFSP5 vs UVB, injury-residualised",
              resid_on(res2$UVB_UIFSP5, inj2)),
  test_anchor("HO6d", "GSE222414", "UVB+TAT-UIFSP6 vs UVB, injury-residualised",
              resid_on(res2$UVB_UIFSP6, inj2)))
H7 <- test_axis("HO7", "GSE329475", "UVA 3 consecutive days vs NC", uva3)
P$BH_FDR <- NA
allp <- c(P$p_matched, H7$p_matched)
bh <- p.adjust(allp, "BH"); P$BH_FDR <- bh[1:4]; H7$BH_FDR <- bh[5]
P$verdict <- ifelse(is.na(P$rho), "NOT EVALUABLE",
             ifelse(P$rho > 0 & P$BH_FDR < 0.05, "SUCCESS",
             ifelse(P$rho < 0 & P$BH_FDR < 0.05, "OPPOSITE", "NOT SUPPORTED")))
H7$verdict <- if (abs(H7$rho_UVA) < 0.15) "PREDICTION HELD" else
              if (H7$rho_UVA > 0.40) "REFUTES THE AXIS INTERPRETATION" else "INTERMEDIATE"
w(P, "HO6_primary_tests.tsv"); w(H7, "HO7_axis_test.tsv")

cat("\n================ HO6 : photoprotection replication ================\n")
print(P[, c("id","dataset","n_shared","rho","p_matched","BH_FDR","verdict")], row.names = FALSE, digits = 3)
n_ok <- sum(P$verdict == "SUCCESS")
both <- length(unique(P$dataset[P$verdict == "SUCCESS"])) == 2
cat(sprintf("\nsuccesses: %d of 4 | both studies represented: %s -> %s\n", n_ok, both,
            if (n_ok >= 2 && both) "SUCCESS" else if (sum(P$verdict == "OPPOSITE") >= 3) "REFUTED" else "NOT SUPPORTED"))
cat("\n================ HO7 : does the UVA axis load accumulated UVA injury? ================\n")
print(H7[, c("id","n_shared","rho_UVA","rho_senescence","D","p_matched","BH_FDR","verdict")],
      row.names = FALSE, digits = 3)
cat("\nfor reference - contrasts that DEFINE the UVA axis: +0.798 and +0.816\n")

## ---------------- SENSITIVITY ----------------
S <- rbind(
  test_anchor("S1", "GSE240486", "0.02% Osmoter vs Control (agent alone, no UV)", res1$noUV_Osmoter_L),
  test_anchor("S2", "GSE240486", "0.06% Osmoter vs Control (agent alone, no UV)", res1$noUV_Osmoter_H),
  test_anchor("S3", "GSE240486", "UV+0.02% vs UV, NOT residualised", res1$UV_Osmoter_L),
  test_anchor("S4", "GSE240486", "UV+0.06% vs UV, NOT residualised", res1$UV_Osmoter_H),
  test_anchor("S5", "GSE222414", "UVB+UIFSP5 vs UVB, NOT residualised", res2$UVB_UIFSP5),
  test_anchor("S6", "GSE222414", "UVB+UIFSP6 vs UVB, NOT residualised", res2$UVB_UIFSP6),
  test_anchor("S7", "GSE240486", "UV vs Control (injury)", inj1),
  test_anchor("S8", "GSE222414", "UVB vs control (injury)", inj2))
w(S, "HO67_sensitivity_tests.tsv")
cat("\n================ SENSITIVITY ================\n")
print(S[, c("id","dataset","contrast","rho","p_matched")], row.names = FALSE, digits = 3)

AXT <- rbind(
  test_axis("ax1", "GSE240486", "UV vs Control (injury)", inj1),
  test_axis("ax2", "GSE222414", "UVB vs control (injury)", inj2),
  test_axis("ax3", "GSE240486", "UV+0.06% Osmoter vs UV", res1$UV_Osmoter_H),
  test_axis("ax4", "GSE222414", "UVB+UIFSP6 vs UVB", res2$UVB_UIFSP6),
  H7[, names(test_axis("x","y","z",uva3))])
w(AXT, "HO67_axis_loadings.tsv")
cat("\n================ axis loadings of the new contrasts ================\n")
print(AXT[, c("dataset","contrast","rho_UVA","rho_senescence","D")], row.names = FALSE, digits = 3)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
