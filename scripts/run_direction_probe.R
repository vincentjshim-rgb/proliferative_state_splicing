#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# run_direction_probe.R
#
# 후속 연구 방향(splicing / proliferation confounding) 결정을 위한 탐색 분석.
# 기존 결과 파일을 수정하지 않고 derived/direction_probe/ 아래에만 기록한다.
#
# P1  553개 Reactome pathway에서 local NES와 decoupling 축(UVA-senescence)의 관계
# P2  module x axis z-map (local / UVA / senescence / MPTR / fibroblast age / lifespan)
# P3  GSE179848(345 samples) 기반 per-gene proliferation loading
# P4  proliferation 보정 후 decoupling index D
# P5  GSE113957(104 donors) module-age 기울기의 proliferation 보정
# P6  로컬 STAR SJ.out.tab 기반 junction PSI 검출력(binomial null 포함)
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db); library(edgeR); library(readxl)
})

root <- "public_data_tierA"
out  <- file.path(root, "derived", "direction_probe")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")

## ---- Reactome sets ---------------------------------------------------------
con <- unz(file.path(root, "network", "ReactomePathways.gmt.zip"), "ReactomePathways.gmt")
gmt <- readLines(con, warn = FALSE); close(con)
sp   <- strsplit(gmt, "\t", fixed = TRUE)
sets <- lapply(sp, function(x) unique(toupper(x[-c(1, 2)])))
names(sets) <- vapply(sp, function(x) x[1], character(1))
pick <- function(nm) unique(unlist(sets[nm]))

SF_REGULON <- toupper(c(
  paste0("SRSF", 1:12), "TRA2A", "TRA2B", "SRRM1", "SRRM2",
  "HNRNPA0","HNRNPA1","HNRNPA2B1","HNRNPA3","HNRNPAB","HNRNPC","HNRNPD","HNRNPDL","HNRNPF",
  "HNRNPH1","HNRNPH2","HNRNPH3","HNRNPK","HNRNPL","HNRNPM","HNRNPR","HNRNPU","HNRNPUL1","HNRNPUL2",
  "PCBP1","PCBP2","SYNCRIP","RALY","DAZAP1","SNRPA","SNRPA1","SNRPB","SNRPB2","SNRPC",
  "SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG","SNRNP70","SNRNP200","SNRNP27","SNRNP40",
  "SF1","SF3A1","SF3A2","SF3A3","SF3B1","SF3B2","SF3B3","SF3B4","SF3B5","SF3B6",
  "U2AF1","U2AF2","U2SURP","PRPF3","PRPF4","PRPF6","PRPF8","PRPF18","PRPF19","PRPF31",
  "PRPF38A","PRPF38B","PRPF39","PRPF40A","EFTUD2","CDC5L","PLRG1","BUD31","AQR","XAB2","ISY1",
  "SART1","SART3","USP39","RBM22","CLK1","CLK2","CLK3","CLK4","SRPK1","SRPK2",
  "ESRP1","ESRP2","MBNL1","MBNL2","MBNL3","QKI","PTBP1","PTBP2","PTBP3","CELF1","CELF2","RBFOX2",
  "KHDRBS1","RBM4","RBM24","RBM38","RBM47","RBM8A","MAGOH","EIF4A3","CASC3","ALYREF","NXF1",
  "DDX39A","DDX39B","THOC1","THOC2","THOC3","THOC5"))

MODULES <- list(
  SF_regulon  = SF_REGULON,
  splicing    = pick(c("mRNA Splicing", "mRNA Splicing - Major Pathway",
                       "mRNA Splicing - Minor Pathway", "snRNP Assembly")),
  premRNA     = pick("Processing of Capped Intron-Containing Pre-mRNA"),
  mRNA_export = pick(c("Transport of Mature mRNA derived from an Intron-Containing Transcript",
                       "Transport of Mature Transcript to Cytoplasm")),
  three_prime = pick(c("mRNA 3'-end processing", "mRNA Polyadenylation")),
  NMD         = pick("Nonsense-Mediated Decay (NMD)"),
  translation = pick(c("Eukaryotic Translation Elongation", "Eukaryotic Translation Initiation")),
  rRNA        = pick(c("rRNA processing", "rRNA modification in the nucleus and cytosol")),
  HSF1_heat_shock = pick(c("Regulation of HSF1-mediated heat shock response", "HSF1 activation",
                           "Cellular response to heat stress", "HSF1-dependent transactivation")),
  UPR_PERK    = pick(c("PERK regulates gene expression", "Unfolded Protein Response (UPR)",
                       "ATF4 activates genes in response to endoplasmic reticulum stress")),
  ROS_detox   = pick(c("Detoxification of Reactive Oxygen Species", "KEAP1-NFE2L2 pathway")),
  autophagy   = pick(c("Autophagy", "Macroautophagy", "Chaperone Mediated Autophagy")),
  proliferation = pick(c("Cell Cycle, Mitotic", "DNA Replication")),
  ECM         = pick("Extracellular matrix organization"),
  IFN         = pick("Interferon Signaling"))

## ===========================================================================
## P1  pathway-level decoupling
## ===========================================================================
pm <- read.delim(file.path(root, "derived/decoupling_validation/pathway_axis_matrix.tsv"),
                 check.names = FALSE)
pm$decoupling_axis_NES <- pm$UVA_meta_NES - pm$senescence_meta_NES
pm$joint_score <- pm$local_NES * pm$decoupling_axis_NES
ct <- suppressWarnings(cor.test(pm$local_NES, pm$decoupling_axis_NES, method = "spearman"))
p1 <- data.frame(
  test = c("local_NES vs (UVA-senescence) decoupling axis",
           "local_NES vs UVA_meta_NES", "local_NES vs senescence_meta_NES"),
  n_pathways = nrow(pm),
  spearman_rho = c(unname(ct$estimate),
                   cor(pm$local_NES, pm$UVA_meta_NES, method = "spearman"),
                   cor(pm$local_NES, pm$senescence_meta_NES, method = "spearman")),
  p_value = c(ct$p.value,
              suppressWarnings(cor.test(pm$local_NES, pm$UVA_meta_NES, method="spearman"))$p.value,
              suppressWarnings(cor.test(pm$local_NES, pm$senescence_meta_NES, method="spearman"))$p.value))
w(p1, "P1_pathway_decoupling_summary.tsv")
w(pm[order(-pm$joint_score), ], "P1_pathway_decoupling_ranked.tsv")
cat("[P1] rho(local_NES, UVA-senescence) =", round(p1$spearman_rho[1], 4),
    " p =", signif(p1$p_value[1], 3), "over", nrow(pm), "pathways\n")

## ===========================================================================
## P2  module x axis z-map
## ===========================================================================
gx <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
# Two of the Tyshkovskiy 2026 sheets carry no expression column, so no
# expression-matched null is definable for them. Those axes fall back to an
# unmatched permutation and are labelled as such; never mix the two silently.
matched_z <- function(val, gene, members, expr, nperm = 5000, seed = 7) {
  has_expr <- !is.null(expr) && any(is.finite(expr))
  keep <- if (has_expr) is.finite(val) & is.finite(expr) else is.finite(val)
  val <- val[keep]; gene <- gene[keep]; if (has_expr) expr <- expr[keep]
  inset <- gene %in% members
  if (sum(inset) < 10) return(c(n = sum(inset), mean = NA, z = NA, p = NA, matched = has_expr))
  obs <- mean(val[inset]); set.seed(seed); n <- sum(inset)
  null <- if (has_expr) {
    dec <- as.integer(cut(rank(expr, ties.method = "first"), breaks = 10, labels = FALSE))
    tab <- table(dec[inset]); idx <- split(seq_along(val), dec)
    replicate(nperm, mean(val[unlist(lapply(names(tab), function(k) sample(idx[[k]], tab[[k]])))]))
  } else replicate(nperm, mean(sample(val, n)))
  c(n = n, mean = obs, z = (obs - mean(null)) / sd(null),
    p = (1 + sum(abs(null - mean(null)) >= abs(obs - mean(null)))) / (nperm + 1),
    matched = has_expr)
}
axes <- list(
  local       = list(gx$gene, gx$Repro_specific_score, gx$mean_log2_expression),
  UVA_meta    = list(gx$gene, gx$UVA_meta,             gx$mean_log2_expression),
  senesc_meta = list(gx$gene, gx$senescence_meta,      gx$mean_log2_expression))
mp <- read.delim(file.path(root, "derived/figure2_mptr/GSE165177_MPTR_paired_DE.tsv"))
mp <- mp[!duplicated(mp$gene), ]
axes$MPTR_reprogramming <- list(toupper(mp$gene), mp$logFC, mp$AveExpr)
ag <- read.delim(file.path(root, "derived/figure2_public_aging/GSE113957_age_per_decade_DE.tsv"))
ag <- ag[!duplicated(ag$gene), ]
axes$fibroblast_age <- list(toupper(ag$gene), ag$logFC, ag$AveExpr)
f26 <- file.path(root, "longevity", "Tyshkovskiy_2026_Table_S2_signatures.xlsx")
for (sh in c("(M) Rodents Max lifespan", "(K) Rodents Mortality rate", "(F) ITP Max lifespan")) {
  tx <- as.data.frame(read_excel(f26, sheet = sh))
  tx <- tx[!duplicated(toupper(tx$Gene.symbol)), ]
  nm <- if (grepl("^\\(F\\)", sh)) "Tyshkovskiy_ITP_max_lifespan" else
        if (grepl("Max", sh)) "Tyshkovskiy_rodent_max_lifespan" else "Tyshkovskiy_rodent_mortality"
  axes[[nm]] <- list(toupper(tx$Gene.symbol), tx$Slope,
                     if ("logCPM" %in% names(tx)) tx$logCPM else NULL)
}
p2 <- do.call(rbind, lapply(names(MODULES), function(mn) do.call(rbind, lapply(names(axes), function(an) {
  a <- axes[[an]]; r <- matched_z(a[[2]], a[[1]], MODULES[[mn]], a[[3]])
  data.frame(module = mn, axis = an, n_genes = unname(r["n"]), module_mean = unname(r["mean"]),
             matched_z = unname(r["z"]), empirical_p = unname(r["p"]),
             null_type = if (isTRUE(as.logical(r["matched"]))) "expression-matched" else "unmatched")
}))))
w(p2, "P2_module_axis_zmap.tsv")
cat("[P2] module x axis map written (", nrow(p2), "rows )\n")

## ===========================================================================
## P3  per-gene proliferation loading from GSE179848
## ===========================================================================
raw <- read.csv(gzfile(file.path(root, "aging",
          "GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz")), check.names = FALSE)
ids <- sub("\\..*$", "", raw[[1]])
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(ids), keytype = "REFSEQ",
                               column = "SYMBOL", multiVals = "first"))[ids]
mat <- as.matrix(raw[, -1]); storage.mode(mat) <- "numeric"
ok  <- !is.na(sym) & nzchar(sym)
mat <- rowsum(mat[ok, ], sym[ok])
dge <- DGEList(mat); dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
Z   <- t(scale(t(cpm(dge, log = TRUE, prior.count = 1))))
prol_idx <- rownames(Z) %in% MODULES$proliferation
P <- colMeans(Z[prol_idx, , drop = FALSE])
loading <- apply(Z, 1, function(v) cor(v, P, method = "spearman"))
w(data.frame(gene = names(loading), proliferation_loading = loading,
             source = "GSE179848_345_fibroblast_samples"),
  "P3_proliferation_loading_GSE179848.tsv")
cat("[P3] proliferation loading for", length(loading), "genes;",
    sum(prol_idx), "genes define the meta-gene\n")

## ===========================================================================
## P4  decoupling D with proliferation removed
## ===========================================================================
a <- gx[gx$comparator_consistent & is.finite(gx$UVA_meta) & is.finite(gx$senescence_meta), ]
a$proliferation_loading <- loading[match(a$gene, names(loading))]
a <- a[is.finite(a$proliferation_loading), ]
Dfun <- function(l, u, s) { ru <- cor(l, u, method = "spearman")
                            rs <- cor(l, s, method = "spearman")
                            c(ru, rs, ru - rs) }
res <- function(y, x) residuals(lm(y ~ x))
neu <- a[abs(a$proliferation_loading) <  0.2, ]
hi  <- a[abs(a$proliferation_loading) >= 0.4, ]
p4 <- rbind(
  data.frame(analysis = "no adjustment (shared genes)", n_genes = nrow(a),
             t(Dfun(a$Repro_specific_score, a$UVA_meta, a$senescence_meta))),
  data.frame(analysis = "all three axes residualized on proliferation loading", n_genes = nrow(a),
             t(Dfun(res(a$Repro_specific_score, a$proliferation_loading),
                    res(a$UVA_meta, a$proliferation_loading),
                    res(a$senescence_meta, a$proliferation_loading)))),
  data.frame(analysis = "proliferation-neutral genes (|loading| < 0.2)", n_genes = nrow(neu),
             t(Dfun(neu$Repro_specific_score, neu$UVA_meta, neu$senescence_meta))),
  data.frame(analysis = "proliferation-linked genes (|loading| >= 0.4)", n_genes = nrow(hi),
             t(Dfun(hi$Repro_specific_score, hi$UVA_meta, hi$senescence_meta))))
names(p4)[3:5] <- c("rho_UVA", "rho_senescence", "decoupling_index_D")
p4$cor_axis_with_proliferation_loading <- c(
  cor(a$Repro_specific_score, a$proliferation_loading, method = "spearman"), NA, NA, NA)
w(p4, "P4_decoupling_proliferation_sensitivity.tsv")
cat("[P4] D unadjusted =", round(p4$decoupling_index_D[1], 3),
    "| residualized =", round(p4$decoupling_index_D[2], 3),
    "| proliferation-neutral genes =", round(p4$decoupling_index_D[3], 3), "\n")

## ===========================================================================
## P5  GSE113957 module-age slope with and without proliferation adjustment
## ===========================================================================
fp  <- read.delim(gzfile(file.path(root, "aging", "GSE113957_fpkm.txt.gz")), check.names = FALSE)
fsym <- toupper(sub("\\|.*", "", fp[["Annotation/Divergence"]]))
fm  <- as.matrix(fp[, 9:ncol(fp)]); storage.mode(fm) <- "numeric"
kp  <- !is.na(fsym) & nzchar(fsym) & fsym != "NA"
fm  <- rowsum(fm[kp, ], fsym[kp])
cl  <- tolower(colnames(fm))
age <- as.numeric(sub("^[^_]+_([0-9]+)(yr|ys)?([0-9]*mos)?_.*$", "\\1", cl))
sex <- sub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", cl)
sex[sex == "f"] <- "female"; sex[sex == "m"] <- "male"
sex[!sex %in% c("male", "female")] <- NA
use <- !grepl("hgps|progeria", cl) & !is.na(age) & !is.na(sex) & age >= 22 & age <= 96
E   <- log2(fm[, use] + 1); E <- E[rowMeans(E > 0.5) > 0.5, ]
ZE  <- t(scale(t(E))); A <- age[use]; S <- factor(sex[use])
msc <- sapply(MODULES, function(g) colMeans(ZE[rownames(ZE) %in% g, , drop = FALSE]))
p5 <- do.call(rbind, lapply(colnames(msc), function(mn) {
  y  <- msc[, mn]
  b1 <- summary(lm(y ~ I(A / 10) + S))$coef["I(A/10)", ]
  if (mn == "proliferation") { b2 <- c(NA, NA, NA, NA) } else {
    b2 <- summary(lm(y ~ I(A / 10) + msc[, "proliferation"] + S))$coef["I(A/10)", ] }
  data.frame(module = mn, n_donors = length(y),
             r_with_proliferation = cor(y, msc[, "proliferation"]),
             beta_per_decade_raw = b1[1], p_raw = b1[4],
             beta_per_decade_proliferation_adjusted = b2[1], p_adjusted = b2[4])
}))
w(p5, "P5_GSE113957_module_age_proliferation_adjusted.tsv")
cat("[P5] splicing module: r(prolif) =",
    round(p5$r_with_proliferation[p5$module == "splicing"], 3),
    "; age slope p_raw =", signif(p5$p_raw[p5$module == "splicing"], 2),
    "-> p_adjusted =", signif(p5$p_adjusted[p5$module == "splicing"], 2), "\n")

## ===========================================================================
## P6  local junction PSI detectability vs binomial null
## ===========================================================================
sjd <- file.path(root, "derived/local_repro_anchor/star_genecounts")
rdsj <- function(s) {
  x <- read.delim(file.path(sjd, s, "SJ.out.tab"), header = FALSE)
  names(x) <- c("chr","start","end","strand","motif","annot","uniq","multi","overhang")
  x <- x[x$annot == 1 & x$strand %in% c(1, 2), ]
  x$id <- paste(x$chr, x$start, x$end, x$strand, sep = ":")
  x[, c("id", "uniq")]
}
SJ <- lapply(c(HDF = "HDF", REP = "REP", IPS = "IPS"), rdsj)
jid <- Reduce(union, lapply(SJ, function(x) x$id))
J <- data.frame(id = jid)
for (n in names(SJ)) J[[n]] <- SJ[[n]]$uniq[match(jid, SJ[[n]]$id)]
J[is.na(J)] <- 0
pp <- do.call(rbind, strsplit(J$id, ":", fixed = TRUE))
J$chr <- pp[, 1]; J$start <- as.integer(pp[, 2]); J$end <- as.integer(pp[, 3]); J$strand <- pp[, 4]

psi_block <- function(keycols, mincov, nsim = 30) {
  key <- do.call(paste, c(J[keycols], sep = "_"))
  tH <- ave(J$HDF, key, FUN = sum); tR <- ave(J$REP, key, FUN = sum); tI <- ave(J$IPS, key, FUN = sum)
  nalt <- ave(rep(1, nrow(J)), key, FUN = sum)
  k <- nalt >= 2 & tH >= mincov & tR >= mincov & tI >= mincov
  o <- data.frame(id = J$id[k], kH = J$HDF[k], kR = J$REP[k], kI = J$IPS[k],
                  nH = tH[k], nR = tR[k], nI = tI[k])
  o$pH <- o$kH/o$nH; o$pR <- o$kR/o$nR; o$pI <- o$kI/o$nI
  o <- o[!(o$pH > 0.98 & o$pR > 0.98 & o$pI > 0.98), ]
  o$dH <- o$pR - o$pH; o$dI <- o$pR - o$pI
  ptrue <- (o$kH + o$kR + o$kI) / (o$nH + o$nR + o$nI)
  sim <- function() { sH <- rbinom(nrow(o), o$nH, ptrue)/o$nH
                      sR <- rbinom(nrow(o), o$nR, ptrue)/o$nR
                      sI <- rbinom(nrow(o), o$nI, ptrue)/o$nI
                      list(a = sR - sH, b = sR - sI) }
  set.seed(11)
  nr <- replicate(nsim, { s <- sim(); cor(s$a, s$b, method = "spearman") })
  nc <- replicate(nsim, { s <- sim(); sum(sign(s$a) == sign(s$b) &
                                          abs(s$a) >= 0.10 & abs(s$b) >= 0.10) })
  obs_n <- sum(sign(o$dH) == sign(o$dI) & abs(o$dH) >= 0.10 & abs(o$dI) >= 0.10)
  list(row = data.frame(
        group_type = paste(keycols, collapse = "+"), min_coverage = mincov,
        n_testable = nrow(o),
        rho_observed = cor(o$dH, o$dI, method = "spearman"),
        rho_binomial_null = mean(nr), rho_null_sd = sd(nr),
        n_candidates_dPSI_0.10 = obs_n,
        expected_by_noise = mean(nc),
        empirical_FDR = ifelse(obs_n > 0, mean(nc)/obs_n, NA)), events = o)
}
p6 <- do.call(rbind, lapply(c(20, 50, 100, 200), function(mc)
  rbind(psi_block(c("chr","start","strand"), mc)$row,
        psi_block(c("chr","end","strand"), mc)$row)))
w(p6, "P6_local_junction_psi_detectability.tsv")
ev <- rbind(psi_block(c("chr","start","strand"), 50)$events,
            psi_block(c("chr","end","strand"), 50)$events)
ev$conservative_dPSI <- ifelse(sign(ev$dH) == sign(ev$dI),
                               sign(ev$dH) * pmin(abs(ev$dH), abs(ev$dI)), 0)
ev <- ev[!duplicated(ev$id), ]
w(ev[order(-abs(ev$conservative_dPSI)), ][1:500, ], "P6_local_junction_top500_candidates.tsv")
cat("[P6] junction-level candidates at |dPSI|>=0.10 carry empirical FDR",
    paste(round(range(p6$empirical_FDR, na.rm = TRUE), 2), collapse = "-"), "\n")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("\nAll probe outputs written to", out, "\n")
