## Fig. 7 (ageing-programme audit) revisited for the npj Aging revision, 2026-09-15.
##
## Two changes and one reviewer request:
##  1. the donor-age half is recomputed in the redefined primary cohort
##     (GSE113957 normal donors aged 20+, HGPS excluded, repository / instrument /
##     sex as covariates) from recount3 gene sums, so that every figure in the
##     manuscript rests on one matrix and one cohort;
##  2. the division-rate half is recomputed with all Reactome cell-cycle genes
##     removed from every programme, to test whether the coupling of the eight
##     "division readout" programmes is carried by cell-cycle genes inside them
##     (reviewer 1, major 7a);
##  3. a second senescence panel (SenMayo) was looked for in the project and is
##     not present; nothing was downloaded.
##
## Writes public_data_tierA/derived/hallmark_revised/ and leaves the published
## outputs in public_data_tierA/derived/hallmark_audit/ untouched.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db); library(fgsea)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"
O <- file.path(D, "hallmark_revised"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")

## ---------------------------------------------------------------- gene sets --
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
P <- P[sapply(P, length) >= 20]; names(P) <- sub("__.*", "", names(P))
## the cell-cycle union used elsewhere in the project (diag_contrast_null.R)
CC <- unique(unlist(P[grep("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(P))]))
CC <- unique(c(CC, P[["Cell Cycle"]]))
say("cell-cycle union: %d genes", length(CC))

HALL <- list(
 `pre-mRNA processing`      = "Processing of Capped Intron-Containing Pre-mRNA",
 `mRNA splicing`            = "mRNA Splicing - Major Pathway",
 `translation`              = "Translation",
 `rRNA processing`          = "rRNA processing in the nucleus and cytosol",
 `proteasome`               = "Antigen processing: Ubiquitination & Proteasome degradation",
 `autophagy`                = "Macroautophagy",
 `chaperone / HSF1`         = "HSF1-dependent transactivation",
 `unfolded protein response`= "Unfolded Protein Response (UPR)",
 `DNA double-strand repair` = "DNA Double-Strand Break Repair",
 `base excision repair`     = "Base Excision Repair",
 `telomere maintenance`     = "Telomere Maintenance",
 `chromatin organisation`   = "Chromatin organization",
 `mitochondrial translation`= "Mitochondrial translation",
 `respiratory electron transport` = "Respiratory electron transport",
 `TCA cycle`                = "Citric acid cycle (TCA cycle)",
 `interferon signalling`    = "Interferon Signaling",
 `NF-kB / TNF`              = "TNFR1-induced NF-kappa-B signaling pathway",
 `extracellular matrix`     = "Extracellular matrix organization",
 `collagen formation`       = "Collagen formation",
 `insulin / IGF signalling` = "Signaling by Insulin receptor",
 `mTOR signalling`          = "mTORC1-mediated signalling",
 `cellular senescence`      = "Cellular Senescence",
 `interleukin signalling`   = "Signaling by Interleukins",
 `cell cycle (positive control)` = "Cell Cycle, Mitotic")
pw <- function(nm) { hit <- if (nm %in% names(P)) nm else grep(nm, names(P), fixed = TRUE, value = TRUE)[1]
                     if (is.na(hit)) NULL else P[[hit]] }

## =============================================================================
## 1. division rate, with and without the cell-cycle genes inside each programme
## =============================================================================
d0 <- read.delim(file.path(D, "division_rate/sample_division_rate.tsv"))
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                check.names = FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", raw[[1]]),
        keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")))
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok, ], toupper(sy[ok]))
cnt <- cnt[, match(d0$sample, colnames(cnt))]
y <- DGEList(cnt); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
lcD <- cpm(y, log = TRUE); zD <- t(scale(t(lcD[apply(lcD, 1, sd) > 0, ])))
say("GSE179848: %d genes x %d counted libraries", nrow(zD), ncol(zD))

div <- do.call(rbind, lapply(names(HALL), function(h) {
  set <- pw(HALL[[h]]); if (is.null(set)) return(NULL)
  g  <- intersect(rownames(zD), set)
  gk <- setdiff(g, CC)                                   # cell-cycle genes removed
  if (length(g) < 15) return(NULL)
  r_all <- cor(colMeans(zD[g, , drop = FALSE]), d0$rate, method = "spearman")
  if (length(gk) >= 15) {
    sc <- colMeans(zD[gk, , drop = FALSE])
    ct <- cor.test(sc, d0$rate, method = "spearman", exact = FALSE)
    r_kept <- ct$estimate; p_kept <- ct$p.value
  } else { r_kept <- NA_real_; p_kept <- NA_real_ }
  data.frame(programme = h, n_genes = length(g), n_cc_removed = length(g) - length(gk),
             n_kept = length(gk), rho_rate = r_all, rho_rate_noCC = r_kept,
             p_noCC = p_kept, row.names = NULL) }))
div$FDR_noCC <- p.adjust(div$p_noCC, "BH")
div <- div[order(-abs(div$rho_rate)), ]
fwrite(div, file.path(O, "division_rate_cellcycle_removed.tsv"), sep = "\t")
say("\n=== division rate: programmes before and after removing cell-cycle genes ===")
print(data.frame(programme = div$programme, genes = div$n_genes, removed = div$n_cc_removed,
                 rho = round(div$rho_rate, 2), rho_noCC = round(div$rho_rate_noCC, 2),
                 FDR = signif(div$FDR_noCC, 2)), row.names = FALSE)

EIGHT <- c("chromatin organisation", "cellular senescence", "pre-mRNA processing", "mRNA splicing",
           "interferon signalling", "DNA double-strand repair", "telomere maintenance", "base excision repair")
say("\nthe eight programmes called division readouts:")
for (h in EIGHT) { r <- div[div$programme == h, ]
  say("  %-26s rho %.2f -> %.2f after removing %d of %d genes", h, r$rho_rate, r$rho_rate_noCC,
      r$n_cc_removed, r$n_genes) }
ccrow <- div[div$programme == "cell cycle (positive control)", ]
say("  reference: the cell-cycle set itself rho = %.2f (%d genes; %d remain after removal)",
    ccrow$rho_rate, ccrow$n_genes, ccrow$n_kept)

## =============================================================================
## 2. donor age in the redefined primary cohort (recount3, normal adults 20+)
## =============================================================================
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))
DON <- merge(DON, MET[, .(srr, log_depth)], by = "srr")
A <- DON[gse113957_keep(DON, "primary"), ]
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
X <- GM[, A$srr, drop = FALSE]
yA <- DGEList(X); yA <- yA[filterByExpr(yA), , keep.lib.sizes = FALSE]; yA <- calcNormFactors(yA)
lcA <- cpm(yA, log = TRUE); zA <- t(scale(t(lcA[apply(lcA, 1, sd) > 0, ])))
A$prolif <- colMeans(zA[intersect(rownames(zA), PROLIF), , drop = FALSE])
say("\nprimary cohort: %d donors, ages %d-%d, %d expressed genes; rho(age, proliferation) = %+.2f",
    nrow(A), min(A$age), max(A$age), nrow(zA), cor(A$age, A$prolif, method = "spearman"))

pubJ <- read.delim(file.path(D, "hallmark_audit/hallmark_division_and_age.tsv"))
ageres <- do.call(rbind, lapply(names(HALL), function(h) {
  set <- pw(HALL[[h]]); if (is.null(set)) return(NULL)
  g <- intersect(rownames(zA), set); if (length(g) < 15) return(NULL)
  A$s <- colMeans(zA[g, , drop = FALSE])
  f0 <- lm(s ~ I(age/10) + log_depth + repo + instr + sex, A)
  f1 <- lm(s ~ I(age/10) + log_depth + repo + instr + sex + prolif, A)
  c0 <- summary(f0)$coef["I(age/10)", ]; c1 <- summary(f1)$coef["I(age/10)", ]
  ci <- c0[1] + c(-1, 1) * qt(.975, f0$df.residual) * c0[2]
  data.frame(programme = h, n_genes = length(g),
             rho_age = cor(A$s, A$age, method = "spearman"),
             beta_decade = c0[1], lo = ci[1], hi = ci[2], p_age = c0[4],
             beta_adj = c1[1], p_adj = c1[4], row.names = NULL) }))
ageres$FDR_age <- p.adjust(ageres$p_age, "BH")
ageres$FDR_adj <- p.adjust(ageres$p_adj, "BH")
J <- merge(ageres, pubJ[, c("programme", "n_genes", "rho_rate", "FDR", "rho_age", "FDR_age")],
           by = "programme", suffixes = c("", "_published"))
J <- J[order(-J$rho_rate), ]
fwrite(J, file.path(O, "donor_age_primary_cohort.tsv"), sep = "\t")
say("\n=== donor age: published (97 donors, FPKM, 22-89) vs primary cohort (%d donors, recount3, 20+) ===", nrow(A))
print(data.frame(programme = J$programme, rho_div = round(J$rho_rate, 2),
                 rho_age_pub = round(J$rho_age_published, 2), FDR_pub = signif(J$FDR_age_published, 2),
                 rho_age_new = round(J$rho_age, 2), beta10 = round(J$beta_decade, 4),
                 FDR_new = signif(J$FDR_age, 2), FDR_adj = signif(J$FDR_adj, 2)), row.names = FALSE)

## =============================================================================
## 3. senescence set dissection (division rate) and the SenMayo question
## =============================================================================
sen <- intersect(rownames(zD), P[["Cellular Senescence"]])
gr  <- ifelse(sen %in% CC, "also in cell cycle", "senescence set only")
rg  <- sapply(sen, function(g) cor(zD[g, ], d0$rate, method = "spearman"))
SG  <- data.frame(gene = sen, group = gr, rho = rg, row.names = NULL)
w <- wilcox.test(rho ~ group, SG)
say("\nsenescence set: %d genes, %d also cell cycle; median rho %.2f vs %.2f, Wilcoxon P = %.2g",
    nrow(SG), sum(gr == "also in cell cycle"),
    median(SG$rho[gr == "also in cell cycle"]), median(SG$rho[gr == "senescence set only"]), w$p.value)
say("aggregate score of the senescence set with cell-cycle members removed: rho = %+.2f",
    div$rho_rate_noCC[div$programme == "cellular senescence"])
fwrite(SG, file.path(O, "senescence_dissection_cellcycle_union.tsv"), sep = "\t")
sm <- length(list.files(".", pattern = "(?i)senmayo|saul_sen", recursive = TRUE))
say("SenMayo / SAUL_SEN_MAYO files present in the project: %d (nothing downloaded)", sm)
sessionInfo_path <- file.path(O, "sessionInfo.txt")
capture.output(sessionInfo(), file = sessionInfo_path)
