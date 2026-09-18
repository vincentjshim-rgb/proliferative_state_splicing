## Which ageing-associated transcriptional programmes in cultured fibroblasts are
## readouts of division rate?  Measured divisions per day, 328 libraries.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db); library(fgsea)})
D <- "public_data_tierA/derived"; O <- file.path(D, "hallmark_audit"); dir.create(O, showWarnings=FALSE)
say <- function(...) cat(sprintf(...), "\n")
d0 <- read.delim(file.path(D, "division_rate/sample_division_rate.tsv"))
gm <- read.delim(file.path(D, "audit/all_geo_samples.tsv"), check.names = FALSE)
gm <- gm[gm$dataset == "GSE179848", ]
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"), check.names=FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys=sub("\\..*$","",raw[[1]]),
        keytype="REFSEQ", column="SYMBOL", multiVals="first")))
cnt <- as.matrix(raw[,-1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok,], toupper(sy[ok]))
gm$mn <- paste0("Sample_", gm$rnaseq_sampleid); gm <- gm[match(colnames(cnt), gm$mn), ]
cnt <- cnt[, match(d0$sample, colnames(cnt))]
y <- DGEList(cnt); y <- y[filterByExpr(y), , keep.lib.sizes=FALSE]; y <- calcNormFactors(y)
lc <- cpm(y, log=TRUE); z <- t(scale(t(lc[apply(lc,1,sd)>0, ])))
say("genes %d, samples %d", nrow(z), ncol(z))

td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir=td)
P <- gmtPathways(list.files(td, pattern="\\.gmt$", full.names=TRUE)[1])
P <- P[sapply(P, length) >= 20]
names(P) <- sub("__.*", "", names(P))

## programmes standing for the recognised hallmarks of ageing
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
res <- do.call(rbind, lapply(names(HALL), function(h) {
  nm <- HALL[[h]]; hit <- names(P)[names(P) == nm]
  if (!length(hit)) { alt <- grep(nm, names(P), fixed=TRUE, value=TRUE); hit <- alt[1] }
  if (is.na(hit) || !length(hit)) return(NULL)
  g <- intersect(rownames(z), P[[hit]]); if (length(g) < 15) return(NULL)
  sc <- colMeans(z[g, , drop=FALSE])
  ct <- cor.test(sc, d0$rate, method="spearman", exact=FALSE)
  cd <- cor.test(sc, d0$days, method="spearman", exact=FALSE)
  data.frame(programme=h, n_genes=length(g), rho_rate=ct$estimate, p_rate=ct$p.value,
             rho_days=cd$estimate, row.names=NULL) }))
res$FDR <- p.adjust(res$p_rate, "BH")
res <- res[order(-abs(res$rho_rate)), ]
say("\n=== ageing-associated programmes against MEASURED division rate (n = %d) ===", ncol(z))
print(data.frame(programme=res$programme, genes=res$n_genes,
      rho_division=round(res$rho_rate,2), rho_days=round(res$rho_days,2),
      FDR=signif(res$FDR,2)), row.names=FALSE)
say("\nprogrammes with |rho| > 0.5 against division rate: %d of %d",
    sum(abs(res$rho_rate) > 0.5), nrow(res))
write.table(res, file.path(O, "hallmark_vs_division.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
sc_all <- t(sapply(names(HALL), function(h) { nm <- HALL[[h]]
  hit <- if (nm %in% names(P)) nm else grep(nm, names(P), fixed=TRUE, value=TRUE)[1]
  if (is.na(hit)) return(rep(NA, ncol(z)))
  g <- intersect(rownames(z), P[[hit]]); if (length(g) < 15) return(rep(NA, ncol(z)))
  colMeans(z[g, , drop=FALSE]) }))
saveRDS(list(scores=sc_all, meta=d0), file.path(O, "programme_scores.rds"))

## ==========================================================================
## 2.  the same programmes against DONOR AGE in 97 adult donors (GSE113957)
##     a programme can be a division readout, an age readout, both or neither
## ==========================================================================
x <- read.delim(gzfile("public_data_tierA/aging/GSE113957_fpkm.txt.gz"), check.names = FALSE)
sym <- toupper(sub("[;|].*$", "", x[["Annotation/Divergence"]]))
ANN <- c("Transcript ID","chr","start","end","strand","Length","Copies","Annotation/Divergence")
scn <- setdiff(names(x), ANN)
age <- suppressWarnings(as.numeric(sub("^[0-9]+_([0-9]+)(yr|YR)?_.*$", "\\1", scn)))
sex <- tolower(sub("^[0-9]+_[0-9]+(yr|YR)?_([A-Za-z]+)_.*$", "\\2", scn))
keep <- is.finite(age) & sex %in% c("male","female") & age >= 22 & age <= 89
scn <- scn[keep]; age <- age[keep]; sex <- sex[keep]
mA <- as.matrix(x[, scn]); rownames(mA) <- sym
mA <- mA[!is.na(rownames(mA)) & rownames(mA) != "", ]
mA <- mA[rowMeans(mA > 1) > 0.5, ]
zA <- t(scale(t(log2(mA + 1))))
say("\nGSE113957: %d adult donors, %d expressed genes", ncol(zA), nrow(zA))
ageres <- do.call(rbind, lapply(names(HALL), function(h) {
  nm <- HALL[[h]]; hit <- if (nm %in% names(P)) nm else grep(nm, names(P), fixed=TRUE, value=TRUE)[1]
  if (is.na(hit)) return(NULL)
  g <- intersect(rownames(zA), P[[hit]]); if (length(g) < 15) return(NULL)
  s <- colMeans(zA[g, , drop=FALSE]); f <- summary(lm(s ~ age + sex))$coef
  data.frame(programme=h, n_age=length(g), beta_age=f["age","Estimate"],
             t_age=f["age","t value"], p_age=f["age","Pr(>|t|)"],
             rho_age=cor(s, age, method="spearman"), row.names=NULL) }))
ageres$FDR_age <- p.adjust(ageres$p_age, "BH")
J <- merge(res, ageres, by="programme")
J <- J[order(-J$rho_rate), ]
say("\n=== division readout vs donor-age readout ===")
print(data.frame(programme=J$programme, rho_division=round(J$rho_rate,2),
      FDR_div=signif(J$FDR,1), rho_age=round(J$rho_age,2), FDR_age=signif(J$FDR_age,1)),
      row.names=FALSE)
write.table(J, file.path(O, "hallmark_division_and_age.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ==========================================================================
## 3.  why does the senescence gene set rise with division rate?
##     per-gene correlation, split by cell-cycle membership
## ==========================================================================
sen <- intersect(rownames(z), P[["Cellular Senescence"]])
cc  <- unique(unlist(P[c("Cell Cycle, Mitotic","Cell Cycle")]))
gr  <- ifelse(sen %in% cc, "also in cell cycle", "senescence set only")
rg  <- sapply(sen, function(g) cor(z[g, ], d0$rate, method="spearman"))
MK  <- c("CDKN1A","CDKN2A","CDKN2B","MDM2","ATM","TP53","LMNB1","RB1","E2F1","CDK4","CDK6","TINF2")
SG  <- data.frame(gene=sen, group=gr, rho=rg, marker=sen %in% MK, row.names=NULL)
say("\n=== Reactome 'Cellular Senescence' (%d genes) dissected ===", nrow(SG))
print(as.data.frame(do.call(rbind, lapply(split(SG, SG$group), function(d)
  data.frame(n=nrow(d), median_rho=round(median(d$rho),2),
             frac_positive=round(mean(d$rho > 0),2))))))
w <- wilcox.test(rho ~ group, SG); say("  Wilcoxon between groups: P = %.2g", w$p.value)
say("  canonical markers (a biologist's senescence genes):")
print(SG[SG$marker, c("gene","group","rho")][order(SG$rho[SG$marker]), ], row.names=FALSE)
say("  aggregate score of the whole 150-gene set vs division rate: rho = %+.2f",
    cor(colMeans(z[sen, ]), d0$rate, method="spearman"))
write.table(SG, file.path(O, "senescence_geneset_dissected.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(z=z, meta=d0, zA=zA, age=age, sex=sex,
             pathways=P[unlist(HALL)[unlist(HALL) %in% names(P)]]),
        file.path(O, "audit_inputs.rds"))
