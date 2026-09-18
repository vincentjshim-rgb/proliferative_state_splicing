## The control the senescence field actually uses: verify arrest directly.
## Sample-level MKI67 (Ki-67) against the pre-mRNA processing core, 345 samples.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
D <- "public_data_tierA/derived"; O <- file.path(D, "core_depth")
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
meta <- read.delim(file.path(D, "audit/all_geo_samples.tsv"), check.names = FALSE)
gm <- meta[meta$dataset == "GSE179848", ]
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                check.names = FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$","",raw[[1]]),
        keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")))
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok, ], toupper(sy[ok]))
gm$matrix_name <- paste0("Sample_", gm$rnaseq_sampleid); gm <- gm[match(colnames(cnt), gm$matrix_name), ]
gm$days <- as.numeric(gm$days_grown_udays)

y <- DGEList(cnt); y <- y[filterByExpr(y, group = factor(gm$treatments)), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
zz <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
cg <- intersect(rownames(zz), CORE)
d <- data.frame(sample = colnames(lc), line = gm$cell_line, cond = gm$treatments,
                oxy = as.character(gm$percent_oxygen), days = gm$days,
                core = colMeans(zz[cg, , drop = FALSE]),
                MKI67 = if ("MKI67" %in% rownames(lc)) lc["MKI67", ] else NA,
                CCNB1 = if ("CCNB1" %in% rownames(lc)) lc["CCNB1", ] else NA,
                CDKN2A = if ("CDKN2A" %in% rownames(lc)) lc["CDKN2A", ] else NA,
                CDKN1A = if ("CDKN1A" %in% rownames(lc)) lc["CDKN1A", ] else NA)
say("samples %d | core genes %d | MKI67 available: %s", nrow(d), length(cg), !all(is.na(d$MKI67)))

ct <- cor.test(d$MKI67, d$core, method = "spearman", exact = FALSE)
say("\nMKI67 vs pre-mRNA processing core, all %d samples: rho = %+.3f, P = %.2e",
    nrow(d), ct$estimate, ct$p.value)
ct2 <- cor.test(d$CDKN2A, d$core, method = "spearman", exact = FALSE)
say("CDKN2A (p16) vs core: rho = %+.3f, P = %.2e", ct2$estimate, ct2$p.value)
ct3 <- cor.test(d$CDKN1A, d$core, method = "spearman", exact = FALSE)
say("CDKN1A (p21) vs core: rho = %+.3f, P = %.2e", ct3$estimate, ct3$p.value)

## within untreated healthy donors only
hc <- d[grepl("^HC", d$line) & d$cond == "Control" & d$oxy == "21", ]
ct4 <- cor.test(hc$MKI67, hc$core, method = "spearman", exact = FALSE)
say("\nhealthy untreated 21%% O2 (n = %d): MKI67 vs core rho = %+.3f, P = %.2e",
    nrow(hc), ct4$estimate, ct4$p.value)

## is the core still tied to days once MKI67 is in the model?
f1 <- lm(core ~ days + line, hc); f2 <- lm(core ~ days + MKI67 + line, hc)
say("days alone      : beta = %+.5f, P = %.2e", coef(f1)["days"], summary(f1)$coef["days","Pr(>|t|)"])
say("days + MKI67    : beta = %+.5f, P = %.3f   | MKI67 beta = %+.4f, P = %.2e",
    coef(f2)["days"], summary(f2)$coef["days","Pr(>|t|)"],
    coef(f2)["MKI67"], summary(f2)$coef["MKI67","Pr(>|t|)"])

## contact inhibition: arrest without senescence
ci <- d[grepl("^HC", d$line) & d$cond %in% c("Control","Contact_Inhibition") & d$oxy == "21", ]
ci$cond <- factor(ci$cond, levels = c("Control","Contact_Inhibition"))
say("\ncontact inhibition (n=%d) vs control (n=%d):", sum(ci$cond!="Control"), sum(ci$cond=="Control"))
say("  core     : %+.3f vs %+.3f, Wilcoxon P = %.2e", median(ci$core[ci$cond!="Control"]),
    median(ci$core[ci$cond=="Control"]), wilcox.test(core ~ cond, ci)$p.value)
say("  MKI67    : %+.2f vs %+.2f, Wilcoxon P = %.2e", median(ci$MKI67[ci$cond!="Control"]),
    median(ci$MKI67[ci$cond=="Control"]), wilcox.test(MKI67 ~ cond, ci)$p.value)
say("  CDKN2A   : %+.2f vs %+.2f, Wilcoxon P = %.3f", median(ci$CDKN2A[ci$cond!="Control"]),
    median(ci$CDKN2A[ci$cond=="Control"]), wilcox.test(CDKN2A ~ cond, ci)$p.value)
write.table(d, file.path(O, "GSE179848_sample_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
