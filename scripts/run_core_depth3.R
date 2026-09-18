.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
D <- "public_data_tierA/derived"; O <- file.path(D, "core_depth")
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
modscore <- function(mat, set) { mat <- mat[apply(mat, 1, function(v) all(is.finite(v)) && sd(v) > 0), , drop=FALSE]
  z <- t(scale(t(mat))); g <- intersect(rownames(z), set)
  list(score = colMeans(z[g,,drop=FALSE]), n = length(g)) }

meta <- read.delim(file.path(D, "audit/all_geo_samples.tsv"), check.names = FALSE)
gm <- meta[meta$dataset == "GSE179848", ]
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                check.names = FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$","",raw[[1]]),
        keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")))
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok, ], toupper(sy[ok]))
gm$matrix_name <- paste0("Sample_", gm$rnaseq_sampleid)
gm <- gm[match(colnames(cnt), gm$matrix_name), ]
gm$days <- as.numeric(gm$days_grown_udays)
say("GSE179848: %d genes x %d samples; days available for %d", nrow(cnt), ncol(cnt), sum(is.finite(gm$days)))

## healthy donors, untreated, standard oxygen - the replicative-ageing series
keep <- gm$clinical_condition %in% c("Control","control") | grepl("^HC", gm$cell_line)
keep <- keep & gm$treatments == "Control" & as.character(gm$percent_oxygen) == "21" &
        is.finite(gm$days) & grepl("^HC", gm$cell_line)
say("healthy untreated 21%% O2 samples: %d across %d donors", sum(keep),
    length(unique(gm$cell_line[keep])))
y <- DGEList(cnt[, keep]); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
ms <- modscore(lc, CORE)
d <- data.frame(sample = colnames(lc), line = gm$cell_line[keep], days = gm$days[keep],
                score = ms$score)
say("core genes measured: %d", ms$n)
f <- lm(score ~ days + line, d)
say("score ~ days + donor:  beta = %+.5f per day, t = %.2f, P = %.2e",
    coef(f)["days"], summary(f)$coef["days","t value"], summary(f)$coef["days","Pr(>|t|)"])
for (l in sort(unique(d$line))) { s <- d[d$line == l, ]
  ct <- cor.test(s$days, s$score, method = "spearman", exact = FALSE)
  say("   %-5s n=%2d  days %3.0f-%3.0f  rho = %+.3f  P = %.3f", l, nrow(s),
      min(s$days), max(s$days), ct$estimate, ct$p.value) }
write.table(d, file.path(O, "GSE179848_longitudinal_core.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## contact inhibition vs control at matched time - arrest without senescence
k2 <- gm$treatments %in% c("Control","Contact_Inhibition") & as.character(gm$percent_oxygen) == "21" &
      grepl("^HC", gm$cell_line) & is.finite(gm$days)
y2 <- DGEList(cnt[, k2]); y2 <- y2[filterByExpr(y2), , keep.lib.sizes = FALSE]
y2 <- calcNormFactors(y2); lc2 <- cpm(y2, log = TRUE)
ms2 <- modscore(lc2, CORE)
d2 <- data.frame(sample = colnames(lc2), line = gm$cell_line[k2], days = gm$days[k2],
                 cond = gm$treatments[k2], score = ms2$score)
d2$cond <- factor(d2$cond, levels = c("Control","Contact_Inhibition"))
f2 <- lm(score ~ cond + days + line, d2)
cf <- summary(f2)$coef; rn <- grep("^cond", rownames(cf), value = TRUE)[1]
say("\ncontact inhibition (n=%d) vs control (n=%d), days and donor adjusted:",
    sum(d2$cond=="Contact_Inhibition"), sum(d2$cond=="Control"))
say("   beta = %+.4f, t = %.2f, P = %.3g", cf[rn,1], cf[rn,3], cf[rn,4])
say("   raw medians: control %+.3f | contact inhibition %+.3f",
    median(d2$score[d2$cond=="Control"]), median(d2$score[d2$cond=="Contact_Inhibition"]))
write.table(d2, file.path(O, "GSE179848_contact_inhibition_core.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
