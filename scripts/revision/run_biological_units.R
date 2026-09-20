## The headline quantities in units a reader can picture.
##
## Almost everything in the manuscript is a mean gene-wise z-score, so effects are reported in
## standard deviations. That is the right unit for comparing datasets whose scales differ, but
## it says nothing about how much RNA changes. This script recomputes the same relationships
## from the set's mean log2 counts per million, which converts directly into a fold change, and
## counts how many of the 177 genes move individually with FDR.
## Output: public_data_tierA/derived/biological_units/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(limma); library(data.table)
                                library(AnnotationDbi); library(org.Hs.eg.db)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "biological_units")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
pct <- function(l2) 100 * (2^l2 - 1)
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
RES <- list()

## =====================================================================================
## 1. The counted resource: how much RNA per division per day
## =====================================================================================
say("=== 1. the counted resource (GSE179848) ===")
gm <- read.delim(file.path(D, "audit/all_geo_samples.tsv"), check.names = FALSE)
gm <- gm[gm$dataset == "GSE179848", ]
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                check.names = FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$", "", raw[[1]]),
        keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")))
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok, ], toupper(sy[ok]))
gm$mn <- paste0("Sample_", gm$rnaseq_sampleid); gm <- gm[match(colnames(cnt), gm$mn), ]
y <- DGEList(cnt); y <- y[filterByExpr(y, group = factor(gm$treatments)), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
d <- data.frame(line = gm$cell_line, cond = gm$treatments,
                pdt = as.numeric(gm$population_doubling_time_uhours_per_division))
d$rate <- 24 / d$pdt
keep <- is.finite(d$rate) & d$rate > 0; d <- d[keep, ]; lc <- lc[, keep]
cg <- intersect(rownames(lc), CORE)
d$set <- colMeans(lc[cg, , drop = FALSE])                       # mean log2 CPM of the set
say("  %d libraries, %d of the 177 genes measured; the set averages %.2f log2 CPM",
    nrow(d), length(cg), mean(d$set))
m <- lm(set ~ rate, d)
say("  per division per day: %+.3f log2 CPM (%+.0f%%), P = %.1e",
    coef(m)[2], pct(coef(m)[2]), coef(summary(m))[2, 4])
say("  per 0.1 divisions per day: %+.3f log2 CPM (%+.1f%%)", coef(m)[2]/10, pct(coef(m)[2]/10))
dec <- cut(d$rate, quantile(d$rate, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
lo <- mean(d$set[dec == 1]); hi <- mean(d$set[dec == 10])
say("  slowest tenth %.2f log2 CPM, fastest tenth %.2f: the set is %.2f-fold (%+.0f%%) higher",
    lo, hi, 2^(hi - lo), pct(hi - lo))
## per gene, against the counted rate
des <- model.matrix(~ rate, d)
fit <- eBayes(lmFit(lc[cg, ], des))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
say("  of the %d genes, %d rise with the counted rate at FDR < 0.05 and %d fall",
    length(cg), sum(tt$adj.P.Val < 0.05 & tt$logFC > 0), sum(tt$adj.P.Val < 0.05 & tt$logFC < 0))
## the set score is a mean of gene-wise z-scores, which weights every gene alike; the median
## gene is the matching summary in fold-change units, and is what the manuscript quotes
say("  the median gene: %+.3f log2 CPM per division per day (%+.0f%%), %+.3f per 0.1 (%+.1f%%)",
    median(tt$logFC), pct(median(tt$logFC)), median(tt$logFC)/10, pct(median(tt$logFC)/10))
gl <- rowMeans(lc[cg, dec == 1, drop = FALSE]); gh <- rowMeans(lc[cg, dec == 10, drop = FALSE])
say("  fastest tenth against slowest, median gene: %.2f-fold (%+.0f%%); genes higher in the fastest tenth: %d of %d",
    2^median(gh - gl), pct(median(gh - gl)), sum(gh > gl), length(cg))
RES$rate <- data.table(
  quantity = c("set mean, per division per day", "set mean, fastest vs slowest tenth",
               "median gene, per division per day", "median gene, per 0.1 divisions per day",
               "median gene, fastest vs slowest tenth"),
  log2FC = c(coef(m)[2], hi - lo, median(tt$logFC), median(tt$logFC)/10, median(gh - gl)),
  percent = pct(c(coef(m)[2], hi - lo, median(tt$logFC), median(tt$logFC)/10, median(gh - gl))),
  P = c(coef(summary(m))[2, 4], NA, NA, NA, NA))
RES$rate_genes <- data.table(gene = rownames(tt), log2FC_per_div = tt$logFC,
                             P = tt$P.Value, FDR = tt$adj.P.Val)

## =====================================================================================
## 2. The donor cohort: how much RNA per decade, before and after proliferation
## =====================================================================================
say("\n=== 2. the donor cohort (GSE113957, 107 adults) ===")
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
M <- fread(file.path(D, "cohort_revised/primary/sample_metrics.tsv"))
## the same identifier path as run_cohort_core.R: recount3 columns are SRR accessions and the
## gene ids are mapped through the deposited GENCODE table rather than through org.Hs.eg.db
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
X <- GM[, M$srr, drop = FALSE]
stopifnot(ncol(X) == nrow(M))
yc <- DGEList(X); yc <- yc[filterByExpr(yc), , keep.lib.sizes = FALSE]
yc <- calcNormFactors(yc); lcd <- cpm(yc, log = TRUE)
cd <- intersect(rownames(lcd), CORE)
M[, set := colMeans(lcd[cd, , drop = FALSE])]
say("  %d donors, %d of the 177 genes measured; the set averages %.2f log2 CPM",
    nrow(M), length(cd), mean(M$set))
M[, sex := factor(sex)][, repo := factor(repo)][, instr := factor(instr)]
m0 <- lm(set ~ I(age/10) + repo + instr + sex + log_depth, M)
m1 <- lm(set ~ I(age/10) + prolif + repo + instr + sex + log_depth, M)
b0 <- coef(m0)["I(age/10)"]; b1 <- coef(m1)["I(age/10)"]
say("  age effect  %+.4f log2 CPM per decade (%+.1f%%), P = %.1e",
    b0, pct(b0), coef(summary(m0))["I(age/10)", 4])
say("  after proliferation %+.4f log2 CPM per decade (%+.1f%%), P = %.3f",
    b1, pct(b1), coef(summary(m1))["I(age/10)", 4])
say("  over the cohort's 76-year span: %+.0f%% before adjustment, %+.0f%% after",
    pct(b0 * 7.6), pct(b1 * 7.6))
dsn <- model.matrix(~ I(age/10) + repo + instr + sex + log_depth, M)
tg <- topTable(eBayes(lmFit(lcd[cd, ], dsn)), coef = 2, number = Inf, sort.by = "none")
dsn1 <- model.matrix(~ I(age/10) + prolif + repo + instr + sex + log_depth, M)
tg1 <- topTable(eBayes(lmFit(lcd[cd, ], dsn1)), coef = 2, number = Inf, sort.by = "none")
say("  per gene: %d of %d fall with age at FDR < 0.05; across all %d the median gene changes",
    sum(tg$adj.P.Val < 0.05 & tg$logFC < 0), length(cd), length(cd))
say("    %+.4f log2 CPM per decade (%+.1f%%), and %+.4f (%+.1f%%) with proliferation in the model",
    median(tg$logFC), pct(median(tg$logFC)), median(tg1$logFC), pct(median(tg1$logFC)))
say("    over the cohort's 76-year span the median gene falls %.0f%%", -pct(median(tg$logFC) * 7.6))
say("  with proliferation in the model, %d keep it", sum(tg1$adj.P.Val < 0.05 & tg1$logFC < 0))
RES$age <- data.table(model = c("set mean, age", "set mean, age + proliferation",
                                "median gene, age", "median gene, age + proliferation"),
                      log2FC_per_decade = c(b0, b1, median(tg$logFC), median(tg1$logFC)),
                      percent_per_decade = pct(c(b0, b1, median(tg$logFC), median(tg1$logFC))),
                      P = c(coef(summary(m0))["I(age/10)", 4], coef(summary(m1))["I(age/10)", 4], NA, NA),
                      genes_FDR05 = c(sum(tg$adj.P.Val < 0.05 & tg$logFC < 0),
                                      sum(tg1$adj.P.Val < 0.05 & tg1$logFC < 0), NA, NA))
RES$age_genes <- data.table(gene = rownames(tg), log2FC_per_decade = tg$logFC,
                            P = tg$P.Value, FDR = tg$adj.P.Val,
                            log2FC_adj = tg1$logFC, P_adj = tg1$P.Value, FDR_adj = tg1$adj.P.Val)

## =====================================================================================
## 3. GTEx: the same slope in log2 units
## =====================================================================================
say("\n=== 3. GTEx, per standard deviation of the proliferation score ===")
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
G <- rbindlist(lapply(c("culture", "legskin", "pubskin"), function(tk) {
  dd <- as.data.frame(TS[[tk]]$d); s <- TS[[tk]]$S[, "splicing 96"]
  o <- complete.cases(dd$prolif, dd$rin, dd$isch, dd$loglib, s); dd <- dd[o, ]; s <- s[o]
  rx <- resid(lm(prolif ~ rin + isch + loglib, dd)); ry <- resid(lm(s ~ rin + isch + loglib, dd))
  data.table(tissue = tk, n = nrow(dd), slope_sd = coef(lm(scale(ry) ~ scale(rx)))[2],
             rho = cor(rx, ry, method = "spearman")) }))
print(G, row.names = FALSE, digits = 3)
say("  the GTEx scores are z-scores of gene-wise z-scores, so they carry no CPM unit;")
say("  the correlation is the quantity reported for those tissues.")
RES$gtex <- G

for (nm in names(RES)) fwrite(RES[[nm]], file.path(O, paste0(nm, ".tsv")), sep = "\t")
say("\nwritten to %s", O)
