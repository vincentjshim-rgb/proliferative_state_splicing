## Gene-level verification of the coupling, and the counted division rate read three ways.
##
## Everything in the manuscript so far compares a SET SCORE with a rate or with another set
## score. Two objections follow from that and neither has been answered in the figures:
##   (i)  a mean z-score can be carried by a handful of genes, or by a global compositional
##        shift, so the set may not behave as a set;
##   (ii) one correlation on one construction of the rate (divisions per day) may be an
##        artefact of that construction or of assuming a straight line.
## This script answers both from the counted resource (GSE179848) without new data:
##   A  every expressed gene against the counted rate, and where the 177 sit in that distribution
##   B  the same against an expression-matched background drawn 1,000 times
##   C  the rate in its native unit (hours per doubling), by decile, and at fixed replicative age
##   D  per-gene values for the named factors, for small-multiple scatters
## Output: public_data_tierA/derived/gene_level/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
D <- "public_data_tierA/derived"; O <- file.path(D, "gene_level")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
set.seed(20260920)

CORE  <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # the 177
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")                  # Holly 2013, 20

## ---- the same preparation as run_division_rate.R ---------------------------------
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

d <- data.frame(sample = colnames(lc), line = gm$cell_line, cond = gm$treatments,
                oxy = as.character(gm$percent_oxygen),
                days  = as.numeric(gm$days_grown_udays),
                pdt   = as.numeric(gm$population_doubling_time_uhours_per_division),
                pdtot = as.numeric(gm$population_doublings_utotal_divisions))
d$rate <- 24 / d$pdt
keep <- is.finite(d$rate) & d$rate > 0
d <- d[keep, ]; lc <- lc[, keep]
z <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
cg <- intersect(rownames(z), CORE)
d$splice <- colMeans(z[cg, , drop = FALSE])
say("libraries with a counted rate: %d | genes after filtering: %d | splicing set measured: %d of %d",
    nrow(d), nrow(lc), length(cg), length(CORE))
say("set score vs counted rate: rho = %.3f  (the value Fig. 2a reports)",
    cor(d$rate, d$splice, method = "spearman"))

## ---- A. every expressed gene against the counted rate ----------------------------
say("\nA. every expressed gene against the counted rate")
G <- lc[rownames(z), , drop = FALSE]
rho <- apply(G, 1, function(v) cor(v, d$rate, method = "spearman"))
expr <- rowMeans(G)
GL <- data.frame(gene = rownames(G), rho = rho, mean_logcpm = expr,
                 in_set = rownames(G) %in% CORE, named = rownames(G) %in% NAMED)
bg <- GL$rho[!GL$in_set]
q95 <- quantile(bg, 0.95); q99 <- quantile(bg, 0.99)
st <- GL$rho[GL$in_set]
say("  background (%d genes not in the set): median rho = %.3f, 95th pct = %.3f, 99th = %.3f",
    length(bg), median(bg), q95, q99)
say("  the 177 set (%d measured):            median rho = %.3f, range %.2f to %.2f",
    length(st), median(st), min(st), max(st))
say("  of the set, %d of %d (%.0f%%) are above the background 95th percentile, %d (%.0f%%) above the 99th",
    sum(st > q95), length(st), 100 * mean(st > q95), sum(st > q99), 100 * mean(st > q99))
say("  positive rho: %.0f%% of the set against %.0f%% of the background",
    100 * mean(st > 0), 100 * mean(bg > 0))
say("  Wilcoxon set vs background: P = %.3g", wilcox.test(st, bg)$p.value)

## how much of the set score is carried by its strongest genes
ord <- sort(st, decreasing = TRUE)
drop_top <- function(k) { g <- setdiff(cg, names(sort(rho[cg], decreasing = TRUE))[seq_len(k)])
                          cor(d$rate, colMeans(z[g, , drop = FALSE]), method = "spearman") }
say("  set score with the top %d genes removed: %.3f | top 20 removed: %.3f | top 50 removed: %.3f",
    10, drop_top(10), drop_top(20), drop_top(50))
say("  single best gene in the set: %s (rho = %.3f); median gene alone: %.3f",
    names(ord)[1], ord[1], median(st))

## ---- B. expression-matched background ---------------------------------------------
say("\nB. 1,000 expression-matched random sets of the same size")
brk <- quantile(expr, seq(0, 1, length.out = 21)); brk[1] <- -Inf; brk[21] <- Inf
bin <- cut(expr, brk, labels = FALSE)
names(bin) <- rownames(G)
pool <- split(setdiff(rownames(G), CORE), bin[setdiff(rownames(G), CORE)])
tab <- table(bin[cg])
draw <- function() unlist(lapply(names(tab), function(b) sample(pool[[b]], tab[[b]], replace = FALSE)))
nullmed <- replicate(1000, median(rho[draw()]))
nullset <- replicate(1000, cor(d$rate, colMeans(z[draw(), , drop = FALSE]), method = "spearman"))
say("  median per-gene rho: set %.3f | matched sets %.3f (95%% range %.3f to %.3f), empirical P %s",
    median(st), median(nullmed), quantile(nullmed, .025), quantile(nullmed, .975),
    ifelse(mean(nullmed >= median(st)) == 0, "< 0.001", sprintf("= %.3f", mean(nullmed >= median(st)))))
say("  set score:           set %.3f | matched sets max %.3f, 95%% range %.3f to %.3f, empirical P %s",
    cor(d$rate, d$splice, method = "spearman"), max(nullset), quantile(nullset, .025),
    quantile(nullset, .975),
    ifelse(mean(nullset >= cor(d$rate, d$splice, method = "spearman")) == 0, "< 0.001",
           sprintf("= %.3f", mean(nullset >= cor(d$rate, d$splice, method = "spearman")))))

## ---- C. the rate read three ways ---------------------------------------------------
say("\nC. the counted rate read three ways")
say("  divisions per day        rho = %+.3f  (P = %.1e)", cor(d$rate, d$splice, method = "spearman"),
    cor.test(d$rate, d$splice, method = "spearman", exact = FALSE)$p.value)
say("  hours per doubling       rho = %+.3f  (the study's own unit; the same measurement inverted)",
    cor(d$pdt, d$splice, method = "spearman"))
say("  cumulative doublings     rho = %+.3f  (replicative age, not speed)",
    cor(d$pdtot, d$splice, method = "spearman"))
dec <- cut(d$rate, quantile(d$rate, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
DEC <- do.call(rbind, lapply(sort(unique(dec)), function(k) data.frame(
  decile = k, n = sum(dec == k), rate_median = median(d$rate[dec == k]),
  splice_median = median(d$splice[dec == k]),
  lo = quantile(d$splice[dec == k], .25), hi = quantile(d$splice[dec == k], .75))))
print(DEC, digits = 3, row.names = FALSE)
say("  monotone across deciles: %s (Spearman of decile medians = %.3f)",
    ifelse(all(diff(DEC$splice_median) > -0.05), "yes", "no"),
    cor(DEC$decile, DEC$splice_median, method = "spearman"))
## at fixed replicative age and time in culture: is speed still what matters?
m0 <- lm(splice ~ rate + pdtot + days + line, d)
say("  partial: rate with cumulative doublings, days and cell line in the model: beta = %+.3f (P = %.1e)",
    coef(summary(m0))["rate", 1], coef(summary(m0))["rate", 4])
rp <- function(v, w, cv) { rx <- resid(lm(v ~ ., data.frame(v, cv))); ry <- resid(lm(w ~ ., data.frame(w, cv)))
                           cor(rx, ry, method = "spearman") }
cv <- data.frame(pdtot = d$pdtot, days = d$days, line = factor(d$line))
say("  partial rho(rate, splicing | cumulative doublings, days, cell line) = %+.3f", rp(d$rate, d$splice, cv))

## ---- D. named factors, per gene -----------------------------------------------------
NM <- GL[GL$gene %in% NAMED, c("gene", "rho", "mean_logcpm")]
NM <- NM[order(-NM$rho), ]
say("\nD. the twenty named factors, individually: median rho = %.3f, %d of %d above the background 95th pct",
    median(NM$rho), sum(NM$rho > q95), nrow(NM))

## ---- deposit -------------------------------------------------------------------------
write.table(GL, file.path(O, "gene_vs_counted_rate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DEC, file.path(O, "rate_decile_dose_response.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(d = d, expr_named = lc[intersect(NAMED, rownames(lc)), , drop = FALSE],
             null_median = nullmed, null_setscore = nullset, q95 = q95, q99 = q99),
        file.path(O, "gene_level_inputs.rds"))
write.table(data.frame(
  quantity = c("background median rho", "background 95th pct", "background 99th pct",
               "set median rho", "set above 95th pct (%)", "set above 99th pct (%)",
               "set score vs rate", "set score, top 10 genes removed",
               "hours per doubling", "cumulative doublings",
               "partial rho given pdtot, days, line", "named factors median rho"),
  value = c(median(bg), q95, q99, median(st), 100 * mean(st > q95), 100 * mean(st > q99),
            cor(d$rate, d$splice, method = "spearman"), drop_top(10),
            cor(d$pdt, d$splice, method = "spearman"), cor(d$pdtot, d$splice, method = "spearman"),
            rp(d$rate, d$splice, cv), median(NM$rho))),
  file.path(O, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\nwritten to %s", O)
