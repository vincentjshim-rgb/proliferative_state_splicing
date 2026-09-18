## Splicing-factor expression against DIRECTLY MEASURED division rate.
## The Cellular Lifespan Study counted cells at every passage, so growth rate is a
## measurement, not a transcriptional proxy.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
D <- "public_data_tierA/derived"; O <- file.path(D, "division_rate"); dir.create(O, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded

gm <- read.delim(file.path(D, "audit/all_geo_samples.tsv"), check.names = FALSE)
gm <- gm[gm$dataset == "GSE179848", ]
raw <- read.csv(gzfile("public_data_tierA/aging/GSE179848_raw_counts_cell_lifespan_RNAseq_data.csv.gz"),
                check.names = FALSE)
sy <- suppressMessages(unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = sub("\\..*$","",raw[[1]]),
        keytype = "REFSEQ", column = "SYMBOL", multiVals = "first")))
cnt <- as.matrix(raw[, -1]); storage.mode(cnt) <- "numeric"
ok <- !is.na(sy) & nzchar(sy); cnt <- rowsum(cnt[ok, ], toupper(sy[ok]))
gm$mn <- paste0("Sample_", gm$rnaseq_sampleid); gm <- gm[match(colnames(cnt), gm$mn), ]

y <- DGEList(cnt); y <- y[filterByExpr(y, group = factor(gm$treatments)), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
z  <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
cg <- intersect(rownames(z), CORE)

d <- data.frame(
  sample = colnames(lc), line = gm$cell_line, cond = gm$treatments,
  oxy = as.character(gm$percent_oxygen), days = as.numeric(gm$days_grown_udays),
  pdt = as.numeric(gm$population_doubling_time_uhours_per_division),
  div = as.numeric(gm$divisions_per_passage_udoublings),
  pdtot = as.numeric(gm$population_doublings_utotal_divisions),
  splice = colMeans(z[cg, , drop = FALSE]),
  MKI67 = lc["MKI67", ])
d$rate <- 24 / d$pdt                                   # divisions per day, measured
d <- d[is.finite(d$rate) & d$rate > 0 & is.finite(d$splice), ]
say("samples with a measured division rate: %d of 345", nrow(d))
say("division rate: median %.3f divisions/day (range %.3f-%.3f)",
    median(d$rate), min(d$rate), max(d$rate))

## ---- the central relationship ----------------------------------------------
ct <- cor.test(d$rate, d$splice, method = "spearman", exact = FALSE)
say("\nsplicing-factor expression vs MEASURED division rate: rho = %+.3f, P = %.2e",
    ct$estimate, ct$p.value)
say("  (for comparison, vs MKI67 transcript: rho = %+.3f)",
    cor(d$MKI67, d$splice, method = "spearman"))
ct2 <- cor.test(d$rate, d$MKI67, method = "spearman", exact = FALSE)
say("  measured rate vs MKI67 transcript: rho = %+.3f", ct2$estimate)

## ---- does every growth-slowing perturbation obey the same rule? ------------
say("\n=== each perturbation, measured rate and splicing-factor expression ===")
ctl <- d[d$cond == "Control" & d$oxy == "21", ]
per <- do.call(rbind, lapply(setdiff(unique(d$cond), "Control"), function(k) {
  s <- d[d$cond == k, ]; if (nrow(s) < 4) return(NULL)
  b <- d[d$cond == "Control" & d$line %in% s$line & d$oxy %in% s$oxy, ]
  if (nrow(b) < 4) b <- ctl
  data.frame(cond = k, n = nrow(s),
             d_rate   = median(s$rate) - median(b$rate),
             d_splice = median(s$splice) - median(b$splice)) }))
o2 <- d[d$oxy == "3", ]; o2b <- d[d$oxy == "21" & d$line %in% o2$line, ]
per <- rbind(per, data.frame(cond = "3% oxygen", n = nrow(o2),
             d_rate = median(o2$rate) - median(o2b$rate),
             d_splice = median(o2$splice) - median(o2b$splice)))
per <- per[order(per$d_rate), ]
print(per, digits = 3, row.names = FALSE)
ctp <- cor.test(per$d_rate, per$d_splice, method = "spearman", exact = FALSE)
say("across %d perturbations with unrelated mechanisms: rho = %+.3f, P = %.4f",
    nrow(per), ctp$estimate, ctp$p.value)
write.table(per, file.path(O, "perturbation_rate_vs_splicing.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- named splicing factors, individually ----------------------------------
## One sourced list is used in every figure: the twenty factors assayed in senescing
## primary human dermal fibroblasts by Holly et al. 2013 (Table 3), all reported there
## as age-associated in this cell type. See scripts/revision/named_splicing_factors_source.md
SF <- readLines("scripts/revision/named_splicing_factors.txt")
sf <- intersect(SF, rownames(lc))
say("\n=== named splicing factors: correlation with measured division rate (n = %d) ===", nrow(d))
R <- do.call(rbind, lapply(sf, function(g) {
  v <- lc[g, d$sample]
  data.frame(gene = g, rho_rate = cor(v, d$rate, method = "spearman"),
             rho_days = cor(v, d$days, method = "spearman"),
             in_core = g %in% CORE) }))
R <- R[order(-R$rho_rate), ]
print(R, digits = 2, row.names = FALSE)
say("median |rho| with division rate across %d named factors: %.2f", nrow(R), median(abs(R$rho_rate)))
write.table(R, file.path(O, "named_splicing_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(d, file.path(O, "sample_division_rate.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
