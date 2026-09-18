## Sensitivity analyses the reviewers asked for that were still open after the
## first pass of the revision.
##   A. technical QC covariates for the junction-outcome index (reviewer 1, major 3)
##   B. a core-spliceosome subset of the 177-gene set (reviewer 1, major 6)
##   C. a negative control for the proliferation adjustment: does adjustment shrink
##      the age effect of ANY expression-matched gene set? (reviewer 2, major 3)
##   D. recount3 gene sums are base-coverage sums, not read counts; read length
##      differs between samples here (51 vs 75 nt), so the scores are recomputed
##      from coverage divided by read length (reviewer 2, major 7)
## Outputs: public_data_tierA/derived/reviewer_sensitivities/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR); library(fgsea)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "reviewer_sensitivities")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")

## ---- donors, metrics, QC ------------------------------------------------------
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))
DON <- merge(DON, MET[, .(srr, depth, unannot_reads, unannot_junc, entropy, log_depth)], by = "srr")
md <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", SRP))))
md[, read_len := as.numeric(sub(".*average:([0-9.]+).*", "\\1", read_info))]
qc <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.recount_qc.%s.MD.gz", SRP))))
qc[, read_len := md$read_len[match(qc$rail_id, md$rail_id)]]
pick <- function(nm) { hit <- grep(nm, names(qc), fixed = TRUE, value = TRUE)[1]; qc[[hit]] }
QC <- data.table(srr = qc$external_id,
                 read_len = qc$read_len,
                 uniq_pct = pick("star.uniquely_mapped_reads_%"),
                 multi_pct = pick("star.%_of_reads_mapped_to_multiple_loci"),
                 intron_pct = pick("intron_sum_%"),
                 exon_pct = pick("exon_fc.unique_%"))
DON <- as.data.table(merge(DON, QC, by = "srr"))
say("QC joined for %d donors | read length: %s", nrow(DON),
    paste(names(table(DON$read_len)), table(DON$read_len), sep = " nt x ", collapse = ", "))
P <- DON[gse113957_keep(DON, "primary")]
say("primary cohort: %d donors", nrow(P))
for (v in c("uniq_pct", "multi_pct", "intron_pct", "exon_pct"))
  say("   %-11s median %.1f (range %.1f-%.1f) | correlation with donor age %+.2f",
      v, median(P[[v]]), min(P[[v]]), max(P[[v]]), cor(P[[v]], P$age, method = "spearman"))

## expression matrix for the primary cohort, as in run_cohort_core.R
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))

score_matrix <- function(X) {
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
  list(lc = lc, z = t(scale(t(lc[apply(lc, 1, sd) > 0, ]))))
}
S <- score_matrix(GM[, P$srr, drop = FALSE])
P[, machinery := colMeans(S$z[intersect(rownames(S$z), CORE), , drop = FALSE])]
P[, prolif    := colMeans(S$z[intersect(rownames(S$z), PROLIF), , drop = FALSE])]

fit1 <- function(d, metric, extra = "") {
  f0 <- as.formula(paste(metric, "~ I(age/10) + log_depth + repo + instr + sex", extra))
  f1 <- as.formula(paste(metric, "~ I(age/10) + log_depth + repo + instr + sex + prolif", extra))
  c0 <- summary(lm(f0, d))$coef["I(age/10)", ]; c1 <- summary(lm(f1, d))$coef["I(age/10)", ]
  data.frame(metric = metric, covariates = if (nzchar(extra)) "cohort + QC" else "cohort",
             beta = c0[1], p = c0[4], beta_adj = c1[1], p_adj = c1[4],
             pct_lost = 100 * (1 - c1[1] / c0[1]), row.names = NULL)
}

## ================= A. QC covariates for the outcome index =====================
say("\nA. technical QC covariates (uniquely mapped %%, multi-mapping %%, intron %%, exon assignment %%)")
QCX <- "+ uniq_pct + multi_pct + intron_pct + exon_pct"
A <- rbind(fit1(P, "machinery"), fit1(P, "machinery", QCX),
           fit1(P, "unannot_reads"), fit1(P, "unannot_reads", QCX),
           fit1(P, "entropy"), fit1(P, "entropy", QCX))
print(A, digits = 3, row.names = FALSE)
write.table(A, file.path(O, "qc_covariates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("   correlation of the outcome index with the QC measures: uniq %+.2f, multi %+.2f, intron %+.2f, exon %+.2f",
    cor(P$unannot_reads, P$uniq_pct, method = "spearman"), cor(P$unannot_reads, P$multi_pct, method = "spearman"),
    cor(P$unannot_reads, P$intron_pct, method = "spearman"), cor(P$unannot_reads, P$exon_pct, method = "spearman"))

## ================= B. core spliceosome subset ==================================
td <- tempdir(); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
RE <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
major <- unique(unlist(RE[grep("^mRNA Splicing - Major Pathway", names(RE))]))
minor <- unique(unlist(RE[grep("^mRNA Splicing - Minor Pathway", names(RE))]))
spliceo <- intersect(CORE, union(major, minor))
rest <- setdiff(CORE, spliceo)
say("\nB. core spliceosome subset: %d of %d genes (remainder %d: capping, 3'-end processing, export)",
    length(spliceo), length(CORE), length(rest))
P[, spliceo := colMeans(S$z[intersect(rownames(S$z), spliceo), , drop = FALSE])]
P[, rest    := colMeans(S$z[intersect(rownames(S$z), rest), , drop = FALSE])]
B <- rbind(fit1(P, "machinery"), fit1(P, "spliceo"), fit1(P, "rest"))
print(B, digits = 3, row.names = FALSE)
write.table(B, file.path(O, "core_spliceosome_subset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## and against the measured division rate in GSE179848
T <- read.delim(file.path(D, "division_rate/sample_division_rate.tsv"))
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z179 <- IN$z; rate <- IN$meta$rate
sc <- function(g) { gg <- intersect(rownames(z179), g); colMeans(z179[gg, , drop = FALSE]) }
for (nm in c("all 177", "core spliceosome", "remainder")) {
  g <- switch(nm, "all 177" = CORE, "core spliceosome" = spliceo, "remainder" = rest)
  s <- sc(g); ct <- cor.test(s, rate, method = "spearman", exact = FALSE)
  say("   %-17s division-rate rho = %+.2f (P = %.1g), %d genes measured", nm, ct$estimate, ct$p.value,
      length(intersect(rownames(z179), g)))
}
writeLines(sort(spliceo), file.path(O, "core_spliceosome_genes.txt"))

## ========== C. negative control for the proliferation adjustment ==============
say("\nC. negative control: the same adjustment applied to 1000 expression-matched random sets")
set.seed(20260916)
lev <- rowMeans(S$lc)
pool <- setdiff(names(lev), c(CORE, PROLIF, unique(unlist(RE[grep("^Cell Cycle|^M Phase|^DNA Replication", names(RE))]))))
bins <- cut(lev, quantile(lev, seq(0, 1, .05)), include.lowest = TRUE, labels = FALSE); names(bins) <- names(lev)
core_b <- table(bins[intersect(CORE, names(lev))])
draw <- function() unlist(lapply(names(core_b), function(b) {
  cand <- pool[bins[pool] == as.integer(b)]; sample(cand, min(core_b[[b]], length(cand))) }))
obs <- fit1(P, "machinery")
NL <- t(replicate(1000, {
  g <- draw(); P[, rnd := colMeans(S$z[intersect(rownames(S$z), g), , drop = FALSE])]
  r <- fit1(P, "rnd"); c(beta = r$beta, p = r$p, beta_adj = r$beta_adj, p_adj = r$p_adj, lost = r$pct_lost) }))
NL <- as.data.frame(NL)
sig <- is.finite(NL$p) & NL$p < 0.05
say("   random sets with an unadjusted age effect at P < 0.05: %d of 1000", sum(sig))
say("   of those, %% of the age effect removed by adjustment: median %.0f%% (95%% range %.0f-%.0f)",
    median(NL$lost[sig]), quantile(NL$lost[sig], .025), quantile(NL$lost[sig], .975))
say("   observed for the splicing set: %.0f%% removed (unadjusted P = %.1g, adjusted P = %.3f)",
    obs$pct_lost, obs$p, obs$p_adj)
say("   random sets keeping an adjusted effect at P < 0.05: %d of %d significant sets",
    sum(sig & NL$p_adj < 0.05), sum(sig))
say("   empirical P for an attenuation as large as observed: %.3f", mean(NL$lost[sig] >= obs$pct_lost))
write.table(NL, file.path(O, "adjustment_negative_control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ========== D. base coverage vs read counts ===================================
say("\nD. gene sums are base-coverage sums; dividing by read length gives read-equivalents")
rl <- setNames(DON$read_len, DON$srr)[P$srr]
GMr <- sweep(GM[, P$srr, drop = FALSE], 2, rl, "/")
S2 <- score_matrix(GMr)
P[, machinery_rl := colMeans(S2$z[intersect(rownames(S2$z), CORE), , drop = FALSE])]
P[, prolif_rl    := colMeans(S2$z[intersect(rownames(S2$z), PROLIF), , drop = FALSE])]
say("   agreement with the coverage-based scores: machinery r = %.4f, proliferation r = %.4f",
    cor(P$machinery, P$machinery_rl), cor(P$prolif, P$prolif_rl))
Pr <- copy(P); Pr[, `:=`(machinery = machinery_rl, prolif = prolif_rl)]
Dt <- rbind(cbind(fit1(P, "machinery"), scale_used = "base coverage"),
            cbind(fit1(Pr, "machinery"), scale_used = "read equivalents"))
print(Dt, digits = 3, row.names = FALSE)
write.table(Dt, file.path(O, "read_length_conversion.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(P, file.path(O, "primary_cohort_with_qc.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
