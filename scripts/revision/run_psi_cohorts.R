## Event-level splicing (PSI) under the redefined cohort.
## Direction (a) takes splicing outcome and COL1A2 out of the core claims; what is
## reported instead is what the data support: the age-associated events of the
## published 142-donor cohort are not recovered once progeria donors and children
## are removed, and the COL1A2 event that led the published list is carried by the
## donors aged 83 and over, all of whom come from one cell repository.
## Outputs: public_data_tierA/derived/psi_cohorts/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(limma); library(data.table)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "psi_cohorts")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

X <- readRDS(file.path(D, "psi_age/psi_matrix.rds")); PSI <- X$PSI; ANN <- X$ANN
M <- as.data.frame(X$M)
DON <- gse113957_donors()
M <- cbind(M, DON[match(M$srr, DON$srr), c("disease", "repo", "instr", "sex", "cell_id")])
M$primary <- DON$primary[match(M$srr, DON$srr)]
stopifnot(ncol(PSI) == nrow(M))
say("PSI matrix: %d events x %d donors", nrow(PSI), ncol(PSI))

lgt <- function(x) { x <- pmin(pmax(x, 0.005), 0.995); log(x / (1 - x)) }
lg <- lgt(PSI)
res0 <- fread(file.path(D, "psi_age/psi_age_events.tsv"))
ev_col <- which(res0$gene == "COL1A2" & res0$FDR < 0.05)

run <- function(keep, extra = "") {
  f <- as.formula(paste("~ I(age/10) + log_depth", extra))
  f1 <- as.formula(paste("~ I(age/10) + prolif + log_depth", extra))
  a0 <- topTable(eBayes(lmFit(lg[, keep], model.matrix(f,  data = M[keep, ]))), coef = 2, number = Inf, sort.by = "none")
  a1 <- topTable(eBayes(lmFit(lg[, keep], model.matrix(f1, data = M[keep, ]))), coef = 2, number = Inf, sort.by = "none")
  list(a0 = a0, a1 = a1)
}
COH <- list(
  "published, all 142"          = list(rep(TRUE, nrow(M)), ""),
  "normal donors (133)"         = list(M$disease == "Normal", "+ repo + instr + sex"),
  "primary: normal adults 20+"  = list(M$primary, "+ repo + instr + sex"),
  "normal adults 20-82"         = list(M$primary & M$age < 83, "+ repo + instr + sex"))
SUM <- do.call(rbind, lapply(names(COH), function(k) {
  keep <- COH[[k]][[1]]; r <- run(keep, COH[[k]][[2]])
  sig <- r$a0$adj.P.Val < .05
  data.frame(cohort = k, n = sum(keep), age_events = sum(sig),
             kept_after_proliferation = sum(sig & r$a1$adj.P.Val < .05),
             col1a2_min_P = min(r$a0$P.Value[ev_col]), col1a2_max_P = max(r$a0$P.Value[ev_col]),
             col1a2_events_sig = sum(r$a0$adj.P.Val[ev_col] < .05)) }))
print(SUM, digits = 3, row.names = FALSE)
write.table(SUM, file.path(O, "psi_events_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## the COL1A2 event that led the published list, by age band and repository
top <- ev_col[which.min(res0$p[ev_col])]
M$psi_top <- PSI[top, ]
M$band <- cut(M$age, c(0, 20, 40, 60, 83, 100), right = FALSE)
tab <- with(M[M$disease == "Normal", ], tapply(psi_top, list(band, repo), median))
say("\nCOL1A2 event %s:%d-%d (%s): median PSI by age band and repository", res0$chromosome[top], res0$start[top], res0$end[top], res0$gene[top])
print(round(tab, 3))
write.table(data.frame(band = rownames(tab), tab), file.path(O, "col1a2_psi_by_band.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(M[, c("srr", "age", "sex", "disease", "repo", "instr", "cell_id", "prolif", "log_depth", "psi_top", "primary")],
       file.path(O, "col1a2_top_event_per_donor.tsv"), sep = "\t")

## age association of that event within and outside the 83+ stratum
for (lab in c("all normal", "normal, age < 83", "normal 83+")) {
  k <- switch(lab, "all normal" = M$disease == "Normal",
              "normal, age < 83" = M$disease == "Normal" & M$age < 83,
              "normal 83+" = M$disease == "Normal" & M$age >= 83)
  ct <- suppressWarnings(cor.test(M$age[k], M$psi_top[k], method = "spearman", exact = FALSE))
  say("  %-18s n = %3d  rho = %+.2f  P = %.2g", lab, sum(k), ct$estimate, ct$p.value)
}

## Reactome enrichment of the age-associated events in the primary cohort
r <- run(M$primary, "+ repo + instr + sex")
sig <- which(r$a0$adj.P.Val < .05)
genes <- unique(na.omit(res0$gene[sig]))
say("\nprimary cohort: %d age-associated events in %d genes", length(sig), length(genes))
writeLines(sort(genes), file.path(O, "primary_cohort_event_genes.txt"))
out <- data.frame(event = with(res0, sprintf("%s:%d-%d%s", chromosome, start, end, strand)),
                  type = res0$type, annotated = res0$annotated, gene = res0$gene,
                  beta_primary = r$a0$logFC, P_primary = r$a0$P.Value, FDR_primary = r$a0$adj.P.Val,
                  beta_primary_adj = r$a1$logFC, FDR_primary_adj = r$a1$adj.P.Val,
                  FDR_published = res0$FDR)
fwrite(out, file.path(O, "psi_age_events_primary.tsv"), sep = "\t")

## local read-level evidence for the COL1A2 junctions (STAR SJ.out.tab, GENCODE v41)
sj <- list.files(file.path(D, "local_bam_v41"), pattern = "SJ.out.tab$", full.names = TRUE)
if (length(sj)) {
  say("\nCOL1A2 junctions leaving the donor site chr7:94,406,304 in the four locally aligned libraries")
  J <- do.call(rbind, lapply(sj, function(f) {
    x <- read.delim(f, header = FALSE)
    names(x)[1:9] <- c("chr", "start", "end", "strand", "motif", "annot", "uniq", "multi", "overhang")
    x <- x[x$chr == "chr7" & x$start == 94406304, ]
    if (!nrow(x)) return(NULL)
    data.frame(library = sub("\\.SJ.*$", "", basename(f)), x[, c("start", "end", "annot", "uniq", "multi")]) }))
  J$intron_kb <- round((J$end - J$start) / 1000, 1)
  print(J[order(J$end, J$library), ], row.names = FALSE)
  write.table(J, file.path(O, "col1a2_local_junctions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
