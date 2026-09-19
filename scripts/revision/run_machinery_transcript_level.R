## The splicing machinery's own transcripts: does the junction-level behaviour of these genes
## move with proliferation the way their abundance does?
##
## Everything reported so far is transcript ABUNDANCE summarised as a set score. The machinery
## also regulates itself at the level of splicing: SR-protein genes carry alternative exons whose
## inclusion sends the transcript to nonsense-mediated decay, so that more protein means more
## inclusion. If abundance is a readout of proliferative state, then in the same donors the
## inclusion of those exons should move with proliferation too. That is a different measurement
## on the same genes -- junction counts rather than gene counts -- and it is tested here.
##
## post hoc: this analysis was designed after the main results were known.
## Output: public_data_tierA/derived/machinery_transcript/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(limma); library(Matrix)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "machinery_transcript")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

CORE  <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # the 177
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")

## ---- events, with the gene each was mapped to ------------------------------------
X <- readRDS(file.path(D, "psi_age/psi_matrix.rds"))
PSI <- X$PSI; ANN <- as.data.table(X$ANN); M <- as.data.frame(X$M)
EV <- fread(file.path(D, "psi_age/psi_age_events.tsv"))
stopifnot(nrow(EV) == nrow(PSI), all(EV$start == ANN$start), all(EV$end == ANN$end))
ANN[, gene := EV$gene]
DON <- gse113957_donors()
M <- cbind(M, DON[match(M$srr, DON$srr), c("disease", "repo", "instr", "sex")])
M$primary <- DON$primary[match(M$srr, DON$srr)]
keep <- M$primary & !is.na(M$primary)
say("events: %s | donors in the primary cohort: %d of %d",
    format(nrow(PSI), big.mark = ","), sum(keep), nrow(M))

lgt <- function(x) { x <- pmin(pmax(x, 0.005), 0.995); log(x / (1 - x)) }
lg <- lgt(PSI[, keep]); Mk <- M[keep, ]

## ---- PSI against proliferation, and against age ----------------------------------
fit <- function(form, coef = 2) {
  tt <- topTable(eBayes(lmFit(lg, model.matrix(form, data = Mk))), coef = coef,
                 number = Inf, sort.by = "none")
  tt }
P  <- fit(~ prolif + log_depth + repo + instr + sex)
A0 <- fit(~ I(age/10) + log_depth + repo + instr + sex)
A1 <- fit(~ I(age/10) + prolif + log_depth + repo + instr + sex)

R <- data.table(gene = ANN$gene, chromosome = ANN$chromosome, start = ANN$start, end = ANN$end,
                strand = ANN$strand, annotated = ANN$annotated, type = ANN$type,
                meanPSI = rowMeans(PSI[, keep]),
                d_prolif = P$logFC, t_prolif = P$t, p_prolif = P$P.Value, FDR_prolif = P$adj.P.Val,
                d_age = A0$logFC, FDR_age = A0$adj.P.Val,
                d_age_adj = A1$logFC, FDR_age_adj = A1$adj.P.Val)
R[, machinery := !is.na(gene) & gene %in% CORE]
R[, named := !is.na(gene) & gene %in% NAMED]

say("\nPSI against the proliferation score (primary cohort, %d donors)", sum(keep))
say("  events at FDR < 0.05: %d of %s (%.1f%%)", sum(R$FDR_prolif < 0.05), format(nrow(R), big.mark = ","),
    100 * mean(R$FDR_prolif < 0.05))
mm <- R[machinery == TRUE]; bg <- R[machinery == FALSE & !is.na(gene) & gene != ""]
say("  events in the 177 machinery genes: %d in %d genes | other genes: %d events",
    nrow(mm), uniqueN(mm$gene), nrow(bg))
say("  FDR < 0.05: machinery %d of %d (%.0f%%) against background %d of %d (%.1f%%)",
    sum(mm$FDR_prolif < 0.05), nrow(mm), 100 * mean(mm$FDR_prolif < 0.05),
    sum(bg$FDR_prolif < 0.05), nrow(bg), 100 * mean(bg$FDR_prolif < 0.05))
say("  |t| median: machinery %.2f against background %.2f (Wilcoxon P = %.3g)",
    median(abs(mm$t_prolif)), median(abs(bg$t_prolif)),
    wilcox.test(abs(mm$t_prolif), abs(bg$t_prolif))$p.value)
say("  direction of the significant machinery events: %d up with proliferation, %d down",
    sum(mm$FDR_prolif < 0.05 & mm$d_prolif > 0), sum(mm$FDR_prolif < 0.05 & mm$d_prolif < 0))
say("  unannotated among the significant machinery events: %d of %d",
    sum(mm$FDR_prolif < 0.05 & mm$annotated == 0), sum(mm$FDR_prolif < 0.05))

say("\nthe same events against donor age")
say("  machinery events at FDR < 0.05 for age: %d; after proliferation is added: %d",
    sum(mm$FDR_age < 0.05), sum(mm$FDR_age < 0.05 & mm$FDR_age_adj < 0.05))

say("\nstrongest machinery events against proliferation:")
print(head(mm[order(p_prolif), .(gene, chromosome, start, end, type, annotated,
                                 meanPSI = round(meanPSI, 3), dPSI = round(d_prolif, 3),
                                 FDR = signif(FDR_prolif, 2))], 14))

## SR-protein genes specifically: the autoregulatory family
SR <- grep("^SRSF", unique(R$gene[!is.na(R$gene)]), value = TRUE)
sr <- R[gene %in% SR]
say("\nSR-protein genes with a testable event: %s", paste(sort(unique(sr$gene)), collapse = ", "))
if (nrow(sr)) print(sr[order(p_prolif), .(gene, start, end, type, annotated, meanPSI = round(meanPSI, 3),
                                          dPSI = round(d_prolif, 3), FDR = signif(FDR_prolif, 2))])

fwrite(R, file.path(O, "psi_vs_proliferation.tsv"), sep = "\t")

## ---- junction-level detail for the leading gene, for a sashimi-style panel ---------
## SRSF3 is the one drawn: of the named factors it gives the strongest event, and its two
## junctions from one donor site are complementary, so the pair reads as a single inclusion
## event -- the short exon at 36,599,820 that sits between exons 3 and 5
TARGET <- if ("SRSF3" %in% mm$gene) "SRSF3" else mm[order(p_prolif)][1]$gene
lead <- mm[gene == TARGET][order(p_prolif)][1]
if (!is.na(lead$gene)) {
  gb <- fread(file.path(D, "ir_annotation_v41/genes.bed"),
              col.names = c("chr", "s", "e", "gene", "score", "strand"))
  loc <- gb[gene == lead$gene][1]
  say("\njunction-level detail for %s (%s:%s-%s)", lead$gene, loc$chr,
      format(loc$s, big.mark = ","), format(loc$e, big.mark = ","))
  rc <- "public_data_tierA/recount3"; SRP <- "SRP144355"
  ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))), what = "",
              sep = "\n", quiet = TRUE)
  ids <- unlist(strsplit(ids, "[\t ]+")); ids <- ids[ids != "" & ids != "rail_id"]
  RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))))
  sel <- which(RR$chromosome == loc$chr & RR$start >= loc$s - 1000 & RR$end <= loc$e + 1000)
  say("  junctions in the locus: %d", length(sel))
  trip <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)),
                            "| tail -n +4"), col.names = c("i", "j", "v"))
  trip <- trip[i %in% sel]
  S <- sparseMatrix(i = match(trip$i, sel), j = trip$j, x = trip$v,
                    dims = c(length(sel), length(ids)))
  colnames(S) <- ids
  J <- as.matrix(S[, match(as.character(M$rail), colnames(S))])
  Jk <- J[, keep, drop = FALSE]
  cov <- rowSums(Jk >= 5)
  use <- cov >= 0.5 * ncol(Jk)
  say("  junctions covered in at least half the donors: %d", sum(use))
  hi <- Mk$prolif >= quantile(Mk$prolif, 2/3); lo <- Mk$prolif <= quantile(Mk$prolif, 1/3)
  JD <- data.table(RR[sel][use], mean_high_prolif = rowMeans(Jk[use, hi, drop = FALSE]),
                   mean_low_prolif  = rowMeans(Jk[use, lo, drop = FALSE]))
  JD[, ratio := (mean_high_prolif + 1) / (mean_low_prolif + 1)]
  print(JD[order(-mean_high_prolif), .(start, end, annotated, mean_high_prolif = round(mean_high_prolif),
                                       mean_low_prolif = round(mean_low_prolif), ratio = round(ratio, 2))])
  ex <- fread(file.path(D, "ir_annotation_v41/exons.bed"),
              col.names = c("chr", "s", "e", "gene", "score", "strand"))[gene == lead$gene]
  fwrite(JD, file.path(O, sprintf("junctions_%s.tsv", lead$gene)), sep = "\t")
  fwrite(ex, file.path(O, sprintf("exons_%s.tsv", lead$gene)), sep = "\t")
  saveRDS(list(J = Jk[use, , drop = FALSE], ann = RR[sel][use], meta = Mk, gene = lead$gene,
               locus = loc, exons = ex, psi = PSI[, keep][which(R$gene == lead$gene &
                 R$start == lead$start & R$end == lead$end)[1], ]),
          file.path(O, sprintf("junction_detail_%s.rds", lead$gene)))
  say("  written: junctions_%s.tsv, exons_%s.tsv, junction_detail.rds", lead$gene, lead$gene)
}
say("\nwritten to %s", O)
