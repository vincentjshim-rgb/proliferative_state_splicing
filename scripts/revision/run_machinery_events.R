## What kind of splicing events move with proliferation in the machinery genes, and do the
## alternative exons belong to productive or to nonsense-mediated-decay isoforms?
##
## run_machinery_transcript_level.R tests every junction-pair event against the proliferation
## score but says nothing about what the event is. Here each event is classified from the
## geometry of its junctions against GENCODE v41: a cassette exon (one junction lands on an exon
## that the other junction skips), an alternative 5' or 3' splice site (both junctions land in the
## same exon), or something else. For every alternative element the transcript types of the
## GENCODE transcripts that carry it are recorded, so an exon that exists only in
## nonsense_mediated_decay transcripts -- the "poison exon" through which SR proteins and hnRNPs
## regulate their own abundance -- is flagged as such.
##
## The event tables group junctions by their shared genomic coordinate ("alternative acceptor" =
## shared start, "alternative donor" = shared end), which is a strand-blind label; the biological
## site is worked out from the strand here.
##
## post hoc. Inputs: machinery_transcript/psi_vs_proliferation.tsv, psi_age/psi_matrix.rds,
## machinery_transcript/gencode_v41_exons_event_genes.tsv (exons of every gene with an event,
## extracted from the GENCODE v41 GTF; the GTF itself is not in the repository).
## Output: public_data_tierA/derived/machinery_transcript/events_*.tsv, named_events_psi.rds
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; O <- file.path(D, "machinery_transcript")
say <- function(...) cat(sprintf(...), "\n")
NMD_TYPES <- c("nonsense_mediated_decay", "non_stop_decay")

R  <- fread(file.path(O, "psi_vs_proliferation.tsv"))
R[, row := .I]
X  <- readRDS(file.path(D, "psi_age/psi_matrix.rds"))
PSI <- X$PSI; M <- as.data.frame(X$M)
stopifnot(nrow(R) == nrow(PSI), all(R$start == X$ANN$start), all(R$end == X$ANN$end))
DON <- gse113957_donors()
keep <- DON$primary[match(M$srr, DON$srr)]; keep[is.na(keep)] <- FALSE
Mk <- M[keep, ]; PSIk <- PSI[, keep, drop = FALSE]
stopifnot(ncol(PSIk) == 107)
lo <- Mk$prolif <= quantile(Mk$prolif, 1/3); hi <- Mk$prolif >= quantile(Mk$prolif, 2/3)
EX <- fread(file.path(O, "gencode_v41_exons_event_genes.tsv"))
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")

## ---- groups: the junctions that share a site --------------------------------------------
R[, grp := ifelse(type == "alternative acceptor", paste(chromosome, strand, start, "S", sep = ":"),
                                                   paste(chromosome, strand, end, "E", sep = ":"))]
R[, shared := ifelse(type == "alternative acceptor", "start", "end")]

## transcript types of the GENCODE transcripts that carry an exon with exactly these bounds
types_of <- function(cc, s, e) sort(unique(EX[chr == cc & start == s & end == e]$transcript_type))
type_label <- function(tt) {
  if (!length(tt)) return("not an annotated exon")
  if (all(tt %in% NMD_TYPES)) return("NMD-only exon")
  if ("protein_coding" %in% tt) return("coding exon")
  "non-coding annotation only" }

## classify one pair of junctions sharing a site; `inc` is the junction that uses the proximal
## site (it includes the alternative element), `skp` the one that uses the distal site
classify_pair <- function(inc, skp) {
  cc <- inc$chromosome; strand <- inc$strand; g <- inc$gene
  exg <- EX[chr == cc & gene == g]
  if (inc$shared == "start") {                           # same start; ends differ
    a1 <- inc$end + 1L; a2 <- skp$end + 1L               # exon starts at the two acceptor sides
    e1 <- exg[start == a1]; e2 <- exg[start == a2]
    cas <- e1[end < skp$end]
    if (nrow(cas) && nrow(e2)) {
      ## the exon the reads use is the one whose downstream junction is in the event table
      cas[, has_jn := mapply(function(s, e) any(R$chromosome == cc & R$start == e + 1L & R$end == skp$end),
                             start, end)]
      cas[, coding := transcript_type == "protein_coding"]
      cas <- cas[order(-has_jn, -coding, -(end - start))][1]
      return(list(type = "cassette exon", el_s = cas$start, el_e = cas$end,
                  tt = types_of(cc, cas$start, cas$end)))
    }
    if (nrow(e1) && nrow(e2) && any(outer(e1$end, e2$end, `==`))) {
      en <- e1$end[e1$end %in% e2$end][1]
      return(list(type = if (strand == "+") "alternative 3' splice site" else "alternative 5' splice site",
                  el_s = a1, el_e = a2 - 1L, tt = types_of(cc, a1, en)))
    }
    if (!nrow(e1) || !nrow(e2)) return(list(type = "unannotated splice site", el_s = a1, el_e = a2 - 1L, tt = character()))
    return(list(type = "other", el_s = a1, el_e = a2 - 1L, tt = character()))
  } else {                                               # same end; starts differ
    d1 <- inc$start - 1L; d2 <- skp$start - 1L           # exon ends at the two donor sides
    e1 <- exg[end == d1]; e2 <- exg[end == d2]
    cas <- e1[start > skp$start]
    if (nrow(cas) && nrow(e2)) {
      cas[, has_jn := mapply(function(s, e) any(R$chromosome == cc & R$start == skp$start & R$end == s - 1L),
                             start, end)]
      cas[, coding := transcript_type == "protein_coding"]
      cas <- cas[order(-has_jn, -coding, -(end - start))][1]
      return(list(type = "cassette exon", el_s = cas$start, el_e = cas$end,
                  tt = types_of(cc, cas$start, cas$end)))
    }
    if (nrow(e1) && nrow(e2) && any(outer(e1$start, e2$start, `==`))) {
      st <- e1$start[e1$start %in% e2$start][1]
      return(list(type = if (strand == "+") "alternative 5' splice site" else "alternative 3' splice site",
                  el_s = d2 + 1L, el_e = d1, tt = types_of(cc, st, d1)))
    }
    if (!nrow(e1) || !nrow(e2)) return(list(type = "unannotated splice site", el_s = d2 + 1L, el_e = d1, tt = character()))
    return(list(type = "other", el_s = d2 + 1L, el_e = d1, tt = character()))
  }
}

## every non-dominant junction of a group is one event against the dominant junction; the
## junction with the proximal site is the "inclusion" junction whatever its share
events_of <- function(G) {
  ref <- G[which.max(ifelse(is.finite(meanPSI), meanPSI, -1))]
  rbindlist(lapply(setdiff(G$row, ref$row), function(r) {
    oth <- G[row == r]
    if (G$shared[1] == "start") { inc <- if (oth$end < ref$end) oth else ref; skp <- if (oth$end < ref$end) ref else oth }
    else                        { inc <- if (oth$start > ref$start) oth else ref; skp <- if (oth$start > ref$start) ref else oth }
    cl <- classify_pair(inc, skp)
    ## the event is measured by the junction that departs from the dominant path, so that a
    ## group of three junctions gives two events with two different measurements; `path` says
    ## whether that junction includes the element or skips it
    path <- if (identical(oth$row, inc$row)) "inclusion" else "skipping"
    data.table(gene = inc$gene, strand = inc$strand, chromosome = inc$chromosome, grp = inc$grp,
               shared = inc$shared, event_type = cl$type, element_start = cl$el_s, element_end = cl$el_e,
               element_length = cl$el_e - cl$el_s + 1L,
               element_transcript_types = paste(cl$tt, collapse = ";"), element_class = type_label(cl$tt),
               inclusion_junction = sprintf("%s:%d-%d", inc$chromosome, inc$start, inc$end),
               skipping_junction  = sprintf("%s:%d-%d", skp$chromosome, skp$start, skp$end),
               path = path, rep_row = oth$row, annotated = min(inc$annotated, skp$annotated),
               meanPSI = oth$meanPSI, d_prolif = oth$d_prolif, t_prolif = oth$t_prolif,
               p_prolif = oth$p_prolif, FDR_prolif = oth$FDR_prolif,
               d_age = oth$d_age, FDR_age = oth$FDR_age, d_age_adj = oth$d_age_adj, FDR_age_adj = oth$FDR_age_adj)
  }))
}

say("classifying %s junction rows in %d groups", format(nrow(R), big.mark = ","), uniqueN(R$grp))
EVall <- rbindlist(lapply(split(R[!is.na(gene) & gene != ""], by = "grp"), events_of))
EVall[, machinery := gene %in% readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))]
EVall[, named := gene %in% NAMED]
EVall[, sig := is.finite(FDR_prolif) & FDR_prolif < 0.05]
## PSI of the inclusion junction in the slowest and fastest thirds of the donors
EVall[, PSI_slow := rowMeans(PSIk[rep_row, lo, drop = FALSE], na.rm = TRUE)]
EVall[, PSI_fast := rowMeans(PSIk[rep_row, hi, drop = FALSE], na.rm = TRUE)]
EVall[, dPSI_pp := 100 * (PSI_fast - PSI_slow)]
## the effect expressed as inclusion of the alternative element, whichever path was measured
EVall[, d_inclusion := ifelse(path == "inclusion", d_prolif, -d_prolif)]
EVall[, d_inclusion_pp := ifelse(path == "inclusion", dPSI_pp, -dPSI_pp)]
## a cassette exon is seen from both of its junctions; the two sides are one event
EVall[, event_key := ifelse(event_type == "cassette exon", sprintf("%s|%d-%d", gene, element_start, element_end),
                      sprintf("%s|%s|%s", gene, grp, inclusion_junction))]
EVall[, n_sides := .N, by = event_key]
EV <- EVall[order(event_key, p_prolif)][, .SD[1], by = event_key]      # the better-supported side represents it
EV[, both_sides_agree := n_sides > 1]
setcolorder(EV, c("gene", "named", "machinery"))

fwrite(EV[order(-machinery, gene, p_prolif)], file.path(O, "events_classified_all_genes.tsv"), sep = "\t")
ME <- EV[machinery == TRUE][order(gene, p_prolif)]
fwrite(ME, file.path(O, "events_classified_machinery.tsv"), sep = "\t")

## ---- what the machinery events are -----------------------------------------------------
say("\nmachinery genes: %d events in %d genes; %d events in %d genes at FDR < 0.05",
    nrow(ME), uniqueN(ME$gene), sum(ME$sig), uniqueN(ME[sig == TRUE]$gene))
MS <- ME[sig == TRUE]
say("  by type: %s", paste(sprintf("%s %d", names(table(MS$event_type)), table(MS$event_type)), collapse = "; "))
say("  alternative element: %s", paste(sprintf("%s %d", names(table(MS$element_class)), table(MS$element_class)), collapse = "; "))
say("  direction (inclusion of the element with proliferation): %d down, %d up", sum(MS$d_inclusion < 0), sum(MS$d_inclusion > 0))
say("  NMD-only elements: inclusion falls with proliferation in %d of %d",
    sum(MS$element_class == "NMD-only exon" & MS$d_inclusion < 0), sum(MS$element_class == "NMD-only exon"))
TAB <- MS[, .(gene, named, event_type, element = sprintf("%d-%d (%d bp)", element_start, element_end, element_length),
              element_class, path, PSI_slow = round(PSI_slow, 3), PSI_fast = round(PSI_fast, 3),
              d_incl_pp = round(d_inclusion_pp, 1), FDR = signif(FDR_prolif, 2), both_sides_agree)]
fwrite(TAB, file.path(O, "events_machinery_significant_table.tsv"), sep = "\t")

## ---- the same events in every other gene, for the composition comparison ----------------
BS <- EV[machinery == FALSE & sig == TRUE]
comp <- function(dt, label) data.table(
  genes = label, n_events = nrow(dt),
  cassette = sum(dt$event_type == "cassette exon"),
  alt5 = sum(dt$event_type == "alternative 5' splice site"),
  alt3 = sum(dt$event_type == "alternative 3' splice site"),
  other = sum(dt$event_type %in% c("other", "unannotated splice site")),
  nmd_only = sum(dt$element_class == "NMD-only exon"),
  nmd_only_down = sum(dt$element_class == "NMD-only exon" & dt$d_inclusion < 0),
  coding = sum(dt$element_class == "coding exon"))
SUMM <- rbind(comp(MS, "177 machinery genes"), comp(BS, "all other genes"),
              comp(ME, "177 machinery genes, all events"), comp(EV[machinery == FALSE], "all other genes, all events"))
print(SUMM)
fwrite(SUMM, file.path(O, "event_type_summary.tsv"), sep = "\t")
ft <- fisher.test(matrix(c(sum(MS$element_class == "NMD-only exon"), nrow(MS) - sum(MS$element_class == "NMD-only exon"),
                           sum(BS$element_class == "NMD-only exon"), nrow(BS) - sum(BS$element_class == "NMD-only exon")), 2))
say("  NMD-only share among significant events: machinery %.0f%% vs other genes %.0f%% (Fisher P = %.3g)",
    100 * mean(MS$element_class == "NMD-only exon"), 100 * mean(BS$element_class == "NMD-only exon"), ft$p.value)

## ---- per-donor inclusion for the named factors' events, for the small multiples ---------
NE <- MS[named == TRUE][order(p_prolif)]
PD <- rbindlist(lapply(seq_len(nrow(NE)), function(i) data.table(
  gene = NE$gene[i], event_key = NE$event_key[i], event_type = NE$event_type[i],
  element_class = NE$element_class[i], element_length = NE$element_length[i], path = NE$path[i],
  FDR = NE$FDR_prolif[i], prolif = Mk$prolif, age = Mk$age, psi = PSIk[NE$rep_row[i], ])))
saveRDS(list(events = NE, psi = PD), file.path(O, "named_events_psi.rds"))
say("\nnamed factors with a significant event: %s", paste(unique(NE$gene), collapse = ", "))
say("written to %s", O)
