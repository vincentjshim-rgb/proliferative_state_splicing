#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 5, rebuilt around splicing OUTCOME.
#   A  senescence changes alternative splice-site usage        (GSE109700)
#   B  part of it survives a matched-PD quiescence control     (GSE93535 2x2)
#   C  but most of it is shared with proliferative arrest
#   D  what the local n=1 data can and cannot say
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2) })
root <- "public_data_tierA"; rc <- file.path(root, "recount3")
ps   <- file.path(root, "derived", "phase2_splicing")
out  <- file.path(root, "derived", "manuscript_figures_v4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) write.table(x, file.path(ps, f), sep = "\t", quote = FALSE, row.names = FALSE)

BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"; GREY <- "#9E9E9E"
INK <- "#222222"; MUTED <- "#5A5A5A"
th <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.3),
        plot.caption = element_text(colour = MUTED, size = 7.1, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8.2))

## ===================== shared helpers ======================================
load_j <- function(SRP, soft, grpfun) {
  L <- readLines(gzfile(file.path(root, "metadata", soft)))
  bl <- split(L, cumsum(grepl("^\\^SAMPLE", L)))[-1]
  gm <- do.call(rbind, lapply(bl, function(b) { g1 <- function(p) { z <- grep(p, b, value = TRUE)[1]
      if (is.na(z)) NA else sub(".*= ", "", z) }
    data.frame(title = g1("!Sample_title"), srx = sub(".*term=", "", g1("!Sample_relation = SRA"))) }))
  md <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", SRP))),
              select = c("rail_id", "experiment_acc"), data.table = FALSE)
  gm$rail_id <- md$rail_id[match(gm$srx, md$experiment_acc)]
  gm$group <- grpfun(gm$title)
  ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))),
              what = character(), quiet = TRUE); ids <- ids[ids != "rail_id"]
  RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))), data.table = FALSE)
  tr <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)), "| tail -n +4"),
              header = FALSE, data.table = FALSE)
  S <- matrix(0L, nrow = nrow(RR), ncol = length(ids), dimnames = list(NULL, ids))
  S[cbind(tr[[1]], tr[[2]])] <- as.integer(tr[[3]]); rm(tr); gc()
  S <- S[, match(as.character(gm$rail_id), colnames(S)), drop = FALSE]; colnames(S) <- gm$title
  list(S = S, RR = RR, gm = gm)
}
build <- function(S, RR, smp, mincov) {
  key <- do.call(paste, c(RR[c("chromosome", "start", "strand")], sep = "_"))
  nalt <- ave(rep(1, nrow(RR)), key, FUN = sum); sel <- which(nalt >= 2)
  Js <- S[sel, smp, drop = FALSE]; ks <- key[sel]
  tot <- apply(Js, 2, function(v) ave(v, ks, FUN = sum))
  ok <- rowMeans(tot >= mincov) == 1
  Js <- Js[ok, , drop = FALSE]; tot <- tot[ok, , drop = FALSE]; P <- Js / pmax(tot, 1)
  P[apply(P, 1, function(r) !(min(r) > 0.98)) & apply(P, 1, sd) > 0, , drop = FALSE]
}
nsig2 <- function(P, a, b) {
  L <- log((P + 0.01) / (1 - P + 0.01))
  g <- factor(ifelse(colnames(L) %in% b, "b", "a"), levels = c("a", "b"))
  tt <- topTable(eBayes(lmFit(L, model.matrix(~ g))), coef = 2, number = Inf, sort.by = "none")
  d <- rowMeans(P[, b, drop = FALSE]) - rowMeans(P[, a, drop = FALSE])
  c(n = nrow(L), sig = sum(tt$adj.P.Val < 0.05), sig10 = sum(tt$adj.P.Val < 0.05 & abs(d) >= 0.10))
}

## ===================== A : GSE109700 =======================================
cat("A: GSE109700\n")
G <- load_j("SRP131506", "GSE109700_family.soft.gz",
            function(t) ifelse(grepl("^Prolif", t), "proliferating",
                        ifelse(grepl("^Early", t), "early", "deep")))
set.seed(109700); tgt <- min(colSums(G$S)); Sd <- G$S
for (j in seq_len(ncol(G$S))) { p <- tgt / sum(G$S[, j]); if (p < 1) Sd[, j] <- rbinom(nrow(G$S), G$S[, j], p) }
PRO <- G$gm$title[G$gm$group == "proliferating"]; EAR <- G$gm$title[G$gm$group == "early"]
DEE <- G$gm$title[G$gm$group == "deep"]
A <- do.call(rbind, lapply(list(list("deep senescent\nvs proliferating", PRO, DEE),
                                list("early senescent\nvs proliferating", PRO, EAR),
                                list("deep vs early\n(both senescent)", EAR, DEE)), function(z) {
  P <- build(Sd, G$RR, c(z[[2]], z[[3]]), 20); r <- nsig2(P, z[[2]], z[[3]])
  data.frame(contrast = z[[1]], n = unname(r["n"]), sig = unname(r["sig"]), sig10 = unname(r["sig10"])) }))
## label-permutation null on the depth-matched early-vs-proliferating comparison
six <- c(PRO, EAR); P6 <- build(Sd, G$RR, six, 20); seen <- c(); nulA <- c(); truA <- NA
for (k in combn(6, 3, simplify = FALSE)) {
  ga <- six[k]; gb <- six[-k]; tg <- paste(sort(ga), collapse = "|")
  if (tg %in% seen) next
  seen <- c(seen, tg, paste(sort(gb), collapse = "|"))
  v <- nsig2(P6, ga, gb)[["sig"]]
  if (setequal(ga, PRO) || setequal(ga, EAR)) truA <- v else nulA <- c(nulA, v)
}
w(A, "Figure5A_GSE109700_contrasts.tsv")
A$contrast <- factor(A$contrast, levels = rev(A$contrast))
pA <- ggplot(A, aes(sig, contrast)) +
  geom_segment(aes(x = 0, xend = sig, yend = contrast), colour = "grey82", linewidth = 1) +
  geom_point(colour = BLUE, size = 3) +
  geom_text(aes(label = sprintf("%d  (%d with |dPSI| >= 0.10)", sig, sig10)),
            hjust = -0.12, size = 2.7, colour = INK) +
  annotate("text", x = 40, y = 0.68, hjust = 0, size = 2.5, colour = MUTED,
           label = sprintf("label-permutation null (early vs proliferating):\ntrue split %d  |  all 9 other 3-vs-3 splits = %s",
                           truA, paste(unique(nulA), collapse = "/"))) +
  scale_x_continuous(limits = c(0, 1450), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.45, 3.5), clip = "off") +
  labs(x = "Alternative junctions with FDR < 0.05", y = NULL,
       title = "A. Senescence changes alternative splice-site usage",
       subtitle = "GSE109700, ribodepleted paired-end, all libraries downsampled to 19.5 M junction reads",
       caption = "Two senescent states barely differ from each other, so the signal is senescence, not batch. ~92-98k junctions testable per contrast.") +
  th
ggsave(file.path(out, "Figure5A_senescence_changes_splicing.png"), pA,
       width = 168, height = 84, units = "mm", dpi = 300)

## ===================== B : GSE93535 2x2 ====================================
cat("B: GSE93535\n")
d <- readRDS(file.path(ps, "H1_inputs.rds")); S0 <- d$S; RR2 <- d$RR; gm2 <- d$gmap
st <- ifelse(grepl("^SIPS", gm2$arm), "senescent", "quiescent")[match(colnames(S0), gm2$title)]
cp <- ifelse(grepl("1201$", gm2$arm), "c1201", "none")[match(colnames(S0), gm2$title)]
keep <- colSums(S0) / 1e6 >= 2
smp <- colnames(S0)[keep]; st2 <- st[keep]; cp2 <- cp[keep]
nsigST <- function(P, s, c) {
  L <- log((P + 0.01) / (1 - P + 0.01))
  dsn <- model.matrix(~ factor(s, levels = c("quiescent", "senescent")) + factor(c))
  if (qr(dsn)$rank < ncol(dsn)) return(NA)
  sum(topTable(eBayes(lmFit(L, dsn)), coef = 2, number = Inf, sort.by = "none")$adj.P.Val < 0.05)
}
PB <- build(S0, RR2, smp, 20); obsB <- nsigST(PB, st2, cp2)
set.seed(2); nulB <- replicate(300, { s <- st2
  for (lv in unique(cp2)) { i <- which(cp2 == lv); s[i] <- sample(s[i]) }
  if (identical(s, st2)) return(NA); nsigST(PB, s, cp2) })
nulB <- nulB[!is.na(nulB)]
pval <- (1 + sum(nulB >= obsB)) / (length(nulB) + 1)
w(data.frame(permuted_sig = nulB), "Figure5B_GSE93535_permutation_null.tsv")
pB <- ggplot(data.frame(x = nulB), aes(x)) +
  geom_histogram(binwidth = 1, fill = GREY, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obsB, colour = ORANGE, linewidth = 0.9) +
  annotate("text", x = obsB, y = Inf, vjust = 1.8, hjust = 1.12, size = 2.9, colour = ORANGE,
           label = sprintf("observed = %d", obsB)) +
  annotate("text", x = obsB, y = Inf, vjust = 3.6, hjust = 1.12, size = 2.6, colour = MUTED,
           label = sprintf("empirical p = %.4f", pval)) +
  labs(x = "Junctions with FDR < 0.05 (senescence main effect)", y = "permutations",
       title = "B. Part of it survives a matched-PD quiescence control",
       subtitle = "GSE93535, all samples at population doubling 15, 2x2 state x compound design",
       caption = paste0("Quiescent and senescent cells are both growth-arrested at the same PD. State labels permuted within compound strata, 298 draws;\n",
                        "the observed count exceeds every permutation, so the empirical p is at the permutation floor.")) +
  th
ggsave(file.path(out, "Figure5B_quiescence_control.png"), pB,
       width = 168, height = 80, units = "mm", dpi = 300)

## ===================== C : magnitude ========================================
Cdat <- data.frame(
  ref = c("vs proliferating\n(GSE109700)", "vs quiescent, matched PD\n(GSE93535)"),
  sig = c(A$sig[A$contrast == "deep senescent\nvs proliferating"], obsB),
  n   = c(A$n[A$contrast == "deep senescent\nvs proliferating"], nrow(PB)))
Cdat$pct <- 100 * Cdat$sig / Cdat$n
Cdat$ref <- factor(Cdat$ref, levels = rev(Cdat$ref))
w(Cdat, "Figure5C_magnitude.tsv")
pC <- ggplot(Cdat, aes(pct, ref)) +
  geom_segment(aes(x = 0, xend = pct, yend = ref), colour = "grey82", linewidth = 1) +
  geom_point(colour = ORANGE, size = 3) +
  geom_text(aes(label = sprintf("%.2f%%   (%d / %s)", pct, sig, format(n, big.mark = ","))),
            hjust = -0.12, size = 2.7, colour = INK) +
  scale_x_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
  labs(x = "Percent of testable junctions reaching FDR < 0.05", y = NULL,
       title = "C. Most of the change is shared with growth arrest",
       subtitle = "Senescence-specific component is roughly five-fold smaller than the total",
       caption = "Different studies, depths and n, so read the ratio as an order of magnitude rather than an exact effect size.") +
  th
ggsave(file.path(out, "Figure5C_magnitude.png"), pC, width = 168, height = 62, units = "mm", dpi = 300)

## ===================== D : local boundary ===================================
b <- read.delim(file.path(root, "derived", "intron_retention", "local_IR_burden.tsv"))
b$label <- factor(b$label, levels = c("UVA 15J, 0 h", "iPSC-CM", "HDF-CM", "Repro-CM"))
b$grp <- ifelse(b$sample == "UVA0", "0 h, no CM (confounded with time)", "24 h conditioned medium")
pD <- ggplot(b, aes(100 * intronic_read_fraction, label, fill = grp)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", 100 * intronic_read_fraction)),
            hjust = -0.25, size = 2.7, colour = INK) +
  scale_fill_manual(values = c(`24 h conditioned medium` = BLUE,
                               `0 h, no CM (confounded with time)` = GREY)) +
  scale_x_continuous(limits = c(0, 7.4), expand = c(0, 0)) +
  labs(x = "Intronic read fraction (%)", y = NULL,
       title = "D. The local data cannot answer this for Repro-CM",
       subtitle = "6,749 genes covered in all four libraries; unambiguous intronic vs exonic blocks, GENCODE v41",
       caption = paste0("DESCRIPTIVE ONLY: n = 1 per condition, no p-values. Repro-CM has the HIGHEST intron retention of the three 24 h conditions,\n",
                        "the opposite of a splicing-efficiency improvement, and per-gene IR tracks the expression anchor (rho = -0.25).\n",
                        "Junction-level dPSI candidates from these libraries carry an empirical FDR of 0.71-0.74. Answering this needs n >= 3 per condition.")) +
  th + theme(panel.grid.major.y = element_blank())
ggsave(file.path(out, "Figure5D_local_limitation.png"), pD, width = 168, height = 84, units = "mm", dpi = 300)
cat("Figure 5 A-D written\n")
