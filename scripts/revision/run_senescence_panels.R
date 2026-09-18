## Which senescence panel reports senescence, and which reports the cell cycle?
##
## Fig. 6c showed that the Reactome senescence set scores higher in faster-dividing
## cultures because 82 of its 150 measured genes are cell-cycle genes, while SenMayo
## does not. This extends that comparison to the panels an intervention study would
## actually choose, so that the choice can be made from a table. Post hoc; definitions
## fixed before any output was read.
##
##  panels    Reactome Cellular Senescence; SenMayo; Fridman & Tainsky up and down;
##            CellAge split by annotated effect (induces / inhibits senescence);
##            Reactome Cell Cycle, Mitotic as the positive control
##  measures  genes measured, share of them in the Reactome cell-cycle union,
##            correlation of the panel score with the measured division rate before
##            and after removing those genes, median per-gene correlation, and the
##            donor-age effect before and after adjustment for proliferation
## Output: public_data_tierA/derived/senescence_panels/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"; P <- "public_data_tierA/senescence"
O <- file.path(D, "senescence_panels"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

gmt <- strsplit(system("unzip -p public_data_tierA/network/ReactomePathways.gmt.zip",
                       intern = TRUE), "\t")
SETS <- setNames(lapply(gmt, function(x) toupper(unique(x[-(1:2)]))), vapply(gmt, `[`, "", 1))
CCU <- unique(unlist(SETS[grep("^Cell Cycle|^Mitotic|^M Phase|^G1|^G2|^S Phase|^Synthesis of DNA",
                               names(SETS))]))          # the Reactome cell-cycle union
grp <- function(f) { x <- readLines(file.path(P, "panels", f)); toupper(x[!grepl("^#", x)][-1]) }
ca <- fread(file.path(P, "panels", "cellage3.tsv"))
setnames(ca, make.names(names(ca)))
PANELS <- list(
  `Reactome Cellular Senescence` = SETS[["Cellular Senescence"]],
  `SenMayo` = toupper(readLines(file.path(P, "senmayo/SAUL_SEN_MAYO.v2026.1.Hs.grp"))[-(1:2)]),
  `Fridman senescence, up` = grp("FRIDMAN_SENESCENCE_UP.grp"),
  `Fridman senescence, down` = grp("FRIDMAN_SENESCENCE_DN.grp"),
  `CellAge, induces senescence` = toupper(unique(ca[Senescence.Effect == "Induces", Gene.symbol])),
  `CellAge, inhibits senescence` = toupper(unique(ca[Senescence.Effect == "Inhibits", Gene.symbol])),
  `Reactome Cell Cycle, Mitotic (control)` = SETS[["Cell Cycle, Mitotic"]])
for (n in names(PANELS)) say("%-40s %4d genes", n, length(PANELS[[n]]))

sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value) }
sc <- function(z, g) { g <- intersect(rownames(z), g)
  if (length(g) < 5) return(rep(NA_real_, ncol(z))); colMeans(z[g, , drop = FALSE]) }

## ---- the counted libraries ---------------------------------------------------
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z; M <- as.data.table(IN$meta)
TAB <- rbindlist(lapply(names(PANELS), function(n) {
  g <- intersect(rownames(z), PANELS[[n]])
  gcc <- intersect(g, CCU); gno <- setdiff(g, CCU)
  s_all <- sc(z, g); s_no <- sc(z, gno)
  per <- apply(z[g, , drop = FALSE], 1, function(v) sp(v, M$rate)["rho"])
  data.table(panel = n, genes_total = length(PANELS[[n]]), genes_measured = length(g),
             cellcycle_genes = length(gcc), pct_cellcycle = 100 * length(gcc) / length(g),
             rho_rate = sp(s_all, M$rate)["rho"], p_rate = sp(s_all, M$rate)["p"],
             rho_rate_noCC = if (length(gno) >= 5) sp(s_no, M$rate)["rho"] else NA_real_,
             genes_after_removal = length(gno), median_gene_rho = median(per, na.rm = TRUE))
}))
say("\n-- 328 counted libraries --")
print(TAB[, .(panel, measured = genes_measured, pct_cc = round(pct_cellcycle),
              rho_rate = round(rho_rate, 2), rho_noCC = round(rho_rate_noCC, 2),
              gene_median = round(median_gene_rho, 2))])
say("across the seven panels, share of cell-cycle genes vs coupling to division rate: %s",
    sprintf("rho = %+.2f", sp(TAB$pct_cellcycle, TAB$rho_rate)["rho"]))

## ---- the donor cohort ---------------------------------------------------------
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))[, .(srr, log_depth)]
DON <- as.data.table(merge(DON, MET, by = "srr"))
d <- DON[gse113957_keep(DON, "primary")]
gs <- fread(cmd = "zcat public_data_tierA/recount3/sra.gene_sums.SRP144355.G026.gz | tail -n +2")
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
y <- DGEList(GM[, d$srr]); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); zz <- t(scale(t(cpm(y, log = TRUE))))
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
d[, prolif := sc(zz, PROLIF)]
AGE <- rbindlist(lapply(names(PANELS), function(n) {
  d[, s := sc(zz, PANELS[[n]])]
  m0 <- lm(s ~ I(age/10) + log_depth + repo + instr + sex, d)
  m1 <- update(m0, . ~ . + prolif)
  g <- function(m) { co <- summary(m)$coef; ci <- confint(m)
    c(b = co["I(age/10)", 1], lo = ci["I(age/10)", 1], hi = ci["I(age/10)", 2], p = co["I(age/10)", 4]) }
  a <- g(m0); b <- g(m1)
  data.table(panel = n, beta = a["b"], lo = a["lo"], hi = a["hi"], p = a["p"],
             beta_adj = b["b"], lo_adj = b["lo"], hi_adj = b["hi"], p_adj = b["p"],
             pct_lost = 100 * (1 - abs(b["b"]) / abs(a["b"]))) }))
AGE[, `:=`(FDR = p.adjust(p, "BH"), FDR_adj = p.adjust(p_adj, "BH"))]
say("\n-- 107 adult donors: age effect per decade --")
print(AGE[, .(panel, beta = round(beta, 3), FDR = signif(FDR, 2),
              beta_adj = round(beta_adj, 3), FDR_adj = signif(FDR_adj, 2),
              pct_lost = round(pct_lost))])
OUT <- merge(TAB, AGE, by = "panel", sort = FALSE)
fwrite(OUT, file.path(O, "senescence_panel_audit.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
