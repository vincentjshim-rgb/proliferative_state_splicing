## Is any of this specific to splicing?
##
## The splicing set follows the cell cycle, but so does the rest of the machinery of
## gene expression (Fig. 5: Metabolism of RNA 0.77, translation 0.58). This script
## asks the question directly, in the three places the paper measures it, and is
## post hoc.
##
## Revised 2026-09-17: selection asymmetry. The 177-gene set was SELECTED for falling
## in three series, one of them late versus early passage in GSE179848 (the resource
## in which division rate falls with passage) and one the GSE113957 donor cohort; the
## comparators were not selected. Every question is therefore run with three splicing
## definitions side by side, and with two lower bounds:
##   selected_177                  the 177-gene set (v2), as in the original analysis
##   reactome_mrna_splicing        Reactome "mRNA Splicing", whole, unselected
##   reactome_premrna_processing   Reactome "Processing of Capped Intron-Containing
##                                 Pre-mRNA", whole, unselected (the pathway the 177
##                                 were drawn from)
##   *_not_selected                the genes of each Reactome pathway that are NOT in
##                                 the 177: selected against, so a lower bound
## These definitions and the rules below were written down before any output of the
## revised script was read.
##
##  scores      mean gene-wise z of each set
##  shared genes  no comparison shares a gene. For the unselected definitions the
##              genes shared with a comparator are removed from BOTH sides. Two
##              exceptions, recorded in the column shared_rule: (i) Metabolism of RNA
##              is the Reactome parent of both splicing pathways, so removing shared
##              genes from both sides would leave nothing to score; the splicing genes
##              are removed from the comparator only. (ii) selected_177 is kept exactly
##              as in the original analysis (the 177 whole, its genes removed from the
##              comparator), so that the published numbers are reproduced. The other
##              rule is reported for Q1 in shared_rule_sensitivity.tsv.
##  Q1 (328 counted libraries) partial Spearman of each score with the measured
##              division rate, holding the other score constant
##  Q2 (107 donors) age effect per decade with the cohort's technical covariates,
##              then adding proliferation, then adding the other machinery score;
##              run in both directions
##  Q3 (511 GTEx cultures) the same partial correlations with the 20-marker
##              proliferation score, from gene-level scores recomputed with the
##              preregistered normalisation so that the same definitions and the same
##              shared-gene rule apply. NOTE: the original Q3 used the programme
##              scores stored in gtex_boundary/tissue_scores.rds ("pre-mRNA
##              processing" and "translation", shared genes not removed) and, as its
##              proliferation variable, the Reactome "Cell Cycle, Mitotic" programme
##              score rather than the 20-marker proliferation score. Those numbers are
##              reproduced in gtex_culture_partials_rds_programmes.tsv with the
##              variable named.
## Output: public_data_tierA/derived/machinery_specificity/
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
source("scripts/revision/cohort_gse113957.R")
D <- "public_data_tierA/derived"
O <- file.path(D, "machinery_specificity"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
NBOOT <- 2000

CORE <- toupper(readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt")))
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
gmt <- strsplit(system("unzip -p public_data_tierA/network/ReactomePathways.gmt.zip",
                       intern = TRUE), "\t")
SETS <- setNames(lapply(gmt, function(x) toupper(unique(x[-(1:2)]))), vapply(gmt, `[`, "", 1))
WANT <- c(translation = "Translation", metabolism_of_rna = "Metabolism of RNA",
          rrna_processing = "rRNA processing")
MRNASPL <- SETS[["mRNA Splicing"]]
PREMRNA <- SETS[["Processing of Capped Intron-Containing Pre-mRNA"]]
DEFS <- list(selected_177 = CORE,
             reactome_mrna_splicing = MRNASPL,
             reactome_premrna_processing = PREMRNA,
             mrna_splicing_not_selected = setdiff(MRNASPL, CORE),
             premrna_processing_not_selected = setdiff(PREMRNA, CORE))
DEF_LABEL <- c(selected_177 = "177-gene set, selected for falling in three series",
               reactome_mrna_splicing = "Reactome mRNA Splicing, whole, unselected",
               reactome_premrna_processing = "Reactome Processing of Capped Intron-Containing Pre-mRNA, whole, unselected",
               mrna_splicing_not_selected = "Reactome mRNA Splicing genes not in the 177 (selected against; lower bound)",
               premrna_processing_not_selected = "Reactome pre-mRNA processing genes not in the 177 (selected against; lower bound)")
SELECTION <- c(selected_177 = "selected", reactome_mrna_splicing = "unselected",
               reactome_premrna_processing = "unselected",
               mrna_splicing_not_selected = "selected against",
               premrna_processing_not_selected = "selected against")
for (k in names(DEFS)) say("%-32s %3d genes, %3d of them in the selected 177", k,
                           length(DEFS[[k]]), length(intersect(DEFS[[k]], CORE)))

## the two sides of one comparison, sharing no gene
pair_sets <- function(def, k, rule = NULL) {
  Dg <- DEFS[[def]]; Kg <- SETS[[WANT[[k]]]]; sh <- intersect(Dg, Kg)
  nested <- all(Dg %in% Kg)
  if (is.null(rule)) rule <- if (nested) "comparator only (definition nested in comparator)"
    else if (def == "selected_177") "comparator only (original analysis)" else "both sides"
  if (nested && rule == "both sides") return(NULL)
  list(D = if (rule == "both sides") setdiff(Dg, sh) else Dg, K = setdiff(Kg, sh),
       n_shared = length(sh), rule = rule)
}
GRID <- CJ(def = names(DEFS), k = names(WANT), sorted = FALSE)

sc <- function(z, g) colMeans(z[intersect(rownames(z), g), , drop = FALSE])
ng <- function(z, g) length(intersect(rownames(z), g))
pcor <- function(x, y, z) {   # partial Spearman via ranks
  rx <- residuals(lm(rank(x) ~ rank(z))); ry <- residuals(lm(rank(y) ~ rank(z)))
  ct <- suppressWarnings(cor.test(rx, ry, method = "pearson"))
  c(rho = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2], p = ct$p.value)
}
## the same, with technical covariates also held constant (ranks regressed on them)
pcor_cov <- function(x, y, z, cov) {
  ok <- complete.cases(cov); cv <- as.data.frame(cov)[ok, , drop = FALSE]
  rx <- residuals(lm(rank(x[ok]) ~ rank(z[ok]) + ., cv)); ry <- residuals(lm(rank(y[ok]) ~ rank(z[ok]) + ., cv))
  r <- cor(rx, ry); n <- sum(ok); df <- n - 3 - ncol(cv)
  t <- r * sqrt(df / (1 - r^2)); c(rho = r, p = 2 * pt(-abs(t), df), n = n)
}
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value) }
prho <- function(R, a, b, c) (R[a, b] - R[a, c] * R[b, c]) / sqrt((1 - R[a, c]^2) * (1 - R[b, c]^2))

## every definition against every comparator, for one expression matrix and one target
partial_table <- function(z, target, rule = NULL, cov = NULL, boot = TRUE) {
  rows <- list(); X <- list(target = target)
  for (i in seq_len(nrow(GRID))) {
    def <- GRID$def[i]; k <- GRID$k[i]; P <- pair_sets(def, k, rule); if (is.null(P)) next
    s <- sc(z, P$D); cmp <- sc(z, P$K); whole <- sc(z, DEFS[[def]])
    a <- sp(s, target); b <- sp(cmp, target); w <- sp(whole, target)
    pa <- pcor(s, target, cmp); pb <- pcor(cmp, target, s)
    r <- data.table(splicing_definition = def, selection = SELECTION[[def]], comparator = k,
      n = length(target), shared_rule = P$rule, n_shared_removed = P$n_shared,
      n_genes_splicing = ng(z, P$D), n_genes_comparator = ng(z, P$K),
      rho_splicing_whole_set = w["rho"],
      rho_splicing = a["rho"], rho_comparator = b["rho"],
      partial_splicing = pa["rho"], partial_splicing_lo = pa["lo"],
      partial_splicing_hi = pa["hi"], partial_splicing_p = pa["p"],
      partial_comparator = pb["rho"], partial_comparator_lo = pb["lo"],
      partial_comparator_hi = pb["hi"], partial_comparator_p = pb["p"],
      r_between = cor(s, cmp), asymmetry = pa["rho"] - pb["rho"])
    if (!is.null(cov)) { ca <- pcor_cov(s, target, cmp, cov); cb <- pcor_cov(cmp, target, s, cov)
      r[, `:=`(n_tech = ca["n"], partial_splicing_tech = ca["rho"], partial_splicing_tech_p = ca["p"],
               partial_comparator_tech = cb["rho"], partial_comparator_tech_p = cb["p"])] }
    rows[[i]] <- r; X[[paste(def, k, "s")]] <- s; X[[paste(def, k, "c")]] <- cmp
  }
  T <- rbindlist(rows)
  if (boot) {   # percentile bootstrap over libraries: the asymmetry, and its shrinkage from selected_177
    X <- do.call(cbind, X); set.seed(1); n <- nrow(X)
    A <- t(replicate(NBOOT, { R <- cor(apply(X[sample(n, replace = TRUE), ], 2, rank))
      vapply(seq_len(nrow(T)), function(i) { s <- paste(T$splicing_definition[i], T$comparator[i], "s")
        cc <- paste(T$splicing_definition[i], T$comparator[i], "c")
        prho(R, s, "target", cc) - prho(R, cc, "target", s) }, 0) }))
    q <- function(v) quantile(v, c(0.025, 0.975))
    T[, c("asymmetry_lo", "asymmetry_hi") := as.data.table(t(apply(A, 2, q)))]
    ref <- match(paste("selected_177", T$comparator), paste(T$splicing_definition, T$comparator))
    T[, shrink_from_selected := asymmetry[ref] - asymmetry]
    ## a percentage of the selected set's asymmetry means something only where that set shows one
    T[, pct_of_selected_asymmetry := ifelse(asymmetry[ref] > 0, 100 * asymmetry / asymmetry[ref], NA_real_)]
    SH <- A[, ref, drop = FALSE] - A
    T[, c("shrink_lo", "shrink_hi") := as.data.table(t(apply(SH, 2, q)))]
    T[splicing_definition == "selected_177", c("shrink_lo", "shrink_hi") := NA_real_]
  }
  T[]
}
show <- function(T) print(T[, .(splicing_definition, comparator, genes = paste(n_genes_splicing, n_genes_comparator, sep = "/"),
  rho_spl = round(rho_splicing, 2), rho_cmp = round(rho_comparator, 2),
  `spl|cmp` = round(partial_splicing, 2), `cmp|spl` = round(partial_comparator, 2),
  asym = round(asymmetry, 2), pct_sel = round(pct_of_selected_asymmetry))], nrows = 40)

## ---- Q1: the counted libraries ----------------------------------------------
IN <- readRDS(file.path(D, "hallmark_audit/audit_inputs.rds"))
z <- IN$z[, as.character(IN$meta$sample)]; M <- as.data.table(IN$meta)
Q1 <- partial_table(z, M$rate)
say("\n-- Q1, %d counted libraries: correlation with the measured division rate --", nrow(M))
show(Q1)
fwrite(Q1, file.path(O, "counted_libraries_partials.tsv"), sep = "\t")
## the other shared-gene rule, where it differs
ALT <- rbind(partial_table(z, M$rate, rule = "both sides", boot = FALSE)[splicing_definition == "selected_177"],
             partial_table(z, M$rate, rule = "comparator only (sensitivity)", boot = FALSE)[splicing_definition != "selected_177"])
ALT <- ALT[n_shared_removed > 0 & comparator != "metabolism_of_rna"]
fwrite(ALT, file.path(O, "shared_rule_sensitivity.tsv"), sep = "\t")
say("\nshared-gene rule sensitivity (Q1), the other rule:")
print(ALT[, .(splicing_definition, comparator, shared_rule, n_shared_removed,
              `spl|cmp` = round(partial_splicing, 3), `cmp|spl` = round(partial_comparator, 3))])

## ---- Q2: the donor cohort ----------------------------------------------------
DON <- gse113957_donors()
MET <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))[, .(srr, log_depth)]
DON <- as.data.table(merge(DON, MET, by = "srr"))
d <- DON[gse113957_keep(DON, "primary")]
rc <- "public_data_tierA/recount3"
gs <- fread(cmd = paste("zcat", file.path(rc, "sra.gene_sums.SRP144355.G026.gz"), "| tail -n +2"))
GM <- as.matrix(gs[, -1, with = FALSE]); rownames(GM) <- sub("\\..*$", "", gs[[1]])
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]; GM <- rowsum(GM, rownames(GM))
y <- DGEList(GM[, d$srr]); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
zz <- t(scale(t(lc[apply(lc, 1, sd) > 0, ])))
## donor_scores.tsv keeps its original columns (run_residual_meta.R reads premrna_reactome);
## translation, metabolism_of_rna and rrna_processing there are the comparators with the
## 177 genes removed, as before
d[, splicing := sc(zz, CORE)][, prolif := sc(zz, PROLIF)]
for (k in names(WANT)) set(d, j = k, value = sc(zz, pair_sets("selected_177", k)$K))
d[, premrna_reactome := sc(zz, PREMRNA)][, mrna_splicing_reactome := sc(zz, MRNASPL)]
d[, mrna_splicing_not_selected := sc(zz, DEFS$mrna_splicing_not_selected)]
d[, premrna_processing_not_selected := sc(zz, DEFS$premrna_processing_not_selected)]
fwrite(d, file.path(O, "donor_scores.tsv"), sep = "\t")
say("\ndonors: %d | genes measured: selected %d, mRNA Splicing %d, pre-mRNA processing %d",
    nrow(d), ng(zz, CORE), ng(zz, MRNASPL), ng(zz, PREMRNA))

base <- "I(age/10) + log_depth + repo + instr + sex"
age_row <- function(def, comparator, side, metric, y, extra, label, n_genes) {
  dd <- cbind(d, .y = y); if (!is.null(extra)) dd <- cbind(dd, extra)
  f <- as.formula(paste(".y ~", base, if (!is.null(extra)) paste("+", paste(names(extra), collapse = " + "))))
  m <- lm(f, dd); s <- summary(m)$coef; ci <- confint(m)
  data.table(splicing_definition = def, selection = SELECTION[[def]], comparator = comparator, side = side,
             metric = metric, model = label, n_donors = nobs(m), n_genes = n_genes,
             beta = s["I(age/10)", 1], lo = ci["I(age/10)", 1], hi = ci["I(age/10)", 2], p = s["I(age/10)", 4])
}
CLAB <- c(translation = "translation", metabolism_of_rna = "Metabolism of RNA", rrna_processing = "rRNA processing")
Q2 <- rbindlist(lapply(names(DEFS), function(def) {
  P <- lapply(setNames(names(WANT), names(WANT)), function(k) pair_sets(def, k))
  whole <- sc(zz, DEFS[[def]]); pr <- data.frame(prolif_score = d$prolif)
  r <- c(
    list(age_row(def, NA, "splicing", "splicing", whole, NULL, "technical covariates only", ng(zz, DEFS[[def]]))),
    lapply(names(WANT), function(k) age_row(def, k, "comparator", k, sc(zz, P[[k]]$K), NULL,
                                           "technical covariates only", ng(zz, P[[k]]$K))),
    list(age_row(def, NA, "splicing", "splicing", whole, pr, "+ proliferation score", ng(zz, DEFS[[def]]))),
    lapply(names(WANT), function(k) age_row(def, k, "comparator", k, sc(zz, P[[k]]$K), pr,
                                           "+ proliferation score", ng(zz, P[[k]]$K))),
    unlist(lapply(names(WANT), function(k) list(
      age_row(def, k, "splicing", "splicing", sc(zz, P[[k]]$D), cbind(pr, other_set = sc(zz, P[[k]]$K)),
              paste("+ proliferation +", CLAB[[k]]), ng(zz, P[[k]]$D)),
      age_row(def, k, "comparator", k, sc(zz, P[[k]]$K), cbind(pr, other_set = sc(zz, P[[k]]$D)),
              "+ proliferation + splicing", ng(zz, P[[k]]$K)))), recursive = FALSE))
  r <- rbindlist(r); r[, FDR := p.adjust(p, "BH")]; r }))   # BH within each splicing definition
say("\n-- Q2, %d adult donors: age effect per decade --", nrow(d))
print(Q2[, .(splicing_definition, metric, model, n_genes, beta = round(beta, 3), lo = round(lo, 3),
             hi = round(hi, 3), p = signif(p, 2))], nrows = 100)
fwrite(Q2, file.path(O, "donor_age_models.tsv"), sep = "\t")

## ---- Q3: GTEx cultures -------------------------------------------------------
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
g <- as.data.table(TS$culture$d); GS <- TS$culture$S
## (i) the stored programme scores, as the original Q3 read them. Definitions are those of
## run_gtex_boundary.R; genes shared between programmes were NOT removed there.
RDS_DEF <- c(`splicing 96` = "preregistered 96-gene selected set",
             `pre-mRNA processing` = "Reactome Processing of Capped Intron-Containing Pre-mRNA (unselected)",
             `mRNA splicing` = "Reactome mRNA Splicing - Major Pathway (unselected)")
RDS_CMP <- c(translation = "Reactome Translation",
             `rRNA processing` = "Reactome rRNA processing in the nucleus and cytosol")
RDS_TGT <- list(`20-marker proliferation score` = g$prolif,
                `Reactome Cell Cycle, Mitotic programme score (the variable the original Q3 used)` =
                  GS[, "cell cycle (positive control)"])
Q3r <- rbindlist(lapply(names(RDS_TGT), function(tg) rbindlist(lapply(names(RDS_DEF), function(a)
  rbindlist(lapply(names(RDS_CMP), function(b) { t <- RDS_TGT[[tg]]
    pa <- pcor(GS[, a], t, GS[, b]); pb <- pcor(GS[, b], t, GS[, a])
    data.table(splicing_definition = a, definition_detail = RDS_DEF[[a]], comparator = b,
               comparator_detail = RDS_CMP[[b]], target = tg, n = nrow(g),
               n_genes_splicing = TS$culture$NG[[a]], n_genes_comparator = TS$culture$NG[[b]],
               shared_genes_removed = FALSE,
               rho_splicing = sp(GS[, a], t)["rho"], rho_comparator = sp(GS[, b], t)["rho"],
               partial_splicing = pa["rho"], partial_splicing_lo = pa["lo"], partial_splicing_hi = pa["hi"],
               partial_comparator = pb["rho"], partial_comparator_lo = pb["lo"], partial_comparator_hi = pb["hi"],
               asymmetry = pa["rho"] - pb["rho"]) }))))))
say("\n-- Q3 (i), %d GTEx cultures, stored programme scores --", nrow(g))
print(Q3r[, .(splicing_definition, comparator, target = substr(target, 1, 22), rho_spl = round(rho_splicing, 3),
              rho_cmp = round(rho_comparator, 3), `spl|cmp` = round(partial_splicing, 3),
              `cmp|spl` = round(partial_comparator, 3))])
fwrite(Q3r, file.path(O, "gtex_culture_partials_rds_programmes.tsv"), sep = "\t")

## (ii) gene-level scores for the same 511 cultures, normalised exactly as in
## run_gtex_boundary.R (same samples, filterByExpr, TMM, log-CPM, gene-wise z), so that the
## selected 177-gene set, the Q1/Q2 Reactome definitions and the shared-gene rule apply
G <- "public_data_tierA/gtex"
hn <- names(fread(cmd = sprintf("zcat %s/gtex.gene_sums.SKIN.G026.gz | grep -v '^##' | head -2", G)))
idcol <- if (sum(as.character(g$rail_id) %in% hn) > sum(g$external_id %in% hn)) "rail_id" else "external_id"
ids <- as.character(g[[idcol]]); stopifnot(all(ids %in% hn))
gx <- fread(cmd = sprintf("zcat %s/gtex.gene_sums.SKIN.G026.gz | grep -v '^##'", G), select = c(hn[1], ids))
sym <- key[sub("\\..*$", "", gx[[1]])]
X <- as.matrix(gx[, -1, with = FALSE])[, ids]; rm(gx); gc()
ok <- !is.na(sym) & sym != ""; X <- rowsum(X[ok, , drop = FALSE], sym[ok])
yg <- DGEList(X); yg <- yg[filterByExpr(yg), , keep.lib.sizes = FALSE]; yg <- calcNormFactors(yg)
lg <- cpm(yg, log = TRUE); rm(X, yg); gc()
zg <- t(scale(t(lg[apply(lg, 1, sd) > 0, ]))); rm(lg); gc()
chk <- c(premrna = max(abs(sc(zg, setdiff(PREMRNA, PROLIF)) - GS[, "pre-mRNA processing"])),
         translation = max(abs(sc(zg, setdiff(SETS[["Translation"]], PROLIF)) - GS[, "translation"])),
         prolif = max(abs(sc(zg, PROLIF) - g$prolif)))
say("\nGTEx gene-level recomputation against the stored scores, max |difference|: %s",
    paste(sprintf("%s %.2e", names(chk), chk), collapse = "; "))
stopifnot(all(chk < 1e-8))
tech <- g[, .(rin, isch, loglib)]
Q3 <- partial_table(zg, g$prolif, cov = tech)
Q3[, target := "20-marker proliferation score"]
say("\n-- Q3 (ii), %d GTEx cultures, gene-level, shared genes removed: correlation with the proliferation score --", nrow(g))
show(Q3)
print(Q3[, .(splicing_definition, comparator, n_tech, `spl|cmp,tech` = round(partial_splicing_tech, 2),
             `cmp|spl,tech` = round(partial_comparator_tech, 2))], nrows = 40)
fwrite(Q3, file.path(O, "gtex_culture_partials.tsv"), sep = "\t")
## sensitivity: the variable the original Q3 used as "proliferation". It is a Reactome
## programme score (Cell Cycle, Mitotic, less the 20 markers) and, unlike the 20-marker
## score, it shares genes with the sets being compared (nucleoporins with pre-mRNA
## processing, proteasome and ubiquitin genes with Translation), so it is not the primary
CCM <- SETS[["Cell Cycle, Mitotic"]]
Q3c <- partial_table(zg, GS[, "cell cycle (positive control)"], cov = tech)
Q3c[, target := "Reactome Cell Cycle, Mitotic programme score (the variable the original Q3 used)"]
Q3c[, n_splicing_genes_in_target := vapply(seq_len(.N), function(i)
  length(intersect(pair_sets(splicing_definition[i], comparator[i])$D, CCM)), 0L)]
Q3c[, n_comparator_genes_in_target := vapply(seq_len(.N), function(i)
  length(intersect(pair_sets(splicing_definition[i], comparator[i])$K, CCM)), 0L)]
say("\n-- Q3 sensitivity, target = Cell Cycle, Mitotic programme score --"); show(Q3c)
fwrite(Q3c, file.path(O, "gtex_culture_partials_cellcycle_target.tsv"), sep = "\t")

fwrite(data.table(splicing_definition = names(DEFS), label = DEF_LABEL[names(DEFS)],
                  selection = SELECTION[names(DEFS)], n_genes = lengths(DEFS),
                  n_in_selected_177 = vapply(DEFS, function(s) length(intersect(s, CORE)), 0L),
                  measured_counted_libraries = vapply(DEFS, function(s) ng(z, s), 0L),
                  measured_donor_cohort = vapply(DEFS, function(s) ng(zz, s), 0L),
                  measured_gtex_cultures = vapply(DEFS, function(s) ng(zg, s), 0L)),
       file.path(O, "splicing_definitions.tsv"), sep = "\t")
sink(file.path(O, "sessionInfo.txt")); print(sessionInfo()); sink()
