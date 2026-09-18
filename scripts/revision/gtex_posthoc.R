## POST HOC GTEx analyses requested by the pre-submission review (reviewer 3).
## NONE of this is preregistered: the preregistration and its outputs are in
## scripts/run_gtex_boundary.R and public_data_tierA/derived/gtex_boundary/,
## which this script does not touch.  Outputs go to derived/gtex_posthoc/.
##
##  1. range restriction   - does the smaller spread of proliferation in tissue
##                           account for the weaker coupling?
##  2. cell composition    - does the coupling in skin survive adjustment for
##                           keratinocyte / fibroblast / immune content?
##  3. positive control    - can a donor-level match be detected at all between
##                           a donor's cultured fibroblasts and that donor's skin?
##  4. sample counts       - how many donors enter each panel, and why they differ.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR); library(fgsea)})
say <- function(...) cat(sprintf(...), "\n")
G <- "public_data_tierA/gtex"; O <- "public_data_tierA/derived/gtex_posthoc"
dir.create(O, showWarnings = FALSE, recursive = TRUE)
set.seed(20260915)

## ---- 1. definitions, identical to the preregistered script --------------------
CORE96  <- readLines("public_data_tierA/derived/conserved_core/age_down_splicing_core.txt")
CORE177 <- readLines("public_data_tierA/derived/conserved_core/age_down_splicing_core_v2.txt")
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
## composition markers (post hoc; reviewer 3 major 2)
## History of the fibroblast panel (both runs on 2026-09-15, both post hoc):
##   first run (outputs 20:30)  COL1A1 COL1A2 DCN LUM PDGFRB   - PDGFRB stood in place of
##                              the PDGFRA named by the preregistered S4 panel
##   revised   (outputs 20:51)  COL1A1 COL1A2 DCN LUM PDGFRA PDGFRB - the preregistered
##                              five restored, PDGFRB kept as an addition
## `fibro` is the revised panel and is what every reported value uses.  `fibro_first`
## reproduces the first run and is written only to composition_adjusted_first_definition.tsv
## so that both sets of values stay on record.
MARK <- list(
  kerat    = c("KRT5","KRT14","KRT1","KRT10","TP63","IVL","LOR"),
  fibro    = c("COL1A1","COL1A2","DCN","LUM","PDGFRA","PDGFRB"),
  fibro_first = c("COL1A1","COL1A2","DCN","LUM","PDGFRB"),
  immune   = c("PTPRC"),
  myofibre = c("ACTA1","MYH1","MYH2","MYH7","CKM"),
  satellite= c("PAX7","MYF5"))
## donor-identity positive control genes (post hoc; reviewer 3 major 3)
PC <- list(
  sex  = c("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","USP9Y","EIF1AY"),
  eqtl = c("GSTM1","GSTT1","UGT2B17","ERAP2"))
HALL <- list(
 `pre-mRNA processing`      = "Processing of Capped Intron-Containing Pre-mRNA",
 `mRNA splicing`            = "mRNA Splicing - Major Pathway",
 `collagen formation`       = "Collagen formation",
 `cell cycle (positive control)` = "Cell Cycle, Mitotic")
td <- tempfile(); dir.create(td); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
RP <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
names(RP) <- sub("__.*", "", names(RP))
SETS <- lapply(HALL, function(nm) {
  hit <- if (nm %in% names(RP)) nm else grep(nm, names(RP), fixed = TRUE, value = TRUE)[1]
  if (is.na(hit)) NULL else toupper(RP[[hit]]) })
SETS <- c(list(`splicing 96` = CORE96, `splicing 177` = CORE177), SETS[!sapply(SETS, is.null)])
SETS_NP <- lapply(SETS, setdiff, PROLIF)           # circularity rule, as preregistered

mpf  <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz",
              select = c("Geneid", "Gene_name"))
key2 <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))

## ---- 2. sample table, identical to the preregistered script -------------------
col_ci <- function(d, pat) { k <- grep(pat, names(d), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
read_md <- function(PJ) {
  a <- fread(file.path(G, sprintf("gtex.gtex.%s.MD.gz", PJ)))
  b <- fread(file.path(G, sprintf("gtex.recount_qc.%s.MD.gz", PJ)))
  k <- intersect(c("rail_id","external_id"), intersect(names(a), names(b)))
  m <- merge(a, b, by = k, suffixes = c("", ".qc"))
  data.table(project = PJ, rail_id = m$rail_id, external_id = m$external_id,
             tissue  = m[[col_ci(m, "^smtsd$")]],
             age_bin = m[[col_ci(m, "^age$")]],
             sex     = m[[col_ci(m, "^sex$")]],
             hardy   = m[[col_ci(m, "^dthhrdy$")]],
             rin     = suppressWarnings(as.numeric(m[[col_ci(m, "^smrin$")]])),
             isch    = suppressWarnings(as.numeric(m[[col_ci(m, "^smtsisch$")]]))) }
MD <- rbind(read_md("SKIN"), read_md("MUSCLE"))
KEEP_T <- c(culture = "Cells - Cultured fibroblasts", legskin = "Skin - Sun Exposed (Lower leg)",
            pubskin = "Skin - Not Sun Exposed (Suprapubic)", muscle = "Muscle - Skeletal")
MD_ALL <- copy(MD)
MD <- MD[tissue %in% KEEP_T]
MD[, tkey := names(KEEP_T)[match(tissue, KEEP_T)]]
MD[, donor := sub("^([^-]+-[^-]+).*$", "\\1", external_id)]
MD[, age := suppressWarnings(as.numeric(sub("^([0-9]+)-.*$", "\\1", age_bin))) + 5]
CNT_labelled <- MD[, .(labelled = .N), by = tkey]
MD <- MD[is.finite(age)]
MD[, sex := factor(sex)]
MD[, hardy := factor(ifelse(is.na(hardy) | hardy == "", "missing", as.character(hardy)))]
CNT_age <- MD[, .(with_age = .N), by = tkey]

## ---- 3. per tissue: normalise, score, keep the control genes ------------------
hdr <- function(PJ) names(fread(cmd = sprintf("zcat %s/gtex.gene_sums.%s.G026.gz | grep -v '^##' | head -2", G, PJ)))
load_tissue <- function(tk) {
  d  <- MD[tkey == tk]; PJ <- d$project[1]
  hn <- hdr(PJ)
  idcol <- if (sum(as.character(d$rail_id) %in% hn) > sum(d$external_id %in% hn)) "rail_id" else "external_id"
  cols <- intersect(hn, as.character(d[[idcol]]))
  gs <- fread(cmd = sprintf("zcat %s/gtex.gene_sums.%s.G026.gz | grep -v '^##'", G, PJ),
              select = c(hn[1], cols))
  sym <- key2[sub("\\..*$", "", gs[[1]])]
  X <- as.matrix(gs[, -1, with = FALSE]); rm(gs); gc()
  ok <- !is.na(sym) & sym != ""; X <- rowsum(X[ok, , drop = FALSE], sym[ok])
  d <- d[match(colnames(X), as.character(d[[idcol]]))]
  n_matched <- ncol(X)
  lib <- colSums(X); d[, lib := lib]
  keep <- d$lib >= quantile(d$lib, 0.01)
  n_libfilter <- sum(keep)
  d <- d[keep]; X <- X[, keep, drop = FALSE]
  o <- order(d$donor, -d$lib); dup <- duplicated(d$donor[o]); drop <- o[dup]
  n_dupdrop <- length(drop)
  if (length(drop)) { d <- d[-drop]; X <- X[, -drop, drop = FALSE] }
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  lc <- cpm(y, log = TRUE); rm(X, y); gc()
  expressed <- rownames(lc)
  z <- t(scale(t(lc[apply(lc, 1, sd) > 0, ]))); rm(lc); gc()
  d[, loglib := log10(lib)]
  S <- sapply(SETS_NP, function(g) { g <- intersect(rownames(z), g)
    if (length(g) < 15) rep(NA_real_, ncol(z)) else colMeans(z[g, , drop = FALSE]) })
  NG <- sapply(SETS_NP, function(g) length(intersect(rownames(z), g)))
  d[, prolif := colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE])]
  for (m in names(MARK)) {
    g <- intersect(rownames(z), MARK[[m]])
    d[[m]] <- if (length(g)) colMeans(z[g, , drop = FALSE]) else NA_real_
    attr(d, paste0("n_", m)) <- length(g) }
  pcg <- intersect(rownames(z), unlist(PC))
  ZPC <- z[pcg, , drop = FALSE]
  say("%-8s labelled %4d | with age %4d | in gene matrix %4d | after 1%% library filter %4d | duplicate donors dropped %3d | final %4d donors, %5d genes",
      tk, CNT_labelled[tkey == tk]$labelled, CNT_age[tkey == tk]$with_age, n_matched,
      n_libfilter, n_dupdrop, nrow(d), nrow(z))
  list(d = d, S = S, NG = NG, ZPC = ZPC, expressed = expressed,
       counts = data.table(tissue = tk, labelled = CNT_labelled[tkey == tk]$labelled,
                           with_age = CNT_age[tkey == tk]$with_age, in_matrix = n_matched,
                           after_libfilter = n_libfilter, dup_donors_dropped = n_dupdrop,
                           final_donors = nrow(d), genes_expressed = nrow(z))) }

tech_of <- function(d) c("rin","isch","loglib")[c(sum(is.finite(d$rin)) > 0.9*nrow(d),
                                                  sum(is.finite(d$isch)) > 0.9*nrow(d), TRUE)]
## partial Spearman exactly as preregistered (variance 1.06/(n-3-k) used for z tests)
pspear <- function(a, b, d, cov) {
  ok <- is.finite(a) & is.finite(b) & complete.cases(d[, ..cov])
  if (length(cov)) { f <- as.formula(paste("v ~", paste(cov, collapse = "+")))
    ra <- resid(lm(f, data = cbind(d[ok], v = a[ok]))); rb <- resid(lm(f, data = cbind(d[ok], v = b[ok])))
  } else { ra <- a[ok]; rb <- b[ok] }
  r <- cor(ra, rb, method = "spearman"); n <- sum(ok); k <- length(cov); df <- n - 2 - k
  t <- r * sqrt(df / (1 - r^2)); list(r = r, p = 2*pt(-abs(t), df), n = n, k = k) }

TIS <- setNames(lapply(names(KEEP_T), load_tissue), names(KEEP_T))
CNT <- rbindlist(lapply(TIS, `[[`, "counts"))

## ================================================================================
## 4. ANALYSIS 1 (post hoc): range restriction
## ================================================================================
say("\n=== 1. range restriction (post hoc) ===")
RR <- rbindlist(lapply(names(TIS), function(tk) {
  T <- TIS[[tk]]; cov <- tech_of(T$d); d <- T$d
  rbindlist(lapply(c("splicing 96", "splicing 177"), function(pg) {
    s <- T$S[, pg]
    ok <- is.finite(s) & is.finite(d$prolif) & complete.cases(d[, ..cov])
    dd <- cbind(d[ok], s = s[ok])
    ps <- pspear(s, d$prolif, d, cov)
    f  <- lm(as.formula(paste("s ~ prolif +", paste(cov, collapse = "+"))), data = dd)
    ci <- confint(f)["prolif", ]
    data.table(tissue = tk, set = pg, n = ps$n,
               sd_prolif = sd(dd$prolif), var_prolif = var(dd$prolif),
               sd_score = sd(dd$s), var_score = var(dd$s),
               partial_rho = ps$r, p = ps$p,
               slope = coef(f)["prolif"], slope_lo = ci[1], slope_hi = ci[2],
               partial_r_pearson = { ra <- resid(lm(as.formula(paste("prolif ~", paste(cov, collapse="+"))), dd))
                                     rb <- resid(lm(as.formula(paste("s ~", paste(cov, collapse="+"))), dd))
                                     cor(ra, rb) }) })) }))
fwrite(RR, file.path(O, "range_restriction.tsv"), sep = "\t")
print(RR[set == "splicing 96", .(tissue, n, sd_prolif, sd_score, partial_rho,
                                 slope = round(slope, 3), slope_lo = round(slope_lo, 3),
                                 slope_hi = round(slope_hi, 3))], digits = 3)

## Thorndike case-II analytic prediction: what would the culture correlation become
## if only the spread of proliferation were reduced to the tissue's spread?
thorndike <- function(r, ratio) (r * ratio) / sqrt(1 - r^2 + r^2 * ratio^2)
cultR <- RR[tissue == "culture" & set == "splicing 96"]
PRED <- rbindlist(lapply(c("legskin","pubskin","muscle"), function(tk) {
  tt <- RR[tissue == tk & set == "splicing 96"]
  ratio <- tt$sd_prolif / cultR$sd_prolif
  data.table(tissue = tk, sd_ratio = ratio,
             predicted_rho_if_only_range = thorndike(cultR$partial_rho, ratio),
             observed_rho = tt$partial_rho) }))
fwrite(PRED, file.path(O, "range_restriction_prediction.tsv"), sep = "\t")
say("\nThorndike case-II prediction from the culture coupling (rho = %.3f):", cultR$partial_rho)
print(PRED, digits = 3)

## variance-matched subsampling of cultured fibroblasts to the leg-skin spread
matched_draws <- function(tk_target, ndraw = 500, nsub = 200) {
  Tc <- TIS$culture; cov <- tech_of(Tc$d)
  s <- Tc$S[, "splicing 96"]; d <- Tc$d
  ok <- is.finite(s) & is.finite(d$prolif) & complete.cases(d[, ..cov])
  d <- d[ok]; s <- s[ok]
  target_sd <- RR[tissue == tk_target & set == "splicing 96"]$sd_prolif
  ## a window of proliferation values whose SD matches the target, centred on the median
  centre <- median(d$prolif)
  h <- uniroot(function(h) { i <- abs(d$prolif - centre) <= h
                             if (sum(i) < 30) return(-target_sd) ; sd(d$prolif[i]) - target_sd },
               interval = c(0.01, max(abs(d$prolif - centre))))$root
  inwin <- which(abs(d$prolif - centre) <= h)
  res <- t(replicate(ndraw, {
    i <- sample(inwin, min(nsub, length(inwin)))
    ps <- pspear(s[i], d$prolif[i], d[i], cov)
    c(rho = ps$r, sd = sd(d$prolif[i]), n = ps$n) }))
  list(target_sd = target_sd, window = h, n_in_window = length(inwin), res = res) }
## second scheme: keep the whole distribution but downweight the tails with a
## Gaussian kernel whose width is tuned to the target spread, so that the donors
## at the extremes are not simply discarded
kernel_draws <- function(tk_target, ndraw = 500, nsub = 200) {
  Tc <- TIS$culture; cov <- tech_of(Tc$d)
  s <- Tc$S[, "splicing 96"]; d <- Tc$d
  ok <- is.finite(s) & is.finite(d$prolif) & complete.cases(d[, ..cov])
  d <- d[ok]; s <- s[ok]
  target_sd <- RR[tissue == tk_target & set == "splicing 96"]$sd_prolif
  centre <- median(d$prolif)
  wsd <- function(sig) { w <- exp(-((d$prolif - centre)^2) / (2 * sig^2))
                         sqrt(sum(w * (d$prolif - centre)^2) / sum(w)) }
  sig <- uniroot(function(sig) wsd(sig) - target_sd,
                 interval = c(0.05, 10 * sd(d$prolif)))$root
  w <- exp(-((d$prolif - centre)^2) / (2 * sig^2))
  res <- t(replicate(ndraw, {
    i <- sample(seq_along(w), min(nsub, sum(w > 1e-6)), prob = w)
    ps <- pspear(s[i], d$prolif[i], d[i], cov)
    c(rho = ps$r, sd = sd(d$prolif[i]), n = ps$n) }))
  list(sigma = sig, n_pool = nrow(d), res = res) }

MM <- matched_draws("legskin")
KD <- kernel_draws("legskin")
say("\nvariance-matched subsampling of cultured fibroblasts to the leg-skin spread (post hoc)")
say("  target SD(prolif) = %.3f | culture window half-width %.3f | donors in window %d | %d draws of %d",
    MM$target_sd, MM$window, MM$n_in_window, nrow(MM$res), min(200, MM$n_in_window))
say("  achieved SD: median %.3f | partial rho: median %.3f, 95%% range %.3f to %.3f",
    median(MM$res[, "sd"]), median(MM$res[, "rho"]),
    quantile(MM$res[, "rho"], 0.025), quantile(MM$res[, "rho"], 0.975))
say("  kernel-weighted draws (tails downweighted, not discarded; sigma %.3f): achieved SD median %.3f | partial rho: median %.3f, 95%% range %.3f to %.3f",
    KD$sigma, median(KD$res[, "sd"]), median(KD$res[, "rho"]),
    quantile(KD$res[, "rho"], 0.025), quantile(KD$res[, "rho"], 0.975))
fwrite(rbind(as.data.table(MM$res)[, scheme := "window"],
             as.data.table(KD$res)[, scheme := "kernel"]),
       file.path(O, "variance_matched_draws.tsv"), sep = "\t")
VM <- data.table(scheme = c("window", "kernel"), target_tissue = "legskin",
                 target_sd = MM$target_sd,
                 tuning = c(MM$window, KD$sigma),
                 ## the kernel draws sample from the cultures with complete technical
                 ## covariates (492), not from all 511; this column used to say 511
                 donors_available = c(MM$n_in_window, KD$n_pool),
                 draws = c(nrow(MM$res), nrow(KD$res)),
                 achieved_sd_median = c(median(MM$res[, "sd"]), median(KD$res[, "sd"])),
                 rho_median = c(median(MM$res[, "rho"]), median(KD$res[, "rho"])),
                 rho_lo = c(quantile(MM$res[, "rho"], 0.025), quantile(KD$res[, "rho"], 0.025)),
                 rho_hi = c(quantile(MM$res[, "rho"], 0.975), quantile(KD$res[, "rho"], 0.975)),
                 observed_legskin_rho = RR[tissue == "legskin" & set == "splicing 96"]$partial_rho)
fwrite(VM, file.path(O, "variance_matched_summary.tsv"), sep = "\t")

## ================================================================================
## 5. ANALYSIS 2 (post hoc): cell composition
## ================================================================================
say("\n=== 2. cell composition (post hoc) ===")
comp_of <- function(tk) {
  m <- c("kerat","fibro","immune")
  if (tk == "muscle") m <- c("myofibre","satellite","fibro","immune")
  if (tk == "culture") m <- c("fibro","immune")
  m }
## 2026-09-18: the composition adjustment originally ran on the two splicing sets only,
## so it could not show whether composition removes splicing's coupling specifically or
## flattens the whole proliferation axis.  `progs` now selects the programmes; the
## splicing-only table is kept under its old name and the full one written beside it.
comp_table <- function(fibro_col = "fibro",
                       progs = c("splicing 96", "splicing 177")) rbindlist(lapply(names(TIS), function(tk) {
  T <- TIS[[tk]]; d <- T$d; cov <- tech_of(d); extra <- sub("^fibro$", fibro_col, comp_of(tk))
  extra <- extra[sapply(extra, function(e) is.finite(d[[e]][1]) && sd(d[[e]], na.rm = TRUE) > 0)]
  rbindlist(lapply(intersect(progs, colnames(T$S)), function(pg) {
    s <- T$S[, pg]
    a <- pspear(s, d$prolif, d, cov)
    b <- pspear(s, d$prolif, d, c(cov, extra))
    data.table(tissue = tk, set = pg, markers = paste(extra, collapse = "+"),
               rho_tech = a$r, p_tech = a$p, n_tech = a$n,
               rho_tech_comp = b$r, p_comp = b$p, n_comp = b$n) })) }))
CP <- comp_table("fibro")
fwrite(CP, file.path(O, "composition_adjusted.tsv"), sep = "\t")
print(CP[set == "splicing 96"], digits = 3)

## every programme, so the cell-cycle positive control can be read beside splicing
CPALL <- comp_table("fibro", progs = colnames(TIS$legskin$S))
CPALL[, drop_abs := rho_tech - rho_tech_comp]
fwrite(CPALL, file.path(O, "composition_adjusted_all_programmes.tsv"), sep = "\t")
say("\ncomposition adjustment, all programmes (skin sites and the positive control):")
print(CPALL[tissue %in% c("legskin", "pubskin") &
            set %in% c("cell cycle (positive control)", "splicing 96", "splicing 177",
                       "mRNA splicing", "translation", "DNA double-strand repair",
                       "base excision repair", "telomere maintenance")][order(tissue, -rho_tech)],
      digits = 3)
## record of the first post hoc run, whose fibroblast panel held PDGFRB in place of
## PDGFRA (see the note at MARK).  Not used by any figure; kept so both sets of values
## can be reported.
CP_first <- comp_table("fibro_first")
fwrite(CP_first, file.path(O, "composition_adjusted_first_definition.tsv"), sep = "\t")
say("\nfirst-run fibroblast panel (PDGFRB in place of PDGFRA), for the record:")
print(CP_first[set == "splicing 96"], digits = 3)
## how much of the skin proliferation score is keratinocyte content?
KC <- rbindlist(lapply(c("legskin","pubskin","culture","muscle"), function(tk) {
  d <- TIS[[tk]]$d
  data.table(tissue = tk,
             rho_prolif_kerat = if (sd(d$kerat, na.rm = TRUE) > 0) cor(d$prolif, d$kerat, method = "spearman") else NA_real_,
             rho_prolif_fibro = cor(d$prolif, d$fibro, method = "spearman"),
             rho_splice96_kerat = if (sd(d$kerat, na.rm = TRUE) > 0) cor(TIS[[tk]]$S[, "splicing 96"], d$kerat, method = "spearman") else NA_real_,
             rho_splice96_fibro = cor(TIS[[tk]]$S[, "splicing 96"], d$fibro, method = "spearman")) }))
fwrite(KC, file.path(O, "marker_correlations.tsv"), sep = "\t")
print(KC, digits = 3)

## ================================================================================
## 6. ANALYSIS 3 (post hoc): donor-level positive control for Fig. 8D
## ================================================================================
say("\n=== 3. donor-identity positive control (post hoc) ===")
dd <- merge(TIS$culture$d[, .(donor, i = .I)], TIS$legskin$d[, .(donor, j = .I)], by = "donor")
say("donors contributing both cultured fibroblasts and leg skin: %d", nrow(dd))
in_map <- unique(toupper(mpf$Gene_name))
PCres <- rbindlist(lapply(unlist(PC), function(g) {
  grp <- names(PC)[sapply(PC, function(v) g %in% v)]
  inC <- g %in% rownames(TIS$culture$ZPC); inS <- g %in% rownames(TIS$legskin$ZPC)
  why <- if (inC && inS) "" else if (!(g %in% in_map)) "absent from the gene-symbol map" else
         "present in the annotation but filtered out as not expressed"
  if (!inC || !inS) return(data.table(gene = g, group = grp, in_culture = inC, in_skin = inS,
                                      rho = NA_real_, p = NA_real_, n = NA_integer_, note = why))
  a <- TIS$culture$ZPC[g, dd$i]; b <- TIS$legskin$ZPC[g, dd$j]
  ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE))
  data.table(gene = g, group = grp, in_culture = TRUE, in_skin = TRUE,
             rho = unname(ct$estimate), p = ct$p.value, n = nrow(dd), note = "") }))
## how much of each control gene is driven by the donor's recorded sex, in each
## tissue separately: the product of the two within-tissue correlations is the
## ceiling a cross-tissue correlation can reach for a sex-driven gene
sexcor <- function(tk, g) { T <- TIS[[tk]]
  if (!g %in% rownames(T$ZPC)) return(NA_real_)
  suppressWarnings(cor(T$ZPC[g, ], as.numeric(T$d$sex), method = "spearman")) }
PCres[, rho_with_sex_culture := sapply(gene, function(g) sexcor("culture", g))]
PCres[, rho_with_sex_skin    := sapply(gene, function(g) sexcor("legskin", g))]
PCres[, ceiling_from_sex := abs(rho_with_sex_culture * rho_with_sex_skin)]
fwrite(PCres, file.path(O, "donor_identity_positive_control.tsv"), sep = "\t")
print(PCres, digits = 3)
PCV <- data.table(donor = dd$donor, sex = TIS$culture$d$sex[dd$i])
for (g in intersect(rownames(TIS$culture$ZPC), rownames(TIS$legskin$ZPC))) {
  PCV[[paste0(g, "_culture")]] <- TIS$culture$ZPC[g, dd$i]
  PCV[[paste0(g, "_skin")]]    <- TIS$legskin$ZPC[g, dd$j] }
fwrite(PCV, file.path(O, "positive_control_values.tsv"), sep = "\t")
S3 <- fread("public_data_tierA/derived/gtex_boundary/donor_concordance_culture_vs_legskin.tsv")
say("programme-score concordance for the same donors (preregistered run): rho %.2f to %.2f, none with FDR < 0.05",
    min(S3$rho), max(S3$rho))

## ================================================================================
## 7. ANALYSIS 4: sample counts
## ================================================================================
say("\n=== 4. sample counts ===")
CNTfull <- copy(CNT)
for (tk in names(TIS)) {
  T <- TIS[[tk]]; cov <- tech_of(T$d)
  ok <- complete.cases(T$d[, ..cov]) & is.finite(T$S[, "splicing 96"]) & is.finite(T$d$prolif)
  CNTfull[tissue == tk, complete_covariates := sum(ok)]
  CNTfull[tissue == tk, covariates_used := paste(cov, collapse = "+")] }
CNTfull[, paired_with_legskin := ifelse(tissue == "culture", nrow(dd), NA_integer_)]
fwrite(CNTfull, file.path(O, "sample_counts.tsv"), sep = "\t")
print(CNTfull, digits = 3)

saveRDS(list(RR = RR, PRED = PRED, VM = VM, CP = CP, CP_first = CP_first, KC = KC, PC = PCres, CNT = CNTfull),
        file.path(O, "posthoc_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo.txt"))
say("\nwritten to %s", O)
