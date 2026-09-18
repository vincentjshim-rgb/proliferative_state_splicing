## Culture-vs-tissue boundary test, exactly as preregistered in
## preregistration_gtex_culture_tissue_boundary_ko.md (SHA-256 locked before data were opened).
## GTEx v8 via recount3: cultured fibroblasts, sun-exposed and non-exposed skin, skeletal muscle.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR); library(fgsea)})
say <- function(...) cat(sprintf(...), "\n")
G <- "public_data_tierA/gtex"; O <- "public_data_tierA/derived/gtex_boundary"
dir.create(O, showWarnings = FALSE, recursive = TRUE)

## ---- 0. the preregistration must be unchanged --------------------------------
PRE <- "preregistration_gtex_culture_tissue_boundary_ko.md"
LOCK <- "61972b991dc6de770478cd7f6df39a5319e11107ef2993aac2faeee8eb8b774e"
h <- sub(" .*", "", system(paste("sha256sum", PRE), intern = TRUE))
if (h != LOCK) stop("preregistration changed since it was locked")
say("preregistration verified: %s", substr(h, 1, 16))

## ---- 1. fixed definitions ------------------------------------------------------
CORE   <- readLines("public_data_tierA/derived/conserved_core/age_down_splicing_core.txt")
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
KERAT  <- c("KRT5","KRT14","KRT1","KRT10")
DFIB   <- c("COL1A1","COL1A2","DCN","LUM","PDGFRA")
HALL <- list(
 `pre-mRNA processing`      = "Processing of Capped Intron-Containing Pre-mRNA",
 `mRNA splicing`            = "mRNA Splicing - Major Pathway",
 `translation`              = "Translation",
 `rRNA processing`          = "rRNA processing in the nucleus and cytosol",
 `proteasome`               = "Antigen processing: Ubiquitination & Proteasome degradation",
 `autophagy`                = "Macroautophagy",
 `chaperone / HSF1`         = "HSF1-dependent transactivation",
 `unfolded protein response`= "Unfolded Protein Response (UPR)",
 `DNA double-strand repair` = "DNA Double-Strand Break Repair",
 `base excision repair`     = "Base Excision Repair",
 `telomere maintenance`     = "Telomere Maintenance",
 `chromatin organisation`   = "Chromatin organization",
 `mitochondrial translation`= "Mitochondrial translation",
 `respiratory electron transport` = "Respiratory electron transport",
 `TCA cycle`                = "Citric acid cycle (TCA cycle)",
 `interferon signalling`    = "Interferon Signaling",
 `NF-kB / TNF`              = "TNFR1-induced NF-kappa-B signaling pathway",
 `extracellular matrix`     = "Extracellular matrix organization",
 `collagen formation`       = "Collagen formation",
 `insulin / IGF signalling` = "Signaling by Insulin receptor",
 `mTOR signalling`          = "mTORC1-mediated signalling",
 `cellular senescence`      = "Cellular Senescence",
 `interleukin signalling`   = "Signaling by Interleukins",
 `cell cycle (positive control)` = "Cell Cycle, Mitotic")
td <- tempfile(); dir.create(td); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
RP <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
names(RP) <- sub("__.*", "", names(RP))
SETS <- lapply(HALL, function(nm) {
  hit <- if (nm %in% names(RP)) nm else grep(nm, names(RP), fixed = TRUE, value = TRUE)[1]
  if (is.na(hit)) NULL else toupper(RP[[hit]]) })
SETS <- SETS[!sapply(SETS, is.null)]
SETS <- c(list(`splicing 96` = CORE), SETS)
## circularity rule: proliferation markers are removed from every programme
SETS_NP <- lapply(SETS, setdiff, PROLIF)
say("programmes: %d (+ splicing 96); markers removed from programmes: %d",
    length(SETS) - 1, sum(sapply(SETS, function(s) length(intersect(s, PROLIF)))))

## ENSG -> symbol, the map already used by run_outcome_vs_expression.R
mpf  <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz",
              select = c("Geneid", "Gene_name"))
key2 <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))

## ---- 2. sample table -----------------------------------------------------------
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
MD <- MD[tissue %in% KEEP_T]
MD[, tkey := names(KEEP_T)[match(tissue, KEEP_T)]]
MD[, donor := sub("^([^-]+-[^-]+).*$", "\\1", external_id)]
MD[, age := suppressWarnings(as.numeric(sub("^([0-9]+)-.*$", "\\1", age_bin))) + 5]
MD <- MD[is.finite(age)]
MD[, sex := factor(sex)]
MD[, hardy := factor(ifelse(is.na(hardy) | hardy == "", "missing", as.character(hardy)))]
say("\nsamples by tissue after label and age filters:"); print(MD[, .N, by = tkey])

## ---- 3. per tissue: normalise, score, test -------------------------------------
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
  lib <- colSums(X)
  ## preregistered exclusions: bottom 1% library size; one sample per donor (largest library)
  d[, lib := lib]
  keep <- d$lib >= quantile(d$lib, 0.01)
  d <- d[keep]; X <- X[, keep, drop = FALSE]
  o <- order(d$donor, -d$lib); dup <- duplicated(d$donor[o]); drop <- o[dup]
  if (length(drop)) { d <- d[-drop]; X <- X[, -drop, drop = FALSE] }
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  lc <- cpm(y, log = TRUE); rm(X, y); gc()
  z <- t(scale(t(lc[apply(lc, 1, sd) > 0, ]))); rm(lc); gc()
  d[, loglib := log10(lib)]
  S <- sapply(SETS_NP, function(g) { g <- intersect(rownames(z), g)
    if (length(g) < 15) rep(NA_real_, ncol(z)) else colMeans(z[g, , drop = FALSE]) })
  NG <- sapply(SETS_NP, function(g) length(intersect(rownames(z), g)))
  d[, prolif := colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE])]
  d[, kerat  := colMeans(z[intersect(rownames(z), KERAT),  , drop = FALSE])]
  d[, dfib   := colMeans(z[intersect(rownames(z), DFIB),   , drop = FALSE])]
  say("%-8s n = %d donors, %d genes expressed, %d proliferation markers present",
      tk, nrow(d), nrow(z), length(intersect(rownames(z), PROLIF)))
  list(d = d, S = S, NG = NG) }

## technical covariates actually available in a tissue (recorded, per the preregistration)
tech_of <- function(d) c("rin","isch","loglib")[c(sum(is.finite(d$rin)) > 0.9*nrow(d),
                                                  sum(is.finite(d$isch)) > 0.9*nrow(d), TRUE)]
pspear <- function(a, b, d, cov) {
  ok <- is.finite(a) & is.finite(b) & complete.cases(d[, ..cov])
  if (length(cov)) { f <- as.formula(paste("v ~", paste(cov, collapse = "+")))
    ra <- resid(lm(f, data = cbind(d[ok], v = a[ok]))); rb <- resid(lm(f, data = cbind(d[ok], v = b[ok])))
  } else { ra <- a[ok]; rb <- b[ok] }
  r <- cor(ra, rb, method = "spearman"); n <- sum(ok); k <- length(cov); df <- n - 2 - k
  t <- r * sqrt(df / (1 - r^2)); list(r = r, p = 2*pt(-abs(t), df), n = n, k = k) }
agefit <- function(s, d, extra = NULL) {
  cov <- c("sex", if (nlevels(droplevels(d$hardy)) > 1) "hardy", tech_of(d), extra)
  dd <- cbind(d, s = s); dd <- dd[is.finite(s) & complete.cases(dd[, ..cov])]
  f  <- summary(lm(as.formula(paste("s ~ age +", paste(cov, collapse = "+"))), data = dd))$coef
  c(beta = f["age","Estimate"]*10, p = f["age","Pr(>|t|)"], n = nrow(dd)) }   # per decade

TIS <- setNames(lapply(names(KEEP_T), load_tissue), names(KEEP_T))

## ---- 4. S1 and S5: every programme in every tissue ----------------------------
S1 <- rbindlist(lapply(names(TIS), function(tk) { T <- TIS[[tk]]; cov <- tech_of(T$d)
  rbindlist(lapply(colnames(T$S), function(pg) { s <- T$S[, pg]; if (all(is.na(s))) return(NULL)
    a <- pspear(s, T$d$prolif, T$d, cov); raw <- suppressWarnings(cor.test(s, T$d$prolif, method = "spearman", exact = FALSE))
    u <- agefit(s, T$d); v <- agefit(s, T$d, "prolif")
    data.table(tissue = tk, programme = pg, n_genes = T$NG[pg], n = a$n,
               rho_prolif = a$r, p_prolif = a$p, rho_prolif_raw = unname(raw$estimate),
               beta_age = u["beta"], p_age = u["p"], beta_age_adj = v["beta"], p_age_adj = v["p"],
               attenuation = 1 - v["beta"]/u["beta"]) })) }))
S1[, FDR_prolif := p.adjust(p_prolif, "BH"), by = tissue]
S1[, FDR_age := p.adjust(p_age, "BH"), by = tissue]
fwrite(S1, file.path(O, "programme_by_tissue.tsv"), sep = "\t")

## ---- 5. primary hypotheses H1-H6 -----------------------------------------------
g <- function(tk, pg, col) S1[tissue == tk & programme == pg][[col]]
fz <- function(r1, n1, k1, r2, n2, k2) {                     # Spearman Fisher z, var 1.06/(n-3)
  se <- sqrt(1.06/(n1 - 3 - k1) + 1.06/(n2 - 3 - k2)); (atanh(r1) - atanh(r2)) / se }
cc <- pspear(TIS$culture$S[, "splicing 96"], TIS$culture$d$prolif, TIS$culture$d, tech_of(TIS$culture$d))
ls <- pspear(TIS$legskin$S[, "splicing 96"], TIS$legskin$d$prolif, TIS$legskin$d, tech_of(TIS$legskin$d))
z2 <- fz(ls$r, ls$n, ls$k, cc$r, cc$n, cc$k)
H <- data.table(
 id = paste0("H", 1:6),
 hypothesis = c("culture: splicing 96 vs proliferation (partial rho)",
                "leg skin coupling weaker than culture (Fisher z, one-sided)",
                "skeletal muscle: splicing 96 vs proliferation (partial rho)",
                "culture: age effect on splicing 96 attenuated by proliferation",
                "culture: collagen formation vs proliferation (partial rho)",
                "leg skin: collagen formation falls with age"),
 estimate = c(cc$r, ls$r - cc$r, g("muscle","splicing 96","rho_prolif"),
              g("culture","splicing 96","attenuation"), g("culture","collagen formation","rho_prolif"),
              g("legskin","collagen formation","beta_age")),
 p = c(cc$p, pnorm(z2), g("muscle","splicing 96","p_prolif"),
       g("culture","splicing 96","p_age"), g("culture","collagen formation","p_prolif"),
       g("legskin","collagen formation","p_age")))
H[, FDR := p.adjust(p, "BH")]
p4 <- g("culture","splicing 96","p_age")
H[, verdict := c(
  ifelse(cc$r >= 0.5, "supported", ifelse(cc$r < 0.3, "not supported", "partial")),
  ifelse(pnorm(z2) < 0.05, "supported", ifelse(ls$r > cc$r, "not supported (opposite)", "not supported")),
  { r <- g("muscle","splicing 96","rho_prolif"); ifelse(r < 0.3, "supported", ifelse(r >= 0.5, "not supported", "indeterminate")) },
  { a <- g("culture","splicing 96","attenuation"); ifelse(p4 >= 0.05, "not evaluable (no unadjusted age effect)",
      ifelse(a >= 0.5, "supported", ifelse(a < 0.25, "not supported", "partial"))) },
  { r <- abs(g("culture","collagen formation","rho_prolif")); ifelse(r < 0.3, "supported", ifelse(r >= 0.5, "not supported", "indeterminate")) },
  { b <- g("legskin","collagen formation","beta_age"); pp <- g("legskin","collagen formation","p_age")
    ifelse(b < 0 & pp < 0.05, "supported", "not supported") })]
H[, p_meaning := c("partial Spearman", "one-sided Fisher z", "partial Spearman",
                    "unadjusted age effect (decides evaluability); attenuation itself has no P",
                    "partial Spearman", "age effect")]
fwrite(H, file.path(O, "primary_hypotheses.tsv"), sep = "\t")
say("\n=== PRIMARY HYPOTHESES (preregistered) ===")
print(H[, .(id, estimate = round(estimate, 3), p = signif(p, 2), FDR = signif(FDR, 2), verdict)], row.names = FALSE)
say("  H2 detail: leg skin rho = %+.3f (n = %d) vs culture rho = %+.3f (n = %d), z = %.2f",
    ls$r, ls$n, cc$r, cc$n, z2)
say("  H4 detail: culture splicing 96 age beta/decade %+.4f (P = %.2g) -> %+.4f (P = %.2g)",
    g("culture","splicing 96","beta_age"), p4, g("culture","splicing 96","beta_age_adj"),
    g("culture","splicing 96","p_age_adj"))

## ---- 6. exploratory S2-S4 --------------------------------------------------------
say("\n=== S1: proliferation coupling by tissue (partial rho) ===")
print(dcast(S1, programme ~ tissue, value.var = "rho_prolif")[order(-culture)], digits = 2)
say("\n=== S5: age effect per decade, unadjusted, by tissue ===")
print(dcast(S1, programme ~ tissue, value.var = "beta_age")[order(culture)], digits = 2)
b1 <- S1[tissue == "legskin" & programme == "collagen formation"]; b2 <- S1[tissue == "pubskin" & programme == "collagen formation"]
say("\nS2 collagen formation age beta: leg (sun-exposed) %+.3f (P = %.2g) vs suprapubic %+.3f (P = %.2g)",
    b1$beta_age, b1$p_age, b2$beta_age, b2$p_age)
dd <- merge(TIS$culture$d[, .(donor, i = .I)], TIS$legskin$d[, .(donor, j = .I)], by = "donor")
S3 <- rbindlist(lapply(colnames(TIS$culture$S), function(pg) {
  a <- TIS$culture$S[dd$i, pg]; b <- TIS$legskin$S[dd$j, pg]
  if (all(is.na(a)) || all(is.na(b))) return(NULL)
  ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE))
  data.table(programme = pg, n_donors = nrow(dd), rho = unname(ct$estimate), p = ct$p.value) }))
S3[, FDR := p.adjust(p, "BH")]
fwrite(S3, file.path(O, "donor_concordance_culture_vs_legskin.tsv"), sep = "\t")
say("\nS3 donors with both culture and leg skin: %d; programmes with FDR < 0.05 concordance: %d of %d",
    nrow(dd), sum(S3$FDR < 0.05), nrow(S3))
for (tk in c("legskin","pubskin")) { T <- TIS[[tk]]
  u <- agefit(T$S[, "splicing 96"], T$d); v <- agefit(T$S[, "splicing 96"], T$d, c("kerat","dfib"))
  say("S4 %s splicing 96 age beta %+.3f (P = %.2g) -> with composition %+.3f (P = %.2g)",
      tk, u["beta"], u["p"], v["beta"], v["p"]) }
saveRDS(lapply(TIS, function(T) list(d = T$d, S = T$S, NG = T$NG)), file.path(O, "tissue_scores.rds"))
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo.txt"))
