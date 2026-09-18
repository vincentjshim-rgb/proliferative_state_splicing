.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
G <- "public_data_tierA/gtex"; S <- commandArgs(TRUE)[1]
PROLIF <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
            "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
CORES <- lapply(list(core96 = readLines("public_data_tierA/derived/conserved_core/age_down_splicing_core.txt"),
                     core177 = readLines(file.path(S, "core_without_GSE226189.txt"))), setdiff, PROLIF)
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz", select = c("Geneid", "Gene_name"))
key2 <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
col_ci <- function(d, pat) { k <- grep(pat, names(d), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
read_md <- function(PJ) {
  a <- fread(file.path(G, sprintf("gtex.gtex.%s.MD.gz", PJ))); b <- fread(file.path(G, sprintf("gtex.recount_qc.%s.MD.gz", PJ)))
  k <- intersect(c("rail_id","external_id"), intersect(names(a), names(b))); m <- merge(a, b, by = k, suffixes = c("", ".qc"))
  data.table(project = PJ, rail_id = m$rail_id, external_id = m$external_id, tissue = m[[col_ci(m, "^smtsd$")]],
             age_bin = m[[col_ci(m, "^age$")]], sex = m[[col_ci(m, "^sex$")]], hardy = m[[col_ci(m, "^dthhrdy$")]],
             rin = suppressWarnings(as.numeric(m[[col_ci(m, "^smrin$")]])), isch = suppressWarnings(as.numeric(m[[col_ci(m, "^smtsisch$")]]))) }
MD <- rbind(read_md("SKIN"), read_md("MUSCLE"))
KEEP_T <- c(culture = "Cells - Cultured fibroblasts", legskin = "Skin - Sun Exposed (Lower leg)", muscle = "Muscle - Skeletal")
MD <- MD[tissue %in% KEEP_T]; MD[, tkey := names(KEEP_T)[match(tissue, KEEP_T)]]
MD[, donor := sub("^([^-]+-[^-]+).*$", "\\1", external_id)]
MD[, age := suppressWarnings(as.numeric(sub("^([0-9]+)-.*$", "\\1", age_bin))) + 5]; MD <- MD[is.finite(age)]
MD[, sex := factor(sex)]; MD[, hardy := factor(ifelse(is.na(hardy) | hardy == "", "missing", as.character(hardy)))]
hdr <- function(PJ) names(fread(cmd = sprintf("zcat %s/gtex.gene_sums.%s.G026.gz | grep -v '^##' | head -2", G, PJ)))
tech_of <- function(d) c("rin","isch","loglib")[c(sum(is.finite(d$rin)) > 0.9*nrow(d), sum(is.finite(d$isch)) > 0.9*nrow(d), TRUE)]
pspear <- function(a, b, d, cov) { ok <- is.finite(a) & is.finite(b) & complete.cases(d[, ..cov])
  f <- as.formula(paste("v ~", paste(cov, collapse = "+")))
  ra <- resid(lm(f, data = cbind(d[ok], v = a[ok]))); rb <- resid(lm(f, data = cbind(d[ok], v = b[ok])))
  r <- cor(ra, rb, method = "spearman"); n <- sum(ok); df <- n - 2 - length(cov); c(r = r, p = 2*pt(-abs(r*sqrt(df/(1-r^2))), df), n = n) }
agefit <- function(s, d, extra = NULL) { cov <- c("sex", if (nlevels(droplevels(d$hardy)) > 1) "hardy", tech_of(d), extra)
  dd <- cbind(d, s = s); dd <- dd[is.finite(s) & complete.cases(dd[, ..cov])]
  f <- summary(lm(as.formula(paste("s ~ age +", paste(cov, collapse = "+"))), data = dd))$coef; c(beta = f["age","Estimate"]*10, p = f["age","Pr(>|t|)"]) }
for (tk in names(KEEP_T)) {
  d <- MD[tkey == tk]; PJ <- d$project[1]; hn <- hdr(PJ); cols <- intersect(hn, d$external_id)
  gs <- fread(cmd = sprintf("zcat %s/gtex.gene_sums.%s.G026.gz | grep -v '^##'", G, PJ), select = c(hn[1], cols))
  sym <- key2[sub("\\..*$", "", gs[[1]])]; X <- as.matrix(gs[, -1, with = FALSE]); rm(gs); gc()
  ok <- !is.na(sym) & sym != ""; X <- rowsum(X[ok, , drop = FALSE], sym[ok]); d <- d[match(colnames(X), d$external_id)]
  d[, lib := colSums(X)]; keep <- d$lib >= quantile(d$lib, 0.01); d <- d[keep]; X <- X[, keep, drop = FALSE]
  o <- order(d$donor, -d$lib); drop <- o[duplicated(d$donor[o])]; if (length(drop)) { d <- d[-drop]; X <- X[, -drop, drop = FALSE] }
  y <- DGEList(X); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y); lc <- cpm(y, log = TRUE); rm(X, y); gc()
  z <- t(scale(t(lc[apply(lc, 1, sd) > 0, ]))); rm(lc); gc(); d[, loglib := log10(lib)]
  d[, prolif := colMeans(z[intersect(rownames(z), PROLIF), , drop = FALSE])]
  for (k in names(CORES)) { g <- intersect(rownames(z), CORES[[k]]); s <- colMeans(z[g, , drop = FALSE])
    a <- pspear(s, d$prolif, d, tech_of(d))
    line <- sprintf("%-8s %-8s genes %3d | partial rho vs proliferation = %+.3f (P = %.1e, n = %d)", tk, k, length(g), a["r"], a["p"], a["n"])
    if (tk == "culture") { u <- agefit(s, d); v <- agefit(s, d, "prolif")
      line <- paste0(line, sprintf(" | age/decade %+.4f (P = %.2f) -> adjusted %+.4f (P = %.4f)", u["beta"], u["p"], v["beta"], v["p"])) }
    cat(line, "\n") } }
