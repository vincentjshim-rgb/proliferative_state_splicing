## Sample-level behaviour of the age-associated pre-mRNA processing core
## in the three richest datasets: two donor-age cohorts and one longitudinal series.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(edgeR))
D <- "public_data_tierA/derived"; O <- file.path(D, "core_depth"); dir.create(O, showWarnings=FALSE)
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
say("splicing core: %d genes", length(CORE))

modscore <- function(mat, set) {                    # mean z across the module
  mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
  keep <- apply(mat, 1, sd) > 0; mat <- mat[keep, , drop = FALSE]
  z <- t(scale(t(mat)))
  g <- intersect(rownames(z), set)
  list(score = colMeans(z[g, , drop = FALSE]), n = length(g)) }

## ============ 1. GSE113957 — 104 donors, age encoded in the column names ====
x <- read.delim(gzfile("public_data_tierA/aging/GSE113957_fpkm.txt.gz"), check.names = FALSE)
ann <- x[["Annotation/Divergence"]]
sym <- toupper(sub("\\|.*$", "", sub("^.*?\\|", "", ifelse(grepl("\\|", ann), ann, ann))))
sym <- toupper(sub("[;|].*$", "", ann))
## two naming conventions are present: "101_19yr_Female_..." and "108_31_female_..."
ANN <- c("Transcript ID","chr","start","end","strand","Length","Copies","Annotation/Divergence")
sc  <- setdiff(names(x), ANN)
age <- suppressWarnings(as.numeric(sub("^[0-9]+_([0-9]+)(yr|YR)?_.*$", "\\1", sc)))
sex <- tolower(sub("^[0-9]+_[0-9]+(yr|YR)?_([A-Za-z]+)_.*$", "\\2", sc))
say("GSE113957 columns parsed: %d, ages %d-%d, NA %d",
    length(sc), min(age, na.rm=TRUE), max(age, na.rm=TRUE), sum(is.na(age)))
ok113 <- is.finite(age) & sex %in% c("male","female")
sc <- sc[ok113]; age <- age[ok113]; sex <- sex[ok113]
m   <- as.matrix(x[, sc]); rownames(m) <- sym
m   <- m[!is.na(rownames(m)) & rownames(m) != "", ]
m   <- m[rowMeans(m > 1) > 0.5, ]
lm113 <- log2(m + 1)
ms <- modscore(lm113, CORE)
d1 <- data.frame(sample = sc, age = age, sex = sex, score = ms$score, study = "GSE113957")
d1 <- d1[d1$age >= 22 & d1$age <= 89, ]                   # locked adult range
f1 <- lm(score ~ age + sex, d1)
say("\nGSE113957: %d adult donors, %d core genes measured", nrow(d1), ms$n)
say("  beta(age) = %+.5f per year, t = %.2f, P = %.2e  | rho = %+.3f",
    coef(f1)["age"], summary(f1)$coef["age","t value"], summary(f1)$coef["age","Pr(>|t|)"],
    cor(d1$age, d1$score, method = "spearman"))

## ============ 2. GSE226189 — 82 donors ======================================
td <- file.path(O, "GSE226189"); dir.create(td, showWarnings = FALSE)
if (!length(list.files(td))) untar("public_data_tierA/aging/GSE226189_RAW.tar", exdir = td)
ff <- list.files(td, pattern = "geneCOUNT", full.names = TRUE)
say("\nGSE226189 files: %d  e.g. %s", length(ff), basename(ff[1]))
rd <- function(f) { d <- read.delim(gzfile(f), header = TRUE, stringsAsFactors = FALSE)
  d <- d[!grepl("^__|^N_", d[[1]]), ]
  v <- suppressWarnings(as.numeric(d[[2]])); names(v) <- sub("\\..*$", "", d[[1]])
  v[is.finite(v)] }
L <- lapply(ff, rd); names(L) <- sub("^(GSM[0-9]+)_.*$", "\\1", basename(ff))
gg <- Reduce(intersect, lapply(L, names))
m2 <- do.call(cbind, lapply(L, function(v) v[gg])); rownames(m2) <- gg
## age and sex are encoded in the deposited file names (SKIN_AGE22_M_...)
bn <- basename(ff)
meta2 <- data.frame(gsm = sub("^(GSM[0-9]+)_.*$", "\\1", bn),
                    age = as.numeric(sub("^.*_AGE([0-9]+)_.*$", "\\1", bn)),
                    sex = sub("^.*_AGE[0-9]+_([MF])_.*$", "\\1", bn))
meta2 <- meta2[match(colnames(m2), meta2$gsm), ]
say("GSE226189 ages parsed: %d, range %d-%d", sum(is.finite(meta2$age)),
    min(meta2$age, na.rm=TRUE), max(meta2$age, na.rm=TRUE))
dge <- DGEList(m2); dge <- calcNormFactors(dge, "TMM")
cp <- cpm(dge); m2 <- log2(cp[rowSums(cp >= 1) >= ncol(cp)*0.5, ] + 0.5)
## ENSG -> symbol via the GSE113957 symbols is not possible; use the local anchor map
rk <- read.delim(file.path(D, "local_repro_anchor/local_Repro_specific_gene_rank.tsv"))
if (all(grepl("^ENSG", rownames(m2)[1:5]))) {
  mp <- read.delim(gzfile("public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz"))
  key <- setNames(mp$Gene_name, sub("\\..*$","",mp$Geneid))
  rownames(m2) <- toupper(key[sub("\\..*$","",rownames(m2))]) }
m2 <- m2[!is.na(rownames(m2)) & rownames(m2) != "", ]
ms2 <- modscore(m2, CORE)
d2 <- data.frame(sample = colnames(m2), age = meta2$age, sex = meta2$sex,
                 score = ms2$score, study = "GSE226189")
f2 <- lm(score ~ age + sex, d2)
say("GSE226189: %d donors, %d core genes measured", nrow(d2), ms2$n)
say("  beta(age) = %+.5f per year, t = %.2f, P = %.2e  | rho = %+.3f",
    coef(f2)["age"], summary(f2)$coef["age","t value"], summary(f2)$coef["age","Pr(>|t|)"],
    cor(d2$age, d2$score, method = "spearman"))

AGE <- rbind(d1[, c("sample","age","sex","score","study")], d2[, c("sample","age","sex","score","study")])
write.table(AGE, file.path(O, "donor_age_core_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(capture.output(list(GSE113957 = summary(f1)$coef, GSE226189 = summary(f2)$coef)),
           file.path(O, "donor_age_models.txt"))
