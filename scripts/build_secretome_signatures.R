## Build the secretome-class contrasts fixed in repro_cm_preregistration_secretome_class_ko.md
## (SHA-256 55bf016d...), plus one extra UV-injury contrast from GSE134533.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(edgeR))
F   <- "public_data_tierA/secretome/files"
OUT <- "public_data_tierA/derived/secretome_class"; dir.create(OUT, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
sigs <- list()

## logFC from raw counts: TMM + logCPM group-mean difference (project convention)
lfc <- function(mat, a, b, min_cpm = 1) {
  mat <- mat[, c(a, b), drop = FALSE]
  d <- DGEList(mat); d <- calcNormFactors(d, method = "TMM")
  cp <- cpm(d); keep <- rowSums(cp >= min_cpm) >= min(2, ncol(mat))
  l  <- log2(cp[keep, , drop = FALSE] + 0.5)
  rowMeans(l[, a, drop = FALSE]) - rowMeans(l[, b, drop = FALSE]) }
add <- function(id, dataset, contrast, v, recipient) {
  v <- v[is.finite(v)]; v <- v[names(v) != "" & !is.na(names(v))]
  v <- tapply(v, names(v), function(x) x[which.max(abs(x))])
  sigs[[id]] <<- data.frame(id, dataset, contrast, recipient,
                            gene = toupper(names(v)), logFC = as.numeric(v))
  say("  %-8s %-11s n_genes = %5d   %s", id, dataset, length(v), contrast) }

## ---- SC1, SC2  GSE139563  IMR90 + senescent-fibroblast CM -----------------
x <- read.delim(gzfile(file.path(F, "GSE139563_Processed_data_Paracrine_Senescence.txt.gz")),
                skip = 1, check.names = FALSE)
cn <- c("ensembl","symbol", paste0(rep(c("D0r1","D0r2","D4r1","D4r2","D7r1","D7r2","D10r1","D10r2"),
        each = 2), c("_counts","_RPM")))
names(x)[seq_along(cn)] <- cn
cm <- as.matrix(x[, grep("_counts$", names(x))]); rownames(cm) <- x$symbol
colnames(cm) <- sub("_counts$", "", colnames(cm))
cm <- cm[!is.na(rownames(cm)) & rownames(cm) != "", ]
add("SC1", "GSE139563", "senescent-cell CM, day 4 vs day 0", lfc(cm, c("D4r1","D4r2"), c("D0r1","D0r2")), "IMR90")
add("SC2", "GSE139563", "senescent-cell CM, day 10 vs day 0", lfc(cm, c("D10r1","D10r2"), c("D0r1","D0r2")), "IMR90")

## ---- SC3 + fractions  GSE251807  NHDF + bone-marrow MSC fractions ---------
y <- read.csv(gzfile(file.path(F, "GSE251807_normalised_counts.csv.gz")), row.names = 1, check.names = FALSE)
ens <- sub("\\..*$", "", rownames(y))
## symbol map harvested from the other files in this batch (no extra download)
mp <- list()
z <- read.delim(gzfile(file.path(F, "GSE282054_raw_counts.txt.gz")), check.names = FALSE)
mp[[1]] <- data.frame(e = sub("\\..*$","",z$Geneid), s = z$Gene_name)
z2 <- read.delim(gzfile(file.path(F, "GSE279804_gene_expression_count.txt.gz")), check.names = FALSE)
mp[[2]] <- data.frame(e = sub("\\..*$","",z2$gene_id), s = z2$gene_name)
z3 <- read.delim(gzfile(file.path(F, "GSE268248_ProcessedDataMatrix.txt.gz")), check.names = FALSE)
mp[[3]] <- data.frame(e = sub("\\..*$","",z3$ENSG), s = z3$GeneSymbol)
MAP <- do.call(rbind, mp); MAP <- MAP[!duplicated(MAP$e) & MAP$s != "" & !is.na(MAP$s), ]
sym <- MAP$s[match(ens, MAP$e)]
say("GSE251807 symbol mapping: %.1f%% of %d genes", 100*mean(!is.na(sym)), length(ens))
ok <- !is.na(sym); ym <- as.matrix(y[ok, ]); rownames(ym) <- sym[ok]
## expression filter matched to the CPM >= 1 rule used for the count-based contrasts
libsz <- median(colSums(ym)); thr <- libsz / 1e6
keepy <- rowSums(ym >= thr) >= 4
say("GSE251807 expression filter: %d of %d genes kept (threshold %.1f normalised counts)",
    sum(keepy), nrow(ym), thr)
ym <- ym[keepy, , drop = FALSE]
gr <- function(p) grep(paste0("^", p, "_"), colnames(ym), value = TRUE)
lfc_norm <- function(a, b) {                       # already normalised counts
  l <- log2(ym + 0.5); rowMeans(l[, a, drop=FALSE]) - rowMeans(l[, b, drop=FALSE]) }
add("SC3",  "GSE251807", "BM-MSC conditioned medium vs DMEM",   lfc_norm(gr("CM"),   gr("DMEM")), "NHDF")
add("SC3b", "GSE251807", "BM-MSC EV-depleted medium vs DMEM",   lfc_norm(gr("DM"),   gr("DMEM")), "NHDF")
add("SC3c", "GSE251807", "BM-MSC small EVs vs DMEM",            lfc_norm(gr("sEV"),  gr("DMEM")), "NHDF")
add("SC3d", "GSE251807", "BM-MSC non-small EVs vs DMEM",        lfc_norm(gr("NsEV"), gr("DMEM")), "NHDF")

## ---- SC4, SC5  GSE266052  HCA2 hTERT fibroblast + adipose MSC-CM ----------
w <- read.delim(gzfile(file.path(F, "GSE266052_Brizio_Fibroblasts_hMSC_CM_raw_count_matrix.txt.gz")),
                row.names = 1, check.names = FALSE)
wm <- as.matrix(w)
add("SC4", "GSE266052", "cytokine-primed hMSC-CM vs untreated (both TGF-beta)",
    lfc(wm, paste0("P", 1:5), paste0("M", 1:5)), "HCA2 hTERT fibroblast")
add("SC5", "GSE266052", "resting hMSC-CM vs untreated (both TGF-beta)",
    lfc(wm, paste0("R", 1:5), paste0("M", 1:5)), "HCA2 hTERT fibroblast")

## ---- SC6  GSE279804  skin fibroblast + ADSC nanovesicles ------------------
v <- read.delim(gzfile(file.path(F, "GSE279804_gene_expression_count.txt.gz")), check.names = FALSE)
vm <- as.matrix(v[, grep("_count$", names(v))]); rownames(vm) <- v$gene_name
vm <- vm[!is.na(rownames(vm)) & rownames(vm) != "", ]
add("SC6", "GSE279804", "ADSC nanovesicles vs control",
    lfc(vm, grep("^Treatment", colnames(vm), value = TRUE),
            grep("^Control",   colnames(vm), value = TRUE)), "skin fibroblast")

## ---- SC7  GSE282054  WI-38 + hTSC-CM --------------------------------------
zm <- as.matrix(z[, grep("CM", names(z))]); rownames(zm) <- z$Gene_name
zm <- zm[!is.na(rownames(zm)) & rownames(zm) != "", ]
add("SC7", "GSE282054", "hTSC secretome vs ES-CM",
    lfc(zm, grep("hTSC", colnames(zm), value = TRUE),
            grep("^ES-CM", colnames(zm), value = TRUE)), "WI-38")

## ---- SC8  GSE293186  fibroblast + endothelial EVs -------------------------
u <- read.csv(gzfile(file.path(F, "GSE293186_gene_count.csv.gz")), check.names = FALSE)
um <- as.matrix(u[, grep("^(CTRL|ECEV)_72h", names(u))]); rownames(um) <- u$gene_name
um <- um[!is.na(rownames(um)) & rownames(um) != "", ]
add("SC8", "GSE293186", "endothelial-cell EVs vs control, 72 h",
    lfc(um, grep("^ECEV", colnames(um), value = TRUE),
            grep("^CTRL", colnames(um), value = TRUE)), "fibroblast")

## ---- SC9, SC10  GSE306748  primary HDF + MSC secretome (DE tables) --------
td <- file.path(OUT, "GSE306748"); dir.create(td, showWarnings = FALSE)
untar(file.path(F, "GSE306748_RAW.tar"), exdir = td)
rd306 <- function(pat) {
  f <- list.files(td, pattern = pat, full.names = TRUE)[1]
  d <- read.csv(gzfile(f), check.names = FALSE)
  lg <- grep("log2FoldChange|logFC", names(d), value = TRUE)[1]
  gn <- grep("^Gene_name$|^gene_symbol$|^symbol$", names(d), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(gn)) gn <- grep("^gene", names(d), value = TRUE, ignore.case = TRUE)[1]
  v <- setNames(d[[lg]], d[[gn]]); v[!is.na(names(v)) & names(v) != ""] }
add("SC9",  "GSE306748", "immortalised MSC secretome vs control", rd306("MSCimmort_vs_C"), "primary HDF")
add("SC10", "GSE306748", "primary MSC secretome vs control",      rd306("MSCwt_vs_C"),     "primary HDF")

## ---- SC11  GSE268248  gingival fibroblast + PRF serum (boundary) ----------
gm <- as.matrix(z3[, grep("^LP", names(z3))]); rownames(gm) <- z3$GeneSymbol
gm <- gm[!is.na(rownames(gm)) & rownames(gm) != "", ]
add("SC11", "GSE268248", "PRF serum vs untreated (gingival, boundary)",
    lfc(gm, c("LP03003","LP03007","LP03011"), c("LP03001","LP03005","LP03009")), "gingival fibroblast")

## ---- extra UV-injury contrast  GSE134533  HDF + UVB -----------------------
tu <- file.path(OUT, "GSE134533"); dir.create(tu, showWarnings = FALSE)
untar(file.path(F, "GSE134533_RAW.tar"), exdir = tu)
hf <- list.files(tu, pattern = "HDF", full.names = TRUE)
say("GSE134533 HDF files: %s", paste(basename(hf), collapse = ", "))
## columns are Id(RefSeq) / Symbol / Count / FPKM / TPM with an in-file header row
rdh <- function(f) { d <- read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE)
  d <- d[d$V1 != "Id", c("V2","V3")]; names(d) <- c("g","c")
  d$c <- suppressWarnings(as.numeric(d$c)); d <- d[is.finite(d$c) & d$g != "", ]
  tapply(d$c, d$g, sum) }
lst <- lapply(hf, rdh); names(lst) <- sub("^GSM[0-9]+_", "", sub("\\.hg19.*$", "", basename(hf)))
gg  <- Reduce(intersect, lapply(lst, names))
hm  <- do.call(cbind, lapply(lst, function(v) v[gg])); rownames(hm) <- gg
add("UVX", "GSE134533", "HDF UVB vs no UV (2 donors)",
    lfc(hm, grep("UVB", colnames(hm), value = TRUE),
            grep("NOUV", colnames(hm), value = TRUE)), "primary HDF")

S <- do.call(rbind, sigs)
write.table(S, file.path(OUT, "secretome_signatures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("\ncontrasts built: %d  (rows %s)", length(sigs), format(nrow(S), big.mark = ","))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
