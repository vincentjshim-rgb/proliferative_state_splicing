## In the same 143 donors: does splicing OUTCOME survive proliferation adjustment
## when splicing MACHINERY EXPRESSION does not?
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(edgeR)})
rc <- "public_data_tierA/recount3"; D <- "public_data_tierA/derived"
O  <- file.path(D, "outcome_vs_expression"); dir.create(O, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
SRP <- "SRP144355"                              # GSE113957, dermal fibroblast donors
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded

## ---- sample table -----------------------------------------------------------
md <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.sra.%s.MD.gz", SRP))))
mp <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.recount_project.%s.MD.gz", SRP))))
say("recount3 samples: %d", nrow(mp))
ttl <- md$sample_title
## titles look like "14_1yr_Male_Asian" or "108_31_female_Caucasian": the age is the
## SECOND number, optionally followed by yr/YR
tl  <- tolower(ttl)
age <- suppressWarnings(as.numeric(sub("^[0-9]+_([0-9]+)(yr)?_.*$", "\\1", tl)))
mos <- grepl("mos", tl)                      # a few paediatric samples are given in months
age[mos] <- suppressWarnings(as.numeric(sub("^[0-9]+_([0-9]+)yr.*$", "\\1", tl[mos])))
say("titles e.g.: %s | ages parsed: %d", paste(head(ttl,2), collapse=" ; "), sum(is.finite(age)))
S <- data.table(rail = as.character(md$rail_id), srr = md$external_id,
                title = ttl, age = age)

## ---- junction matrix --------------------------------------------------------
ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))),
            what = "", sep = "\n", quiet = TRUE)
ids <- unlist(strsplit(ids, "[\t ]+")); ids <- ids[ids != "" & ids != "rail_id"]
RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))))
trip <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)),
                          "| tail -n +4"), col.names = c("i","j","v"))
say("junctions %s | samples %d | nonzero %s", format(nrow(RR), big.mark=","), length(ids),
    format(nrow(trip), big.mark=","))
J <- matrix(0, nrow(RR), length(ids)); J[cbind(trip$i, trip$j)] <- trip$v
colnames(J) <- ids
ann <- RR[[which(names(RR) == "annotated")[1]]] == 1
say("annotated junctions: %d (%.1f%%)", sum(ann), 100*mean(ann))

## ---- splicing OUTCOME metrics ----------------------------------------------
depth <- colSums(J)
nz    <- colSums(J > 0)
## entropy of junction usage within 5' splice-site groups: higher = less precise
cn <- names(RR); scol <- cn[grep("^start$|^chromosome_start$|^left", cn)][1]
ccol <- cn[grep("^chromosome$|^chr", cn)][1]
key <- paste(RR[[ccol]], RR[[scol]])
grp <- split(seq_len(nrow(J)), key)
grp <- grp[lengths(grp) >= 2]
say("5' splice sites with >= 2 junctions: %s", format(length(grp), big.mark=","))
ent <- sapply(seq_len(ncol(J)), function(s) {
  v <- vapply(grp, function(ix) { x <- J[ix, s]; t <- sum(x)
    if (t < 20) return(NA_real_); p <- x[x > 0]/t; -sum(p*log(p)) }, numeric(1))
  mean(v, na.rm = TRUE) })
M <- data.table(rail = colnames(J), depth = depth,
                unannot_reads = colSums(J[!ann, ])/depth,
                unannot_junc  = colSums(J[!ann, ] > 0)/nz,
                entropy = ent)
M <- merge(M, S[, .(rail, srr, age)], by = "rail")
M <- M[is.finite(age)]
say("donors with age and junctions: %d (range %.0f-%.0f)", nrow(M), min(M$age), max(M$age))

## ---- machinery EXPRESSION and proliferation, same samples -------------------
gs <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.gene_sums.%s.G026.gz", SRP)),
                        "| tail -n +2"))
gid <- gs[[1]]; GM <- as.matrix(gs[, -1, with = FALSE])
rownames(GM) <- sub("\\..*$", "", gid)
say("gene sums: %d x %d; columns e.g. %s", nrow(GM), ncol(GM), paste(head(colnames(GM),2), collapse=", "))
## map ENSG -> symbol with the map already in the project
mpf <- fread(cmd = "zcat public_data_tierA/secretome/files/GSE282054_raw_counts.txt.gz")
key2 <- setNames(toupper(mpf$Gene_name), sub("\\..*$", "", mpf$Geneid))
rownames(GM) <- key2[rownames(GM)]
GM <- GM[!is.na(rownames(GM)) & rownames(GM) != "", ]
GM <- rowsum(GM, rownames(GM))
idcol <- if (length(intersect(colnames(GM), M$srr)) > 10) "srr" else "rail"
common <- intersect(colnames(GM), M[[idcol]])
say("joining gene sums on %s", idcol)
say("samples shared between junctions and gene sums: %d", length(common))
GM <- GM[, common, drop = FALSE]; M <- M[match(common, M[[idcol]])]
y <- DGEList(GM); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); lc <- cpm(y, log = TRUE)
z <- t(scale(t(lc[apply(lc,1,sd) > 0, ])))
M[, machinery := colMeans(z[intersect(rownames(z), CORE), , drop = FALSE])]
PR <- readLines(textConnection("")); rm(PR)
prolif_genes <- c("MKI67","CCNB1","CCNA2","CDK1","TOP2A","BIRC5","BUB1","PLK1","AURKA",
                  "TYMS","RRM2","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TK1","UBE2C")
M[, prolif := colMeans(z[intersect(rownames(z), prolif_genes), , drop = FALSE])]
M[, log_depth := log10(depth)]
say("proliferation meta-gene from %d markers; cor(prolif, age) = %+.3f",
    length(intersect(rownames(z), prolif_genes)), cor(M$prolif, M$age))
say("cor(machinery, prolif) = %+.3f", cor(M$machinery, M$prolif))

## ---- the decisive comparison ------------------------------------------------
test <- function(v) {
  f0 <- lm(M[[v]] ~ I(age/10) + log_depth, data = M)
  f1 <- lm(M[[v]] ~ I(age/10) + prolif + log_depth, data = M)
  c0 <- summary(f0)$coef["I(age/10)", ]; c1 <- summary(f1)$coef["I(age/10)", ]
  data.frame(metric = v, r_prolif = cor(M[[v]], M$prolif),
             beta = c0[1], p = c0[4], beta_adj = c1[1], p_adj = c1[4],
             pct_lost = 100*(1 - c1[1]/c0[1])) }
RES <- do.call(rbind, lapply(c("machinery","unannot_reads","unannot_junc","entropy"), test))
say("\n=== same %d donors: age effect before and after proliferation adjustment ===", nrow(M))
print(RES, digits = 3, row.names = FALSE)
write.table(RES, file.path(O, "outcome_vs_expression_models.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
fwrite(M, file.path(O, "GSE113957_sample_metrics.tsv"), sep = "\t")
