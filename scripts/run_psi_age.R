## Event-level alternative splicing across 143 dermal fibroblast donors:
## PSI at shared splice sites, tested against donor age with and without proliferation.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(limma); library(Matrix)})
rc <- "public_data_tierA/recount3"; D <- "public_data_tierA/derived"
O  <- file.path(D, "psi_age"); dir.create(O, showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n"); SRP <- "SRP144355"
M <- fread(file.path(D, "outcome_vs_expression/GSE113957_sample_metrics.tsv"))
say("donors with age, proliferation and depth: %d", nrow(M))

ids <- scan(gzfile(file.path(rc, sprintf("sra.junctions.%s.ALL.ID.gz", SRP))), what="", sep="\n", quiet=TRUE)
ids <- unlist(strsplit(ids, "[\t ]+")); ids <- ids[ids != "" & ids != "rail_id"]
RR <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.RR.gz", SRP))))
trip <- fread(cmd = paste("zcat", file.path(rc, sprintf("sra.junctions.%s.ALL.MM.gz", SRP)), "| tail -n +4"),
              col.names = c("i","j","v"))
S <- sparseMatrix(i = trip$i, j = trip$j, x = trip$v, dims = c(nrow(RR), length(ids)))
colnames(S) <- ids; rm(trip); gc()
S <- S[, match(as.character(M$rail), colnames(S))]
say("junction matrix (sparse): %s x %d", format(nrow(S), big.mark=","), ncol(S))
ok <- Matrix::rowSums(S >= 5) >= 0.5*ncol(S)
J <- as.matrix(S[ok, ]); RRk <- RR[ok]; rm(S); gc()
say("junctions passing coverage: %s (%.1f%% annotated)", format(nrow(J), big.mark=","),
    100*mean(RRk$annotated == 1))

psi_events <- function(keycols, label) {
  key <- do.call(paste, c(lapply(keycols, function(k) RRk[[k]]), sep = "_"))
  n <- ave(rep(1L, length(key)), key, FUN = length)
  idx <- which(n >= 2); if (!length(idx)) return(NULL)
  kk  <- key[idx]
  tot <- rowsum(J[idx, , drop = FALSE], kk)
  mm  <- match(kk, rownames(tot))
  psi <- J[idx, , drop = FALSE] / tot[mm, , drop = FALSE]
  good <- rowSums(tot[mm, , drop = FALSE] >= 20) >= 0.8*ncol(J)
  say("  %-22s %s junctions in %s groups -> %s testable events", label,
      format(length(idx), big.mark=","), format(length(unique(kk)), big.mark=","),
      format(sum(good), big.mark=","))
  list(psi = psi[good, , drop = FALSE], ann = RRk[idx][good], type = label) }
say("\nbuilding events:")
E5 <- psi_events(c("chromosome","start","strand"), "alternative acceptor")
E3 <- psi_events(c("chromosome","end","strand"),   "alternative donor")
PSI <- rbind(E5$psi, E3$psi)
ANN <- rbind(E5$ann, E3$ann)
ANN$type <- c(rep(E5$type, nrow(E5$psi)), rep(E3$type, nrow(E3$psi)))
say("total testable events: %s", format(nrow(PSI), big.mark=","))

lgt <- function(x) { x <- pmin(pmax(x, 0.005), 0.995); log(x/(1-x)) }
lg <- lgt(PSI)
run <- function(form, lab) {
  fit <- eBayes(lmFit(lg, model.matrix(form, data = M)))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  say("  %-30s FDR<0.05: %5d of %s", lab, sum(tt$adj.P.Val < 0.05), format(nrow(tt), big.mark=","))
  tt }
say("\nPSI against donor age:")
A0 <- run(~ I(age/10) + log_depth,          "age + depth")
A1 <- run(~ I(age/10) + prolif + log_depth, "age + proliferation + depth")
res <- data.table(ANN[, .(chromosome, start, end, strand, annotated)], type = ANN$type,
                  dlogit = A0$logFC, p = A0$P.Value, FDR = A0$adj.P.Val,
                  dlogit_adj = A1$logFC, p_adj = A1$P.Value, FDR_adj = A1$adj.P.Val,
                  meanPSI = rowMeans(PSI))
say("events retained after proliferation adjustment: %d of %d (%.0f%%)",
    sum(res$FDR_adj < 0.05), sum(res$FDR < 0.05),
    100*sum(res$FDR_adj < 0.05)/max(1, sum(res$FDR < 0.05)))

gb <- fread("public_data_tierA/derived/ir_annotation_v41/genes.bed",
            col.names = c("chr","s","e","gene","score","strand2"))
setkey(gb, chr, s, e)
res[, mid := (start + end) %/% 2L]
q <- res[, .(chr = chromosome, s = mid, e = mid)]
ov <- foverlaps(q, gb, by.x = c("chr","s","e"), type = "within", mult = "first")
res[, gene := ov$gene]
say("events mapped to a gene: %.1f%%", 100*mean(!is.na(res$gene)))
fwrite(res, file.path(O, "psi_age_events.tsv"), sep = "\t")

sig <- res[FDR < 0.05 & !is.na(gene)]
say("\ngenes with an age-associated splicing change: %d", uniqueN(sig$gene))
say("   after proliferation adjustment: %d", uniqueN(res[FDR_adj < 0.05 & !is.na(gene)]$gene))
SF <- c("SRSF1","SRSF2","SRSF3","SRSF5","SRSF6","SRSF7","HNRNPA1","HNRNPA2B1","HNRNPD",
        "HNRNPK","HNRNPM","SF1","SF3B1","SF3A3","SF3B5","U2AF1","U2AF2","PRPF8","PRPF3",
        "SNRPA","SNRPB","SNRPD2","SNRPF","RBM39","TRA2B","PTBP1")
say("   of which are splicing factors themselves: %d", sum(unique(sig$gene) %in% SF))
writeLines(unique(sig$gene), file.path(O, "genes_age_spliced.txt"))
writeLines(unique(res$gene[!is.na(res$gene)]), file.path(O, "genes_universe.txt"))
saveRDS(list(PSI = PSI, ANN = ANN, M = M), file.path(O, "psi_matrix.rds"))
say("\ntop 12 age-associated events:")
print(head(sig[order(p), .(gene, chromosome, start, end, type, annotated,
                           dlogitPSI = round(dlogit, 3), FDR = signif(FDR, 2))], 12))
