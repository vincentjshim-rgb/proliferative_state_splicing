## Reprogramming and secretome interventions: do they restore splicing-factor
## expression, and is that restoration separable from restored proliferation?
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
D <- "public_data_tierA/derived"; O <- file.path(D, "reprog_splicing"); dir.create(O, showWarnings=FALSE)
say <- function(...) cat(sprintf(...), "\n")
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))   # v2: GSE226189 excluded
SF <- c("SRSF1","SRSF2","SRSF3","SRSF5","SRSF6","SRSF7","HNRNPA1","HNRNPA2B1","HNRNPD",
        "HNRNPK","HNRNPM","SF1","SF3B1","SF3A3","SF3B5","U2AF1","U2AF2","PRPF8","PRPF3",
        "SNRPA","SNRPB","SNRPD2","SNRPF","RBM39","TRA2B","PTBP1")
R <- read.delim(file.path(D, "conserved_core/splicing_vs_proliferation.tsv"))
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
R$class[is.na(R$class)] <- "metabolic / culture"
R$name <- N$contrast[match(R$id, N$id)]
R$name[is.na(R$name)] <- CP$contrast[match(R$id[is.na(R$name)], CP$id)]
R$name[R$id == "LOCAL_ReproCM"] <- "reprogramming-phase conditioned medium"

f <- lm(spl ~ cc, R); R$resid <- resid(f); R$pred <- predict(f)
say("=== reprogramming and secretome interventions ===")
S <- R[R$class %in% c("reprogramming","secretome","Repro-CM"), ]
S <- S[order(-S$spl), ]
print(data.frame(class = substr(S$class,1,13), cellcycle = round(S$cc,3),
                 splicing = round(S$spl,3), residual = round(S$resid,3),
                 name = substr(S$name, 1, 42)), row.names = FALSE)

for (k in c("reprogramming","secretome")) {
  s <- R[R$class == k, ]
  ct <- cor.test(s$cc, s$spl, method = "spearman", exact = FALSE)
  say("\n%s (n = %d): rho(cell cycle, splicing factors) = %+.2f, P = %.4f",
      k, nrow(s), ct$estimate, ct$p.value)
  say("   median splicing-factor change %+.3f, median residual %+.3f",
      median(s$spl), median(s$resid)) }

## how large is the reprogramming effect relative to the ageing loss?
sen <- R[R$class == "senescence", ]; age <- R[R$class == "donor age", ]
say("\nfor scale: senescence %+.3f, donor age %+.3f, reprogramming %+.3f, secretome %+.3f",
    median(sen$spl), median(age$spl), median(R$spl[R$class=="reprogramming"]),
    median(R$spl[R$class=="secretome"]))
write.table(S, file.path(O, "reprog_secretome_splicing.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- named splicing factors in the reprogramming contrasts -----------------
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (fl in list.files(file.path(D,"compendium/signatures"), full.names=TRUE)) {
  d <- read.delim(fl); vec[[sub("\\.tsv$","",basename(fl))]] <- setNames(d$logFC, toupper(d$gene)) }
AX <- read.delim(file.path(D,"decoupling_validation/primary_gene_axis_matrix.tsv"))
loc <- AX[AX$comparator_consistent %in% c(TRUE,"TRUE") & AX$Repro_specific_score != 0, ]
vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))
for (fl in c(GSE165177_MPTR="figure2_mptr/GSE165177_MPTR_paired_DE.tsv")) {
  p <- file.path(D, fl); if (file.exists(p)) { d <- read.delim(p)
    vec[["GSE165177_MPTR"]] <- setNames(d$logFC, toupper(d$gene)) } }
rid <- S$id[S$class %in% c("reprogramming","Repro-CM")]
M <- sapply(rid, function(i) { v <- vec[[i]]; if (is.null(v)) return(rep(NA, length(SF)))
  v[match(SF, names(v))] })
rownames(M) <- SF
say("\n=== named splicing factors across reprogramming contrasts (log2FC) ===")
keep <- rowSums(is.finite(M)) >= ncol(M)*0.6
print(round(M[keep, , drop=FALSE], 2))
write.table(data.frame(gene=rownames(M), M, check.names=FALSE),
            file.path(O, "named_SF_reprogramming.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
