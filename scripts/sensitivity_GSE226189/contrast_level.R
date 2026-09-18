.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages(library(fgsea))
D <- "public_data_tierA/derived"; S <- commandArgs(TRUE)[1]
CORES <- list(core96 = readLines(file.path(D, "conserved_core/age_down_splicing_core.txt")),
              core177 = readLines(file.path(S, "core_without_GSE226189.txt")))
td <- tempfile(); dir.create(td); unzip("public_data_tierA/network/ReactomePathways.gmt.zip", exdir = td)
P <- gmtPathways(list.files(td, pattern = "\\.gmt$", full.names = TRUE)[1])
CC <- unique(unlist(P[grep("^Cell Cycle, Mitotic|^DNA Replication|^M Phase", names(P))]))
Z <- readRDS(file.path(D, "benchmark_figs/figdata.rds")); M <- Z$meta
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
AX <- read.delim(file.path(D, "decoupling_validation/primary_gene_axis_matrix.tsv"))
vec <- lapply(split(N, N$id), function(d) setNames(d$logFC, d$gene))
for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
  d <- read.delim(f); vec[[sub("\\.tsv$", "", basename(f))]] <- setNames(d$logFC, toupper(d$gene)) }
lgf <- list(GSE109700_deep = "figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
 GSE179848_late = "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
 GSE113957_age = "figure2_public_aging/GSE113957_age22_89_per_decade_DE.tsv",
 GSE226189_age = "figure2_public_aging/GSE226189_age_per_decade_DE.tsv",
 GSE93535_SIPS = "figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
 GSE191055_P27 = "figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
 GSE109700_early = "figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv")
for (id in names(lgf)) { f <- file.path(D, lgf[[id]]); if (!file.exists(f)) next
  d <- read.delim(f); vec[[id]] <- setNames(d$logFC, toupper(d$gene)) }
loc <- AX[AX$comparator_consistent %in% c(TRUE, "TRUE") & AX$Repro_specific_score != 0, ]
vec[["LOCAL_ReproCM"]] <- setNames(loc$Repro_specific_score, toupper(loc$gene))
dlt <- function(v, set) { g <- intersect(names(v), set); if (length(g) < 15) return(NA)
  bg <- v[is.finite(v)]; mean(bg[names(bg) %in% set]) - mean(bg) }
FIX <- c(GSE116968_uv_vs_control_1h = "UV injury", GSE116968_uv_vs_control_4h = "UV injury", GSE179848_late = "senescence")
run <- function(core, label, drop = character(0)) {
  R <- do.call(rbind, lapply(names(vec), function(i) data.frame(id = i, spl = dlt(vec[[i]], core), cc = dlt(vec[[i]], CC))))
  R$class <- M$class[match(R$id, M$id)]; hit <- R$id %in% names(FIX); R$class[hit] <- FIX[R$id[hit]]
  R$class[is.na(R$class)] <- "metabolic / culture"
  R <- R[is.finite(R$spl) & is.finite(R$cc) & !(R$id %in% drop), ]
  f <- lm(spl ~ cc, R); R$resid <- resid(f); ct <- cor.test(R$cc, R$spl, method = "spearman", exact = FALSE)
  cls <- do.call(rbind, lapply(sort(unique(R$class)), function(k) { v <- R$resid[R$class == k]
    if (length(v) < 3) return(data.frame(class = k, n = length(v), p = NA)); data.frame(class = k, n = length(v), p = t.test(v)$p.value) }))
  cls$FDR <- NA; ok <- !is.na(cls$p); cls$FDR[ok] <- p.adjust(cls$p[ok], "BH")
  wc <- function(k) { s <- R[R$class == k, ]; c <- cor.test(s$cc, s$spl, method = "spearman", exact = FALSE); sprintf("%.2f (P = %.3f, n = %d)", c$estimate, c$p.value, nrow(s)) }
  cat(sprintf("\n== %s: %d contrasts | rho = %.3f | R2 = %.3f | slope %.3f\n", label, nrow(R), ct$estimate, summary(f)$r.squared, coef(f)[2]))
  cat("   class tests (n, BH-FDR):", paste(sprintf("%s %d %s", cls$class, cls$n, ifelse(is.na(cls$FDR), "n<3", sprintf("%.2f", cls$FDR))), collapse = " | "), "\n")
  cat("   reprogramming:", wc("reprogramming"), "| secretome:", wc("secretome"), "\n") }
for (k in names(CORES)) { run(CORES[[k]], paste(k, "all")); run(CORES[[k]], paste(k, "without GSE226189 contrast"), "GSE226189_age") }
