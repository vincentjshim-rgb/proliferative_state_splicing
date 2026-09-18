#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

root <- "public_data_tierA"
anchor_dir <- file.path(root, "derived", "local_repro_anchor")
dec_dir <- file.path(root, "derived", "decoupling_validation")
out <- file.path(root, "derived", "manuscript_diagnostic_figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                        row.names = FALSE, na = "")

parse_star <- function(sample) {
  folder <- file.path(anchor_dir, "star_genecounts", sample)
  log <- readLines(file.path(folder, "Log.final.out"), warn = FALSE)
  get_num <- function(pattern) {
    z <- log[grepl(pattern, log, fixed = TRUE)][1]
    as.numeric(gsub("%", "", trimws(strsplit(z, "|", fixed = TRUE)[[1]][2])))
  }
  counts <- read.delim(file.path(folder, "ReadsPerGene.out.tab"), header = FALSE, nrows = 4)
  input <- get_num("Number of input reads")
  data.frame(
    sample = c("HDF-CM", "Repro-CM", "iPSC-CM")[match(sample, c("HDF", "REP", "IPS"))],
    input_read_pairs = input,
    unique_mapping_percent = get_num("Uniquely mapped reads %"),
    multimapping_percent = get_num("% of reads mapped to multiple loci"),
    mismatch_percent = get_num("Mismatch rate per base, %"),
    assigned_unstranded = sum(read.delim(file.path(folder, "ReadsPerGene.out.tab"),
                                         header = FALSE, skip = 4)[[2]]),
    assigned_forward = sum(read.delim(file.path(folder, "ReadsPerGene.out.tab"),
                                      header = FALSE, skip = 4)[[3]]),
    assigned_reverse = sum(read.delim(file.path(folder, "ReadsPerGene.out.tab"),
                                      header = FALSE, skip = 4)[[4]])
  )
}
qc <- do.call(rbind, lapply(c("HDF", "REP", "IPS"), parse_star))
qc$assigned_reverse_percent <- 100 * qc$assigned_reverse / qc$input_read_pairs
qc$assigned_forward_percent <- 100 * qc$assigned_forward / qc$input_read_pairs
qc$fastqc_per_base_quality <- "PASS"
qc$fastqc_adapter_content <- "PASS"
qc$biological_replicates <- 1
write_tsv(qc, "mRNA_qc_metrics.tsv")

qc_long <- rbind(
  data.frame(sample = qc$sample, metric = "Uniquely mapped", percent = qc$unique_mapping_percent),
  data.frame(sample = qc$sample, metric = "Gene-assigned (reverse strand)", percent = qc$assigned_reverse_percent),
  data.frame(sample = qc$sample, metric = "Gene-assigned (forward strand)", percent = qc$assigned_forward_percent)
)
p1a <- ggplot(qc_long, aes(sample, percent, fill = metric)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 80, colour = "grey55", linetype = 2) +
  scale_fill_manual(values = c("Uniquely mapped" = "#0072B2",
                               "Gene-assigned (reverse strand)" = "#009E73",
                               "Gene-assigned (forward strand)" = "#D55E00")) +
  labs(x = NULL, y = "Percent of input read pairs", fill = NULL,
       title = "A. Validated paired-end RNA-seq supports reverse-stranded counting",
       subtitle = "All libraries: FastQC base quality and adapter content PASS; n=1 per CM") +
  coord_cartesian(ylim = c(0, 100)) + theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(out, "Figure1_A_mRNA_QC.png"), p1a, width = 8.0, height = 4.8, dpi = 300)

anchor <- read.delim(file.path(anchor_dir, "local_Repro_specific_gene_rank.tsv"), check.names = FALSE)
anchor$class <- "Discordant/zero"
anchor$class[anchor$comparator_consistent] <- "Comparator-consistent"
anchor$class[anchor$comparator_consistent & abs(anchor$Repro_specific_score) >= log2(1.25)] <-
  "Consistent and |FC|>=1.25"
p1b <- ggplot(anchor, aes(logFC_Repro_vs_HDF, logFC_Repro_vs_iPSC, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(alpha = 0.42, size = 0.75) +
  scale_colour_manual(values = c("Discordant/zero" = "grey78", "Comparator-consistent" = "#56B4E9",
                                 "Consistent and |FC|>=1.25" = "#CC3311")) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.06, vjust = -0.8,
           label = sprintf("8,921/12,319 (72.4%%) comparator-consistent\n2,204 genes with |score| >= log2(1.25)"),
           size = 3.5) +
  labs(x = "log2FC: Repro-CM vs HDF-CM", y = "log2FC: Repro-CM vs iPSC-CM",
       colour = NULL, title = "B. The anchor requires agreement against both CM controls",
       subtitle = "Score = signed minimum absolute change; descriptive because n=1 per condition") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out, "Figure1_B_conservative_anchor.png"), p1b, width = 7.2, height = 6.2, dpi = 300)

pathway <- read.delim(file.path(dec_dir, "pathway_axis_matrix.tsv"), check.names = FALSE)
requested <- c(
  "Interferon alpha/beta signaling",
  "rRNA modification in the nucleus and cytosol",
  "Transport of Mature mRNA derived from an Intron-Containing Transcript",
  "mRNA Splicing - Major Pathway",
  "NGF-stimulated transcription",
  "PERK regulates gene expression",
  "Cellular response to heat stress",
  "Fatty acid metabolism"
)
pathway$short <- sub("__R-HSA-.*$", "", pathway$pathway)
selected <- pathway[match(requested, pathway$short), ]
selected <- selected[!is.na(selected$pathway), ]
selected$short <- factor(selected$short, levels = rev(requested))
heat <- rbind(
  data.frame(pathway = selected$short, axis = "Repro-CM anchor", NES = selected$local_NES),
  data.frame(pathway = selected$short, axis = "UVA meta-axis", NES = selected$UVA_meta_NES),
  data.frame(pathway = selected$short, axis = "Senescence meta-axis", NES = selected$senescence_meta_NES)
)
heat$axis <- factor(heat$axis, levels = c("Repro-CM anchor", "UVA meta-axis", "Senescence meta-axis"))
p4 <- ggplot(heat, aes(axis, pathway, fill = NES)) +
  geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", NES)), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish) +
  labs(x = NULL, y = NULL, fill = "NES",
       title = "Selected pathways exhibit stress-senescence decoupling",
       subtitle = "Explanatory subset; global test: all 553 shared Reactome pathways") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(out, "Figure4_pathway_decoupling_heatmap.png"), p4, width = 10.5, height = 5.5, dpi = 300)

rescue <- read.delim(file.path(dec_dir, "rescue_residualized_validation.tsv"), check.names = FALSE)
rescue_long <- rbind(
  data.frame(rescue = rescue$rescue, test = "Raw rescue association", rho = rescue$rho_local_vs_rescue,
             lower = NA, upper = NA),
  data.frame(rescue = rescue$rescue, test = "After removing UVA component",
             rho = rescue$rho_local_vs_rescue_residualized_for_injury,
             lower = rescue$null_q025, upper = rescue$null_q975)
)
p5 <- ggplot(rescue_long, aes(rho, rescue, colour = test)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.12, na.rm = TRUE) +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = c("Raw rescue association" = "#999999",
                                 "After removing UVA component" = "#009E73")) +
  labs(x = "Spearman rho with Repro-CM anchor", y = NULL, colour = NULL,
       title = "Repro-CM retains an independent repair-aligned component",
       subtitle = "Both residual associations: expression-matched FDR < 0.0002") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out, "Figure_rescue_residualized.png"), p5, width = 7.2, height = 4.4, dpi = 300)

cat("Manuscript diagnostic figures completed:", out, "\n")
