#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

base <- "/home/shim/Downloads/project/project_shim/rejuvenation"
out <- "public_data_tierA/derived/local_mirna_audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, filename) write.table(x, file.path(out, filename), sep = "\t",
                                                quote = FALSE, row.names = FALSE, na = "")

files <- c(
  HDF = file.path(base, "mirDeep2/HDF_miRDeep2/miRNAs_expressed_all_samples_HDF.csv"),
  iPSC = file.path(base, "mirDeep2/IPS_miRDeep2/miRNAs_expressed_all_samples_IPS.csv"),
  Repro = file.path(base, "mirDeep2/REP_miRDeep2/miRNAs_expressed_all_samples_REP.csv")
)

# miRDeep2 repeats the mature count for each compatible precursor. Taking the
# maximum per mature miRNA avoids counting the same reads multiple times.
counts <- lapply(files, function(f) {
  x <- read.csv(f, check.names = FALSE)
  aggregate(read_count ~ miRNA, x, max)
})
merged <- Reduce(function(x, y) merge(x, y, by = "miRNA", all = TRUE),
                 lapply(names(counts), function(nm) {
                   x <- counts[[nm]]
                   names(x)[2] <- nm
                   x
                 }))
merged[is.na(merged)] <- 0
for (nm in names(counts)) merged[[paste0(nm, "_CPM")]] <- 1e6 * merged[[nm]] / sum(merged[[nm]])
d_hdf <- log2((merged$Repro_CPM + 0.5) / (merged$HDF_CPM + 0.5))
d_ips <- log2((merged$Repro_CPM + 0.5) / (merged$iPSC_CPM + 0.5))
merged$Repro_specific_score <- ifelse(sign(d_hdf) == sign(d_ips),
                                      sign(d_hdf) * pmin(abs(d_hdf), abs(d_ips)), 0)
write_tsv(merged[order(-abs(merged$Repro_specific_score)), ], "mature_mirna_composition.tsv")

family <- function(x) {
  ifelse(grepl("^hsa-miR-(302[a-z]?|367)(-|$)", x, ignore.case = TRUE), "miR-302/367",
  ifelse(grepl("^hsa-miR-(371[ab]?|372|373)(-|$)", x, ignore.case = TRUE), "miR-371/372/373",
  ifelse(grepl("^hsa-let-7", x, ignore.case = TRUE), "let-7", "other")))
}
merged$family <- family(merged$miRNA)
family_rows <- do.call(rbind, lapply(c("HDF", "iPSC", "Repro"), function(sample) {
  value <- tapply(merged[[paste0(sample, "_CPM")]], merged$family, sum)
  data.frame(sample = sample, family = names(value), CPM = as.numeric(value))
}))
write_tsv(family_rows, "mirna_family_composition.tsv")

qc <- data.frame(
  sample = c("HDF", "iPSC", "Repro"),
  raw_reads = c(17636826, 14499092, 20571973),
  adapter_fraction = c(0.995, 0.981, 0.967),
  post_trim_base_fraction = c(0.661, 0.431, 0.339),
  mirdeep_quantifier_reads = c(17104024, 10072780, 10865219),
  mature_miRNA_mapped_reads = c(124439, 322279, 1166775),
  mature_miRNA_mapping_fraction = c(0.00728, 0.03200, 0.10739),
  biological_replicates = 1
)
write_tsv(qc, "mirna_qc_metrics.tsv")

focus <- family_rows[family_rows$family != "other", ]
p <- ggplot(focus, aes(x = sample, y = CPM / 1e4, fill = family)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("let-7" = "#0072B2", "miR-302/367" = "#D55E00",
                               "miR-371/372/373" = "#CC79A7")) +
  labs(x = NULL, y = "Family abundance (10,000 CPM)", fill = NULL,
       title = "Reprogramming-associated miRNA families mark Repro-CM",
       subtitle = "Composition only; n=1 per CM and unequal mature-miRNA mapping preclude inference") +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out, "Figure_miRNA_family_composition.png"), p, width = 7.2, height = 4.6, dpi = 300)

cat("Local miRNA audit completed:", out, "\n")
