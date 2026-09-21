## DRAFT A  The 177 genes across the 63 contrasts, one cell per gene and contrast.
##
## The meta-analysis form (gene x study heatmap of effect sizes, de Magalhaes 2009, Peters 2015,
## Palmer 2021): every gene in every study, so that the reader sees the machinery move as a
## block with the cell cycle rather than a summary score. Rows are the 177 genes ordered by
## their correlation with the counted division rate; columns are the 63 contrasts ordered by
## the change in cell-cycle genes; the cell is the gene's log2 fold change in that contrast.
## Output: figures_sciadv/draft_A_gene_by_contrast.png (not wired into the manuscript)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
GV   <- fread(file.path(D, "gene_level/gene_vs_counted_rate.tsv"))
R63  <- fread(file.path(D, "revision_stats/fig6_contrasts_revised.tsv"))
NAMED <- readLines("scripts/revision/named_splicing_factors.txt")

## ---- the contrast vectors, assembled exactly as run_stats_supplements.R (C2) does ----
N <- fread(file.path(D, "secretome_class/secretome_signatures.tsv"))
vec <- lapply(split(N, N$id), function(x) setNames(x$logFC, toupper(x$gene)))
for (f in list.files(file.path(D, "compendium/signatures"), full.names = TRUE)) {
  x <- fread(f); vec[[sub("\\.tsv$", "", basename(f))]] <- setNames(x$logFC, toupper(x$gene)) }
lgf <- list(GSE109700_deep = "figure2_public_aging/GSE109700_deep_vs_proliferating_DE.tsv",
            GSE179848_late = "figure2_public_aging/GSE179848_late_vs_early_donor_level_DE.tsv",
            GSE93535_SIPS  = "figure2_public_aging/GSE93535_SIPS_vs_Q_DE.tsv",
            GSE191055_P27  = "figure2_public_aging/GSE191055_P27_vs_P4_DE.tsv",
            GSE109700_early = "figure2_public_aging/GSE109700_early_vs_proliferating_DE.tsv")
for (id in names(lgf)) { f <- file.path(D, lgf[[id]]); if (file.exists(f)) {
  x <- fread(f); vec[[id]] <- setNames(x$logFC, toupper(x$gene)) } }
DE <- fread(file.path(D, "revision_stats/GSE113957_primary_cohort_per_decade_DE.tsv"))
vec[["GSE113957_age"]] <- setNames(DE$logFC, toupper(DE$gene))
stopifnot(all(R63$id %in% names(vec)))
vec <- vec[R63$id]
## the set mean recomputed here must equal the stored contrast score (same definition)
dlt <- function(v, set) { bg <- v[is.finite(v)]; mean(bg[names(bg) %in% set]) - mean(bg) }
chk <- sapply(vec, dlt, set = CORE)
stopifnot(max(abs(chk - R63$spl)) < 1e-6)

## ---- the matrix: gene x contrast, log2 fold change against the contrast's own background ----
## each contrast is centred on its background mean so that a cell reads as "this gene, relative
## to the transcriptome, in this contrast" -- the same centring the set score uses
M <- sapply(vec, function(v) { bg <- v[is.finite(v)]; (v - mean(bg))[CORE] })
rownames(M) <- CORE
rord <- order(GV$rho[match(CORE, GV$gene)], decreasing = TRUE)       # most coupled at the top
M <- M[rord, ]
cord <- order(R63$cc)                                                 # cell cycle falling -> rising
M <- M[, cord]; C <- R63[cord]
cat(sprintf("matrix %d genes x %d contrasts; %.1f%% cells measured\n", nrow(M), ncol(M), 100 * mean(is.finite(M))))
cat(sprintf("column order: cc from %.2f to %.2f; class counts: %s\n", min(C$cc), max(C$cc),
            paste(sprintf("%s %d", names(table(C$class)), table(C$class)), collapse = "; ")))

CC <- c(senescence = RED, `donor age` = BROWN, secretome = TEAL, `photoprotection / rescue` = ORANGE,
        `UV injury` = PURPLE, `metabolic / culture` = GREY, reprogramming = "#2C4A6E")
LIM <- 1.5
H <- as.data.table(as.table(M)); setnames(H, c("gene", "id", "lfc"))
H[, gi := match(gene, rownames(M))]; H[, ci := match(id, colnames(M))]
H[, v := pmax(pmin(lfc, LIM), -LIM)]
NG <- nrow(M); NC <- ncol(M)

## the accession of each contrast; the full contrast names are Supplementary Data 1
acc <- sub("^(GSE[0-9]+).*$", "\\1", C$id); acc[!grepl("^GSE", acc)] <- "secretome"
CT <- data.table(ci = seq_len(NC), cc = C$cc, spl = C$spl, class = C$class, acc = acc)

pH <- ggplot(H, aes(ci, gi)) +
  geom_tile(aes(fill = v), colour = NA) +
  scale_fill_gradientn(colours = c(BLUE, "#F7F7F7", RED), limits = c(-LIM, LIM),
                       breaks = c(-LIM, 0, LIM), labels = c(paste0(MINUS, LIM), "0", LIM), na.value = "#E3E6EA",
                       name = "log₂ fold change of the gene,\nrelative to the transcriptome",
                       guide = guide_colourbar(barwidth = unit(48, "pt"), barheight = unit(4, "pt"),
                                               title.position = "top", ticks.colour = "white")) +
  scale_x_continuous(NULL, breaks = seq_len(NC), labels = CT$acc, expand = c(0, 0)) +
  scale_y_reverse(sprintf("the 177 genes, ordered by %s with the counted division rate", RHO),
                  breaks = NULL, expand = c(0, 0)) +
  theme_sa(8) + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5, colour = INK2),
                      legend.position = "bottom", legend.title = element_text(size = 7, lineheight = 1.0),
                      legend.text = element_text(size = 7), legend.margin = margin(t = 2),
                      plot.margin = margin(2, 3, 3, 3))
## the named factors, labelled at the right with leaders (11 of the 20 belong to the set)
LB <- data.table(gene = intersect(NAMED, rownames(M)))
LB[, gi := match(gene, rownames(M))]
LB <- LB[order(gi)]; LB[, ly := seq(4, NG - 3, length.out = .N)]
pH <- pH + annotate("segment", x = NC + 0.6, xend = NC + 3.2, y = LB$gi, yend = LB$ly, colour = GREY_L, linewidth = 0.25) +
  annotate("text", x = NC + 3.8, y = LB$ly, label = LB$gene, hjust = 0, family = FONT, size = pt(7),
           fontface = "italic", colour = INK2) +
  coord_cartesian(xlim = c(0.5, NC + 14), clip = "off")

## top tracks: the contrast's cell-cycle change (bars) and its class (strip)
pT <- ggplot(CT, aes(ci, cc)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_col(aes(fill = class), width = 0.85) +
  geom_point(aes(y = spl), size = 0.7, colour = INK) +
  scale_fill_manual(values = CC, name = NULL) +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous("change in\ncell-cycle genes", breaks = c(-1, 0), labels = num_axis(0)) +
  coord_cartesian(xlim = c(0.5, NC + 14)) +
  annotate("text", x = NC + 1, y = max(CT$cc), hjust = 0, vjust = 1, family = FONT, size = pt(7), colour = INK2,
           lineheight = 1.05, label = "bars: cell-cycle genes\npoints: the 177 (set mean)") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
  theme_sa(8) + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                      axis.title.y = element_text(size = 7.5, lineheight = 0.95),
                      legend.position = "top", legend.justification = "left", legend.key.size = unit(6, "pt"),
                      legend.text = element_text(size = 7), legend.margin = margin(b = -3),
                      plot.margin = margin(3, 3, 0, 3))
save_fig(plot_grid(pT, pH, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.22, 1)),
         "draft_A_gene_by_contrast.png", 183, 150)
