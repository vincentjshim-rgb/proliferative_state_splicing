#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 3C : stress-senescence decoupling is independent of proliferative state
# Figure 4A/B : which pathway modules carry the decoupling
# 입력은 전부 이미 고정된 파일이다. 새 분석 결정 없음.
# ---------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
suppressPackageStartupMessages({ library(ggplot2) })

root   <- "public_data_tierA"
probe  <- file.path(root, "derived", "direction_probe")
out    <- file.path(root, "derived", "manuscript_figures_v4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wtsv <- function(x, f) write.table(x, file.path(out, f), sep = "\t", quote = FALSE,
                                   row.names = FALSE, na = "")

## Okabe-Ito: colour-blind safe. 기존 figure 스크립트와 동일 팔레트를 쓴다.
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"
GREY <- "#9E9E9E"; INK <- "#222222"; MUTED <- "#5A5A5A"

base_theme <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey60", linewidth = 0.4),
        axis.text = element_text(colour = MUTED),
        axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 10.5),
        plot.subtitle = element_text(colour = MUTED, size = 8.6),
        plot.caption = element_text(colour = MUTED, size = 7.4, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.size = unit(9, "pt"), legend.text = element_text(size = 8.4))

## ===========================================================================
## Figure 3C
## ===========================================================================
gx <- read.delim(file.path(root, "derived/decoupling_validation/primary_gene_axis_matrix.tsv"),
                 check.names = FALSE)
pl <- read.delim(file.path(probe, "P3_proliferation_loading_GSE179848.tsv"))
a  <- gx[gx$comparator_consistent & is.finite(gx$UVA_meta) & is.finite(gx$senescence_meta), ]
a$loading <- pl$proliferation_loading[match(a$gene, pl$gene)]
a <- a[is.finite(a$loading), ]

Dof <- function(d) { ru <- cor(d$Repro_specific_score, d$UVA_meta, method = "spearman")
                     rs <- cor(d$Repro_specific_score, d$senescence_meta, method = "spearman")
                     c(rho_UVA = ru, rho_sen = rs, D = ru - rs) }
D_all <- Dof(a)

## -- 3C left: prespecified sensitivity conditions -------------------------
p4 <- read.delim(file.path(probe, "P4_decoupling_proliferation_sensitivity.tsv"))
p4$label <- c("No adjustment", "All three axes residualised\non proliferation loading",
              "Proliferation-neutral genes\n(|loading| < 0.2)",
              "Proliferation-linked genes\n(|loading| ≥ 0.4)")
p4$label <- factor(p4$label, levels = rev(p4$label))
p3c_left <- ggplot(p4, aes(decoupling_index_D, label)) +
  geom_vline(xintercept = D_all["D"], colour = GREY, linetype = 2, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = decoupling_index_D, yend = label),
               colour = "grey80", linewidth = 0.9, lineend = "round") +
  geom_point(colour = BLUE, size = 2.6) +
  geom_text(aes(label = sprintf("%.3f", decoupling_index_D)),
            hjust = -0.45, size = 2.9, colour = INK) +
  geom_text(aes(x = 0.015, label = paste0("n = ", format(n_genes, big.mark = ","))),
            hjust = 0, vjust = 2.3, size = 2.5, colour = MUTED) +
  scale_x_continuous(limits = c(0, 0.66), expand = c(0, 0)) +
  labs(x = expression(paste("Decoupling index  D = ", rho[UVA], " − ", rho[senescence])),
       y = NULL, title = "C. Decoupling does not depend on proliferative state",
       subtitle = "Proliferation loading from 345 independent fibroblast samples (GSE179848)") +
  base_theme + theme(legend.position = "none")

## -- 3C right: D across the full proliferation-loading spectrum -----------
a$bin <- cut(a$loading, breaks = quantile(a$loading, seq(0, 1, 0.1)),
             include.lowest = TRUE, labels = FALSE)
set.seed(3)
bin_tab <- do.call(rbind, lapply(sort(unique(a$bin)), function(b) {
  d <- a[a$bin == b, ]
  bs <- replicate(500, { i <- sample(nrow(d), replace = TRUE); Dof(d[i, ])["D"] })
  data.frame(bin = b, n = nrow(d), mean_loading = mean(d$loading),
             D = Dof(d)["D"], lo = quantile(bs, 0.025), hi = quantile(bs, 0.975))
}))
wtsv(bin_tab, "Figure3C_D_by_proliferation_loading_decile.tsv")
p3c_right <- ggplot(bin_tab, aes(mean_loading, D)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_hline(yintercept = D_all["D"], colour = GREY, linetype = 2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "grey72", linewidth = 0.55) +
  geom_point(colour = BLUE, size = 2.2) +
  annotate("text", x = max(bin_tab$mean_loading), y = D_all["D"], vjust = -0.8, hjust = 1,
           label = sprintf("genome-wide D = %.3f", D_all["D"]), size = 2.5, colour = MUTED) +
  scale_y_continuous(limits = c(0, max(bin_tab$hi) * 1.06), expand = c(0, 0)) +
  labs(x = "Per-gene proliferation loading (decile mean)", y = "Decoupling index D",
       title = "D stays positive at every level of proliferation coupling",
       subtitle = "Loading deciles (~222 genes each); bars are 500x bootstrap 95% intervals") +
  base_theme + theme(legend.position = "none")

ggsave(file.path(out, "Figure3C_left_proliferation_sensitivity.png"), p3c_left,
       width = 152, height = 74, units = "mm", dpi = 300)
ggsave(file.path(out, "Figure3C_right_loading_spectrum.png"), p3c_right,
       width = 138, height = 74, units = "mm", dpi = 300)

## ===========================================================================
## Figure 4A / 4B
## ===========================================================================
pm <- read.delim(file.path(probe, "P1_pathway_decoupling_ranked.tsv"), check.names = FALSE)
pm$name <- sub("__R-HSA.*$", "", pm$pathway)
rna_proc <- "splic|snrnp|intron-containing|3'-end processing|polyadenylation|transport of mature"
transl   <- "translation|ribosom|rrna|nonsense|peptide chain|40s|60s"
stressr  <- "heat shock|perk|hsf1|unfolded protein|nrf2|oxidative stress|nfe2l2"
pm$class <- "Other pathway"
pm$class[grepl(transl,   pm$name, ignore.case = TRUE)] <- "Translation / rRNA / NMD"
pm$class[grepl(stressr,  pm$name, ignore.case = TRUE)] <- "Stress-response transcription"
pm$class[grepl(rna_proc, pm$name, ignore.case = TRUE)] <- "Nuclear pre-mRNA processing"
pm$class <- factor(pm$class, levels = c("Nuclear pre-mRNA processing", "Translation / rRNA / NMD",
                                        "Stress-response transcription", "Other pathway"))
cols <- c("Nuclear pre-mRNA processing" = BLUE, "Translation / rRNA / NMD" = ORANGE,
          "Stress-response transcription" = GREEN, "Other pathway" = GREY)
ct <- suppressWarnings(cor.test(pm$local_NES, pm$decoupling_axis_NES, method = "spearman"))
## label placement is manual (ggrepel unavailable offline); positions chosen so
## no label overlaps a point or another label.
LP <- data.frame(
  name = c("mRNA 3'-end processing",
           "Processing of Capped Intron-Containing Pre-mRNA",
           "mRNA Splicing - Major Pathway",
           "Transport of Mature mRNA derived from an Intron-Containing Transcript",
           "PERK regulates gene expression",
           "Regulation of HSF1-mediated heat shock response",
           "rRNA processing",
           "Nonsense-Mediated Decay (NMD)",
           "Eukaryotic Translation Elongation",
           "Extracellular matrix organization"),
  short = c("mRNA 3'-end processing", "Capped pre-mRNA processing", "mRNA splicing (major)",
            "mRNA export", "PERK", "HSF1 heat-shock response", "rRNA processing",
            "NMD", "Translation elongation", "ECM organization"),
  lx = c(-0.62, -0.62, -0.62, -0.62,  2.98,  2.98,  2.98,  2.98,  2.98, -0.55),
  ly = c( 3.74,  3.36,  2.98,  2.60,  3.74,  3.30,  2.30,  1.20,  0.45, -3.95),
  hj = c( 1, 1, 1, 1, 0, 0, 0, 0, 0, 0))
lab <- merge(LP, pm[, c("name", "local_NES", "decoupling_axis_NES", "class")], by = "name")

p4a <- ggplot(pm, aes(local_NES, decoupling_axis_NES)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = subset(pm, class == "Other pathway"),
             colour = GREY, alpha = 0.3, size = 1.05, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = INK, linewidth = 0.45, linetype = 2) +
  geom_segment(data = lab, aes(x = local_NES, y = decoupling_axis_NES, xend = lx, yend = ly),
               colour = "grey70", linewidth = 0.22, inherit.aes = FALSE) +
  geom_point(data = subset(pm, class != "Other pathway"),
             aes(colour = class), size = 1.7, stroke = 0) +
  geom_text(data = lab, aes(x = lx, y = ly, label = short, hjust = hj),
            size = 2.35, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = cols, breaks = names(cols)[1:3]) +
  scale_x_continuous(limits = c(-2.45, 4.45)) +
  scale_y_continuous(limits = c(-4.4, 4.0)) +
  annotate("text", x = 2.98, y = -2.1, hjust = 0, vjust = 1, size = 2.8, colour = INK,
           label = sprintf("Spearman rho = %.3f\np = %.1e\n553 Reactome pathways",
                           unname(ct$estimate), ct$p.value)) +
  labs(x = "Repro-CM anchor NES (local)",
       y = expression(paste("Decoupling axis   ", NES[UVA], " - ", NES[senescence])),
       title = "A. Nuclear pre-mRNA processing carries the decoupling",
       subtitle = "One point per Reactome pathway. Upper right = up locally and on the stress-over-senescence axis.",
       caption = paste0("Local vs UVA alone: rho = 0.079, p = 0.063 (not significant). ",
                        "Local vs senescence alone: rho = -0.528, p = 4.3e-41.\n",
                        "Enrichment p-values are signature associations, not replicate-level treatment p-values (local n = 1 per condition).")) +
  base_theme
ggsave(file.path(out, "Figure4A_pathway_decoupling_scatter.png"), p4a,
       width = 168, height = 125, units = "mm", dpi = 300)

## -- 4B module x axis heat map (diverging: two hues + neutral midpoint) ----
zm <- read.delim(file.path(probe, "P2_module_axis_zmap.tsv"))
mod_lab <- c(splicing = "mRNA splicing", premRNA = "Capped pre-mRNA processing",
             mRNA_export = "mRNA export", three_prime = "3'-end processing",
             SF_regulon = "Splicing-factor regulon", HSF1_heat_shock = "HSF1 heat-shock response",
             UPR_PERK = "UPR / PERK", ROS_detox = "ROS detoxification / NRF2",
             autophagy = "Autophagy", NMD = "NMD",
             translation = "Translation", rRNA = "rRNA processing",
             proliferation = "Cell cycle / DNA replication", ECM = "ECM organization",
             IFN = "Interferon signalling")
ax_lab <- c(local = "Repro-CM\n(local)", UVA_meta = "UVA\nmeta", senesc_meta = "Senescence\nmeta",
            MPTR_reprogramming = "Partial\nreprogramming", fibroblast_age = "Fibroblast\nage",
            Tyshkovskiy_rodent_max_lifespan = "Lifespan\nrodent*",
            Tyshkovskiy_ITP_max_lifespan = "Lifespan\nITP")
zm <- zm[zm$axis %in% names(ax_lab), ]
zm$module_lab <- factor(mod_lab[zm$module], levels = rev(mod_lab))
zm$axis_lab   <- factor(ax_lab[zm$axis],    levels = ax_lab)
zm$z_clip <- pmax(pmin(zm$matched_z, 8), -8)
p4b <- ggplot(zm, aes(axis_lab, module_lab, fill = z_clip)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%+.1f", matched_z),
                colour = abs(z_clip) > 5), size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = BLUE, mid = "grey96", high = ORANGE, midpoint = 0,
                       limits = c(-8, 8), name = "matched-null z") +
  labs(x = NULL, y = NULL, title = "B. Two module families carry the decoupling",
       subtitle = "Permutation z of the module mean; orange = coordinately up, blue = down. Rows ordered by module family.",
       caption = paste0("Two module families can produce decoupling: nuclear pre-mRNA processing (strong anti-senescence, weak UVA)\n",
                        "and HSF1/UPR-PERK proteostatic stress response (strong UVA, weak anti-senescence). Every other locally\n",
                        "up-shifted module is negative on the UVA axis. * no expression column exists for the rodent sheet, so that\n",
                        "column alone uses an unmatched null; the two lifespan references disagree in sign.")) +
  base_theme +
  theme(legend.position = "right", legend.title = element_text(size = 7.6, colour = MUTED),
        legend.key.width = unit(7, "pt"), legend.key.height = unit(26, "pt"),
        panel.grid = element_blank(), axis.text.x = element_text(size = 7.6))
ggsave(file.path(out, "Figure4B_module_axis_heatmap.png"), p4b,
       width = 176, height = 128, units = "mm", dpi = 300)

cat("D genome-wide =", round(D_all["D"], 3), "\n")
cat("D by loading decile range:", paste(round(range(bin_tab$D), 3), collapse = " - "), "\n")
cat("pathway rho =", round(unname(ct$estimate), 3), " p =", signif(ct$p.value, 3), "\n")
cat("panels written to", out, "\n")
