## Supplementary Fig. S7: is any of it specific to splicing, and does the answer survive
## when the splicing definition is not the selected one? The 177-gene set was selected for
## falling in three series (one of them late versus early passage in the counted resource),
## the comparators were not. Each question is therefore shown for the selected set, for two
## whole, unselected Reactome definitions, and for the pathway genes left outside the 177.
## Outputs from scripts/revision/run_machinery_specificity.R.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
D <- "public_data_tierA/derived/machinery_specificity"
Q1 <- read.delim(file.path(D, "counted_libraries_partials.tsv"))
Q2 <- read.delim(file.path(D, "donor_age_models.tsv"))
Q3 <- read.delim(file.path(D, "gtex_culture_partials.tsv"))

CMP <- c(translation = "translation", metabolism_of_rna = "Metabolism of RNA",
         rrna_processing = "rRNA processing")
SHOW <- c("selected_177", "reactome_mrna_splicing", "reactome_premrna_processing",
          "premrna_processing_not_selected")
STRIP <- c(selected_177 = "selected\n177-gene set",
           reactome_mrna_splicing = "unselected\nReactome mRNA Splicing",
           reactome_premrna_processing = "unselected, Reactome\npre-mRNA processing",
           premrna_processing_not_selected = "pre-mRNA processing genes\noutside the 177")
SIDE <- c("splicing definition, other set held constant",
          "other set, splicing definition held constant")
LEG <- theme(legend.position = "top", legend.direction = "horizontal",
             legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
             legend.key.size = unit(7, "pt"), legend.text = element_text(size = 7.5),
             axis.text.y = element_text(size = 7.5),
             strip.text = element_text(size = 7.5, lineheight = 0.95, margin = margin(b = 3)))

## ---- a, b. who keeps the coupling, by splicing definition ---------------------
partial_panel <- function(Q, xlab, xlim, breaks, legend) {
  Q <- Q[Q$splicing_definition %in% SHOW, ]
  A <- rbind(
    data.frame(def = Q$splicing_definition, comparator = Q$comparator, side = SIDE[1],
               rho = Q$partial_splicing, lo = Q$partial_splicing_lo, hi = Q$partial_splicing_hi),
    data.frame(def = Q$splicing_definition, comparator = Q$comparator, side = SIDE[2],
               rho = Q$partial_comparator, lo = Q$partial_comparator_lo, hi = Q$partial_comparator_hi))
  U <- unique(data.frame(def = Q$splicing_definition, rho = Q$rho_splicing_whole_set))
  lab <- setNames(sprintf("%s\nunadjusted %s = %s", STRIP[U$def], RHO, num(U$rho)), U$def)
  A$def <- factor(A$def, levels = SHOW); U$def <- factor(U$def, levels = SHOW)
  A$comparator <- factor(CMP[A$comparator], levels = rev(CMP))
  A$side <- factor(A$side, levels = SIDE)
  ggplot(A, aes(rho, comparator, colour = side)) +
    geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
    geom_vline(data = U, aes(xintercept = rho), colour = INK, linewidth = 0.3, linetype = "22") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5,
                   position = position_dodge(width = 0.55)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.55)) +
    facet_wrap(~def, nrow = 1, labeller = as_labeller(lab)) +
    scale_colour_manual(values = c(BLUE, ORANGE), name = NULL) +
    scale_x_continuous(limits = xlim, breaks = breaks, labels = num_axis(2)) +
    labs(x = xlab, y = NULL) + theme_sa() + LEG +
    theme(legend.position = if (legend) "top" else "none")
}
pa <- partial_panel(Q1, sprintf("partial %s with the measured division rate, %d counted libraries (dashed: the splicing definition unadjusted)",
                                RHO, Q1$n[1]), c(-0.36, 0.80), seq(-0.25, 0.75, 0.25), TRUE)
pb <- partial_panel(Q3, sprintf("partial %s with the proliferation score, %d GTEx fibroblast cultures, one per donor (dashed: the splicing definition unadjusted)",
                                RHO, Q3$n[1]), c(-0.66, 0.86), c(-0.5, 0, 0.5), FALSE)

## ---- c. the donor cohort: which set carries the residual age effect ----------
CDEF <- c(selected_177 = "selected 177-gene set",
          reactome_mrna_splicing = "unselected Reactome mRNA Splicing")
ROW <- c(splicing = "splicing definition", translation = "translation",
         metabolism_of_rna = "Metabolism of RNA")
C <- Q2[Q2$splicing_definition %in% names(CDEF) & Q2$metric %in% names(ROW), ]
## for the splicing definition the other set is translation
C <- C[!(C$side == "splicing" & grepl("Metabolism of RNA|rRNA processing", C$model)), ]
MOD <- c("age, technical covariates", "+ proliferation score", "+ proliferation + the other set")
C$model2 <- factor(ifelse(C$model == "technical covariates only", MOD[1],
                   ifelse(C$model == "+ proliferation score", MOD[2], MOD[3])), levels = MOD)
C$metric <- factor(ROW[C$metric], levels = rev(ROW))
C$def <- factor(CDEF[C$splicing_definition], levels = CDEF)
pc <- ggplot(C, aes(beta, metric, colour = model2)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.66)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.66)) +
  ## P values in one column at the right, clear of the zero line and of the other intervals
  geom_text(aes(x = 0.043, label = pfmt_v(p)), hjust = 0, family = FONT, size = pt(7),
            position = position_dodge(width = 0.66), show.legend = FALSE) +
  facet_wrap(~def, nrow = 1) +
  scale_colour_manual(values = c(GREY, BLUE, TEAL), name = NULL) +
  scale_x_continuous(limits = c(-0.19, 0.098), breaks = seq(-0.15, 0.05, 0.05), labels = num_axis(2)) +
  labs(x = sprintf("age effect per decade (s.d. of the score), %d adult donors", C$n_donors[1]), y = NULL) +
  theme_sa() + LEG

save_fig(lab_grid(pa, pb, pc, labels = c("a", "b", "c"), ncol = 1, rel_heights = c(1.08, 1, 1.22)),
         "FigS5.png", 183, 158)
