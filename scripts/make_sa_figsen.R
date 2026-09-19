## Supplementary Fig. S8 (six-figure restructure 2026-09-17; was Fig. S7): what a
## senescence panel reports in cultured fibroblasts is predicted by how much of it is
## cell-cycle genes. Outputs from run_senescence_panels.R.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
T <- read.delim("public_data_tierA/derived/senescence_panels/senescence_panel_audit.tsv",
                check.names = FALSE)
SHORT <- c("Reactome Cellular Senescence" = "Reactome senescence", "SenMayo" = "SenMayo",
           "Fridman senescence, up" = "Fridman up", "Fridman senescence, down" = "Fridman down",
           "CellAge, induces senescence" = "CellAge induces",
           "CellAge, inhibits senescence" = "CellAge inhibits",
           "Reactome Cell Cycle, Mitotic (control)" = "cell cycle (control)")
T$lab <- SHORT[T$panel]
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value) }

## ---- a. composition predicts what the panel reports --------------------------
## 2026-09-17: "SenMayo" used to be printed over the Fridman-up point and "Fridman up"
## beside the open symbol of Fridman down; both now sit next to their own points
## across the seven panels: Spearman with the exact permutation P (n = 7, no ties)
ctp <- cor.test(T$pct_cellcycle, T$rho_rate, method = "spearman", exact = TRUE)
s <- list(rho = unname(ctp$estimate), p = ctp$p.value, n = nrow(T))
cat(sprintf("  across %d panels: rho = %.4f, exact two-sided P = %.4f\n", s$n, s$rho, s$p))
## grouped bars (2026-09-19) in place of the dumbbell plot; panels ordered by their
## cell-cycle content, which is printed in the label
TA <- rbind(data.frame(lab = T$lab, k = "all genes", rho = T$rho_rate),
            data.frame(lab = T$lab, k = "without cell-cycle genes", rho = T$rho_rate_noCC))
TA$k <- factor(TA$k, levels = c("all genes", "without cell-cycle genes"))
ordl <- T$lab[order(T$pct_cellcycle)]
lab_of <- function(l) sprintf("%s (%s%%)", l, round(T$pct_cellcycle[match(l, T$lab)]))
TA$labp <- factor(lab_of(TA$lab), levels = lab_of(ordl))
pa <- ggplot(TA, aes(rho, labp, fill = k)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = NA, na.rm = TRUE) +
  scale_fill_manual(values = c(`all genes` = BLUE, `without cell-cycle genes` = "#C9CFD6"), name = NULL) +
  scale_x_continuous(limits = c(-0.28, 0.86), labels = num_axis(1)) +
  labs(x = sprintf("%s with the measured division rate", RHO), y = NULL,
       subtitle = sprintf("label: cell-cycle share of the panel\nacross seven panels %s = %s, %s",
                          RHO, num(s$rho), pfmt(s$p))) +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -4), legend.key.size = unit(6, "pt"),
                     legend.text = element_text(size = 7), axis.text.y = element_text(size = 7.5),
                     plot.subtitle = element_text(size = 7, colour = INK2, lineheight = 1.05))

## ---- b. and whether it looks age-associated ----------------------------------
B <- rbind(
  data.frame(lab = T$lab, model = "unadjusted", b = T$beta, lo = T$lo, hi = T$hi, p = T$FDR),
  data.frame(lab = T$lab, model = "proliferation-adjusted", b = T$beta_adj,
             lo = T$lo_adj, hi = T$hi_adj, p = T$FDR_adj))
B$lab <- factor(B$lab, levels = rev(T$lab[order(-T$pct_cellcycle)]))
B$model <- factor(B$model, levels = c("unadjusted", "proliferation-adjusted"))
B$sig <- B$p < 0.05
pb <- ggplot(B, aes(b, lab, fill = model, group = model)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(aes(alpha = sig), position = position_dodge(width = 0.7), width = 0.62, colour = NA) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.4, colour = INK2,
                 position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c(GREY, BLUE), name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), name = NULL,
                     labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR \u2265 0.05")) +
  scale_x_continuous(labels = num_axis(2)) +
  guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
  labs(x = "age effect per decade, 107 adult donors", y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.direction = "horizontal", legend.box = "vertical",
                     legend.box.just = "left", legend.spacing.y = unit(0, "pt"),
                     legend.margin = margin(0, 0, 1, 0), legend.key.size = unit(6, "pt"),
                     legend.text = element_text(size = 7),
                     axis.text.y = element_text(size = 7.5))

save_fig(lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.06)),
         "FigS6.png", 183, 84)
