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
NUD <- list(`Reactome senescence` = c(0, 0.075), `SenMayo` = c(2.5, -0.105),
            `Fridman up` = c(3.5, -0.005), `Fridman down` = c(0, 0.075),
            `CellAge induces` = c(0, -0.20), `CellAge inhibits` = c(-5, 0.075),
            ## was c(-5, 0.075): the label ended left of its point on the row of "Reactome
            ## senescence" and read as part of it; it now sits under its own point at 100%
            `cell cycle (control)` = c(6, -0.085))
T$nx <- vapply(T$lab, function(k) NUD[[k]][1], 0)
T$ny <- vapply(T$lab, function(k) NUD[[k]][2], 0)
T$hj <- ifelse(T$lab == "cell cycle (control)", 1, ifelse(T$lab == "Fridman up", 0, 0.5))
## across the seven panels: Spearman with the exact permutation P (n = 7, no ties)
ctp <- cor.test(T$pct_cellcycle, T$rho_rate, method = "spearman", exact = TRUE)
s <- list(rho = unname(ctp$estimate), p = ctp$p.value, n = nrow(T))
cat(sprintf("  across %d panels: rho = %.4f, exact two-sided P = %.4f\n", s$n, s$rho, s$p))
pa <- ggplot(T, aes(pct_cellcycle, rho_rate)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_segment(aes(xend = pct_cellcycle, yend = rho_rate_noCC), colour = GREY,
               linewidth = 0.3, na.rm = TRUE) +
  geom_point(aes(y = rho_rate_noCC), shape = 21, size = 1.9, fill = "white",
             colour = INK2, stroke = 0.5, na.rm = TRUE) +
  geom_point(size = 2.2, colour = BLUE) +
  geom_text(aes(pct_cellcycle + nx, rho_rate + ny, label = lab, hjust = hj),
            family = FONT, size = pt(7), colour = INK) +      # 7 pt is the floor
  annotate("label", x = Inf, y = -Inf, hjust = 1.02, vjust = -0.05, family = FONT,
           size = pt(7), colour = INK, fill = "white", label.size = 0,
           label.r = unit(0, "pt"), lineheight = 1.1,
           label = sprintf("open symbols: cell-cycle genes removed\nacross the seven panels %s = %s\n%s (exact), n = %d panels",
                           RHO, num(s$rho), pfmt(s$p), s$n)) +
  scale_x_continuous(limits = c(-6, 108), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(-0.28, 0.86), labels = num_axis(1)) +
  labs(x = "genes of the panel that are also cell-cycle genes (%)",
       y = sprintf("%s with the measured division rate", RHO)) + theme_sa()

## ---- b. and whether it looks age-associated ----------------------------------
B <- rbind(
  data.frame(lab = T$lab, model = "unadjusted", b = T$beta, lo = T$lo, hi = T$hi, p = T$FDR),
  data.frame(lab = T$lab, model = "proliferation-adjusted", b = T$beta_adj,
             lo = T$lo_adj, hi = T$hi_adj, p = T$FDR_adj))
B$lab <- factor(B$lab, levels = rev(T$lab[order(-T$pct_cellcycle)]))
B$model <- factor(B$model, levels = c("unadjusted", "proliferation-adjusted"))
B$sig <- B$p < 0.05
pb <- ggplot(B, aes(b, lab, colour = model)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5,
                 position = position_dodge(width = 0.6)) +
  geom_point(aes(fill = ifelse(sig, ifelse(model == "unadjusted", GREY, BLUE), "white")),
             shape = 21, size = 2, stroke = 0.5, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(GREY, BLUE), name = NULL) +
  scale_fill_identity(guide = "none") +
  scale_x_continuous(labels = num_axis(2)) +
  labs(x = "age effect per decade, 107 adult donors", y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.direction = "horizontal",
                     legend.margin = margin(0, 0, 1, 0), legend.key.size = unit(6, "pt"),
                     legend.text = element_text(size = 7),
                     axis.text.y = element_text(size = 7.5))

save_fig(lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.06)),
         "FigS8.png", 183, 76)
