FIGURES_WRITTEN <- c("Fig5.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/journal_theme.R")
D <- "public_data_tierA/derived"
## the label column of these tables contains an embedded newline, so a single
## record is spread over two rows: rejoin the text before keeping the data rows.
rdclean <- function(f, num) {
  x <- read.delim(f, stringsAsFactors = FALSE)
  keep <- rowSums(is.na(x[, num, drop = FALSE])) == 0
  lab  <- x[[1]]
  for (i in which(keep)) if (i > 1 && !keep[i-1]) lab[i] <- paste(lab[i-1], lab[i])
  x[[1]] <- trimws(lab); x[keep, , drop = FALSE] }

## ============= A. depth-matched differential splicing, GSE109700 ===========
A <- rdclean(file.path(D, "phase2_splicing/Figure5A_GSE109700_contrasts.tsv"), c("n","sig","sig10"))
A$contrast <- sub(" vs ", "\nvs ", A$contrast)
A$lab <- factor(A$contrast, levels = rev(A$contrast))
AL <- data.frame(lab = rep(A$lab, 2),
                 v = c(A$sig - A$sig10, A$sig10),
                 k = rep(c("FDR < 0.05", "and |ΔPSI| ≥ 0.10"), each = nrow(A)))
AL$k <- factor(AL$k, levels = c("and |ΔPSI| ≥ 0.10", "FDR < 0.05"))
p5a <- ggplot(AL, aes(v, lab, fill = k)) +
  geom_col(width = 0.56, colour = NA) +
  geom_text(data = A, aes(sig, lab, label = format(sig, big.mark = ",")), inherit.aes = FALSE,
            hjust = -0.22, size = pt(6.2), family = FONT, colour = INK) +
  scale_fill_manual(values = c(`FDR < 0.05` = "#A9C7DF", `and |ΔPSI| ≥ 0.10` = NAVY),
                    name = NULL, breaks = c("FDR < 0.05","and |ΔPSI| ≥ 0.10")) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE)) +
  scale_x_continuous(limits = c(0, 980), breaks = c(0, 300, 600, 900)) +
  labs(x = "alternative junctions", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_text(size = 6.3, lineheight = 0.9),
        legend.position = "bottom", legend.justification = "right",
        legend.box.background = element_blank(), legend.margin = margin(t = -4))

## ================ B. label-permutation null, GSE93535 (2x2) ================
NB <- read.delim(file.path(D, "phase2_splicing/Figure5B_GSE93535_permutation_null.tsv"))
obs <- 21
p5b <- ggplot(NB, aes(permuted_sig)) +
  geom_histogram(binwidth = 1, fill = "#C3CBD4", colour = "white", linewidth = 0.15) +
  annotate("segment", x = obs, xend = obs, y = 46, yend = 8, colour = SEN, linewidth = 0.45,
           arrow = arrow(length = unit(2.4, "pt"), type = "closed")) +
  annotate("text", x = obs, y = 50, label = "observed\n21", colour = SEN, family = FONT,
           size = pt(6.2), vjust = 0, lineheight = 0.95) +
  annotate("text", x = 13.5, y = 108, hjust = 1, family = FONT, size = pt(6), colour = INK2,
           lineheight = 1.1, label = "298 permutations of the\nsenescence label within\ncompound strata\nmax = 11   P = 0.0033") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_x_continuous(limits = c(-1, 24)) +
  labs(x = "junctions at FDR < 0.05 under permutation", y = "permutations") +
  theme_j(7)

## =========== C. senescence-specific fraction of testable junctions =========
Cc <- rdclean(file.path(D, "phase2_splicing/Figure5C_magnitude.tsv"), c("sig","n","pct"))
Cc$ref <- sub(" \\(", "\n(", Cc$ref)
Cc$lab <- factor(Cc$ref, levels = rev(Cc$ref))
p5c <- ggplot(Cc, aes(pct, lab)) +
  geom_col(width = 0.5, fill = NAVY) +
  geom_text(aes(label = sprintf("%.2f%%  (%s / %s)", pct, format(sig, big.mark=","),
                                format(n, big.mark=","))),
            hjust = -0.06, size = pt(6), family = FONT, colour = INK) +
  scale_x_continuous(limits = c(0, 1.35), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = "% of testable junctions differentially spliced", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.3, lineheight = 0.9))

## ================== D. local intron retention (descriptive) ================
IR <- read.delim(file.path(D, "intron_retention/local_IR_burden.tsv"))
IR$label[IR$sample == "UVA0"] <- "UVA 0 h"
IR$lab <- factor(IR$label, levels = c("UVA 0 h","iPSC-CM","HDF-CM","Repro-CM"))
IR$pc  <- 100 * IR$intronic_read_fraction
IR$col <- c(`Repro-CM` = REPRO, `HDF-CM` = GREY, `iPSC-CM` = NAVY, `UVA 0 h` = UV)[as.character(IR$lab)]
p5d <- ggplot(IR, aes(lab, pc)) +
  geom_col(width = 0.5, fill = IR$col) +
  geom_text(aes(label = sprintf("%.2f", pc)), vjust = -0.6, size = pt(6.2), family = FONT, colour = INK) +
  annotate("text", x = 2.5, y = 6.55, family = FONT, size = pt(5.9), colour = INK2,
           label = "n = 1 per condition; descriptive only, no test") +
  scale_y_continuous(limits = c(0, 7.0), breaks = c(0, 2, 4, 6)) +
  labs(x = NULL, y = "intronic read fraction (%)") +
  theme_j(7)

top <- lab_grid(p5a, p5b, labels = c("A","B"), ncol = 2, rel_widths = c(1.15, 1))
bot <- lab_grid(p5c, p5d, labels = c("C","D"), ncol = 2, rel_widths = c(1.3, 1))
save_fig(lab_grid(top, bot, labels = c("",""), ncol = 1, rel_heights = c(1.1, 1)),
         "Fig5.png", 183, 104)
