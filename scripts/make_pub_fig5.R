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
PUBDIR <- "public_data_tierA/derived/figures_publication"
source("scripts/pub_theme.R"); D <- "public_data_tierA/derived"

## a  senescence vs proliferating, depth matched
A <- read.delim(file.path(D, "phase2_splicing/Figure5A_GSE109700_contrasts.tsv"))
A <- A[!is.na(A$sig) & !is.na(A$n), ]   # the contrast field contains newlines; keep complete rows
A$lab <- c("deep senescent\nvs proliferating", "early senescent\nvs proliferating",
           "deep vs early\n(both senescent)")
A$lab <- factor(A$lab, levels = rev(A$lab))
p5a <- ggplot(A, aes(sig, lab)) +
  geom_segment(aes(x = 0, xend = sig, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(size = 1.5, colour = BLUE) +
  geom_text(aes(label = sprintf("%d  (%d at |ΔPSI| ≥ 0.10)", sig, sig10)),
            hjust = -0.12, size = pt(5.2), family = FONT, colour = INK) +
  scale_x_continuous(limits = c(0, 1500), breaks = c(0, 400, 800)) +
  labs(x = "junctions at FDR < 0.05", y = NULL) +
  theme_pub(7) + theme_cat_y() + theme(axis.text.y = element_text(lineheight = 0.92))

## b  quiescence control, permutation null
nb <- read.delim(file.path(D, "phase2_splicing/Figure5B_GSE93535_permutation_null.tsv"))$permuted_sig
obs <- 21
p5b <- ggplot(data.frame(x = nb), aes(x)) +
  geom_histogram(binwidth = 1, fill = "#D6DCE2", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs, colour = VERM, linewidth = 0.45) +
  annotate("text", x = obs, y = Inf, vjust = 1.5, hjust = 1.1, size = pt(5.4), family = FONT,
           colour = VERM, label = "observed 21") +
  annotate("text", x = obs, y = Inf, vjust = 3.1, hjust = 1.1, size = pt(5.2), family = FONT,
           colour = GREY_D, label = "P = 0.0033") +
  scale_x_continuous(breaks = seq(0, 20, 5)) +
  labs(x = "junctions at FDR < 0.05 (senescence main effect)", y = "permutations") +
  theme_pub(7)

## c  proportion
Cd <- read.delim(file.path(D, "phase2_splicing/Figure5C_magnitude.tsv"))
Cd <- Cd[!is.na(Cd$sig) & !is.na(Cd$pct), ]
Cd$lab <- c("vs proliferating", "vs quiescent,\nmatched PD")
Cd$lab <- factor(Cd$lab, levels = rev(Cd$lab))
p5c <- ggplot(Cd, aes(pct, lab)) +
  geom_segment(aes(x = 0, xend = pct, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(size = 1.5, colour = VERM) +
  geom_text(aes(label = sprintf("%.2f%%  (%d / %s)", pct, sig, format(n, big.mark = ","))),
            hjust = -0.12, size = pt(5.2), family = FONT, colour = INK) +
  scale_x_continuous(limits = c(0, 1.35), breaks = c(0, 0.3, 0.6)) +
  labs(x = "% of testable junctions at FDR < 0.05", y = NULL) +
  theme_pub(7) + theme_cat_y() + theme(axis.text.y = element_text(lineheight = 0.92))

## d  local intron retention
b <- read.delim(file.path(D, "intron_retention/local_IR_burden.tsv"))
b$label <- factor(b$label, levels = c("UVA 15J, 0 h", "iPSC-CM", "HDF-CM", "Repro-CM"))
b$g <- ifelse(b$sample == "UVA0", "0 h, no CM", "24 h CM")
p5d <- ggplot(b, aes(100 * intronic_read_fraction, label, colour = g)) +
  geom_segment(aes(x = 0, xend = 100 * intronic_read_fraction, yend = label),
               colour = GREY_L, linewidth = 0.5) +
  geom_point(size = 1.5) +
  geom_text(aes(label = sprintf("%.2f", 100 * intronic_read_fraction)), hjust = -0.4,
            size = pt(5.2), family = FONT, colour = INK) +
  scale_colour_manual(values = c(`24 h CM` = SLATE, `0 h, no CM` = GREY), guide = "none") +
  scale_x_continuous(limits = c(0, 7.2), breaks = c(0, 2, 4, 6)) +
  labs(x = "intronic read fraction (%)", y = NULL) +
  theme_pub(7) + theme_cat_y()

save_fig(compose(p5a, p5b, p5c, p5d, labels = c("a", "b", "c", "d"), ncol = 2,
                 rel_heights = c(1, 0.94)), "Fig5.png", 183, 88)
cat("Fig5 saved\n")
