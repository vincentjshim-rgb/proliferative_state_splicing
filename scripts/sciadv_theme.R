## Figure system for the npj Aging submission (revised 2026-09-15).
## Layout follows Kerepesi et al. 2021 Sci. Adv. (statistic printed on the panel,
## small multiples, simple schematics); typography follows the npj Aging figure
## rules: lower-case bold panel letters, Greek letters as symbols rather than
## spelled out, a real minus sign, scientific P values as 9 x 10^-4, and nothing
## printed below about 7 pt at the reproduction width.
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(grid)})
FONT <- "Nimbus Sans"                      # metrically Helvetica, as the journal asks
INK <- "#111417"; INK2 <- "#454C54"; GREY <- "#8A9099"; GREY_L <- "#CFD4D9"
BLUE <- "#2C6FAF"; RED <- "#C0392B"; TEAL <- "#0E7C7B"; ORANGE <- "#E08214"
PURPLE <- "#7B52AB"; GREEN <- "#3E8E41"; BROWN <- "#8C6D4F"
DIVERGE <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B")
## verdict colours: blue/orange rather than green/red (red-green contrast is not allowed)
PASS <- "#2C6FAF"; FAIL <- "#E08214"; NEUTRAL <- "#8A9099"
pt <- function(x) x/.pt; mm <- function(x) x/25.4
MINUS <- "−"; RHO <- "ρ"; R2 <- "R²"; TIMES <- "×"

## ---- number and P-value formatting -----------------------------------------
sup <- function(s) { d <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴",
                            "5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹",
                            "-"="⁻","+"="")
  paste(d[strsplit(s, "")[[1]]], collapse = "") }
## -0.12 with a real minus; used for hand-written labels
num <- function(x, digits = 2) sub("-", MINUS, formatC(x, format = "f", digits = digits))
## axis labels with a real minus; rounds first so that ggplot's floating-point
## zero (1.11e-16) prints as 0 rather than in scientific notation
num_axis <- function(digits = 2) function(x) {
  v <- formatC(round(x, digits + 2), format = "f", digits = digits, drop0trailing = TRUE)
  sub("-", MINUS, trimws(v)) }
sci <- function(x, digits = 0) {
  e <- floor(log10(abs(x))); m <- signif(x / 10^e, digits + 1)
  if (isTRUE(all.equal(m, 1))) sprintf("10%s", sup(as.character(e)))
  else sprintf("%s %s 10%s", formatC(m, format = "g"), TIMES, sup(as.character(e))) }
## P as the manuscript prints it: "P = 0.004", "P = 9 x 10^-4", "P < 2 x 10^-16"
pfmt <- function(p, italic = TRUE) {
  lab <- if (italic) "P" else "P"
  if (!is.finite(p)) return("")
  if (p < 2.2e-16) return(sprintf("%s < %s", lab, sci(2e-16)))
  if (p >= 0.001)  return(sprintf("%s = %s", lab, sub("-", MINUS, formatC(signif(p, 2), format = "g"))))
  sprintf("%s = %s", lab, sci(p)) }
## vectorised, for aes(label = ...) inside a data frame
pfmt_v <- function(p, italic = TRUE) vapply(p, pfmt, character(1), italic = italic)
## two lines: statistic and P, as printed on the panels
rp <- function(r, p, lab = RHO, digits = 2)
  sprintf("%s = %s\n%s", lab, num(r, digits), pfmt(p))
## with n
rpn <- function(r, p, n, lab = RHO, digits = 2)
  sprintf("%s = %s\n%s\nn = %d", lab, num(r, digits), pfmt(p), n)
stars <- function(p) ifelse(is.na(p), "", ifelse(p<1e-4,"****", ifelse(p<1e-3,"***",
                    ifelse(p<1e-2,"**", ifelse(p<0.05,"*","ns")))))

## ---- theme ------------------------------------------------------------------
## base 8 pt: axis text 7.5, nothing in a panel should be set below 7
theme_sa <- function(base = 8) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(colour = INK, linewidth = 0.35),
          axis.ticks = element_line(colour = INK, linewidth = 0.35),
          axis.text = element_text(colour = INK, size = base - 0.5),
          axis.title = element_text(colour = INK, size = base),
          strip.background = element_blank(),
          strip.text = element_text(colour = INK, size = base, margin = margin(b = 2)),
          legend.background = element_blank(), legend.key = element_blank(),
          legend.key.size = unit(7, "pt"),
          legend.text = element_text(size = base - 0.5),
          legend.title = element_text(size = base - 0.5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_blank(), plot.margin = margin(3, 4, 3, 3),
          panel.spacing = unit(7, "pt")) }

## significance bracket for bar/box comparisons
sigbar <- function(x1, x2, y, label, tick, size = 7) list(
  annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.3, colour = INK),
  annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, linewidth = 0.3, colour = INK),
  annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, linewidth = 0.3, colour = INK),
  annotate("text", x = (x1+x2)/2, y = y, label = label, vjust = -0.25,
           size = pt(size), family = FONT, colour = INK))

PUBDIR <- "public_data_tierA/derived/figures_sciadv"
dir.create(PUBDIR, showWarnings = FALSE, recursive = TRUE)
## panel letters are lower-case bold, as the journal requires
lab_grid <- function(..., labels, ncol = 2, rel_heights = 1, rel_widths = 1,
                     align = "none", axis = "none")
  plot_grid(..., labels = tolower(labels), ncol = ncol, label_fontfamily = FONT,
            label_fontface = "bold", label_size = 10, label_colour = INK,
            hjust = 0, vjust = 1, label_x = 0.004, label_y = 0.998,
            rel_heights = rel_heights, rel_widths = rel_widths, align = align, axis = axis)
save_fig <- function(plot, file, w, h) { ggsave(file.path(PUBDIR, file), plot,
  width = mm(w), height = mm(h), units = "in", dpi = 600, bg = "white")
  cat(sprintf("  %s  %.0f x %.0f mm\n", file, w, h)) }
