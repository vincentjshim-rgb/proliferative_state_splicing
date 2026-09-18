## Figure system for a conventional ageing / molecular-biology journal
## (npj Aging, Aging Cell): Arial/Helvetica metrics, axis line + outward ticks,
## uppercase bold panel letters, boxed legends, significance brackets with
## asterisks, RdBu for expression, subject-grounded categorical hues.
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(grid)})

FONT <- "Nimbus Sans"                      # Helvetica metric-compatible

## ---- palette ---------------------------------------------------------------
INK   <- "#17191C"; INK2 <- "#4B5058"; GREY <- "#8A8F98"; GREY_L <- "#C9CDD3"
PALE  <- "#EDEFF2"
UV    <- "#7B52AB"   # ultraviolet waveband
SEN   <- "#C0392B"   # senescence
REPRO <- "#0E7C7B"   # the intervention
AMBER <- "#D98C1F"   # photoprotection / rescue
NAVY  <- "#28527A"   # proliferation, structural series
RED   <- "#B2182B"; BLUE <- "#2166AC"      # RdBu poles, for expression only
DIVERGE <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B")

sig_stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
  ifelse(p < 1e-2, "**",  ifelse(p < 0.05, "*", "n.s.")))))
fmtp <- function(p) ifelse(p < 1e-16, "italic(P)<10^-16",
  sprintf("italic(P)=='%s'", formatC(p, format = "g", digits = 2)))

pt <- function(x) x / .pt
mm <- function(x) x / 25.4

## ---- theme -----------------------------------------------------------------
theme_j <- function(base = 7) {
  theme_bw(base_size = base, base_family = FONT) +
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.line        = element_line(colour = INK, linewidth = 0.3),
      axis.ticks       = element_line(colour = INK, linewidth = 0.3),
      axis.ticks.length = unit(-1.6, "pt"),          # outward ticks
      axis.text.x      = element_text(colour = INK, size = base - 0.5,
                                      margin = margin(t = 4)),
      axis.text.y      = element_text(colour = INK, size = base - 0.5,
                                      margin = margin(r = 4)),
      axis.title       = element_text(colour = INK, size = base),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.key.size  = unit(7, "pt"),
      legend.text      = element_text(size = base - 0.7, colour = INK),
      legend.title     = element_text(size = base - 0.7, colour = INK),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.background = element_rect(fill = NA, colour = GREY_L, linewidth = 0.25),
      legend.box.margin= margin(1, 2, 1, 2),
      strip.background = element_blank(),
      strip.text       = element_text(colour = INK, size = base, face = "plain",
                                      margin = margin(b = 3)),
      plot.margin      = margin(2, 3, 2, 2),
      plot.title       = element_blank(), plot.subtitle = element_blank())
}
## categorical y axis: drop the y line, keep labels
theme_caty <- function() theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank())

## ---- significance bracket --------------------------------------------------
## horizontal comparison bracket between x1 and x2 at height y
brk <- function(x1, x2, y, label, tick = 0.012, size = 5.6, colour = INK) {
  list(annotate("segment", x = x1, xend = x2, y = y, yend = y,
                linewidth = 0.25, colour = colour),
       annotate("segment", x = x1, xend = x1, y = y, yend = y - tick,
                linewidth = 0.25, colour = colour),
       annotate("segment", x = x2, xend = x2, y = y, yend = y - tick,
                linewidth = 0.25, colour = colour),
       annotate("text", x = (x1 + x2)/2, y = y, label = label, vjust = -0.25,
                size = pt(size), family = FONT, colour = colour))
}
## vertical version, for categorical-y plots
brk_v <- function(y1, y2, x, label, tick = 0.012, size = 5.6, colour = INK) {
  list(annotate("segment", y = y1, yend = y2, x = x, xend = x,
                linewidth = 0.25, colour = colour),
       annotate("segment", y = y1, yend = y1, x = x, xend = x - tick,
                linewidth = 0.25, colour = colour),
       annotate("segment", y = y2, yend = y2, x = x, xend = x - tick,
                linewidth = 0.25, colour = colour),
       annotate("text", y = (y1 + y2)/2, x = x, label = label, hjust = -0.25,
                angle = 90, size = pt(size), family = FONT, colour = colour))
}

## ---- dendrogram helper (no extra packages) ---------------------------------
## returns a data.frame of segments in [0,1] x [1,n] coordinates
dendro_segments <- function(hc) {
  n <- length(hc$order); pos <- order(hc$order)
  xh <- numeric(n - 1); yc <- numeric(n - 1); seg <- list()
  coord <- function(k) if (k < 0) c(x = 0, y = pos[-k]) else c(x = xh[k], y = yc[k])
  for (i in seq_len(n - 1)) {
    a <- coord(hc$merge[i, 1]); b <- coord(hc$merge[i, 2]); h <- hc$height[i]
    xh[i] <- h; yc[i] <- (a["y"] + b["y"]) / 2
    seg[[i]] <- data.frame(x = c(a["x"], h, b["x"]), xend = c(h, h, h),
                           y = c(a["y"], a["y"], b["y"]), yend = c(a["y"], b["y"], b["y"]))
  }
  do.call(rbind, seg)
}

## ---- composition -----------------------------------------------------------
PUBDIR <- "public_data_tierA/derived/figures_journal"
dir.create(PUBDIR, showWarnings = FALSE, recursive = TRUE)

lab_grid <- function(..., labels, ncol = 2, rel_heights = 1, rel_widths = 1,
                     align = "none", axis = "none") {
  plot_grid(..., labels = labels, ncol = ncol, label_fontfamily = FONT,
            label_fontface = "bold", label_size = 9.5, label_colour = INK,
            hjust = 0, vjust = 1, label_x = 0.004, label_y = 0.997,
            rel_heights = rel_heights, rel_widths = rel_widths,
            align = align, axis = axis)
}
save_fig <- function(plot, file, width_mm, height_mm) {
  ggsave(file.path(PUBDIR, file), plot, width = mm(width_mm), height = mm(height_mm),
         units = "in", dpi = 600, bg = "white")
  cat(sprintf("  %s  %.0f x %.0f mm\n", file, width_mm, height_mm))
}
