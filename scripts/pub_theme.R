# ---------------------------------------------------------------------------
# Publication figure conventions (Nature-family):
#   - no title, subtitle or caption inside the artwork; those live in the legend
#   - no panel border, no grid; left and bottom axis lines only
#   - Helvetica-metric sans (Nimbus Sans) at 6-7 pt
#   - direct labelling in preference to legends
#   - lower-case bold panel letters placed outside the plotting area
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(cowplot); library(grid) })

FONT <- "Nimbus Sans"
INK  <- "#1A1D21"; GREY <- "#9AA3AD"; GREY_L <- "#D3D8DE"; GREY_D <- "#5A626B"
BLUE <- "#1F6FB2"; VERM <- "#D55E00"; GREEN <- "#10836A"; GOLD <- "#D9A21B"
PURP <- "#9A6FB0"; SLATE <- "#44505E"

theme_pub <- function(base = 7, grid = "none") {
  g <- element_blank()
  th <- theme_void(base_size = base, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = INK, size = base),
      axis.text = element_text(colour = GREY_D, size = base - 0.5),
      axis.title = element_text(colour = INK, size = base),
      axis.title.x = element_text(margin = margin(t = 3)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 3)),
      axis.line = element_line(colour = INK, linewidth = 0.25, lineend = "square"),
      axis.ticks = element_line(colour = INK, linewidth = 0.25),
      axis.ticks.length = unit(1.4, "pt"),
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.key.size = unit(6, "pt"),
      legend.text = element_text(size = base - 0.5, colour = GREY_D),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      strip.text = element_text(size = base, colour = INK, hjust = 0,
                                margin = margin(b = 2)),
      strip.background = element_blank(),
      plot.margin = margin(2, 3, 2, 2)
    )
  if (grid %in% c("x", "xy"))
    th <- th + theme(panel.grid.major.x = element_line(colour = "#EDF0F2", linewidth = 0.25))
  if (grid %in% c("y", "xy"))
    th <- th + theme(panel.grid.major.y = element_line(colour = "#EDF0F2", linewidth = 0.25))
  th
}
## y-axis category labels: no ticks, no axis line
theme_cat_y <- function(...) theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                                   ...)
lab <- function(p, ...) geom_text(data = p, family = FONT, ...)
pt  <- function(x) x / .pt                       # pt -> ggplot size for text
mm  <- function(x) x / 25.4                      # mm -> inches

save_fig <- function(plot, file, width_mm, height_mm) {
  ggsave(file.path(PUBDIR, file), plot, width = mm(width_mm), height = mm(height_mm),
         units = "in", dpi = 600, bg = "white")
  invisible(NULL)
}
compose <- function(..., labels, ncol = 2, rel_heights = 1, rel_widths = 1,
                    align = "none", axis = "none") {
  plot_grid(..., labels = labels, ncol = ncol, label_fontfamily = FONT,
            label_fontface = "bold", label_size = 9, label_colour = INK,
            hjust = 0, vjust = 1, label_x = 0.002, label_y = 0.998,
            rel_heights = rel_heights, rel_widths = rel_widths,
            align = align, axis = axis)
}
