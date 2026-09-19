## Fig. 7. Graphical abstract (2026-09-20; replaces the text-box schematic of 2026-09-19).
## One panel drawn to scale on a 180 x 112 mm canvas. Three columns, each an object of the
## study drawn as a picture with the measurement made on it plotted underneath from the
## deposited tables:
##   1  a culture whose division rate was counted at every passage  -> splicing score vs counted rate
##   2  cultures from donors of different ages                      -> splicing score vs donor age, and
##                                                                     the share of that effect proliferation takes
##   3  skin from the same donor pool, mostly non-cycling cells     -> the splicing set and the cell-cycle
##                                                                     programme against the same score
## Every point, line and printed number comes from the tables; the figure stops before it is
## written if a printed value disagrees with the table it came from.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(ggplot2); library(cowplot); library(data.table)})
D <- "public_data_tierA/derived"
## text extents are measured on a null pdf device before the PNG device exists; register the
## family there with the Helvetica metrics it is built on, so the measurement is exact and silent
if (!FONT %in% names(grDevices::pdfFonts()))
  do.call(grDevices::pdfFonts, setNames(list(grDevices::Type1Font(FONT,
    c("Helvetica.afm", "Helvetica-Bold.afm", "Helvetica-Oblique.afm",
      "Helvetica-BoldOblique.afm", "Symbol.afm"))), FONT))

## ---------------- data -------------------------------------------------------------
DR <- fread(file.path(D, "division_rate/sample_division_rate.tsv"))                   # Fig. 2a
DS <- fread(file.path(D, "machinery_specificity/donor_scores.tsv"))[primary == TRUE]  # Fig. 3a
CR <- fread(file.path(D, "cohort_revised/sensitivity_metrics.tsv"))
GT <- fread(file.path(D, "gtex_boundary/programme_by_tissue.tsv"))
TS <- readRDS(file.path(D, "gtex_boundary/tissue_scores.rds"))
CN <- fread(file.path(D, "gtex_posthoc/sample_counts.tsv"))
stopifnot(nrow(DR) == 328, nrow(DS) == 107)
CC <- "cell cycle (positive control)"
g <- function(tk, pg, col) GT[tissue == tk & programme == pg][[col]]
r_rate <- cor(DR$rate, DR$splice, method = "spearman")
r_age  <- cor(DS$age, DS$splicing)
prim   <- CR[metric == "machinery" & cohort == "primary" & covariates == TRUE]
stopifnot(nrow(prim) == 1)
chk <- function(x, y, tol = 0.006) stopifnot(abs(x - y) < tol)
chk(0.71, r_rate); chk(-0.50, r_age); chk(70, prim$pct_lost, 1)
chk(0.74, g("culture", "splicing 96", "rho_prolif")); chk(0.23, g("legskin", "splicing 96", "rho_prolif"))
chk(0.67, g("legskin", CC, "rho_prolif")); chk(0.94, g("culture", CC, "rho_prolif"))
cat(sprintf("Fig. 7: counted rate %s = %.3f (n = %d) | donor age r = %.3f (n = %d) | %.0f%% via proliferation\n",
            RHO, r_rate, nrow(DR), r_age, nrow(DS), prim$pct_lost))
cat(sprintf("        leg skin: splicing %.3f, cell cycle %.3f (n = %d donors); culture splicing %.3f\n",
            g("legskin", "splicing 96", "rho_prolif"), g("legskin", CC, "rho_prolif"),
            g("legskin", "splicing 96", "n"), g("culture", "splicing 96", "rho_prolif")))

## ---------------- drawing helpers ---------------------------------------------------
POLY <- list(); SEG <- list(); TXT <- list(); PT <- list()
nid <- local({ i <- 0; function() { i <<- i + 1; i } })
poly <- function(x, y, fill, col = NA, lwd = 0.3)
  POLY[[length(POLY) + 1]] <<- data.table(x = x, y = y, id = nid(),
                                          fill = if (is.na(fill)) NA_character_ else fill,
                                          col = if (is.na(col)) NA_character_ else col, lwd = lwd)
ell <- function(cx, cy, rx, ry, ang = 0, fill = "white", col = INK, lwd = 0.3, n = 44) {
  t <- seq(0, 2 * pi, length.out = n); a <- ang * pi / 180
  px <- rx * cos(t); py <- ry * sin(t)
  poly(cx + px * cos(a) - py * sin(a), cy + px * sin(a) + py * cos(a), fill, col, lwd)
}
rect <- function(x0, y0, x1, y1, fill = "white", col = NA, lwd = 0.3)
  poly(c(x0, x1, x1, x0), c(y0, y0, y1, y1), fill, col, lwd)
tri <- function(x, y, dir = 1, s = 1.4, col = INK)
  poly(c(x, x + dir * s, x), c(y + s * 0.62, y, y - s * 0.62), col, NA, 0)
seg <- function(x, y, xend, yend, col = INK, lwd = 0.3, lty = "solid")
  SEG[[length(SEG) + 1]] <<- data.table(x = x, y = y, xend = xend, yend = yend,
                                        col = col, lwd = lwd, lty = lty)
tx <- function(x, y, s, size = 7, hjust = 0.5, vjust = 0.5, face = "plain", col = INK,
               lh = 0.94, ang = 0)
  TXT[[length(TXT) + 1]] <<- data.table(x = x, y = y, lab = s, size = size, hjust = hjust,
                                        vjust = vjust, face = face, col = col, lh = lh, ang = ang)
pts <- function(x, y, col, size = 0.26, alpha = 0.32)
  PT[[length(PT) + 1]] <<- data.table(x = x, y = y, col = col, size = size, alpha = alpha)

## a fibroblast: spindle body, nucleus, two processes
fibro <- function(cx, cy, ang = 0, s = 1, col = INK) {
  a <- ang * pi / 180
  ell(cx, cy, 3.1 * s, 0.95 * s, ang, fill = "white", col = col, lwd = 0.25)
  ell(cx, cy, 0.8 * s, 0.58 * s, ang, fill = col, col = NA, lwd = 0)
  seg(cx + 3.0 * s * cos(a), cy + 3.0 * s * sin(a), cx + 4.4 * s * cos(a), cy + 4.4 * s * sin(a),
      col = col, lwd = 0.25)
  seg(cx - 3.0 * s * cos(a), cy - 3.0 * s * sin(a), cx - 4.4 * s * cos(a), cy - 4.4 * s * sin(a),
      col = col, lwd = 0.25)
}
## a cell caught dividing
mitotic <- function(cx, cy, s = 1, col = RED) {
  ell(cx - 1.0 * s, cy, 1.15 * s, 1.1 * s, fill = "white", col = col, lwd = 0.4)
  ell(cx + 1.0 * s, cy, 1.15 * s, 1.1 * s, fill = "white", col = col, lwd = 0.4)
  ell(cx - 1.0 * s, cy, 0.42 * s, 0.42 * s, fill = col, col = NA, lwd = 0)
  ell(cx + 1.0 * s, cy, 0.42 * s, 0.42 * s, fill = col, col = NA, lwd = 0)
}
dish <- function(cx, cy, r, fill = "#EAF1F8") {
  ell(cx, cy, r, r, fill = fill, col = GREY, lwd = 0.5)
  ell(cx, cy, r - 0.85, r - 0.85, fill = NA, col = GREY_L, lwd = 0.3)
}
scatter_in_dish <- function(cx, cy, ang, s = 0.6) for (a in ang)
  fibro(cx + a[2] * cos(a[1] * pi / 180), cy + a[2] * sin(a[1] * pi / 180), ang = a[3], s = s)

## an inset plot: axis corner, points, least-squares line, end-of-axis values
inset <- function(x0, y0, w, h, dx, dy, col, xlab, xr = NULL, yr = NULL, alpha = 0.32,
                  size = 0.26, lwd = 0.6, xann = TRUE, digits = 1) {
  if (is.null(xr)) xr <- quantile(dx, c(0.005, 0.995), na.rm = TRUE)
  if (is.null(yr)) yr <- quantile(dy, c(0.005, 0.995), na.rm = TRUE)
  sx <- function(v) x0 + (pmin(pmax(v, xr[1]), xr[2]) - xr[1]) / diff(xr) * w
  sy <- function(v) y0 + (pmin(pmax(v, yr[1]), yr[2]) - yr[1]) / diff(yr) * h
  pts(sx(dx), sy(dy), col, size = size, alpha = alpha)
  m <- lm(dy ~ dx); xs <- c(xr[1], xr[2]); ys <- coef(m)[1] + coef(m)[2] * xs
  seg(sx(xs[1]), sy(ys[1]), sx(xs[2]), sy(ys[2]), col = col, lwd = lwd)
  seg(x0, y0, x0 + w, y0, col = INK2, lwd = 0.3)
  seg(x0, y0, x0, y0 + h, col = INK2, lwd = 0.3)
  if (xann) {
    tx(x0, y0 - 2.4, num(xr[1], digits), 7, col = INK2)
    tx(x0 + w, y0 - 2.4, num(xr[2], digits), 7, col = INK2)
    tx(x0 + w / 2, y0 - 6.0, xlab, 7, col = INK2)
  }
  invisible(list(sx = sx, sy = sy, fit = function(v) coef(m)[1] + coef(m)[2] * v))
}

## ---------------- canvas -------------------------------------------------------------
W <- 180; H <- 112
c1 <- 32; c2 <- 90; c3 <- 148; hw <- 26

hd <- function(cx, a, b) { tx(cx, 108.6, a, 8, face = "bold"); tx(cx, 104.7, b, 8, face = "bold")
                           seg(cx - hw, 101.8, cx + hw, 101.8, col = GREY_L, lwd = 0.5) }
hd(c1, "Fibroblasts in a dish,", "division counted at every passage")
hd(c2, "Cultures from 107 donors,", "20 to 96 years old")
hd(c3, "Skin from the same donor pool,", "where most cells do not cycle")

## ---- column 1: one culture, early and late ----------------------------------------
dish(c1 - 13, 89, 9.2)
scatter_in_dish(c1 - 13, 89, list(c(20, 4.0, -20), c(75, 5.6, 40), c(135, 4.3, -35),
                                  c(170, 5.9, 15), c(300, 3.6, -50), c(245, 5.7, 25)))
mitotic(c1 - 15.6, 93.2, s = 0.7); mitotic(c1 - 9.0, 84.8, s = 0.7)
tx(c1 - 13, 77.4, "early passage", 7)
dish(c1 + 13, 89, 9.2)
scatter_in_dish(c1 + 13, 89, list(c(40, 4.6, -30), c(155, 5.2, 20), c(285, 4.4, -10)))
tx(c1 + 13, 77.4, "late passage", 7)
seg(c1 - 2.8, 89, c1 + 2.2, 89, col = INK, lwd = 0.5); tri(c1 + 2.0, 89, 1, 1.5, INK)
tx(c1, 72.8, "every transcriptome carries a division\nrate measured by counting, not inferred", 7, col = INK2)
tx(c1, 66.6, sprintf("%s = %s", RHO, num(r_rate)), 8.5, face = "bold", col = BLUE, hjust = 1)
tx(c1 + 1.5, 66.6, "  328 libraries, 7 cell lines", 7, col = INK2, hjust = 0)
inset(c1 - 20, 40, 40, 23, DR$rate, DR$splice, BLUE, "divisions per day, counted")
tx(c1 - 23.8, 51.5, "splicing-machinery score", 7, col = INK2, ang = 90)
tx(c1, 28.6, "the machinery is transcribed\nin step with how fast the cells grow", 7, col = INK2)

## ---- column 2: young and old donors ------------------------------------------------
dish(c2 - 13, 89, 9.2)
scatter_in_dish(c2 - 13, 89, list(c(25, 5.2, -25), c(95, 5.4, 35), c(165, 5.0, -15),
                                  c(225, 5.3, 20), c(295, 5.1, -40)))
mitotic(c2 - 13, 89, s = 0.7)
tx(c2 - 13, 77.4, "young donor", 7)
dish(c2 + 13, 89, 9.2)
scatter_in_dish(c2 + 13, 89, list(c(60, 5.0, -20), c(200, 5.2, 30), c(320, 4.8, -5)))
tx(c2 + 13, 77.4, "old donor", 7)
seg(c2 - 2.8, 89, c2 + 2.2, 89, col = INK, lwd = 0.5); tri(c2 + 2.0, 89, 1, 1.5, INK)
tx(c2, 72.8, "cells from older donors divide more\nslowly, and score lower for splicing", 7, col = INK2)
tx(c2, 66.6, sprintf("r = %s", num(r_age)), 8.5, face = "bold", col = BLUE, hjust = 1)
tx(c2 + 1.5, 66.6, "  107 adult donors", 7, col = INK2, hjust = 0)
inset(c2 - 20, 47, 40, 16, DS$age, DS$splicing, BLUE, "donor age (years)", digits = 0)
tx(c2 - 23.8, 55, "splicing-machinery score", 7, col = INK2, ang = 90)
f <- prim$pct_lost / 100
rect(c2 - 20, 34.4, c2 - 20 + 40 * f, 39.4, fill = BLUE)
rect(c2 - 20 + 40 * f, 34.4, c2 + 20, 39.4, fill = GREY_L)
tx(c2 - 20 + 20 * f, 36.9, sprintf("proliferation, %.0f%%", prim$pct_lost), 7, col = "white")
tx(c2 + 20 - 20 * (1 - f), 36.9, "residual", 7)
tx(c2, 28.6, "of the age effect on these genes\nwhat is left is not robust, and no age effect\nindependent of proliferation is reported", 7, col = INK2)

## ---- column 3: skin from the same donors -------------------------------------------
SKIN_E <- "#F7E6D2"; SKIN_D <- "#FBF3E8"
rect(c3 - 20, 92.6, c3 + 20, 98.6, fill = SKIN_E, col = GREY_L, lwd = 0.3)
for (i in 0:9) for (j in 0:1)
  ell(c3 - 18.8 + i * 3.05 + j * 1.5, 94.3 + j * 2.6, 1.22, 0.98, fill = "white", col = ORANGE, lwd = 0.25)
rect(c3 - 20, 79.4, c3 + 20, 92.6, fill = SKIN_D, col = GREY_L, lwd = 0.3)
for (i in seq(-18.5, 17, by = 4.4)) seg(c3 + i, 80.4, c3 + i + 3.2, 91.6, col = "#E7D9C3", lwd = 0.4)
FB <- list(c(-14.5, 89.2, -25), c(-4.5, 84.6, 15), c(6.5, 89.6, -10),
           c(15.0, 83.6, 30), c(-9.5, 81.4, -35), c(3.0, 91.2, 5))
for (p in FB) fibro(c3 + p[1], p[2], ang = p[3], s = 0.58)
mitotic(c3 - 18.0, 95.4, s = 0.55)
tx(c3 + 15.5, 96.9, "epidermis", 7, col = INK2)
tx(c3 + 15.0, 81.0, "dermis", 7, col = INK2)
tx(c3, 75.6, sprintf("%d donors gave both a culture and skin",
                     CN[tissue == "culture"]$paired_with_legskin), 7, col = INK2)
tx(c3, 70.6, "what varies between donors is how many\ncells cycle, not how fast each one divides", 7, col = INK2)
## both programmes as standardised residuals of the technical covariates, so that the slope
## drawn for each is the partial correlation quoted for it
sk <- as.data.frame(TS$legskin$d)
sp <- TS$legskin$S[, "splicing 96"]; cy <- TS$legskin$S[, CC]
ok <- complete.cases(sk$prolif, sk$rin, sk$isch, sk$loglib, sp, cy)
sk <- sk[ok, ]; sp <- sp[ok]; cy <- cy[ok]
rz <- function(v) as.numeric(scale(resid(lm(v ~ rin + isch + loglib, sk))))
rx <- rz(sk$prolif); rs <- rz(sp); rc <- rz(cy)
chk(g("legskin", "splicing 96", "rho_prolif"), cor(rx, rs, method = "spearman"), 0.002)
chk(g("legskin", CC, "rho_prolif"), cor(rx, rc, method = "spearman"), 0.002)
XR <- c(-2.6, 2.6); YR <- c(-2.6, 2.6)
## the two programmes in separate boxes of the same size and scale, so that the line each
## draws is its partial correlation and the two can be read against one another
tx(c3, 65.2, sprintf("cell-cycle programme, partial %s = %s", RHO,
                     num(g("legskin", CC, "rho_prolif"))), 7, col = INK2)
inset(c3 - 20, 53.0, 40, 10.5, rx, rc, GREY, "", xr = XR, yr = YR,
      alpha = 0.22, size = 0.2, lwd = 0.7, xann = FALSE)
tx(c3, 50.4, sprintf("splicing machinery, partial %s = %s", RHO,
                     num(g("legskin", "splicing 96", "rho_prolif"))), 7, col = INK2, face = "bold")
inset(c3 - 20, 38.2, 40, 10.5, rx, rs, ORANGE, "proliferation score, standardised residual",
      xr = XR, yr = YR, alpha = 0.26, size = 0.2, lwd = 0.9)
tx(c3 - 23.8, 51, "standardised residual", 7, col = INK2, ang = 90)
tx(c3, 25.9, sprintf("the same score still reads proliferation;\nthe splicing genes no longer follow it\n(%s = %s in cultures from the same donors)",
                     RHO, num(g("culture", "splicing 96", "rho_prolif"))), 7, col = INK2)

## ---- arrows between columns ---------------------------------------------------------
for (x in c(c1 + hw + 0.8, c2 + hw + 0.8)) { seg(x, 89, x + 2.6, 89, col = GREY, lwd = 0.6)
                                             tri(x + 2.4, 89, 1, 1.6, GREY) }

## ---- reading strip -------------------------------------------------------------------
rect(3, 3.5, W - 3, 21, fill = "#EEF1F4")
tx(8, 12.2, "Reading", 8, hjust = 0, face = "bold")
tx(25, 16.6, "In a cultured fibroblast, a change in splicing-factor transcripts is first of all a change in proliferative state.", 7.6, hjust = 0)
tx(25, 12.2, "In skin the same score still recovers the replication programmes, while the splicing genes no longer follow it.", 7.6, hjust = 0)
tx(25, 7.8, "An ageing signature defined in culture is a statement about proliferative state until it is shown in tissue.", 7.6, hjust = 0, col = INK2)

## ---------------- assemble --------------------------------------------------------------
PG <- rbindlist(POLY); SG <- rbindlist(SEG); TT <- rbindlist(TXT); PP <- rbindlist(PT)
PG[, fill2 := fifelse(is.na(fill), NA_character_, fill)]
p <- ggplot() +
  geom_polygon(data = PG, aes(x, y, group = id), fill = PG$fill2, colour = PG$col,
               linewidth = PG$lwd) +
  geom_point(data = PP, aes(x, y), colour = PP$col, size = PP$size, alpha = PP$alpha, shape = 16) +
  geom_segment(data = SG, aes(x = x, y = y, xend = xend, yend = yend), colour = SG$col,
               linewidth = SG$lwd, linetype = SG$lty) +
  geom_text(data = TT, aes(x, y, label = lab), size = pt(TT$size), hjust = TT$hjust,
            vjust = TT$vjust, fontface = TT$face, colour = TT$col, family = FONT,
            lineheight = TT$lh, angle = TT$ang) +
  coord_fixed(ratio = 1, xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(1, 1, 1, 1),
                       plot.background = element_rect(fill = "white", colour = NA))
save_fig(p, "Fig7.png", W + 2, H + 2)
