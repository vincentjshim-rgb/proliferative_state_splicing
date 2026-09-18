FIGURES_WRITTEN <- c("FigS3.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

## Coverage-based sashimi from the four locally aligned libraries.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages(library(data.table))
D <- "public_data_tierA/derived"; O <- file.path(D, "sashimi_local")
GENE <- "COL1A2"; CHR <- "chr7"
FULL0 <- 94394894L; FULL1 <- 94431227L          # whole locus, used for the coverage read
S0 <- 94403500L;  S1 <- 94420500L               # window shown; widened so that the distal
DISTAL <- 94418498L                             # acceptor discussed in Fig. S1c is inside it
SAMP <- c(UVA0 = "UVA 0 h", HDF = "HDF-CM", REP = "Repro-CM", IPS = "iPSC-CM")
COLS <- c(`UVA 0 h` = "#7B52AB", `HDF-CM` = "#8A9099", `Repro-CM` = "#0E7C7B", `iPSC-CM` = "#28527A")
LIB <- c(UVA0 = 13955400, HDF = 15756799, REP = 14372377, IPS = 15945042)

## ---- coverage, binned for rendering ---------------------------------------
BIN <- 10L
cov <- rbindlist(lapply(names(SAMP), function(s) {
  d <- fread(file.path(O, sprintf("cov_%s.tsv", s)), col.names = c("s","pos","dep"))
  d[, bin := (pos %/% BIN) * BIN]
  d <- d[, .(dep = mean(dep)), by = bin][bin >= S0 & bin <= S1]
  d[, cpm := dep / LIB[s] * 1e6][, samp := SAMP[s]][] }))
cov[, samp := factor(samp, levels = SAMP)]
say <- function(...) cat(sprintf(...), "\n")
say("coverage bins per sample: %d", nrow(cov)/4)

## ---- junctions from the STAR SJ.out.tab -----------------------------------
jl <- rbindlist(lapply(names(SAMP), function(s) {
  j <- fread(file.path(D, "local_bam_v41", paste0(s, ".SJ.out.tab")),
             col.names = c("chr","start","end","strand","motif","annot","uniq","multi","overhang"))
  j <- j[chr == CHR & start >= S0 & end <= S1 & uniq >= 20]
  j[, `:=`(samp = SAMP[s], cpm = uniq / LIB[s] * 1e6)][] }))
jl[, samp := factor(samp, levels = SAMP)]
say("junctions with >=20 unique reads: %s", paste(jl[, .N, by = samp]$N, collapse = ", "))

## keep the junctions present in all four libraries plus the age-associated donor site
AGE_DONOR <- 94406304L
jl[, key := paste(start, end)]
common <- jl[, .N, by = key][N == 4]$key
jl <- jl[key %in% common | start == AGE_DONOR | end == AGE_DONOR]
say("junctions drawn: %d distinct", uniqueN(jl$key))

mx <- max(cov$cpm)
jl[, h := mx * (0.22 + 0.60 * (rank(-cpm, ties.method = "first") %% 3) / 3), by = samp]
arc <- jl[, { t <- seq(0, pi, length.out = 70)
  .(x = start + (end - start) * (1 - cos(t)) / 2, y = sin(t) * h) },
  by = .(samp, key, start, end, uniq, cpm, h)]

## ---- gene model ------------------------------------------------------------
EX <- fread(file.path(D, "ir_annotation_v41/exons_merged_all.bed"),
            col.names = c("chr","s","e"))[chr == CHR & e >= S0 & s <= S1]
eh <- mx * 0.075
p <- ggplot() +
  geom_area(data = cov, aes(bin, cpm, fill = samp), alpha = 0.85, colour = NA) +
  geom_line(data = arc, aes(x, y, group = key, linewidth = cpm), colour = INK, alpha = 0.55,
            lineend = "round") +
  geom_label(data = jl[, .SD[rank(-uniq, ties.method="first") <= 4], by = samp],
             aes(x = (start + end)/2, y = h, label = format(uniq, big.mark = ",")),
             size = pt(5.4), family = FONT, colour = INK, fill = "white",
             label.size = 0, label.padding = unit(0.4, "pt")) +
  geom_rect(data = EX, aes(xmin = s, xmax = e, ymin = -eh*1.9, ymax = -eh*0.5),
            fill = INK, colour = NA, inherit.aes = FALSE) +
  annotate("segment", x = S0, xend = S1, y = -eh*1.2, yend = -eh*1.2,
           colour = INK, linewidth = 0.3) +
  geom_vline(xintercept = AGE_DONOR, colour = "#B2182B", linewidth = 0.4, linetype = "22") +
  geom_vline(xintercept = DISTAL, colour = ORANGE, linewidth = 0.4, linetype = "22") +
  facet_wrap(~ samp, ncol = 1, strip.position = "right") +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_linewidth_continuous(range = c(0.2, 1.9), guide = "none") +
  scale_x_continuous(labels = function(x) sprintf("%.1f kb", (x - S0)/1e3),
                     breaks = seq(S0, S1, by = 2000), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = bquote(italic(.(GENE)) ~ "  " ~ .(sprintf("%s:%s-%s   (position within the window)", CHR,
                   format(S0, big.mark = ","), format(S1, big.mark = ",")))),
       y = "coverage (reads per million)") +
  theme_sa() + theme(strip.text.y.right = element_text(angle = 0, size = 7.5, hjust = 0),
                     panel.spacing = unit(3, "pt"),
                     axis.text.x = element_text(size = 7),
                     plot.margin = margin(16, 4, 3, 3))
## the two sites the donor-cohort analysis singled out, named above the top panel
lab <- data.table(samp = factor(SAMP[1], levels = SAMP),
                  x = c(AGE_DONOR, DISTAL), y = c(mx * 1.30, mx * 1.30),
                  txt = c("donor site chr7:94,406,304",
                          "distal acceptor: no uniquely mapped read\nin any of the four libraries"),
                  col = c("#B2182B", ORANGE), hj = c(0, 1))
p <- p + geom_text(data = lab, aes(x, y, label = txt, hjust = hj), colour = lab$col,
                   family = FONT, size = pt(7), lineheight = 1.05, vjust = 0) +
  coord_cartesian(clip = "off")
save_fig(p, "FigS3.png", 183, 112)

## report the junction usage at the age-associated donor site
say("\njunction usage at the age-associated donor site chr7:%s", format(AGE_DONOR, big.mark=","))
tb <- dcast(jl[start == AGE_DONOR], end ~ samp, value.var = "uniq", fill = 0)
print(tb)
fwrite(jl, file.path(O, "junctions_drawn.tsv"), sep = "\t")
