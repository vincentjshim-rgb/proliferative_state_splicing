## Sashimi-style junction arc plot drawn directly from the recount3 junction counts:
## exon structure from GENCODE, arcs scaled by mean junction reads, young vs old.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages(library(data.table))
D <- "public_data_tierA/derived"; O <- file.path(D, "psi_age")
X <- readRDS(file.path(O, "psi_matrix.rds")); PSI <- X$PSI; ANN <- X$ANN; M <- X$M
res <- fread(file.path(O, "psi_age_events.tsv"))
EX <- fread("public_data_tierA/derived/intron_retention/exonic_unique.gsorted.bed",
            col.names = c("chr","s","e","gene"))

sashimi <- function(gene, ttl = NULL) {
  ex <- EX[gene == gene][order(s)]
  if (!nrow(ex)) return(NULL)
  rg <- range(c(ex$s, ex$e))
  jn <- res[gene == gene & chromosome == ex$chr[1] & start >= rg[1] - 2000 & end <= rg[2] + 2000]
  if (!nrow(jn)) return(NULL)
  ## young and old thirds
  yo <- M$age <= quantile(M$age, 1/3); ol <- M$age >= quantile(M$age, 2/3)
  key <- paste(ANN$chromosome, ANN$start, ANN$end, sep = "_")
  jn[, k := paste(chromosome, start, end, sep = "_")]
  mi <- match(jn$k, key)
  jn[, psi_y := rowMeans(PSI[mi, yo, drop = FALSE], na.rm = TRUE)]
  jn[, psi_o := rowMeans(PSI[mi, ol, drop = FALSE], na.rm = TRUE)]
  jn <- jn[is.finite(psi_y) & is.finite(psi_o)]
  if (!nrow(jn)) return(NULL)
  arcs <- rbind(data.table(jn[, .(start, end, psi = psi_y)], grp = "young tertile"),
                data.table(jn[, .(start, end, psi = psi_o)], grp = "old tertile"))
  arcs[, grp := factor(grp, levels = c("young tertile","old tertile"))]
  arcs[, id := .I]
  curve <- arcs[, { t <- seq(0, pi, length.out = 60)
      .(x = start + (end - start)*(1 - cos(t))/2, y = sin(t)*psi) }, by = .(id, grp, psi)]
  ggplot() +
    geom_rect(data = ex, aes(xmin = s, xmax = e, ymin = -0.13, ymax = -0.02),
              fill = INK, colour = NA) +
    geom_segment(data = data.frame(x = rg[1], xend = rg[2], y = -0.075, yend = -0.075),
                 aes(x, y, xend = xend, yend = yend), colour = INK, linewidth = 0.3) +
    geom_line(data = curve, aes(x, y, group = id, linewidth = psi, colour = grp),
              lineend = "round", alpha = 0.85) +
    geom_text(data = arcs[psi > 0.06], aes(x = (start+end)/2, y = psi + 0.03,
              label = sprintf("%.2f", psi), colour = grp), size = pt(5.2), family = FONT,
              show.legend = FALSE) +
    facet_wrap(~grp, ncol = 1) +
    scale_linewidth_continuous(range = c(0.2, 2.1), guide = "none") +
    scale_colour_manual(values = c(`young tertile` = "#2166AC", `old tertile` = "#B2182B"),
                        guide = "none") +
    scale_y_continuous(limits = c(-0.16, 1.12), breaks = c(0, 0.5, 1)) +
    labs(x = sprintf("%s  (%s)", gene, ex$chr[1]), y = "PSI") +
    theme_sa(7.5) + theme(strip.text = element_text(size = 6.8),
                          axis.text.x = element_text(size = 6)) }
saveRDS(sashimi, file.path(O, "sashimi_fn.rds"))
cat("sashimi function ready\n")
