FIGURES_WRITTEN <- c("Fig6.png")
## RETIRED 2026-09-18. This script writes a figure filename the current manuscript uses
## (CLAUDE.md section 9). Running it would overwrite a live figure with output from a
## superseded analysis -- in some cases one that still contains the withdrawn local
## contrast. Set ALLOW_RETIRED_FIGURE_SCRIPT=1 only if you mean to do that.
if (!nzchar(Sys.getenv("ALLOW_RETIRED_FIGURE_SCRIPT"))) {
  stop("retired script: would overwrite a live manuscript figure (", 
       paste(FIGURES_WRITTEN, collapse = ", "), "). See CLAUDE.md section 9.", call. = FALSE)
}

.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
R <- read.delim(file.path(D, "conserved_core/splicing_vs_proliferation.tsv"))
R <- R[is.finite(R$spl) & is.finite(R$cc), ]
R$class[is.na(R$class)] <- "metabolic / culture"
f <- lm(spl ~ cc, R); R$pred <- predict(f); R$resid <- resid(f)
N <- read.delim(file.path(D, "secretome_class/secretome_signatures.tsv"))
S <- R[R$class %in% c("secretome","Repro-CM"), ]
S$name <- N$contrast[match(S$id, N$id)]
S$name[S$id == "LOCAL_ReproCM"] <- "reprogramming-phase secretome"
S$src  <- N$dataset[match(S$id, N$id)]; S$src[S$id == "LOCAL_ReproCM"] <- "this study"
S$loc  <- S$id == "LOCAL_ReproCM"

## ========== A. observed change lies on the line predicted by proliferation =
lim <- range(c(S$spl, S$pred)) + c(-0.03, 0.06)
p4a <- ggplot(S, aes(pred, spl)) +
  geom_abline(slope = 1, intercept = 0, colour = GREY, linewidth = 0.4, linetype = "22") +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_segment(aes(xend = pred, yend = pred), colour = GREY_L, linewidth = 0.35) +
  geom_point(aes(colour = loc, size = loc)) +
  geom_text(data = S[S$loc | abs(S$resid) > 0.12, ],
            aes(label = sub(" vs .*", "", substr(name, 1, 24)),
                hjust = ifelse(resid > 0, -0.12, 1.12)),
            vjust = -0.55, size = pt(5.8), family = FONT, colour = INK) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = pt(6.4), family = FONT,
           colour = INK, lineheight = 1.05,
           label = "points on the dashed line are\nfully explained by proliferation") +
  scale_colour_manual(values = c(`FALSE` = TEAL, `TRUE` = "#0B3B3A"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 2.6), guide = "none") +
  coord_cartesian(xlim = lim, ylim = lim) +
  labs(x = "change predicted from the cell-cycle module",
       y = "observed pre-mRNA processing change") +
  theme_sa(7.5)

## ========== B. residual per secretome, with the null band =================
sdr <- sd(R$resid)
S2 <- S[order(S$resid), ]
S2$lab <- factor(sprintf("%s   %s", substr(S2$name, 1, 40), S2$src),
                 levels = sprintf("%s   %s", substr(S2$name, 1, 40), S2$src))
p4b <- ggplot(S2, aes(resid, lab)) +
  annotate("rect", xmin = -1.96*sdr, xmax = 1.96*sdr, ymin = -Inf, ymax = Inf,
           fill = "#E8EDF1", colour = NA) +
  geom_vline(xintercept = c(-1.96*sdr, 1.96*sdr), colour = GREY, linewidth = 0.32,
             linetype = "22") +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = resid, yend = lab), colour = GREY_L, linewidth = 0.5) +
  geom_point(aes(colour = loc, size = loc)) +
  scale_colour_manual(values = c(`FALSE` = TEAL, `TRUE` = "#0B3B3A"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.7, `TRUE` = 2.6), guide = "none") +
  annotate("text", x = 0, y = 16.0, label = "95% range of all 65 contrasts", size = pt(5.8),
           family = FONT, colour = GREY, vjust = 0) +
  scale_x_continuous(limits = c(-0.30, 0.17)) +
  coord_cartesian(ylim = c(0.5, 16.6), clip = "off") +
  labs(x = "residual after cell-cycle adjustment", y = NULL) +
  theme_sa(7.5) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.text.y = element_text(size = 5.9))

## ========== C. what would have been concluded, before and after ===========
top <- S[order(-S$spl), ][1:6, ]
Cd <- rbind(transform(top, v = spl,   k = "apparent effect"),
            transform(top, v = resid, k = "after adjustment"))
Cd$k <- factor(Cd$k, levels = c("apparent effect","after adjustment"))
Cd$lab <- factor(sub(" vs .*","",substr(Cd$name, 1, 28)),
                 levels = rev(sub(" vs .*","",substr(top$name, 1, 28))))
p4c <- ggplot(Cd, aes(v, lab, fill = k)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = NA) +
  geom_text(aes(label = sprintf("%+.2f", v)), position = position_dodge(width = 0.72),
            hjust = -0.15, size = pt(5.7), family = FONT, colour = INK) +
  scale_fill_manual(values = c(`apparent effect` = TEAL, `after adjustment` = "#C3CBD4"),
                    name = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(-0.09, 0.37)) +
  labs(x = "pre-mRNA processing change", y = NULL) +
  theme_sa(7.5) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                        axis.text.y = element_text(size = 6.0),
                        legend.position = "bottom", legend.margin = margin(t = -5))

save_fig(lab_grid(p4a, lab_grid(p4b, p4c, labels = c("B","C"), ncol = 1, rel_heights = c(1.25, 1)),
                  labels = c("A",""), ncol = 2, rel_widths = c(1, 1.12)),
         "Fig6.png", 183, 118)
