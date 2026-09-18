FIGURES_WRITTEN <- c("Fig2.png")
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
D <- "public_data_tierA/derived"; CS <- file.path(D, "classical_stats")

## ============ A. random-effects meta-analysis forest (Fisher-z, DL) =========
per <- read.delim(file.path(CS, "per_study_rho_fisherCI.tsv"))
met <- read.delim(file.path(CS, "meta_analysis_axes.tsv"))
per$w <- (1/per$se^2); per$w <- per$w/max(per$w)
ax   <- c(UVA = "UVA reference axis", senescence = "Senescence reference axis")
rows <- list()
y <- 0
for (a in c("UVA","senescence")) {
  y <- y - 1
  rows[[length(rows)+1]] <- data.frame(y, kind = "header", lab = ax[a], est = NA, lo = NA, hi = NA,
                                       w = NA, axis = a, txt = "")
  s <- per[per$axis == a, ]
  for (i in seq_len(nrow(s))) { y <- y - 1
    rows[[length(rows)+1]] <- data.frame(y, kind = "study", lab = s$label[i], est = s$rho[i],
      lo = s$ci_lo[i], hi = s$ci_hi[i], w = s$w[i], axis = a,
      txt = sprintf("%+.3f (%+.3f, %+.3f)", s$rho[i], s$ci_lo[i], s$ci_hi[i])) }
  m <- met[met$axis == a, ]; y <- y - 1
  rows[[length(rows)+1]] <- data.frame(y, kind = "pooled", lab = "Random-effects estimate",
    est = m$pooled_rho, lo = m$ci_lo, hi = m$ci_hi, w = NA, axis = a,
    txt = sprintf("%+.3f (%+.3f, %+.3f)", m$pooled_rho, m$ci_lo, m$ci_hi))
  y <- y - 0.75
  rows[[length(rows)+1]] <- data.frame(y, kind = "het", axis = a, est = NA, lo = NA, hi = NA, w = NA,
    lab = sprintf("Q = %.1f, df = %d, %s;  I² = %.0f%%", m$Q, m$Q_df,
                  ifelse(m$Q_p < 1e-4, "P < 0.0001", sprintf("P = %.3f", m$Q_p)), m$I2), txt = "")
  y <- y - 0.5
}
F <- do.call(rbind, rows)
dia <- do.call(rbind, lapply(which(F$kind == "pooled"), function(i)
  data.frame(x = c(F$lo[i], F$est[i], F$hi[i], F$est[i]),
             y = F$y[i] + c(0, .30, 0, -.30), g = i, axis = F$axis[i])))
AXC <- c(UVA = UV, senescence = SEN)
XL <- c(-1.30, 1.00)
AXB <- c(-0.2, 0, 0.2, 0.4); ybot <- min(F$y) - 1.15
p2a <- ggplot() +
  annotate("segment", x = 0, xend = 0, y = ybot, yend = 0.55, colour = INK, linewidth = 0.3) +
  geom_segment(data = F[F$kind=="study",], aes(x = lo, xend = hi, y = y, yend = y),
               colour = GREY, linewidth = 0.35) +
  geom_point(data = F[F$kind=="study",], aes(x = est, y = y, size = w, colour = axis),
             shape = 15) +
  geom_polygon(data = dia, aes(x, y, group = g, fill = axis), colour = NA) +
  geom_text(data = F[F$kind %in% c("study","pooled"),],
            aes(x = XL[1], y = y, label = lab), hjust = 0, size = pt(6), family = FONT,
            colour = ifelse(F$kind[F$kind %in% c("study","pooled")]=="pooled", INK, INK2),
            fontface = ifelse(F$kind[F$kind %in% c("study","pooled")]=="pooled","bold","plain")) +
  geom_text(data = F[F$kind == "header",], aes(x = XL[1], y = y, label = lab),
            hjust = 0, size = pt(6.4), family = FONT, fontface = "bold", colour = AXC[F$axis[F$kind=="header"]]) +
  geom_text(data = F[F$kind == "het",], aes(x = XL[1]+0.005, y = y, label = lab),
            hjust = 0, size = pt(5.5), family = FONT, colour = GREY, fontface = "italic") +
  geom_text(data = F[F$txt != "",], aes(x = XL[2], y = y, label = txt), hjust = 1,
            size = pt(5.8), family = FONT, colour = INK2) +
  scale_size_continuous(range = c(1.1, 2.6), guide = "none") +
  scale_colour_manual(values = AXC, guide = "none") +
  scale_fill_manual(values = AXC, guide = "none") +
  annotate("text", x = XL[2], y = 0.05, label = "ρ (95% CI)", hjust = 1, size = pt(5.8),
           family = FONT, colour = INK, fontface = "bold") +
  annotate("segment", x = min(AXB)-0.06, xend = max(AXB)+0.06, y = ybot, yend = ybot,
           colour = INK, linewidth = 0.3) +
  annotate("segment", x = AXB, xend = AXB, y = ybot, yend = ybot - 0.16,
           colour = INK, linewidth = 0.3) +
  annotate("text", x = AXB, y = ybot - 0.30, label = sprintf("%.1f", AXB), vjust = 1,
           size = pt(6.5), family = FONT, colour = INK) +
  annotate("text", x = mean(range(AXB)), y = ybot - 1.15,
           label = "Spearman ρ   (anchor vs reference contrast)",
           size = pt(7), family = FONT, colour = INK) +
  scale_x_continuous(limits = XL) +
  coord_cartesian(ylim = c(ybot - 1.5, 0.75), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_j(7) + theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_blank(), axis.text.y = element_blank(),
                     plot.margin = margin(2, 3, 2, 2))

## ================= B. GSEA running enrichment against both axes =============
rd <- function(ax, set) {
  c <- read.delim(file.path(CS, sprintf("gsea_curve_%s_%s.tsv", ax, set)))
  t <- read.delim(file.path(CS, sprintf("gsea_ticks_%s_%s.tsv", ax, set)))
  c$axis <- ax; c$set <- set; t$axis <- ax; t$set <- set; list(c = c, t = t) }
cu <- do.call(rbind, lapply(c("UVA_meta","senescence_meta"), function(a)
        do.call(rbind, lapply(c("ReproCMup","ReproCMdown"), function(s) rd(a, s)$c))))
tk <- do.call(rbind, lapply(c("UVA_meta","senescence_meta"), function(a)
        do.call(rbind, lapply(c("ReproCMup","ReproCMdown"), function(s) rd(a, s)$t))))
cu$axl <- factor(c(UVA_meta = "ranked by UVA axis", senescence_meta = "ranked by senescence axis")[cu$axis],
                 levels = c("ranked by UVA axis","ranked by senescence axis"))
tk$axl <- factor(c(UVA_meta = "ranked by UVA axis", senescence_meta = "ranked by senescence axis")[tk$axis],
                 levels = c("ranked by UVA axis","ranked by senescence axis"))
SETC <- c(ReproCMup = REPRO, ReproCMdown = AMBER)
tk$yy <- ifelse(tk$set == "ReproCMup", -0.30, -0.40)
gs <- read.delim(file.path(CS, "fgsea_anchor_vs_references.tsv"))
p2b <- ggplot() +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.25) +
  geom_segment(data = tk, aes(x = rank, xend = rank, y = yy, yend = yy + 0.055, colour = set),
               linewidth = 0.12, alpha = 0.75) +
  geom_line(data = cu, aes(rank, ES, colour = set), linewidth = 0.5) +
  facet_wrap(~ axl, nrow = 1) +
  scale_colour_manual(values = SETC, labels = c(ReproCMdown = "Repro-CM down (287)",
                      ReproCMup = "Repro-CM up (236)"), name = NULL,
                      breaks = c("ReproCMup","ReproCMdown")) +
  scale_y_continuous(limits = c(-0.46, 0.60), breaks = c(-0.25, 0, 0.25, 0.5)) +
  scale_x_continuous(breaks = c(0, 1500, 3000), labels = c("0","1,500","3,000")) +
  labs(x = "gene rank in the reference signature", y = "enrichment score") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 0.9))) +
  theme_j(7) + theme(legend.position = "top", legend.justification = "center",
                     legend.box.background = element_blank(), legend.margin = margin(b=-5),
                     legend.key.width = unit(9, "pt"), panel.spacing = unit(9, "pt"))

## ============== C. NES per reference study, both gene sets =================
gs$lab <- factor(gs$label, levels = rev(unique(gs$label)))
gs$star <- sig_stars(gs$FDR)
p2c <- ggplot(gs, aes(NES, lab, fill = pathway)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = NA) +
  geom_text(aes(label = star, hjust = ifelse(NES > 0, -0.2, 1.2)),
            position = position_dodge(width = 0.72), size = pt(6), family = FONT, colour = INK2) +
  scale_fill_manual(values = c(`Repro-CM up` = REPRO, `Repro-CM down` = AMBER), name = NULL,
                    breaks = c("Repro-CM up","Repro-CM down")) +
  scale_x_continuous(limits = c(-2.6, 4.0), breaks = c(-2, 0, 2, 4)) +
  labs(x = "normalised enrichment score", y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_j(7) + theme_caty() +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.box.background = element_blank(),
        legend.margin = margin(t = -4))

## ================= D. transferability of each axis ==========================
ug  <- read.delim(file.path(D, "axis_characterisation/UVA_axis_generalisation.tsv"))
sg2 <- read.delim(file.path(D, "axis_characterisation/senescence_axis_generalisation.tsv"))
GG  <- rbind(transform(ug[, c("kind","rho_axis")], axis = "UVA"),
             transform(sg2[, c("kind","rho_axis")], axis = "senescence"))
GG  <- GG[GG$kind != "kind" & is.finite(suppressWarnings(as.numeric(GG$rho_axis))), ]
GG$rho <- as.numeric(GG$rho_axis)
GG$grp <- ifelse(grepl("\\(axis\\)", GG$kind), "defines the axis",
          ifelse(grepl("same study|within the axis", GG$kind), "same study",
          ifelse(grepl("photoprotection", GG$kind), "independent\nphotoprotection",
          ifelse(grepl("^injury", GG$kind), "independent\nUV injury",
          ifelse(grepl("independent", GG$kind), "fully\nindependent", NA)))))
GG <- GG[!is.na(GG$grp), ]
GG$glab <- factor(GG$grp, levels = c("defines the axis","same study",
   "independent\nphotoprotection","independent\nUV injury","fully\nindependent"))
GG$axl <- factor(c(UVA = "UVA axis", senescence = "Senescence axis")[GG$axis],
                 levels = c("UVA axis","Senescence axis"))
inj <- GG$rho[GG$glab == "independent\nUV injury"]; pro <- GG$rho[GG$glab == "independent\nphotoprotection"]
wt  <- wilcox.test(inj, pro)
p2d <- ggplot(GG, aes(glab, rho)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(aes(colour = axl), outlier.shape = NA, width = 0.5, linewidth = 0.32,
               fill = NA, show.legend = FALSE,
               position = position_dodge2(width = 0.62, preserve = "single")) +
  geom_point(aes(colour = axl), size = 0.95, alpha = 0.9,
             position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.62, seed = 7)) +
  geom_hline(yintercept = 0.289, colour = REPRO, linewidth = 0.35, linetype = "22") +
  annotate("text", x = 0.62, y = 0.289, label = "Repro-CM", colour = REPRO, size = pt(5.6),
           family = FONT, vjust = -0.55, hjust = 0) +
  brk(3, 4, 0.44, sprintf("%s   Wilcoxon P = %.4f", sig_stars(wt$p.value), wt$p.value),
      tick = 0.035, size = 5.6) +
  scale_colour_manual(values = c(`UVA axis` = UV, `Senescence axis` = SEN), name = NULL) +
  scale_y_continuous(limits = c(-0.30, 1.0), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = "ρ with its own reference axis") +
  theme_j(7) + theme(axis.text.x = element_text(size = 6, lineheight = 0.92),
                     legend.position = "inside", legend.position.inside = c(0.20, 0.90),
                     legend.box.background = element_blank())
top <- lab_grid(p2a, p2b, labels = c("A","B"), ncol = 2, rel_widths = c(1.05, 1))
bot <- lab_grid(p2c, p2d, labels = c("C","D"), ncol = 2, rel_widths = c(1, 1.18))
save_fig(lab_grid(top, bot, labels = c("",""), ncol = 1, rel_heights = c(1.15, 1)),
         "Fig2.png", 183, 124)
