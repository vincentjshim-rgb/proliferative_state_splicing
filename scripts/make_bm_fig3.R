FIGURES_WRITTEN <- c("Fig3.png")
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
PUBDIR <- "public_data_tierA/derived/figures_benchmark"
D <- "public_data_tierA/derived/benchmark_figs"; O <- "public_data_tierA/derived/secretome_class"
WB <- read.delim(file.path(D, "within_between_by_class.tsv"))
P  <- read.delim(file.path(O, "pairwise_magnitude_vs_similarity.tsv"))
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", secretome = REPRO,
        `photoprotection / rescue` = AMBER, `UV injury` = UV,
        `metabolic / culture` = GREY, reprogramming = NAVY)
WB <- WB[WB$class %in% names(CC), ]

## ============ A. class level: magnitude does not predict reproducibility ====
ct <- cor.test(WB$mag, WB$between, method = "spearman", exact = FALSE)
WB$hj <- ifelse(WB$mag > 0.45, 1.1, -0.10)
WB$vj <- c(senescence = -0.9, `donor age` = -0.9, secretome = -0.9,
           `photoprotection / rescue` = -0.9, `UV injury` = 1.9,
           `metabolic / culture` = 2.9, reprogramming = -0.9)[WB$class]
WB$lab <- sub(" / rescue", "", sub(" / culture", "", WB$class))
p3a <- ggplot(WB, aes(mag, between)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_point(aes(colour = class), size = 2.6) +
  geom_text(aes(label = lab, hjust = hj, vjust = vj, colour = class), size = pt(5.9),
            family = FONT, show.legend = FALSE) +
  annotate("text", x = 0.66, y = 0.50, hjust = 1, family = FONT, size = pt(6.3), colour = INK,
           lineheight = 1.05,
           label = sprintf("ρ = %+.3f\nP = %.2f  (n.s.)\n7 classes", ct$estimate, ct$p.value)) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(limits = c(0.02, 0.70)) +
  scale_y_continuous(limits = c(-0.08, 0.55)) +
  labs(x = expression(paste("effect magnitude   median |", log[2], "FC| of the class")),
       y = "median between-study ρ") +
  theme_j(7)

## ============ B. pair level, with the nested model ==========================
P <- P[P$class %in% names(CC), ]
m1 <- lm(r ~ gmag, P); m2 <- lm(r ~ gmag + factor(class), P); an <- anova(m1, m2)
nd <- data.frame(gmag = seq(min(P$gmag), max(P$gmag), length.out = 80))
nd <- cbind(nd, predict(m1, nd, interval = "confidence"))
p3b <- ggplot(P, aes(gmag, r)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_ribbon(data = nd, aes(gmag, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = GREY, alpha = 0.20) +
  geom_line(data = nd, aes(gmag, fit), inherit.aes = FALSE, colour = INK2,
            linewidth = 0.45, linetype = "22") +
  geom_point(aes(colour = class), size = 0.95, alpha = 0.8) +
  scale_colour_manual(values = CC, name = NULL) +
  annotate("text", x = 0.62, y = -0.28, hjust = 1, family = FONT, size = pt(6), colour = INK,
           lineheight = 1.05,
           label = sprintf("magnitude alone  R² = %.3f\n+ perturbation class  R² = %.3f\nclass effect  F = %.1f,  P = %.0e",
                           summary(m1)$r.squared, summary(m2)$r.squared, an$F[2], an$`Pr(>F)`[2])) +
  labs(x = expression(paste("geometric-mean magnitude of the pair   |", log[2], "FC|")),
       y = "Spearman ρ of the pair") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.5))) +
  theme_j(7) + theme(legend.position = "bottom", legend.justification = "center",
                     legend.box.background = element_blank(), legend.key.size = unit(6, "pt"),
                     legend.margin = margin(t = -4))

## ====== C. magnitude-matched restriction keeps the ordering ================
rng <- range(P$gmag[P$class == "secretome"])
sub <- P[P$gmag >= rng[1] & P$gmag <= rng[2], ]
ag  <- do.call(rbind, lapply(unique(sub$class), function(k) {
  v <- sub$r[sub$class == k]; if (length(v) < 3) return(NULL)
  data.frame(class = k, n = length(v), med = median(v)) }))
ag <- ag[order(ag$med), ]; sub <- sub[sub$class %in% ag$class, ]
sub$class <- factor(sub$class, levels = ag$class)
p3c <- ggplot(sub, aes(r, class, colour = class)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.34, fill = NA, show.legend = FALSE) +
  geom_point(size = 0.9, alpha = 0.75, position = position_jitter(height = 0.14, seed = 9)) +
  geom_text(data = ag, aes(x = 0.50, y = class, label = sprintf("%+.3f  (n = %d)", med, n)),
            hjust = 1, size = pt(5.7), family = FONT, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(limits = c(-0.30, 0.52), breaks = c(-0.25, 0, 0.25, 0.5)) +
  labs(x = sprintf("Spearman ρ of pairs restricted to |log2FC| %.2f-%.2f", rng[1], rng[2]),
       y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.4))

top <- lab_grid(p3a, p3b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1.1))
save_fig(lab_grid(top, p3c, labels = c("","C"), ncol = 1, rel_heights = c(1, 0.70)),
         "Fig3.png", 183, 108)
