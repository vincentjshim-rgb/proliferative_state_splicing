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
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
S <- read.delim(file.path(D, "core_depth/GSE179848_sample_scores.tsv"))
S$grp <- ifelse(S$cond == "Control", "untreated", ifelse(S$cond == "Contact_Inhibition",
                "contact inhibition", "other treatment"))

## ============= A. the module tracks Ki-67 across all 345 samples ===========
ctA <- cor.test(S$MKI67, S$core, method = "spearman", exact = FALSE)
p2a <- ggplot(S, aes(MKI67, core)) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = GREY, alpha = 0.22, linewidth = 0.6) +
  geom_point(aes(colour = grp), size = 1.25, alpha = 0.85) +
  annotate("text", x = min(S$MKI67), y = max(S$core), hjust = 0, vjust = 1, family = FONT,
           size = pt(6.8), colour = INK, lineheight = 1.05,
           label = sprintf("rho = %.2f\nP = %s", ctA$estimate,
                           formatC(ctA$p.value, format = "e", digits = 0))) +
  annotate("text", x = max(S$MKI67), y = min(S$core), hjust = 1, vjust = 0, family = FONT,
           size = pt(6), colour = GREY, label = sprintf("n = %d samples", nrow(S))) +
  scale_colour_manual(values = c(untreated = RED, `contact inhibition` = BLUE,
                                 `other treatment` = GREY), name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
  labs(x = expression(paste("MKI67  ", log[2], " CPM")), y = "pre-mRNA processing score") +
  theme_sa(7.5) + theme(legend.position = "bottom", legend.margin = margin(t = -5))

## ============ B. and in untreated healthy cells alone ======================
H <- S[grepl("^HC", S$line) & S$cond == "Control" & S$oxy == 21, ]
ctB <- cor.test(H$MKI67, H$core, method = "spearman", exact = FALSE)
p2b <- ggplot(H, aes(MKI67, core)) +
  geom_smooth(method = "lm", se = TRUE, colour = RED, fill = RED, alpha = 0.16, linewidth = 0.6) +
  geom_point(aes(shape = line), size = 1.6, colour = RED) +
  annotate("text", x = min(H$MKI67), y = max(H$core), hjust = 0, vjust = 1, family = FONT,
           size = pt(6.8), colour = INK, lineheight = 1.05,
           label = sprintf("rho = %.2f\nP = %s", ctB$estimate,
                           formatC(ctB$p.value, format = "e", digits = 0))) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  guides(shape = guide_legend(nrow = 1)) +
  labs(x = expression(paste("MKI67  ", log[2], " CPM")), y = "pre-mRNA processing score") +
  theme_sa(7.5) + theme(legend.position = "bottom", legend.margin = margin(t = -5))

## ====== C. reversible arrest lowers the module as much as ageing does ======
CI <- S[grepl("^HC", S$line) & S$cond %in% c("Control","Contact_Inhibition") & S$oxy == 21, ]
CI$cond <- factor(ifelse(CI$cond == "Control", "untreated", "contact\ninhibition"),
                  levels = c("untreated","contact\ninhibition"))
mk <- function(var, ylab, cols) {
  d <- data.frame(cond = CI$cond, v = CI[[var]])
  wt <- wilcox.test(v ~ cond, d); yr <- range(d$v); pad <- diff(yr)*0.18
  ggplot(d, aes(cond, v)) +
    geom_boxplot(aes(colour = cond), outlier.shape = NA, width = 0.5, linewidth = 0.38, fill = NA) +
    geom_point(aes(colour = cond), size = 1.1, alpha = 0.8,
               position = position_jitter(width = 0.13, seed = 1)) +
    sigbar(1, 2, yr[2] + pad*0.55, stars(wt$p.value), tick = diff(yr)*0.035, size = 7) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_y_continuous(limits = c(yr[1] - pad*0.15, yr[2] + pad*1.25)) +
    labs(x = NULL, y = ylab) + theme_sa(7.5) +
    theme(axis.text.x = element_text(size = 6.6, lineheight = 0.85)) }
p2c1 <- mk("core",  "pre-mRNA processing", c(untreated = RED, `contact\ninhibition` = BLUE))
p2c2 <- mk("MKI67", "MKI67",               c(untreated = RED, `contact\ninhibition` = BLUE))
p2c3 <- mk("CDKN2A","CDKN2A (p16)",        c(untreated = RED, `contact\ninhibition` = BLUE))
p2c <- plot_grid(p2c1, p2c2, p2c3, nrow = 1)

## ========= D. Ki-67 absorbs the time-in-culture effect =====================
f1 <- lm(core ~ days + line, H); f2 <- lm(core ~ days + MKI67 + line, H)
Md <- data.frame(model = factor(c("time in culture","time in culture\n+ MKI67"),
                                levels = c("time in culture","time in culture\n+ MKI67")),
                 beta = c(coef(f1)["days"], coef(f2)["days"]),
                 se = c(summary(f1)$coef["days","Std. Error"], summary(f2)$coef["days","Std. Error"]),
                 p = c(summary(f1)$coef["days","Pr(>|t|)"], summary(f2)$coef["days","Pr(>|t|)"]))
p2d <- ggplot(Md, aes(model, beta)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.35) +
  geom_col(width = 0.5, fill = c(RED, "#B9C4CF")) +
  geom_errorbar(aes(ymin = beta - 1.96*se, ymax = beta + 1.96*se), width = 0.14, linewidth = 0.35) +
  geom_text(aes(y = beta - 1.96*se, label = sprintf("%s   P = %s", stars(p),
            ifelse(p < 1e-4, formatC(p, format="e", digits=0), sprintf("%.3f", p)))),
            vjust = 1.8, size = pt(6.2), family = FONT, colour = INK) +
  annotate("text", x = 1.5, y = 0.0018, hjust = 0.5, vjust = 0, family = FONT, size = pt(6.4),
           colour = INK, fontface = "bold",
           label = sprintf("%.0f%% of the effect absorbed by MKI67",
                           100*(1 - coef(f2)["days"]/coef(f1)["days"]))) +
  scale_y_continuous(limits = c(-0.0152, 0.0038)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "slope per day in culture") +
  theme_sa(7.5) + theme(axis.text.x = element_text(size = 6.6, lineheight = 0.85))

top <- lab_grid(p2a, p2b, labels = c("A","B"), ncol = 2)
bot <- lab_grid(p2c, p2d, labels = c("C","D"), ncol = 2, rel_widths = c(1.9, 1))
save_fig(lab_grid(top, bot, labels = c("",""), ncol = 1, rel_heights = c(1, 0.85)),
         "Fig2.png", 183, 112)
