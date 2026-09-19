## Supplementary Fig. S3  The splicing outcome index (six-figure restructure, 2026-09-17).
## Panels a, c and e of the former main Fig. 4 (make_sa_fig4new.R):
##   a  the two quantities taken from the same libraries (schematic)
##   b  share of junction reads on unannotated junctions against donor age in the
##      133 normal donors, the 83-and-over stratum (one repository) marked
##   c  age effect per decade on that share, before and after adjustment for
##      proliferation, under the four cohort definitions
## Direction (a): the age association of unannotated-junction use is a negative
## result. It is present only in the published 142-sample cohort. The index was
## never prespecified; it is exploratory (CLAUDE.md section 10).
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"; CR <- file.path(D, "cohort_revised")
M  <- read.delim(file.path(CR, "primary/sample_metrics.tsv"))     # 107 normal adults
N  <- read.delim(file.path(CR, "normal/sample_metrics.tsv"))      # 133 normal donors
S  <- read.delim(file.path(CR, "sensitivity_metrics.tsv"))
NCORE <- length(readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt")))
AR <- arrow(length = unit(1.8, "pt"), type = "closed")
EXPC <- "#B2182B"; OUTC <- "#2166AC"
say <- function(...) cat(sprintf(...), "\n")

## =========== b. the outcome measure against donor age, 83+ marked ==========
GRP <- c("aged 83+ (all from one repository)", "under 83, same repository (AG)", "under 83, other repository")
N$grp <- factor(ifelse(N$age >= 83, GRP[1], ifelse(N$repo == "AG", GRP[2], GRP[3])), levels = GRP)
stopifnot(all(N$repo[N$age >= 83] == "AG"))
ctN <- cor.test(N$age, N$unannot_reads, method = "spearman", exact = FALSE)
Y <- N[N$age < 83, ]
ctY <- cor.test(Y$age, Y$unannot_reads, method = "spearman", exact = FALSE)
say("b  all normal donors: rho = %.4f, P = %.3g, n = %d donors | under 83: rho = %.4f, P = %.3g, n = %d donors",
    ctN$estimate, ctN$p.value, nrow(N), ctY$estimate, ctY$p.value, nrow(Y))
say("b  groups: %s", paste(sprintf("%s %d", levels(N$grp), as.integer(table(N$grp))), collapse = "; "))
say("b  median unannotated read share (%%): %s",
    paste(sprintf("%s %.3f", levels(N$grp), tapply(N$unannot_reads * 100, N$grp, median)), collapse = "; "))
ymax <- max(N$unannot_reads * 1e2)
pb <- ggplot(N, aes(age, unannot_reads * 1e2)) +
  geom_point(aes(colour = grp, shape = grp), size = 1.5, alpha = 0.9, stroke = 0.45) +
  geom_smooth(data = Y, method = "lm", formula = y ~ x, se = FALSE, colour = INK2, linewidth = 0.5,
              linetype = "22") +
  scale_colour_manual(values = setNames(c(ORANGE, GREY, INK2), GRP), name = NULL) +
  scale_shape_manual(values = setNames(c(17, 1, 16), GRP), name = NULL) +
  annotate("text", x = 0, y = ymax * 1.13, hjust = 0, vjust = 1, family = FONT,
           size = pt(7), colour = INK, lineheight = 1.15,
           label = sprintf("all normal donors (n = %d): %s = %s, %s\nunder 83 years (n = %d): %s = %s, %s",
                           nrow(N), RHO, num(ctN$estimate), pfmt(ctN$p.value),
                           nrow(Y), RHO, num(ctY$estimate), pfmt(ctY$p.value))) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(NA, ymax * 1.13), labels = num_axis()) +
  labs(x = "donor age (years)", y = "unannotated junction reads (%)") +
  theme_sa() + theme(legend.position = "bottom", legend.direction = "vertical",
                     legend.key.size = unit(7, "pt"), legend.margin = margin(t = -4),
                     legend.spacing.y = unit(0, "pt"), legend.text = element_text(size = 7),
                     plot.margin = margin(10, 8, 3, 3))

## =========== c. age effect per decade, before and after adjustment ==========
COHORTS <- c("published", "normal", "primary", "adult2082")
COHLAB  <- c(published = "deposited set\n(n = 142 samples)",
             normal    = "all normal donors\n(n = 133)",
             primary   = "normal adults 20+\n(primary, n = 107)",
             adult2082 = "normal adults 20–82\n(n = 76)")
K <- S[S$metric == "unannot_reads", ]; K <- K[match(COHORTS, K$cohort), ]
stopifnot(all(K$n == c(142, 133, 107, 76)))
K$row <- rev(seq_along(COHORTS)); SC <- 1e4
FB <- rbind(
  data.frame(row = K$row, k = "unadjusted", b = K$beta, lo = K$lo, hi = K$hi, p = K$p),
  data.frame(row = K$row, k = "adjusted for proliferation", b = K$beta_adj, lo = K$lo_adj,
             hi = K$hi_adj, p = K$p_adj))
FB[, c("b", "lo", "hi")] <- FB[, c("b", "lo", "hi")] * SC
FB$k <- factor(FB$k, levels = c("unadjusted", "adjusted for proliferation"))
FB$y <- FB$row + ifelse(FB$k == "unadjusted", 0.17, -0.17)
FB$plab <- pfmt_v(FB$p)
for (i in seq_len(nrow(K)))
  say("c  %-10s n = %3d  unadjusted %.3e [%.3e, %.3e] P = %.3g | adjusted %.3e [%.3e, %.3e] P = %.3g",
      K$cohort[i], K$n[i], K$beta[i], K$lo[i], K$hi[i], K$p[i], K$beta_adj[i], K$lo_adj[i],
      K$hi_adj[i], K$p_adj[i])
## grouped bars with 95% CI (2026-09-19); every P is in Supplementary Table 1
FB$cohort <- factor(rep(K$cohort, 2), levels = K$cohort[order(K$row, decreasing = TRUE)])
## short axis labels, the same wording as panel c
SHORT <- c(published = "deposited\n142", normal = "all normal\n133",
           primary = "adults 20+\n107 (primary)", adult2082 = "adults 20 to 82\n76")
SHORT <- SHORT[names(SHORT) %in% K$cohort]
if (length(SHORT) < nrow(K)) SHORT <- COHLAB   # unknown cohort keys: fall back to the long labels
pc <- ggplot(FB, aes(cohort, b, fill = k)) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.72), width = 0.22,
                colour = INK2, linewidth = 0.4) +
  scale_fill_manual(values = c(unadjusted = OUTC, `adjusted for proliferation` = GREY), name = NULL) +
  scale_x_discrete(labels = SHORT) +
  scale_y_continuous(breaks = c(-1, 0, 1), labels = num_axis()) +
  labs(x = NULL, y = sprintf("age effect per decade,\nunannotated read fraction (%s 10%s)", TIMES, sup("-4"))) +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -6), legend.key.size = unit(6, "pt"),
                     axis.text.x = element_text(size = 7, lineheight = 0.95),
                     plot.margin = margin(3, 8, 3, 3))

## ================================================================================
## 2026-09-18: the event-level panels, previously a supplementary figure of their own,
## were folded in here.  Both halves report the same negative result -- neither the
## outcome index nor the events it is built from carries an age effect in the adult
## cohort -- so they belong in one figure.  Code moved from make_splicing_fig.R, which
## is now retired.
## ================================================================================
suppressPackageStartupMessages({library(data.table)})
P  <- file.path(D, "psi_cohorts")
CB <- fread(file.path(P, "psi_events_by_cohort.tsv"))
DN <- fread(file.path(P, "col1a2_top_event_per_donor.tsv"))
pv <- function(p) vapply(p, pfmt, character(1))

## ---------------- a. how many events change with age, by cohort -------------
CB[, lab := c("deposited\n142", "all normal\n133", "adults 20+\n107 (primary)", "adults 20 to 82\n76")]
CB[, lab := factor(lab, levels = lab)]
EA <- melt(CB[, .(lab, `age-associated` = age_events,
                 `also after proliferation adjustment` = kept_after_proliferation)],
          id.vars = "lab", variable.name = "k", value.name = "n")
pd <- ggplot(EA, aes(lab, n, fill = k)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(width = 0.68), vjust = -0.45,
            size = pt(7), family = FONT, colour = INK) +
  scale_fill_manual(values = setNames(c(BLUE, GREY_L), levels(EA$k)), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "alternative splicing events\nassociated with donor age (of 3,401)") +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -4), legend.key.size = unit(6, "pt"),
                     axis.text.x = element_text(size = 7, lineheight = 1.05))

## ---------------- b. the COL1A2 event that led the published list ------------
EV <- DN[disease == "Normal"]
EV[, grp := ifelse(age >= 83, "aged 83+ (all repository AG)",
                  ifelse(repo == "AG", "under 83, repository AG", "under 83, other repository"))]
ct_all <- cor.test(EV$age, EV$psi_top, method = "spearman", exact = FALSE)
ct_u83 <- cor.test(EV$age[EV$age < 83], EV$psi_top[EV$age < 83], method = "spearman", exact = FALSE)
pe <- ggplot(EV, aes(age, psi_top)) +
  geom_point(aes(colour = grp, shape = grp), size = 1.5, alpha = 0.9) +
  geom_smooth(data = EV[age < 83], method = "lm", se = FALSE, colour = INK2, linewidth = 0.5, linetype = "22") +
  scale_colour_manual(values = setNames(c(ORANGE, INK2, GREY), sort(unique(EV$grp))), name = NULL) +
  scale_shape_manual(values = setNames(c(17, 16, 1), sort(unique(EV$grp))), name = NULL) +
  annotate("text", x = 0, y = max(EV$psi_top), hjust = 0, vjust = 1, family = FONT, size = pt(7),
           colour = INK, lineheight = 1.2,
           label = sprintf("all normal donors  %s = %s, %s\nunder 83 years  %s = %s, %s",
                           RHO, num(ct_all$estimate), pfmt(ct_all$p.value),
                           RHO, num(ct_u83$estimate), pfmt(ct_u83$p.value))) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(labels = num_axis(), expand = expansion(mult = c(0.06, 0.16))) +
  labs(x = "donor age (years)",
       y = expression(paste("PSI, distal acceptor of ", italic("COL1A2"), " chr7:94,406,304"))) +
  theme_sa() + theme(legend.position = "bottom", legend.direction = "vertical",
                     legend.key.size = unit(6, "pt"), legend.margin = margin(t = -4),
                     legend.spacing.y = unit(0, "pt"),
                     axis.title.y = element_text(size = 7.5))


top <- lab_grid(pb, pc, labels = c("a", "b"), ncol = 2, rel_widths = c(1.08, 1))
bot <- lab_grid(pd, pe, labels = c("c", "d"), ncol = 2, rel_widths = c(1, 1.02))
save_fig(lab_grid(top, bot, labels = c("", ""), ncol = 1, rel_heights = c(1, 1)),
         "FigS3.png", 183, 168)
