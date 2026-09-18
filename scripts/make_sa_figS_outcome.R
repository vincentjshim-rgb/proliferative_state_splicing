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

## ================== a. the two measurements, same libraries =================
bx <- function(x, y, w, h, lab, fill = "white", col = INK, fc = INK)
  data.frame(x, y, w, h, lab, fill, col, fc, stringsAsFactors = FALSE)
B <- rbind(
  bx(50, 93, 98, 10, sprintf("%d normal dermal fibroblast donors, 20 to 96 years\nuniformly reprocessed RNA-seq\n(of 143 deposited libraries: 10 progeria, 26 under 20 excluded)", nrow(M)),
     "#F2F4F6", INK, INK),
  bx(25, 64, 47, 26, "", "#FBECEA", EXPC, EXPC),
  bx(75, 64, 47, 26, "", "#EAF0F8", OUTC, OUTC),
  bx(50, 37.5, 80, 8, "age effect per decade,\nbefore and after adjusting for proliferation", "white", INK, INK))
B$xmin <- B$x - B$w/2; B$xmax <- B$x + B$w/2; B$ymin <- B$y - B$h/2; B$ymax <- B$y + B$h/2
GY <- 59.6                                              # centre line of the gene model
pa <- ggplot() +
  ## donors -> the two measurements
  annotate("segment", x = 50, xend = 50, y = 88, yend = 83, colour = GREY, linewidth = 0.32) +
  annotate("segment", x = 25, xend = 75, y = 83, yend = 83, colour = GREY, linewidth = 0.32) +
  annotate("segment", x = c(25, 75), xend = c(25, 75), y = 83, yend = 77.5, colour = GREY,
           linewidth = 0.32, arrow = AR) +
  ## the two measurements -> the age model
  annotate("segment", x = c(25, 75), xend = c(25, 75), y = 51, yend = 46.5, colour = GREY, linewidth = 0.32) +
  annotate("segment", x = 25, xend = 75, y = 46.5, yend = 46.5, colour = GREY, linewidth = 0.32) +
  annotate("segment", x = 50, xend = 50, y = 46.5, yend = 42, colour = GREY, linewidth = 0.32, arrow = AR) +
  geom_rect(data = B, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = B$fill, colour = B$col, linewidth = 0.4) +
  geom_text(data = B, aes(x, y, label = lab), family = FONT, colour = B$fc, size = pt(7), lineheight = 1.0) +
  ## machinery
  annotate("text", x = 25, y = 74.4, family = FONT, colour = EXPC, size = pt(7.5), fontface = "bold",
           label = "machinery") +
  annotate("text", x = 25, y = 69.6, family = FONT, colour = EXPC, size = pt(7), lineheight = 1.0,
           label = "how much spliceosome\nis transcribed") +
  annotate("text", x = 25, y = 61.4, family = FONT, colour = EXPC, size = pt(7), fontface = "italic",
           lineheight = 1.05, label = "SRSF1   SRSF2\nHNRNPD   HNRNPA1\nHNRNPK   TRA2B") +
  annotate("text", x = 25, y = 53.8, family = FONT, colour = EXPC, size = pt(7),
           label = sprintf("%d genes, mean z", NCORE)) +
  ## outcome
  annotate("text", x = 75, y = 74.4, family = FONT, colour = OUTC, size = pt(7.5), fontface = "bold",
           label = "outcome") +
  annotate("text", x = 75, y = 69.6, family = FONT, colour = OUTC, size = pt(7), lineheight = 1.0,
           label = "where the spliceosome\nactually cuts") +
  ## gene model: annotated junctions grey, an unannotated junction blue
  annotate("segment", x = 57, xend = 93, y = GY, yend = GY, colour = INK2, linewidth = 0.3) +
  annotate("rect", xmin = c(57, 71, 85), xmax = c(65, 79, 93), ymin = GY - 1, ymax = GY + 1,
           fill = "white", colour = INK2, linewidth = 0.35) +
  annotate("curve", x = 65, xend = 71, y = GY + 1, yend = GY + 1, curvature = -0.75, ncp = 12,
           colour = GREY, linewidth = 0.4) +
  annotate("curve", x = 79, xend = 85, y = GY + 1, yend = GY + 1, curvature = -0.75, ncp = 12,
           colour = GREY, linewidth = 0.4) +
  annotate("curve", x = 65, xend = 89, y = GY - 1, yend = GY - 1, curvature = 0.16, ncp = 16,
           colour = OUTC, linewidth = 0.55) +
  annotate("segment", x = 89, xend = 89, y = GY - 1, yend = GY + 1, colour = OUTC, linewidth = 0.45) +
  annotate("text", x = 75, y = 64.4, family = FONT, colour = GREY, size = pt(7), label = "annotated") +
  annotate("text", x = 75, y = 54.4, family = FONT, colour = OUTC, size = pt(7), lineheight = 1.0,
           label = "unannotated reads /\nall junction reads") +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(32, 99), expand = c(0, 0)) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(12, 4, 6, 4))

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
XP <- max(FB$hi) + 0.12
pc <- ggplot(FB, aes(b, y, colour = k)) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.4) +
  geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = 0.55) +
  geom_point(size = 1.7) +
  geom_text(aes(x = XP, label = plab), hjust = 0, size = pt(7), family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(unadjusted = OUTC, `adjusted for proliferation` = GREY), name = NULL) +
  scale_y_continuous(breaks = K$row, labels = COHLAB[K$cohort], limits = c(0.5, max(K$row) + 0.5),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(min(FB$lo) - 0.05, XP + 0.62), breaks = c(-1, 0, 1), labels = num_axis()) +
  labs(x = sprintf("age effect per decade, unannotated read fraction (%s 10%s)", TIMES, sup("-4")), y = NULL) +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -6), legend.key.size = unit(6, "pt"),
                     axis.text.y = element_text(size = 7, lineheight = 0.95),
                     axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.title.x = element_text(hjust = 1),
                     plot.margin = margin(10, 8, 3, 3))

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


right <- lab_grid(pb, pc, labels = c("b", "c"), ncol = 1, rel_heights = c(1.12, 0.88))
top   <- lab_grid(pa, right, labels = c("a", ""), ncol = 2, rel_widths = c(0.86, 1.14))
bot   <- lab_grid(pd, pe, labels = c("d", "e"), ncol = 2, rel_widths = c(1, 1.02))
save_fig(lab_grid(top, bot, labels = c("", ""), ncol = 1, rel_heights = c(120, 88)),
         "FigS3.png", 183, 200)
