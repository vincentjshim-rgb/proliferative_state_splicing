## Supplementary Fig. S4 (six-figure restructure 2026-09-17; was Fig. S1). Event-level
## splicing in the donor cohort, and why it is not interpreted. Direction (a) removes
## the COL1A2 example and the event-level age signal from the claims: the counts
## depend on how the cohort is defined, the flagship event is a step confined to the donors aged 83 and over (all from one
## cell repository), and the junction carrying it is absent from GENCODE v41.
## 2026-09-17: the panel showing that junction's read support in four locally
## sequenced libraries was removed with the rest of the local data; the junction
## counts remain in psi_cohorts/col1a2_local_junctions.tsv.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
D <- "public_data_tierA/derived"; P <- file.path(D, "psi_cohorts")
CB <- fread(file.path(P, "psi_events_by_cohort.tsv"))
DN <- fread(file.path(P, "col1a2_top_event_per_donor.tsv"))
pv <- function(p) vapply(p, pfmt, character(1))

## ---------------- a. how many events change with age, by cohort -------------
CB[, lab := c("deposited\n142", "all normal\n133", "adults 20+\n107 (primary)", "adults 20 to 82\n76")]
CB[, lab := factor(lab, levels = lab)]
A <- melt(CB[, .(lab, `age-associated` = age_events,
                 `also after proliferation adjustment` = kept_after_proliferation)],
          id.vars = "lab", variable.name = "k", value.name = "n")
pa <- ggplot(A, aes(lab, n, fill = k)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(width = 0.68), vjust = -0.45,
            size = pt(7), family = FONT, colour = INK) +
  scale_fill_manual(values = setNames(c(BLUE, GREY_L), levels(A$k)), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "alternative splicing events\nassociated with donor age (of 3,401)") +
  theme_sa() + theme(legend.position = "top", legend.justification = "left",
                     legend.margin = margin(b = -4), legend.key.size = unit(6, "pt"),
                     axis.text.x = element_text(size = 7, lineheight = 1.05))

## ---------------- b. the COL1A2 event that led the published list ------------
N <- DN[disease == "Normal"]
N[, grp := ifelse(age >= 83, "aged 83+ (all repository AG)",
                  ifelse(repo == "AG", "under 83, repository AG", "under 83, other repository"))]
ct_all <- cor.test(N$age, N$psi_top, method = "spearman", exact = FALSE)
ct_u83 <- cor.test(N$age[N$age < 83], N$psi_top[N$age < 83], method = "spearman", exact = FALSE)
pb <- ggplot(N, aes(age, psi_top)) +
  geom_point(aes(colour = grp, shape = grp), size = 1.5, alpha = 0.9) +
  geom_smooth(data = N[age < 83], method = "lm", se = FALSE, colour = INK2, linewidth = 0.5, linetype = "22") +
  scale_colour_manual(values = setNames(c(ORANGE, INK2, GREY), sort(unique(N$grp))), name = NULL) +
  scale_shape_manual(values = setNames(c(17, 16, 1), sort(unique(N$grp))), name = NULL) +
  annotate("text", x = 0, y = max(N$psi_top), hjust = 0, vjust = 1, family = FONT, size = pt(7),
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

save_fig(lab_grid(pa, pb, labels = c("a", "b"), ncol = 2, rel_widths = c(1, 1.02)),
         "FigS4.png", 183, 84)
