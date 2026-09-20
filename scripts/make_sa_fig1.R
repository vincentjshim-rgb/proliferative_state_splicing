## Fig. 1  Do splicing-machinery transcripts fall with donor age and time in culture?
## Revised 2026-09-15: panel c is rebuilt on the redefined cohort (107 normal donors
## aged 20-96, recount3 matrix) so that Figs. 1, 3 and 4 share one matrix and one
## cohort; panel d prints exact P values for the short series.
## Revised 2026-09-20 (user request): the figure opens with the question. Panel a is a
## schematic drawn from the fitted lines of the real data -- the splicing set and the
## cell-cycle markers fall together, along donor age and along time in culture -- and the
## donor-age and culture-time scatters that were two panels are one row of small multiples.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
suppressPackageStartupMessages(library(cowplot))
D <- "public_data_tierA/derived"
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))

## ============ inputs shared by a and c ======================================
A <- read.delim(file.path(D, "cohort_revised/primary/sample_metrics.tsv"))        # 107 donors
L  <- read.delim(file.path(D, "core_depth/GSE179848_longitudinal_core.tsv"))       # 4 lines, days
CT <- read.delim(file.path(D, "revision_stats/fig1d_culture_time.tsv"))            # exact P where n < 10
LS <- read.delim(file.path(D, "treatment_perturbation/library_scores_vs_counted_rate.tsv"))
stopifnot(nrow(A) == 107)

## ============ a. the question, drawn from the fitted lines =================
## Two programmes, standardised within each dataset, fitted against the two time axes the
## paper starts from. Lines and 95% bands are least-squares fits to the real samples; no
## point is drawn because the same samples appear as points in panel c and in Fig. 2.
##   donor age: the 107-donor cohort (splicing set score, 20-marker proliferation score)
##   time in culture: the untreated cultures of the four healthy lines at 21% oxygen, the
##   libraries Fig. 1c follows, scored in the counted resource
U <- LS[LS$cond == "Control" & LS$oxy == 21 & LS$line %in% c("HC1", "HC2", "HC3", "HC4"), ]
zz <- function(v) (v - mean(v)) / sd(v)
SA <- rbind(data.frame(axis = "donor age", x = A$age, y = zz(A$machinery), prog = "splicing set (177 genes)"),
            data.frame(axis = "donor age", x = A$age, y = zz(A$prolif),     prog = "cell-cycle markers (20 genes)"),
            data.frame(axis = "time in culture", x = U$days, y = zz(U$splice), prog = "splicing set (177 genes)"),
            data.frame(axis = "time in culture", x = U$days, y = zz(U$prolif), prog = "cell-cycle markers (20 genes)"))
SA$axis <- factor(SA$axis, levels = c("donor age", "time in culture"),
                  labels = c(sprintf("donor age\n%d to %d years, %d donors", min(A$age), max(A$age), nrow(A)),
                             sprintf("time in culture\n%d to %d days, %d cultures", round(min(U$days)), round(max(U$days)), nrow(U))))
SA$prog <- factor(SA$prog, levels = c("splicing set (177 genes)", "cell-cycle markers (20 genes)"))
FT <- do.call(rbind, lapply(split(SA, list(SA$axis, SA$prog), drop = TRUE), function(s) {
  f <- lm(y ~ x, s); nd <- data.frame(x = seq(min(s$x), max(s$x), length.out = 60))
  cbind(axis = s$axis[1], prog = s$prog[1], nd, predict(f, nd, interval = "confidence")) }))
SLOPE <- do.call(rbind, lapply(split(SA, list(SA$axis, SA$prog), drop = TRUE), function(s)
  data.frame(axis = s$axis[1], prog = s$prog[1], r = cor(s$x, s$y))))
cat("a  correlations behind the schematic:\n"); print(SLOPE, row.names = FALSE)
PCOL <- c(`splicing set (177 genes)` = BLUE, `cell-cycle markers (20 genes)` = INK2)
p1a <- ggplot(FT, aes(x, fit, colour = prog, fill = prog)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.14) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ axis, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = PCOL, name = NULL) + scale_fill_manual(values = PCOL, name = NULL) +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous("transcript level\n(standardised fit)", breaks = NULL) +
  labs(title = "Two programmes fall together in cultured fibroblasts",
       caption = "How much of the splicing signal reports how fast the cells divide?") +
  theme_sa(8) + theme(strip.text = element_text(size = 7.2, hjust = 0, lineheight = 0.95, margin = margin(b = 2)),
                      axis.title.y = element_text(size = 7.5, lineheight = 0.95),
                      plot.title = element_text(size = 8, colour = INK, hjust = 0, margin = margin(b = 4)),
                      plot.caption = element_text(size = 7.5, colour = INK, hjust = 0, margin = margin(t = 5)),
                      panel.spacing.x = unit(8, "pt"),
                      legend.position = "top", legend.justification = "left", legend.margin = margin(t = -2, b = -2),
                      legend.key.size = unit(7, "pt"), legend.text = element_text(size = 7.2),
                      plot.margin = margin(3, 8, 3, 3))

## ============ b. what the three series agree on ============================
ENT <- read.delim(file.path(D, "conserved_core/age_core_v2_reactome.tsv"))
SHOW <- c("Cell Cycle, Mitotic" = "Cell cycle, mitotic",
          "Processing of Capped Intron-Containing Pre-mRNA" = "Processing of capped intron-\ncontaining pre-mRNA",
          "M Phase" = "M phase",
          "Cell Cycle Checkpoints" = "Cell cycle checkpoints",
          "mRNA Splicing - Major Pathway" = "mRNA splicing,\nmajor pathway",
          "mRNA 3'-end processing" = "mRNA 3'-end processing",
          "DNA Repair" = "DNA repair")
FAM  <- c("cell cycle", "RNA processing", "cell cycle", "cell cycle", "RNA processing", "RNA processing", "DNA repair")
EN <- ENT[match(names(SHOW), ENT$pathway), ]
stopifnot(!anyNA(EN$k))
EN$label <- unname(SHOW); EN$fam <- FAM; EN$pct <- 100 * EN$k / EN$size
EN <- EN[order(EN$rank), ]; EN$label <- factor(EN$label, levels = rev(EN$label))
EN$txt <- sprintf("%d/%d   FDR %s", EN$k, EN$size, vapply(EN$FDR, sci, character(1)))
p1b <- ggplot(EN, aes(pct, label, fill = fam)) +
  geom_col(width=0.62) +
  geom_text(aes(label=txt), hjust=-0.06, size=pt(7), family=FONT, colour=INK) +
  scale_fill_manual(values=c(`cell cycle`=RED, `RNA processing`=ORANGE, `DNA repair`=PURPLE), name=NULL,
                    breaks=c("cell cycle","RNA processing","DNA repair")) +
  scale_x_continuous(limits=c(0,150), breaks=seq(0,100,25), expand=expansion(mult=c(0,0))) +
  labs(x="% of pathway falling in all three series", y=NULL) +
  theme_sa(8) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(),
                      axis.text.y=element_text(size=7, lineheight=0.95),
        ## the key sits within the panel width: small keys, tight spacing, 7 pt text
        legend.position="top", legend.justification="left", legend.key.size=unit(5,"pt"),
        legend.spacing.x=unit(2,"pt"), legend.margin=margin(b=-4),
        legend.text=element_text(size=7, margin=margin(r=3)),
        plot.margin = margin(3, 6, 3, 3))

## ============ c. donor age and time in culture, one row =====================
## left: the redefined cohort; right: the four healthy lines followed in culture. One y-axis
## title serves the row; the two x-axes keep their own units.
ct1 <- cor.test(A$age, A$machinery)                       # Pearson: the panel shows a linear fit
f1 <- lm(machinery ~ age, A)
nd <- data.frame(age = seq(min(A$age), max(A$age), length.out = 60))
nd <- cbind(nd, predict(f1, nd, interval = "confidence"))
p1c1 <- ggplot(A, aes(age, machinery)) +
  geom_ribbon(data=nd, aes(age, ymin=lwr, ymax=upr), inherit.aes=FALSE, fill=BROWN, alpha=0.18) +
  geom_line(data=nd, aes(age, fit), inherit.aes=FALSE, colour=BROWN, linewidth=0.6) +
  geom_point(size=1.3, colour=BROWN, alpha=0.85) +
  annotate("text", x=21, y=min(A$machinery), hjust=0, vjust=0, family=FONT, size=pt(7),
           colour=INK, lineheight=1.1, label=rp(ct1$estimate, ct1$p.value, lab="r")) +
  scale_y_continuous(labels = num_axis()) +
  labs(x="donor age (years)", y="pre-mRNA processing score",
       subtitle = sprintf("%d donors", nrow(A))) +
  theme_sa(8) + theme(plot.subtitle = element_text(size = 7.5, hjust = 0.5, margin = margin(b = 1)),
                      plot.margin = margin(3, 6, 3, 3))

st <- do.call(rbind, lapply(sort(unique(L$line)), function(l) {
  s <- L[L$line == l, ]; r <- CT[CT$line == l, ]
  data.frame(line = l, lab = sprintf("%s,  n = %d,  %d days", l, r$n, r$days),
             txt = rp(r$rho, r$p), x = max(s$days), y = max(L$score)) }))
L$lab <- st$lab[match(L$line, st$line)]
p1c2 <- ggplot(L, aes(days, score)) +
  geom_smooth(method="lm", se=TRUE, colour=RED, fill=RED, alpha=0.16, linewidth=0.6) +
  geom_point(size=1.2, colour=RED, alpha=0.85) +
  geom_text(data=st, aes(x=x, y=y, label=txt), hjust=1, vjust=1, size=pt(7),
            family=FONT, colour=INK, lineheight=1.1) +
  facet_wrap(~lab, nrow=1, scales="free_x") +
  scale_y_continuous(labels = num_axis()) +
  labs(x="days in culture", y=NULL) +
  theme_sa(8) + theme(strip.text=element_text(size=7.5, lineheight=1.05),
                      plot.margin = margin(3, 3, 3, 4))
p1c <- plot_grid(p1c1, p1c2, ncol = 2, rel_widths = c(1, 3.15), align = "h", axis = "tb")

top <- lab_grid(p1a, p1b, labels=c("a","b"), ncol=2, rel_widths=c(1, 1.02))
save_fig(lab_grid(top, lab_grid(p1c, labels="c", ncol=1),
                  labels=c("",""), ncol=1, rel_heights=c(1, 0.92)), "Fig1.png", 183, 118)
