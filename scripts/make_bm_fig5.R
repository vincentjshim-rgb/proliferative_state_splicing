FIGURES_WRITTEN <- c("Fig5.png")
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
Z <- readRDS(file.path(D, "figdata.rds")); R <- Z$R; M <- Z$meta
N <- read.delim(file.path(O, "secretome_signatures.tsv"))
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", secretome = REPRO,
        `photoprotection / rescue` = AMBER, `UV injury` = UV,
        `metabolic / culture` = GREY, reprogramming = NAVY, `Repro-CM` = "#0B3B3A")
fzi <- function(r, n) { se <- 1/sqrt(n-3); c(tanh(atanh(r)-1.96*se), tanh(atanh(r)+1.96*se)) }

## ====== A. every secretome contrast on the one reproducible axis ===========
S <- M[M$class %in% c("secretome","Repro-CM") & is.finite(M$rho_sen), ]
S$lab <- N$contrast[match(S$id, N$id)]
S$lab[S$id == "LOCAL_ReproCM"] <- "reprogramming-phase secretome (this study)"
S$src <- N$dataset[match(S$id, N$id)]; S$src[S$id == "LOCAL_ReproCM"] <- "this study"
S$lab <- sprintf("%s   %s", substr(S$lab, 1, 46), S$src)
ci <- t(mapply(fzi, S$rho_sen, S$n_sen)); S$lo <- ci[,1]; S$hi <- ci[,2]
S <- S[order(S$rho_sen), ]; S$lab <- factor(S$lab, levels = S$lab)
S$isloc <- S$id == "LOCAL_ReproCM"
p5a <- ggplot(S, aes(rho_sen, lab)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F1F5F4", colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.36, colour = GREY) +
  geom_point(aes(colour = isloc, size = isloc)) +
  geom_text(aes(label = sprintf("%+.3f", rho_sen)), hjust = -0.22, vjust = -0.7,
            size = pt(5.8), family = FONT, colour = INK) +
  scale_colour_manual(values = c(`FALSE` = REPRO, `TRUE` = "#0B3B3A"), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 2.4), guide = "none") +
  annotate("text", x = -0.02, y = 16.7, label = "opposes senescence", hjust = 1,
           size = pt(5.8), family = FONT, colour = INK2) +
  annotate("text", x =  0.02, y = 16.7, label = "resembles", hjust = 0,
           size = pt(5.8), family = FONT, colour = INK2) +
  scale_x_continuous(limits = c(-0.40, 0.30), breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2)) +
  coord_cartesian(ylim = c(0.5, 17.0), clip = "off") +
  labs(x = "ρ with the senescence reference axis", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 5.9))

## ====== B. where a claim can and cannot be anchored ========================
WB <- read.delim(file.path(D, "within_between_by_class.tsv"))
WB <- WB[WB$class %in% names(CC), ]
WB$use <- ifelse(WB$b_lo > 0.10, "usable as a reference axis",
          ifelse(WB$b_hi < 0.10, "not usable", "indeterminate"))
WB$class <- factor(WB$class, levels = WB$class[order(WB$between)])
p5b <- ggplot(WB, aes(between, class, colour = use)) +
  annotate("rect", xmin = 0.10, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#EEF4F3", colour = NA) +
  geom_vline(xintercept = 0.10, colour = REPRO, linewidth = 0.4, linetype = "22") +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = b_lo, xmax = b_hi), height = 0, linewidth = 0.45) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`usable as a reference axis` = REPRO,
                                 `not usable` = GREY, indeterminate = AMBER), name = NULL) +
  guides(colour = guide_legend(nrow = 3)) +
  annotate("text", x = 0.105, y = 7.7, label = "ρ = 0.10", hjust = 0, size = pt(5.6),
           family = FONT, colour = REPRO) +
  scale_x_continuous(limits = c(-0.09, 0.60), breaks = c(0, 0.2, 0.4, 0.6)) +
  coord_cartesian(ylim = c(0.6, 8.0), clip = "off") +
  labs(x = "median between-study ρ  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(axis.text.y = element_text(size = 6.4),
        legend.position = "inside", legend.position.inside = c(0.60, 0.24),
        legend.box.background = element_blank(), legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = 6))

save_fig(lab_grid(p5a, p5b, labels = c("A","B"), ncol = 2, rel_widths = c(1.60, 1)),
         "Fig5.png", 183, 86)
