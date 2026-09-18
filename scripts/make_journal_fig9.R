.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/journal_theme.R")
D <- "public_data_tierA/derived"; OUT <- file.path(D, "compendium_structure")
S <- readRDS(file.path(OUT, "structure.rds")); R <- S$R; meta <- S$meta
CP <- read.delim(file.path(D, "compendium/compendium_D_scores.tsv"))
meta$circ <- CP$circularity[match(meta$id, CP$id)]
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", `UV injury` = UV,
        `photoprotection / rescue` = AMBER, reprogramming = NAVY,
        `metabolic / culture` = GREY, `Repro-CM` = REPRO)

## =============== A. 60 x 60 similarity matrix, clustered ==================
Rf <- R; Rf[is.na(Rf)] <- 0
hc <- hclust(as.dist(1 - Rf), method = "average"); ord <- hc$order
ids <- rownames(R)[ord]; nn <- length(ids)
long <- expand.grid(xi = seq_len(nn), yi = seq_len(nn))
long$v <- mapply(function(a, b) R[ids[a], ids[b]], long$xi, long$yi)
long$yi <- nn + 1 - long$yi
ann <- data.frame(xi = seq_len(nn), cl = meta$class[match(ids, meta$id)])
p9a <- ggplot(long, aes(xi, yi, fill = v)) +
  geom_raster() +
  scale_fill_gradientn(colours = DIVERGE, limits = c(-0.55, 0.55), na.value = "#F2F4F6",
    breaks = c(-0.5, 0, 0.5), name = "Spearman ρ between signatures",
    guide = guide_colourbar(barwidth = unit(52,"pt"), barheight = unit(3.6,"pt"),
      ticks.colour = "white", frame.colour = NA, title.position = "top")) +
  annotate("rect", xmin = ann$xi - 0.5, xmax = ann$xi + 0.5,
           ymin = -2.6, ymax = -0.4, fill = CC[ann$cl], colour = NA) +
  annotate("text", x = -1.2, y = -1.5, label = "class", hjust = 1, size = pt(5.6),
           family = FONT, colour = INK2) +
  scale_x_continuous(expand = c(0,0), limits = c(-9, nn + 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(-3.2, nn + 0.5)) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_j(7) + theme(axis.line.x = element_blank(), axis.line.y = element_blank(),
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.x = element_blank(), axis.text.y = element_blank(),
                     legend.position = "bottom", legend.justification = "center",
                     legend.box.background = element_blank(),
                     legend.margin = margin(t = 2), plot.margin = margin(3, 3, 2, 3))

## ================== B. within-class coherence with CIs ====================
CH <- read.delim(file.path(OUT, "class_coherence_CI.tsv"))
CH$class <- factor(CH$class, levels = CH$class[order(CH$within)])
CH$star  <- sig_stars(CH$FDR)
p9b <- ggplot(CH, aes(within, class, colour = class)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(aes(x = between, xend = within, yend = class), linewidth = 0.4, colour = GREY_L) +
  geom_point(aes(x = between), size = 1.2, shape = 1, colour = GREY) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.38) +
  geom_point(size = 1.9) +
  geom_text(aes(label = sprintf("%+.3f %s", within, star)), hjust = -0.18, vjust = -0.9,
            size = pt(5.9), family = FONT, colour = INK) +
  geom_text(aes(label = sprintf("n = %d", n_sig), x = -0.10), hjust = 0,
            size = pt(5.5), family = FONT, colour = GREY) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(limits = c(-0.12, 0.80), breaks = c(0, 0.25, 0.5, 0.75)) +
  annotate("text", x = 0.78, y = 1.4, hjust = 1, family = FONT, size = pt(5.6), colour = GREY,
           label = "open circles: median similarity to all other classes") +
  labs(x = "median within-class similarity  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y = element_text(size = 6.4))

## ========== C. classical MDS of the 60 signatures, all points shown =======
meta$cl <- factor(meta$class, levels = names(CC))
meta$ind <- meta$circ == "independent"
lab_ids <- c("LOCAL_ReproCM","GSE109700_deep","GSE179848_ContactInhibition",
             "GSE326951_sauchinone","GSE240226_UVA","GSE165177_MPTR","GSE113957_age")
LB <- c(LOCAL_ReproCM="Repro-CM", GSE109700_deep="deep senescence",
        GSE179848_ContactInhibition="contact inhibition", GSE326951_sauchinone="sauchinone",
        GSE240226_UVA="UVA (defines axis)", GSE165177_MPTR="MPTR", GSE113957_age="donor age")
meta$lab <- LB[meta$id]
p9c <- ggplot(meta, aes(MDS1, MDS2)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_point(aes(colour = cl, shape = ind), size = 1.6) +
  geom_point(data = meta[meta$id %in% lab_ids, ], size = 2.5, shape = 21,
             colour = INK, fill = NA, stroke = 0.45) +
  geom_text(data = meta[meta$id %in% lab_ids, ], aes(label = lab),
            size = pt(5.8), family = FONT, colour = INK, vjust = -1.0) +
  scale_colour_manual(values = CC, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.5))) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.7, family = FONT,
           size = pt(5.8), colour = GREY,
           label = "axis 1 correlates with the senescence axis at r = +0.82,\nwith the UVA axis at r = +0.25") +
  labs(x = "MDS axis 1  (12.8% of eigenvalue mass)", y = "MDS axis 2  (9.6%)") +
  theme_j(7) + theme(legend.position = "inside", legend.position.inside = c(0.30, 0.90),
                     legend.box.background = element_blank(), legend.key.size = unit(6, "pt"))

## ===== D. Repro-CM's neighbourhood, independent contrasts only ============
nb <- read.delim(file.path(OUT, "reproCM_neighbours.tsv"))
nb$circ <- CP$circularity[match(nb$id, CP$id)]
nb <- nb[is.finite(nb$r) & nb$circ == "independent", ]
nb$cl <- factor(nb$class, levels = names(sort(tapply(nb$r, nb$class, median))))
ph <- nb$r[nb$class == "photoprotection / rescue"]; se <- nb$r[nb$class == "senescence"]
wt <- wilcox.test(ph, se)
p9d <- ggplot(nb, aes(cl, r)) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.3) +
  geom_boxplot(aes(colour = cl), outlier.shape = NA, width = 0.52, linewidth = 0.34, fill = NA) +
  geom_point(aes(colour = cl), size = 1.3, position = position_jitter(width = 0.15, seed = 5)) +
  brk(1, 6, 0.46, sprintf("%s   Wilcoxon P = %.4f", sig_stars(wt$p.value), wt$p.value),
      tick = 0.022, size = 5.6) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_y_continuous(limits = c(-0.22, 0.54)) +
  labs(x = NULL, y = "ρ with the Repro-CM signature") +
  theme_j(7) + theme(axis.text.x = element_text(angle = 24, hjust = 1, size = 6.2))

top <- lab_grid(p9a, p9b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1.04))
bot <- lab_grid(p9c, p9d, labels = c("C","D"), ncol = 2, rel_widths = c(1.25, 1))
save_fig(lab_grid(top, bot, labels = c("",""), ncol = 1, rel_heights = c(1, 0.92)),
         "Fig9_structure.png", 183, 150)
