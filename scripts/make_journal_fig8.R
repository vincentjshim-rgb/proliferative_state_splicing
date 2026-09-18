FIGURES_WRITTEN <- c("Fig8.png")
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
D <- "public_data_tierA/derived"
MG <- read.delim(file.path(D, "magnitude_vs_alignment/magnitude_vs_alignment.tsv"))
MG <- MG[MG$id != "LOCAL_ReproCM", ]   # the anchor cannot be correlated with itself
MO <- read.delim(file.path(D, "what_drives_alignment/module_scores_vs_alignment.tsv"))
GC <- c(senescence = SEN, `donor age` = "#8C6D4F", `UV injury` = UV, injury = UV,
        photoprotection = AMBER, rescue = AMBER, reprogramming = NAVY, other = GREY)

## ============== A. perturbation magnitude does not predict alignment =======
mx <- grep("median_abs|magnitude", names(MG), value = TRUE)[1]
ry <- grep("rho_anchor|alignment", names(MG), value = TRUE)[1]
MG$cl <- if ("class" %in% names(MG)) MG$class else "other"
MG$col <- ifelse(MG$cl %in% names(GC), GC[MG$cl], GREY)
ct <- suppressWarnings(cor.test(MG[[mx]], MG[[ry]], method = "spearman", exact = FALSE))
dd  <- data.frame(x = MG[[mx]], y = MG[[ry]])
dd  <- dd[is.finite(dd$x) & is.finite(dd$y), ]
fit <- lm(y ~ x, dd)
nd  <- data.frame(x = seq(min(dd$x), max(dd$x), length.out = 80))
nd  <- cbind(nd, predict(fit, nd, interval = "confidence"))
sau <- MG[grepl("sauchinone", MG$id, ignore.case = TRUE), ]
p8a <- ggplot(MG, aes(.data[[mx]], .data[[ry]])) +
  geom_hline(yintercept = 0, colour = GREY_L, linewidth = 0.25) +
  geom_ribbon(data = nd, aes(x, ymin = lwr, ymax = upr), inherit.aes = FALSE,
              fill = GREY, alpha = 0.20) +
  geom_line(data = nd, aes(x, fit), inherit.aes = FALSE, colour = INK2, linewidth = 0.45,
            linetype = "22") +
  geom_point(size = 1.5, colour = MG$col) +
  geom_point(data = sau, size = 2.4, shape = 21, colour = INK, fill = NA, stroke = 0.45) +
  geom_text(data = sau, aes(label = "sauchinone"), hjust = 1.15, size = pt(5.8),
            family = FONT, colour = INK) +
  annotate("text", x = max(MG[[mx]]), y = 0.39, hjust = 1, vjust = 1, family = FONT, size = pt(6.3),
           colour = INK, lineheight = 1.05,
           label = sprintf("ρ = %+.3f\nP = %.2f  (n.s.)\n59 signatures", ct$estimate, ct$p.value)) +
  coord_cartesian(ylim = c(-0.22, 0.42)) +
  labs(x = expression(paste("perturbation magnitude   median |", log[2], "FC|")),
       y = "ρ with the Repro-CM anchor") +
  theme_j(7)

## ================ B. module content does predict alignment ================
mods <- c(HSF1 = "HSF1 heat-shock response", premRNA = "Capped pre-mRNA processing",
          splicing = "mRNA splicing", UPR_PERK = "UPR / PERK", prolif = "Cell cycle",
          ECM = "Extracellular matrix")
MO$cl <- if ("class" %in% names(MO)) MO$class else "other"
res <- do.call(rbind, lapply(names(mods), function(m) {
  ct <- suppressWarnings(cor.test(MO[[m]], MO$rho_anchor, method = "spearman", exact = FALSE))
  n  <- sum(is.finite(MO[[m]]) & is.finite(MO$rho_anchor)); se <- 1/sqrt(n - 3)
  data.frame(module = mods[[m]], rho = unname(ct$estimate), p = ct$p.value,
             lo = tanh(atanh(ct$estimate) - 1.96*se), hi = tanh(atanh(ct$estimate) + 1.96*se)) }))
ctm <- suppressWarnings(cor.test(MO$median_abs_logFC, MO$rho_anchor, method = "spearman", exact = FALSE))
nm <- sum(is.finite(MO$median_abs_logFC)); sem <- 1/sqrt(nm - 3)
res <- rbind(res, data.frame(module = "Perturbation magnitude", rho = unname(ctm$estimate),
             p = ctm$p.value, lo = tanh(atanh(ctm$estimate) - 1.96*sem),
             hi = tanh(atanh(ctm$estimate) + 1.96*sem)))
res <- res[order(res$rho), ]; res$module <- factor(res$module, levels = res$module)
res$fdr <- p.adjust(res$p, "BH"); res$star <- sig_stars(res$fdr)
res$fam <- ifelse(grepl("HSF1|UPR", res$module), "proteostatic stress",
           ifelse(grepl("pre-mRNA|splicing", res$module), "nuclear pre-mRNA processing", "other"))
p8b <- ggplot(res, aes(rho, module, colour = fam)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.35) +
  geom_point(size = 1.8, shape = 15) +
  geom_text(aes(label = sprintf("%+.3f %s", rho, star),
                hjust = ifelse(rho > 0, -0.12, 1.12)), size = pt(5.9),
            family = FONT, colour = INK) +
  scale_colour_manual(values = c(`proteostatic stress` = REPRO,
                     `nuclear pre-mRNA processing` = NAVY, other = GREY), name = NULL) +
  guides(colour = guide_legend(nrow = 3)) +
  scale_x_continuous(limits = c(-0.90, 1.02), breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75)) +
  labs(x = "ρ  module content vs alignment with the anchor  (95% CI)", y = NULL) +
  theme_j(7) + theme_caty() +
  theme(legend.position = "inside", legend.position.inside = c(0.23, 0.34),
        legend.box.background = element_blank(), legend.key.size = unit(6, "pt"))

save_fig(lab_grid(p8a, p8b, labels = c("A","B"), ncol = 2, rel_widths = c(1, 1.18)),
         "Fig8.png", 183, 76)
