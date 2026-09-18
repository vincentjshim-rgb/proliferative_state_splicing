FIGURES_WRITTEN <- c("Fig1.png")
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
PUBDIR <- "public_data_tierA/derived/figures_benchmark"; dir.create(PUBDIR, showWarnings = FALSE)
D <- "public_data_tierA/derived/benchmark_figs"
Z <- readRDS(file.path(D, "figdata.rds")); R <- Z$R; M <- Z$meta
CC <- c(senescence = SEN, `donor age` = "#8C6D4F", secretome = REPRO,
        `photoprotection / rescue` = AMBER, `UV injury` = UV,
        `metabolic / culture` = GREY, reprogramming = NAVY, `Repro-CM` = "#0B3B3A")
AR <- arrow(length = unit(1.7, "pt"), type = "closed")

## ============================= A. design schematic ==========================
bx <- function(x,y,w,h,lab,fill="white",col=INK,fc=INK,tsz=6.1,lw=0.32)
  data.frame(x,y,w,h,lab,fill,col,fc,tsz,lw,stringsAsFactors=FALSE)
B <- rbind(
  bx(50, 95, 92, 8, "29 GEO series + this study   ·   human fibroblasts only"),
  bx(50, 79, 92, 9, "76 perturbation contrasts, each a genome-wide log₂FC vector"),
  bx(50, 62, 92, 9, "pairwise Spearman ρ over shared genes  (≥ 1,000 required)"),
  bx(22, 43, 40, 13, "WITHIN a study\nsame experiment,\ndifferent arm", "#EFF3F5", INK2, INK2, 6.0),
  bx(72, 43, 40, 13, "BETWEEN studies\ndifferent laboratory,\nsame perturbation class", "#E4F1F0", REPRO, REPRO, 6.0),
  bx(50, 23, 92, 9, "reproducibility  =  median between-study ρ within a class"),
  bx(50,  6, 92, 10, "controls:  effect magnitude · sequencing platform · gene coverage", "white", INK, INK, 6.0, 0.5))
B$xmin <- B$x-B$w/2; B$xmax <- B$x+B$w/2; B$ymin <- B$y-B$h/2; B$ymax <- B$y+B$h/2
sg <- function(x,xend,y,yend) data.frame(x,xend,y,yend)
S <- rbind(sg(50,50,91,83.5), sg(50,50,74.5,66.5), sg(50,50,57.5,53), sg(22,72,53,53),
           sg(22,22,53,49.5), sg(72,72,53,49.5), sg(22,22,36.5,32), sg(72,72,36.5,32),
           sg(22,72,32,32), sg(50,50,32,27.5), sg(50,50,18.5,11))
p1a <- ggplot() +
  geom_segment(data=S, aes(x,y,xend=xend,yend=yend), colour=GREY, linewidth=0.3, arrow=AR) +
  geom_rect(data=B, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=B$fill, colour=B$col, linewidth=B$lw) +
  geom_text(data=B, aes(x,y,label=lab), family=FONT, colour=B$fc, size=pt(B$tsz), lineheight=0.95) +
  scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(0,100)) +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA),
                       plot.margin=margin(3,3,3,3))

## ====================== B. 76 x 76 clustered similarity =====================
Rf <- R; Rf[is.na(Rf)] <- 0
hc <- hclust(as.dist(1-Rf), method="average"); ord <- hc$order
ids <- rownames(R)[ord]; nn <- length(ids)
lg <- expand.grid(xi=seq_len(nn), yi=seq_len(nn))
lg$v <- mapply(function(a,b) R[ids[a], ids[b]], lg$xi, lg$yi)
lg$yi <- nn + 1 - lg$yi
ann <- data.frame(xi=seq_len(nn), cl=M$class[match(ids, M$id)])
p1b <- ggplot(lg, aes(xi, yi, fill=v)) + geom_raster() +
  scale_fill_gradientn(colours=DIVERGE, limits=c(-0.6,0.6), na.value="#F2F4F6",
    breaks=c(-0.5,0,0.5), name="Spearman ρ between contrasts",
    guide=guide_colourbar(barwidth=unit(56,"pt"), barheight=unit(3.6,"pt"),
      ticks.colour="white", frame.colour=NA, title.position="top")) +
  annotate("rect", xmin=ann$xi-0.5, xmax=ann$xi+0.5, ymin=-3.0, ymax=-0.5,
           fill=CC[ann$cl], colour=NA) +
  annotate("text", x=-1.5, y=-1.75, label="class", hjust=1, size=pt(5.6), family=FONT, colour=INK2) +
  scale_x_continuous(expand=c(0,0), limits=c(-11, nn+0.5)) +
  scale_y_continuous(expand=c(0,0), limits=c(-3.6, nn+0.5)) +
  coord_fixed(clip="off") + labs(x=NULL, y=NULL) +
  theme_j(7) + theme(axis.line.x=element_blank(), axis.line.y=element_blank(),
    axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank(),
    legend.position="bottom", legend.justification="center",
    legend.box.background=element_blank(), legend.margin=margin(t=2),
    plot.margin=margin(3,3,2,3))

## ========================= C. composition of the panel ======================
cmp <- do.call(rbind, lapply(names(CC), function(k) {
  i <- which(M$class == k); if (!length(i)) return(NULL)
  data.frame(class=k, n_sig=length(i), n_study=length(unique(M$dataset[i]))) }))
cmp <- cmp[order(cmp$n_sig), ]; cmp$class <- factor(cmp$class, levels=cmp$class)
cl <- rbind(transform(cmp, v=n_sig,   k="contrasts"),
            transform(cmp, v=n_study, k="independent studies"))
cl$k <- factor(cl$k, levels=c("contrasts","independent studies"))
cl$fillv <- ifelse(cl$k == "contrasts", as.character(cl$class), "studies")
p1c <- ggplot(cl, aes(v, class, fill=fillv)) +
  geom_col(position=position_dodge(width=0.7), width=0.62, colour=NA) +
  geom_text(aes(label=v), position=position_dodge(width=0.7), hjust=-0.3,
            size=pt(5.8), family=FONT, colour=INK) +
  scale_fill_manual(values=c(CC, studies="#D2D8DE"), guide="none") +
  annotate("point", x=15.5, y=1.2, size=1.6, shape=15, colour="#D2D8DE") +
  annotate("text", x=16.0, y=1.2, label="independent studies", hjust=0, size=pt(5.6),
           family=FONT, colour=INK2) +
  annotate("text", x=16.0, y=2.2, label="contrasts (coloured as in B)", hjust=0, size=pt(5.6),
           family=FONT, colour=INK2) +
  scale_x_continuous(limits=c(0,26), breaks=c(0,5,10,15)) +
  labs(x="number", y=NULL) +
  theme_j(7) + theme_caty() + theme(axis.text.y=element_text(size=6.3))

top <- lab_grid(p1a, p1b, labels=c("A","B"), ncol=2, rel_widths=c(1, 1.12))
save_fig(lab_grid(top, p1c, labels=c("","C"), ncol=1, rel_heights=c(1, 0.60)),
         "Fig1.png", 183, 136)
