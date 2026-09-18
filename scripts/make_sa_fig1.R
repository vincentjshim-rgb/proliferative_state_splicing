## Fig. 1  Do splicing-machinery transcripts fall with donor age and time in culture?
## Revised 2026-09-15: panel c is rebuilt on the redefined cohort (107 normal donors
## aged 20-96, recount3 matrix) so that Figs. 1, 3 and 4 share one matrix and one
## cohort; panel d prints exact P values for the short series.
.libPaths(c(normalizePath("analysis_r_lib"), .libPaths()))
source("scripts/sciadv_theme.R")
D <- "public_data_tierA/derived"
CORE <- readLines(file.path(D, "conserved_core/age_down_splicing_core_v2.txt"))
AR <- arrow(length = unit(1.8, "pt"), type = "closed")

## ================== a. model schematic =====================================
bx <- function(x,y,w,h,lab,fill="white",col=INK,fc=INK,tsz=7)
  data.frame(x,y,w,h,lab,fill,col,fc,tsz,stringsAsFactors=FALSE)
B <- rbind(
  bx(50, 92, 84, 9,  "Human dermal fibroblast in culture"),
  bx(22, 62, 38, 13, "Replicative ageing\ntime in culture", "#FBECEA", RED, RED),
  bx(78, 62, 38, 13, "Donor age\n20 to 96 years", "#F6F1EA", BROWN, BROWN),
  bx(50, 34, 86, 11, "splicing-factor transcript abundance falls", "#F2F4F6", INK, INK, 7.2),
  bx(50, 11, 86, 12, "a property of ageing,\nor a readout of how fast the cells divide?", "white", INK, INK, 7.2))
B$xmin<-B$x-B$w/2; B$xmax<-B$x+B$w/2; B$ymin<-B$y-B$h/2; B$ymax<-B$y+B$h/2
SG <- rbind(data.frame(x=50,xend=50,y=87.5,yend=80), data.frame(x=22,xend=78,y=80,yend=80),
            data.frame(x=22,xend=22,y=80,yend=68.8), data.frame(x=78,xend=78,y=80,yend=68.8),
            data.frame(x=22,xend=22,y=55.5,yend=39.5), data.frame(x=78,xend=78,y=55.5,yend=39.5),
            data.frame(x=22,xend=78,y=39.5,yend=39.5), data.frame(x=50,xend=50,y=28.5,yend=17.4))
p1a <- ggplot() +
  geom_segment(data=SG, aes(x,y,xend=xend,yend=yend), colour=GREY, linewidth=0.32) +
  geom_segment(data=data.frame(x=c(22,78,50), xend=c(22,78,50),
               y=c(70.5,70.5,30.5), yend=c(68.8,68.8,17.4)),
               aes(x,y,xend=xend,yend=yend), colour=GREY, linewidth=0.32, arrow=AR) +
  geom_rect(data=B, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=B$fill, colour=B$col, linewidth=0.35) +
  geom_text(data=B, aes(x,y,label=lab), family=FONT, colour=B$fc, size=pt(B$tsz), lineheight=0.95) +
  annotate("text", x=50, y=47, label="both reduce", family=FONT, size=pt(7), colour=GREY,
           fontface="italic") +
  scale_x_continuous(limits=c(0,100)) + scale_y_continuous(limits=c(3,98)) +
  theme_void() + theme(plot.background=element_rect(fill="white",colour=NA),
                       plot.margin=margin(3,3,3,3))

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
  scale_x_continuous(limits=c(0,142), breaks=seq(0,100,25), expand=expansion(mult=c(0,0))) +
  labs(x="% of pathway falling in all three series", y=NULL) +
  theme_sa(8) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(),
                      axis.text.y=element_text(size=7, lineheight=0.95),
        legend.position="top", legend.justification="left", legend.key.size=unit(6,"pt"),
        legend.margin=margin(b=-4), legend.text=element_text(size=7.5))

## ================ c. donor age, redefined cohort, 107 donors ===============
A <- read.delim(file.path(D, "cohort_revised/primary/sample_metrics.tsv"))
ct1 <- cor.test(A$age, A$machinery)                       # Pearson: the panel shows a linear fit
f1 <- lm(machinery ~ age, A)
nd <- data.frame(age = seq(min(A$age), max(A$age), length.out = 60))
nd <- cbind(nd, predict(f1, nd, interval = "confidence"))
p1c <- ggplot(A, aes(age, machinery)) +
  geom_ribbon(data=nd, aes(age, ymin=lwr, ymax=upr), inherit.aes=FALSE, fill=BROWN, alpha=0.18) +
  geom_line(data=nd, aes(age, fit), inherit.aes=FALSE, colour=BROWN, linewidth=0.6) +
  geom_point(size=1.5, colour=BROWN, alpha=0.85) +
  annotate("text", x=21, y=min(A$machinery), hjust=0, vjust=0, family=FONT, size=pt(7.5),
           colour=INK, lineheight=1.1, label=rp(ct1$estimate, ct1$p.value, lab="r")) +
  annotate("text", x=96, y=max(A$machinery)*1.03, hjust=1, vjust=1, family=FONT, size=pt(7),
           colour=GREY, label=sprintf("n = %d donors", nrow(A))) +
  scale_y_continuous(labels = num_axis()) +
  labs(x="donor age (years)", y="pre-mRNA processing score") +
  theme_sa(8)

## ======= d. time in culture, GSE179848, four donor lines ===================
L  <- read.delim(file.path(D, "core_depth/GSE179848_longitudinal_core.tsv"))
CT <- read.delim(file.path(D, "revision_stats/fig1d_culture_time.tsv"))   # exact P where n < 10
st <- do.call(rbind, lapply(sort(unique(L$line)), function(l) {
  s <- L[L$line == l, ]; r <- CT[CT$line == l, ]
  data.frame(line = l, lab = sprintf("%s\nn = %d,  %d days", l, r$n, r$days),
             txt = rp(r$rho, r$p), x = max(s$days), y = max(L$score)) }))
L$lab <- st$lab[match(L$line, st$line)]
p1d <- ggplot(L, aes(days, score)) +
  geom_smooth(method="lm", se=TRUE, colour=RED, fill=RED, alpha=0.16, linewidth=0.6) +
  geom_point(size=1.4, colour=RED, alpha=0.85) +
  geom_text(data=st, aes(x=x, y=y, label=txt), hjust=1, vjust=1, size=pt(7),
            family=FONT, colour=INK, lineheight=1.1) +
  facet_wrap(~lab, nrow=1, scales="free_x") +
  scale_y_continuous(labels = num_axis()) +
  labs(x="days in culture", y="pre-mRNA processing score") +
  theme_sa(8) + theme(strip.text=element_text(size=7.5, lineheight=1.05))

top <- lab_grid(p1a, p1b, labels=c("a","b"), ncol=2, rel_widths=c(0.85,1))
save_fig(lab_grid(top, lab_grid(p1c, p1d, labels=c("c","d"), ncol=2, rel_widths=c(0.62,1)),
                  labels=c("",""), ncol=1, rel_heights=c(1,0.86)), "Fig1.png", 183, 124)
