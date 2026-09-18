---
title: "Splicing machinery expression reports proliferative state, not ageing, in cultured human fibroblasts"
subtitle: "Draft v2 · target journal npj Aging (Article)"
date: "12 September 2026"
---

# Abstract

Loss of splicing homeostasis is emerging as a hallmark of ageing, established most convincingly by outcome-based measures in tissues. A large parallel literature instead uses the expression of spliceosome and pre-mRNA processing genes in cultured cells as a proxy for that loss. Because these genes are among the most tightly cell-cycle-regulated in the genome, and because proliferation rate is the principal variable that changes as cultured fibroblasts age, we tested whether the proxy is valid. We defined a 96-gene pre-mRNA processing module from genes whose direction of change was unanimous across four independent human fibroblast ageing cohorts. The module behaves exactly as reported: it falls with donor age across 97 donors (*r* = −0.35, *P* = 4 × 10⁻⁴) and declines steadily over 200 days of culture in three of four healthy donor lines (ρ = −0.78 to −0.96). It also tracks the proliferation marker MKI67 at ρ = 0.90 in untreated cells, and contact inhibition — a reversible arrest that is not senescence — lowers it as much as ageing does (*P* = 5 × 10⁻³). Including MKI67 in the model absorbs 72% of the effect of time in culture, and at gene level 79 of 96 module genes show an age association that only 27 retain after adjustment, with survival predicted by each gene's own coupling to proliferation (ρ = −0.39, *P* = 7 × 10⁻⁵). Decisively, in 142 donors for whom both measures could be computed from the same libraries, an outcome-based measure of splicing precision rises with age (ρ = 0.38, *P* = 3 × 10⁻⁶) and is entirely unaffected by the same adjustment, while machinery expression loses 67% of its age effect. Across 65 contrasts from 29 studies, a cell-cycle module explains 87% of the variance in the module, and no perturbation class retains a residual; the largest apparent restoration by a cell-free secretome, +0.27, becomes −0.04 once proliferation is accounted for. Splicing machinery expression in cultured fibroblasts is therefore a proliferation readout and should not be used as evidence for age-associated splicing dysfunction; outcome-based measures in the same cells are not subject to this confounding.

# Introduction

Changes in RNA splicing with age are attracting sustained attention, and loss of splicing homeostasis has recently been formalised as a hallmark of ageing¹. The strongest evidence is outcome-based. Integrative analysis of mouse and human tissues shows a systematic deterioration in the fidelity of splicing with age — a rising proportion of isoforms carrying premature stop codons, frameshifts or disrupted protein domains — which increases with age, is alleviated by calorie restriction and rapamycin, and is mechanistically associated with age-related changes in specific splicing factors¹⁹.

Alongside that work sits a much larger literature that uses a different and more convenient measurement: the expression level of spliceosome components, splicing factors and other pre-mRNA processing genes, usually in cultured cells²⁻⁴. These transcripts fall in old donors, in late-passage cultures and in senescent cells with great consistency, and the decline is commonly read as a loss of splicing capacity.

There is a structural reason to check whether the two measurements can be used interchangeably. Spliceosome assembly, capping and 3′-end processing are required in proportion to the nascent transcription a cell performs, and the genes encoding this machinery are transcriptionally coupled to the cell cycle, rising in S phase and falling when cells stop dividing⁵. Cultured fibroblasts slow and then stop dividing as they age, and donor age, replicative exhaustion, senescence induction and most interventions all change proliferation rate. A transcript set that is cell-cycle-regulated, measured in a system whose proliferation rate is the principal thing that changes with age, will decline with age whether or not splicing biology is involved.

The same hazard has been met elsewhere and dealt with. In breast cancer, proliferation is the dominant axis of tumour expression and the reason that most randomly chosen signatures predict outcome⁶,⁷. In blood epigenetic clocks, apparent age acceleration is decomposed against cell composition before interpretation⁸. In the senescence field the accepted control is experimental: a quiescence arm at matched population doubling, or direct verification of arrest, separates a state-specific change from a consequence of having stopped dividing⁹⁻¹¹. The tissue-level splicing work notes explicitly that shifts in cell composition were not controlled for¹⁹, a caveat of the same family.

Here we apply that control to the machinery-expression proxy in the system where it is most used. We show that the proxy behaves as reported, that it is nonetheless a readout of proliferative state, and — importantly — that an outcome-based measure computed from the very same libraries is not, so that the tissue-level conclusions of the splicing-degeneration literature are unaffected by what we report.

# Results

## A pre-mRNA processing module defined by agreement across four ageing cohorts

To avoid defining the module from a single study we assembled 76 perturbation contrasts in human fibroblasts from 29 GEO series, grouped them by what was done to the cells, and asked which genes changed in a consistent direction across the independent studies of each class (Fig. 1A). Four independent donor-age cohorts agreed in direction on 3,444 genes, 2.2-fold more than expected by chance. Genes falling with age were enriched overwhelmingly for one process: 96 of 277 genes in "Processing of capped intron-containing pre-mRNA" (BH-FDR 2 × 10⁻¹⁸), 76 of 193 in mRNA splicing (3 × 10⁻¹⁸) and 60 of 148 in mRNA 3′-end processing (4 × 10⁻¹⁵) (Fig. 1B). We took the 96 pre-mRNA processing genes as the module and scored it in individual samples as the mean *z*-score across genes.

The module reproduces the published behaviour. Across 97 primary dermal fibroblast donors aged 22 to 89 years¹⁷ it falls with age (*r* = −0.35, *P* = 4 × 10⁻⁴; Fig. 1C). In a longitudinal resource in which four healthy donor lines were cultured under identical conditions and sampled repeatedly for up to 200 days¹⁸, it declines steadily in three lines (ρ = −0.89, −0.78 and −0.96; all *P* < 0.001) and shows no trend in the fourth, which was sampled only to day 55 (Fig. 1D). Across all 71 samples the slope is −0.0099 units per day (*P* = 1 × 10⁻¹⁵).

![](<public_data_tierA/derived/figures_sciadv/Fig1.png>){width=6.3in}

**Fig. 1. A pre-mRNA processing module declines with donor age and with time in culture.** (**A**) The question. (**B**) Reactome pathways enriched among genes falling with donor age in all four independent cohorts; orange, pre-mRNA processing terms. (**C**) Module score against donor age in 97 primary dermal fibroblast donors; line and band are a least-squares fit with its 95% confidence interval. (**D**) Module score against days in culture in four healthy donor lines profiled longitudinally under identical conditions; Spearman ρ and *P* are given per line.

## The module tracks Ki-67 and falls with reversible arrest

The longitudinal resource allows the control that the senescence field uses. Across all 345 samples, spanning nine treatments, two oxygen tensions and a mitochondrial mutation, the module correlates with MKI67 expression at ρ = 0.80 (*P* = 1 × 10⁻⁷⁶; Fig. 2A). Restricted to untreated cells from healthy donors, where the only systematic variable is time in culture, the correlation is ρ = 0.90 (*P* = 4 × 10⁻²⁷; Fig. 2B). For comparison, the module correlates with CDKN2A (p16) at ρ = −0.38 and with CDKN1A (p21) at ρ = −0.36.

Correlation with a proliferation marker does not establish direction, so we used a manipulation that arrests cells without making them senescent. Contact inhibition, in the same donors and the same laboratory, reduced MKI67 as expected (*P* = 3 × 10⁻⁸) and reduced the module to the level seen in aged cultures (median −0.42 against +0.18 in untreated cells, *P* = 5 × 10⁻³; Fig. 2C). These cells are not senescent by the usual transcriptional criterion: CDKN2A is only marginally higher than in untreated cells, and the arrest is reversible by design. Adding MKI67 to the model of module score against time in culture absorbs 72% of the time effect, reducing the slope from −0.0099 per day (*P* = 1 × 10⁻¹⁵) to −0.0028 (*P* = 0.011), while MKI67 itself remains strongly associated (*P* = 2 × 10⁻¹²; Fig. 2D).

![](<public_data_tierA/derived/figures_sciadv/Fig2.png>){width=6.3in}

**Fig. 2. The module is a readout of proliferative state.** (**A**) Module score against MKI67 across all 345 samples of the longitudinal resource. (**B**) The same restricted to untreated cells from four healthy donors; symbols distinguish donors. (**C**) Contact inhibition compared with untreated cells in the same donors: module score, MKI67 and CDKN2A. Boxes show median and interquartile range, points are individual samples; **, *P* < 0.01; ****, *P* < 10⁻⁴ (Wilcoxon rank-sum). (**D**) Slope of module score per day in culture, before and after including MKI67; error bars are 95% confidence intervals.

## Which genes survive adjustment is set by their own coupling to proliferation

If the module is a proliferation readout, individual genes should differ in how much of their age association survives, in proportion to how tightly each is tied to the cell cycle. In a 142-donor cohort, 79 of the 96 module genes show an age association at BH-FDR < 0.05 and 27 retain one after a proliferation score is included — a loss of 66% (Fig. 3D). The genes are individually well coupled to proliferation, with a median correlation of 0.68 and a range from −0.63 to 0.92 (Fig. 3B), and the fraction of the age effect a gene retains falls as that coupling rises (ρ = −0.39, *P* = 7 × 10⁻⁵; Fig. 3C).

The 27 survivors are not a random subset: their median coupling to proliferation is 0.60 against 0.69 for the rest (*P* = 0.022). They include core spliceosomal and 3′-processing components — SF1, SF3A3, SF3B5, SNRPA, SNRPB, SNRPD2, SNRPF, PRPF3, PPIE, PPIH, RBM39, SRSF9, BUD13, GPKOW, ISY1 and CTNNBL1 among them — and we provide the list as a proliferation-independent subset for studies that require a machinery-expression readout in this system.

![](<public_data_tierA/derived/figures_sciadv/Fig3.png>){width=6.3in}

**Fig. 3. Gene-level dissection of the module.** (**A**) Expression of the 96 module genes across 142 donors ordered young to old; genes are ordered by their own correlation with proliferation, and annotation bars show donor age and the proliferation score. (**B**) Distribution of each gene's correlation with proliferation. (**C**) Fraction of the age effect retained after adjustment against that coupling; blue, genes retaining an association at BH-FDR < 0.05. (**D**) Number of module genes with an age association before and after adjustment.

## Splicing outcome in the same donors is unaffected by the adjustment

Machinery expression and the outcome of splicing are different measurements, and the tissue-level evidence for age-associated splicing degeneration rests on outcome¹⁹. Because the 142-donor cohort was uniformly reprocessed with junction-level output, both could be computed from the same libraries and subjected to the same adjustment (Fig. 4A).

They behave in opposite ways. Machinery expression is nearly collinear with proliferation in these donors (*r* = 0.86, *P* = 1 × 10⁻⁴³; Fig. 4B), whereas the fraction of junction-spanning reads falling on unannotated junctions — a direct readout of splicing precision — is only weakly related to it (*r* = −0.23, *P* = 0.007; Fig. 4C). Unannotated read usage rises with donor age (ρ = 0.38, *P* = 3 × 10⁻⁶; Fig. 4E), the direction expected if splicing precision deteriorates.

Adjustment separates them decisively. The age effect on machinery expression falls by 67%, from *P* = 8 × 10⁻⁷ to *P* = 0.007, whereas the age effect on splicing outcome is unchanged, from *P* = 0.002 to *P* = 0.004, with a coefficient 2% larger rather than smaller (Fig. 4D). A third metric, the entropy of junction usage at 5′ splice sites, is intermediate: significant unadjusted (*P* = 0.001) and not after adjustment (*P* = 0.12).

Splicing outcome in these cells therefore carries an age association that proliferation does not explain. Machinery expression largely does not.

![](<public_data_tierA/derived/figures_sciadv/Fig4.png>){width=6.3in}

**Fig. 4. Splicing outcome and machinery expression dissociate in the same donors.** (**A**) Both measures were computed from the same uniformly reprocessed libraries of 142 donors. (**B**) Machinery expression against the proliferation score. (**C**) Splicing outcome, the fraction of junction reads on unannotated junctions, against the same score. (**D**) Age effect for each measure before and after proliferation adjustment, scaled to its unadjusted estimate. (**E**) Splicing outcome against donor age.

## Across 29 studies the cell cycle accounts for 87% of the variation

If the module is a proliferation readout, any perturbation that changes proliferation should move it by a predictable amount irrespective of what the perturbation is. We scored each of 65 contrasts for its effect on the module and, independently, on a 533-gene cell-cycle module with which it shares only eight genes.

The two move together almost deterministically: ρ = 0.75 with a linear fit of *R*² = 0.87 (*P* = 5 × 10⁻¹³; Fig. 5A). Senescence contrasts sit at the extreme negative end and ageing cohorts beside them, but on the same line as everything else, and the relationship holds within classes examined separately (within-class ρ 0.19 to 0.95; Fig. 5B). Consequently nothing is left over: after adjustment no perturbation class — senescence, donor age, ultraviolet injury, photoprotection, reprogramming media, metabolic perturbation or secretome — retains a significant residual (Fig. 5C). The median residual for senescence is −0.04 against an apparent effect of −0.78.

![](<public_data_tierA/derived/figures_sciadv/Fig5.png>){width=6.3in}

**Fig. 5. The relationship generalises across 29 studies and leaves no residual.** (**A**) Change in the pre-mRNA processing module against change in a cell-cycle module for 65 contrasts; line and band are a least-squares fit with its 95% confidence interval. (**B**) The same within each perturbation class on common axes; the dashed line is the fit from (**A**). (**C**) Residuals after adjustment, by class; ns, not significant (Wilcoxon signed-rank against zero, BH-adjusted).

## Interventions that appear to restore splicing capacity are restoring proliferation

The practical consequence concerns interventions. We scored 15 cell-free secretome preparations — conditioned media, extracellular-vesicle fractions and extruded nanovesicles from bone-marrow and adipose mesenchymal stromal cells, trophoblast stem cells, endothelial cells, senescent fibroblasts, platelet-rich fibrin and a transient reprogramming phase — applied to recipient fibroblasts in eight independent studies.

Ten of the fifteen raise the module, and the apparent effects are sizeable, but plotted against the change predicted from their effect on the cell cycle alone they fall on the identity line (Fig. 6A). The largest apparent effect, +0.27 for endothelial-cell extracellular vesicles, becomes −0.04 after adjustment; the second largest, +0.10 for a primary mesenchymal stromal cell secretome, becomes −0.13 (Fig. 6C). No secretome has a residual outside the 95% range of all 65 contrasts (Fig. 6B). Two preparations sit at the edges of that range in opposite directions and are reported for completeness rather than as findings: a platelet-rich fibrin serum applied to gingival fibroblasts (+0.10) and a secretome collected during the transient reprogramming phase of iPSC induction (+0.07), the latter from a single library per condition in our own data¹², against a trophoblast stem cell secretome at −0.23.

![](<public_data_tierA/derived/figures_sciadv/Fig6.png>){width=6.3in}

**Fig. 6. Apparent restoration of pre-mRNA processing by secretome interventions is explained by proliferation.** (**A**) Observed change against the change predicted from the cell-cycle module for 15 secretome contrasts; vertical lines are residuals. (**B**) Residuals ranked; shading and dashed lines mark the 95% range of residuals across all 65 contrasts. (**C**) The six largest apparent effects before and after adjustment.

# Discussion

The decline of pre-mRNA processing transcripts with age in cultured fibroblasts is real and reproducible, and on the evidence here it is not informative about splicing. It is present in donor cohorts and in longitudinal culture; it is equally present when cells are simply prevented from dividing by contact with their neighbours; it tracks Ki-67 at ρ = 0.90; and across 29 studies it is predicted by a cell-cycle module with *R*² = 0.87, after which no class of perturbation retains an effect.

This does not conflict with the case for splicing as a hallmark of ageing. That case rests on outcome — the proportion of isoforms with disrupted reading frames and domains, measured across mouse and human tissues¹⁹ — and our own outcome measure, computed from the same libraries as the expression measure, behaves in the opposite way: it rises with age and is untouched by the adjustment that removes two-thirds of the expression effect. The distinction we draw is between two measurements that are often treated as interchangeable, and the recommendation follows directly: in cultured proliferating cells, report outcome.

Two features make the confounding hard to escape in this system. The association is not weak or conditional: at ρ = 0.90 between the module and Ki-67 in untreated cells there is little variance left in which a splicing-specific signal could live. And the interventions studied here almost all change proliferation — conditioned media and vesicle preparations are developed and selected for their ability to promote growth and repair — so it is expected, not surprising, that they raise a cell-cycle-linked module, with apparent restoration of splicing capacity following arithmetically.

The remedy is the one the senescence field already applies to other claims. A quiescence or contact-inhibition arm at matched population doubling separates arrest from the state of interest at the design stage, and a proliferation readout reported alongside the module lets a reader judge the claim. Where only public data are available, the cell-cycle adjustment used here is an approximate substitute, and we provide the module, the 27-gene proliferation-independent subset and the reference relationship for that purpose.

The study has clear limits. It concerns cultured human fibroblasts and does not speak to tissues, to post-mitotic cells, or to systems in which proliferation is held constant; in a non-dividing tissue the confounding we describe cannot arise and machinery expression may well be interpretable. Our cell-cycle and pre-mRNA processing modules are both expression measurements of a shared transcriptional programme rather than independent variables, which is why the argument is anchored on the contact-inhibition experiment rather than on the regression. The outcome metric we could compute from deposited junctions is a measure of splicing precision and not the damaging-isoform proportion used in tissue work, so the two are related but not identical. Contrasts were assembled from processed data deposited by others, with the heterogeneity that entails. Finally, one donor-age cohort in the compendium has processing batches that coincide almost exactly with age strata and was excluded from sample-level analysis for that reason — a reminder that proliferation is not the only confounder in this field.

# Methods

## Contrast assembly

Seventy-six perturbation contrasts in human fibroblasts were assembled from 29 GEO series and one local experiment (Supplementary Table 1). Where processed count matrices were deposited, counts were TMM-normalised¹³,¹⁴ and a contrast was computed as the difference in mean log₂(CPM + 0.5) between arms, retaining genes with CPM ≥ 1 in at least two samples; where only a differential-expression table was deposited the published log₂ fold change was used. Identifiers were mapped to HGNC symbols and sample-to-arm assignments were verified manually against series metadata. Classes were assigned by experimental manipulation, never by transcriptional content.

## Module definition

Within each class, contrasts from the same study were averaged, each study vector was scaled to unit median absolute deviation, and a gene was called unanimous if every study measuring it agreed in direction. Enrichment among unanimous genes was assessed by hypergeometric test against Reactome¹⁵ sets of 15 to 500 genes restricted to the genes tested, with Benjamini–Hochberg correction¹⁶. The module is the 96 genes of "Processing of capped intron-containing pre-mRNA" that fell with age in all four cohorts. The cell-cycle module is the union of Reactome "Cell Cycle, Mitotic", "DNA Replication" and "M Phase" (533 genes, eight shared with the module).

## Sample-level scoring and adjustment

Within a dataset, expression was log-transformed, each gene standardised across samples, and the module score taken as the mean standardised value across module genes. Donor-age cohorts were modelled as score ~ age + sex, or score ~ age + proliferation + sequencing depth where a proliferation score was available. The longitudinal resource was modelled as score ~ days + donor, with MKI67 added where stated. Contact inhibition was compared with untreated cells from the same donors at the same oxygen tension by Wilcoxon rank-sum test. The proliferation score in the 142-donor cohort is the mean standardised expression of 20 canonical proliferation markers.

## Splicing outcome

Junction counts for the 142-donor cohort were taken from recount3²⁰. Splicing outcome is the fraction of junction-spanning reads assigned to junctions absent from the reference annotation; two alternatives, the fraction of distinct unannotated junctions and the mean Shannon entropy of junction usage within 5′ splice-site groups with at least 20 reads, are reported alongside. Sequencing depth was included as a covariate in every outcome model.

## Contrast-level scoring

For a contrast, the module effect is the mean log₂ fold change of module genes minus the mean over all measured genes. The cell-cycle adjustment is the residual from a least-squares regression of the module effect on the cell-cycle effect across all 65 contrasts. Residuals were tested against zero within class by Wilcoxon signed-rank test with Benjamini–Hochberg correction.

## Reporting

Spearman correlations are used for monotone associations and Pearson where a linear fit is shown; *P* values are two-sided. Every *P* value quantifies an association computed over genes or over samples as stated. None is a test of a treatment effect in the local unreplicated dataset, which contributes one contrast and appears only in Fig. 6.

## Data and code availability

Accessions are listed in Supplementary Table 1. The module gene list, the 27-gene proliferation-independent subset, all contrast vectors, class assignments and analysis code will be deposited on acceptance.

# References

1. Donega, S. *et al.* Loss of splicing homeostasis as a hallmark of aging. *Mol. Cell. Biol.* **46**, 673–691 (2026).
2. Harries, L. W. *et al.* Human aging is characterized by focused changes in gene expression and deregulation of alternative splicing. *Aging Cell* **10**, 868–878 (2011).
3. Deschênes, M. & Chabot, B. The emerging role of alternative splicing in senescence and aging. *Aging Cell* **16**, 918–933 (2017).
4. Casella, G. *et al.* Transcriptome signature of cellular senescence. *Nucleic Acids Res.* **47**, 7294–7305 (2019).
5. Whitfield, M. L. *et al.* Identification of genes periodically expressed in the human cell cycle and their expression in tumors. *Mol. Biol. Cell* **13**, 1977–2000 (2002).
6. Venet, D., Dumont, J. E. & Detours, V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput. Biol.* **7**, e1002240 (2011).
7. Ein-Dor, L., Kela, I., Getz, G., Givol, D. & Domany, E. Outcome signature genes in breast cancer: is there a unique set? *Bioinformatics* **21**, 171–178 (2005).
8. Zhang, Z. *et al.* Deciphering the role of immune cell composition in epigenetic age acceleration: insights from cell-type deconvolution applied to human blood epigenetic clocks. *Aging Cell* **23**, e14071 (2024).
9. Coller, H. A., Sang, L. & Roberts, J. M. A new description of cellular quiescence. *PLoS Biol.* **4**, e83 (2006).
10. Fujimaki, K. *et al.* Graded regulation of cellular quiescence depth between proliferation and senescence by a lysosomal dimmer switch. *Proc. Natl Acad. Sci. USA* **116**, 22624–22634 (2019).
11. Gorgoulis, V. *et al.* Cellular senescence: defining a path forward. *Cell* **179**, 813–827 (2019).
12. Shim, V. *et al.* Secretomes from the transient iPSC reprogramming phase reverse cellular photoaging and extend *Caenorhabditis elegans* longevity. *J. Tissue Eng.* (submitted; manuscript JTE-Aug-26-0256).
13. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* **11**, R25 (2010).
14. Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* **26**, 139–140 (2010).
15. Ragueneau, E. *et al.* The Reactome Knowledgebase 2026. *Nucleic Acids Res.* **54**, D673–D681 (2026).
16. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).
17. Fleischer, J. G. *et al.* Predicting age from the transcriptome of human dermal fibroblasts. *Genome Biol.* **19**, 221 (2018).
18. Sturm, G. *et al.* OxPhos defects cause hypermetabolism and reduce lifespan in cells and in patients with mitochondrial diseases. *Commun. Biol.* **6**, 22 (2023).
19. Zhang, S., Tyshkovskiy, A., Ying, K., Wang, S. & Gladyshev, V. N. Mammalian aging involves genome-wide splicing degeneration leading to functional decline. *bioRxiv* https://doi.org/10.64898/2026.06.26.734787 (2026). Preprint.
20. Wilks, C. *et al.* recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol.* **22**, 323 (2021).

# Supplementary Table 1. Datasets contributing to the compendium

Twenty-nine GEO series and one local experiment contributed the 65 contrasts scored in Fig. 5 and Fig. 6. Class is assigned by experimental manipulation, not by transcriptional content.

| Accession | Class | *n* | Contrast(s) |
|---|---|---:|---|
| GSE113957 | donor age | 1 | donor age, per decade |
| GSE179848 | donor age | 1 | late vs early passage |
| GSE226189 | donor age | 1 | donor age, per decade |
| GSE307377 | donor age | 1 | old vs young dermal fibroblast |
| GSE116968 | metabolic / culture | 2 | uv_vs_control 1h; uv_vs_control 4h |
| GSE179848 | metabolic / culture | 11 | 2-Deoxyglucose vs Control (21% O2); betahydroxybutyrate vs Control (21% O2); Contact_Inhibition vs Control (21% O2); DEX vs Contro |
| GSE116968 | photoprotection / rescue | 6 | prered_vs_uv 1h; prered_vs_uv 4h; Pre-Red vs UV, injury-residualised, 1h; Pre-Red vs UV, injury-residualised, 4h; red_vs_control 1 |
| GSE222414 | photoprotection / rescue | 2 | UVB + TAT-UIFSP5 vs UVB (rescue); UVB + TAT-UIFSP6 vs UVB (rescue) |
| GSE240486 | photoprotection / rescue | 4 | 0.06% Osmoter alone, no UV; 0.02% Osmoter alone, no UV; UV + 0.06% Osmoter vs UV (rescue); UV + 0.02% Osmoter vs UV (rescue) |
| GSE326951 | photoprotection / rescue | 1 | UVB + sauchinone vs UVB (rescue) |
| this study | Repro-CM | 1 | reprogramming-phase secretome |
| GSE149694 | reprogramming | 6 | 5iLAF-D13 vs Fibroblast-D3; Fibroblast-D7 vs Fibroblast-D3; NHSM-D13 vs Fibroblast-D3; Primed-D13 vs Fibroblast-D3; RSeT-D13 vs Fi |
| GSE297233 | reprogramming | 2 | O4YRSK +dox vs -dox, day 4; OSK +dox vs -dox, day 4 |
| GSE139563 | secretome | 2 | senescent-cell CM, day 4 vs day 0; senescent-cell CM, day 10 vs day 0 |
| GSE251807 | secretome | 4 | BM-MSC conditioned medium vs DMEM; BM-MSC EV-depleted medium vs DMEM; BM-MSC small EVs vs DMEM; BM-MSC non-small EVs vs DMEM |
| GSE266052 | secretome | 2 | cytokine-primed hMSC-CM vs untreated (both TGF-b; resting hMSC-CM vs untreated (both TGF-beta) |
| GSE268248 | secretome | 1 | PRF serum vs untreated (gingival, boundary) |
| GSE279804 | secretome | 1 | ADSC nanovesicles vs control |
| GSE282054 | secretome | 1 | hTSC secretome vs ES-CM |
| GSE293186 | secretome | 1 | endothelial-cell EVs vs control, 72 h |
| GSE306748 | secretome | 2 | primary MSC secretome vs control; immortalised MSC secretome vs control |
| GSE109700 | senescence | 2 | deep replicative senescence; early replicative senescence |
| GSE191055 | senescence | 1 | P27 vs P4 |
| GSE306957 | senescence | 4 | DDX58_KO senescent vs proliferating; IFIH1_KO senescent vs proliferating; pLenti senescent vs proliferating; siScramble senescent  |
| GSE93535 | senescence | 1 | SIPS vs quiescent |
| GSE134533 | UV injury | 1 | HDF UVB vs no UV (2 donors) |
| GSE222414 | UV injury | 1 | UVB vs control (injury) |
| GSE240486 | UV injury | 1 | UV vs control (injury) |
| GSE329475 | UV injury | 1 | UVA, 3 consecutive days vs NC (injury) |
