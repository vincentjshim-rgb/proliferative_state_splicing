# Decoupling senescence from acute stress: a proliferation-independent transcriptional axis across 60 human fibroblast perturbations

**Draft v1 · 2026-09-12 · target: npj Aging (Article)**

---

## Abstract

Transcriptional signatures of cellular senescence are dominated by proliferative arrest, so an intervention that merely restores growth can be mistaken for one that opposes senescence. Across 1,169 pathways in 104 fibroblast donors, whether an age-associated pathway survives adjustment for proliferation is set by its collinearity with proliferation (77% below |r| = 0.5; 5–9% above 0.8), not by its biology. We therefore built two reference axes from public human fibroblast data and defined a decoupling index, D, as the difference between a signature's loading on each. D resists the confound: residualising on a proliferation loading estimated in 345 independent samples leaves it unchanged (0.461 to 0.468), and it is larger among proliferation-neutral genes (0.546). Scoring 60 perturbations places nine senescence contrasts at the negative pole (−0.50 to −0.73) while contact inhibition — arrest without senescence — sits positive (+0.168). Fourteen preregistered held-out tests in four correction families gave eight successes and no refutations. Two module families carry the axis, proteostatic stress signalling and nuclear pre-mRNA processing, and the same two predict which perturbations resemble a signature, whereas perturbation magnitude does not. A published transcriptome of UVA-damaged fibroblasts given a reprogramming-phase secretome ranks first of 39 perturbations independent of both axes; that intervention is unreplicated, and these data do not reverse chronological ageing.

## Introduction

Cellular senescence is defined by a stable proliferative arrest together with a distinct secretory and chromatin state. In transcriptomic practice the arrest dominates: genes of the cell cycle, ribosome biogenesis and general biosynthesis move together and account for much of the variance that separates senescent from proliferating cells. The consequence is a recurring inferential hazard. Any intervention that restores proliferative capacity will anti-correlate with a senescence signature, whether or not it engages senescence biology. Reports that a compound, a secretome or a partial reprogramming protocol "opposes senescence" are therefore difficult to separate from reports that it makes cells grow.

The problem is not hypothetical. Here we show that among 1,169 Reactome pathways whose expression tracks donor age in cultured dermal fibroblasts, only 19% retain that association after adjustment for a proliferation module score, and that survival is set almost entirely by how collinear the pathway is with proliferation rather than by its biology. Splicing-factor expression, which is widely reported to decline with age, sits at r ≈ 0.89 with proliferation in this setting and cannot be separated from it — a limitation of the data rather than a refutation of the literature, but one that is rarely stated.

A second hazard concerns reference signatures themselves. Meta-signatures assembled from a small number of studies are routinely reused as if they generalised, and their transfer to independent data is seldom tested. We show that the two axes used here behave very differently in this respect, and that the difference changes what each can be claimed to measure.

We set out to construct a quantitative axis that separates the senescence programme from acute stress programmes, to test whether that separation is independent of proliferative state, and to establish what the axis does and does not capture. We then apply it to a specific question in skin ageing: what a secretome collected during the transient reprogramming phase does to the transcriptome of ultraviolet-A (UVA)-damaged human dermal fibroblasts. That transcriptome, together with the phenotypic and lifespan results that accompany it, has been reported separately¹; the present work is a reanalysis and does not re-derive those findings.

Two design commitments constrain what follows. First, every generalisation test was preregistered: contrasts, directional predictions, success and refutation criteria, exclusion rules and multiplicity plans were written and hash-locked before the corresponding data were downloaded, and were not revised afterwards. Four such families were run, and each was corrected within itself; no closed family was re-corrected when a later one was added. Second, negative and unsupported results are reported in full, because the boundary of the claim is set by them.

---

## Results

### The two reference axes transfer to independent data very differently

We assembled a senescence axis by rank-averaging four independent contrasts — longitudinal replicative senescence, deep replicative senescence, late-versus-early passage, and oxidative stress-induced premature senescence — and a UVA axis from two UVA-versus-control studies (Methods). Scoring every contrast in our perturbation compendium on the axis to which it belongs shows that the two axes are not equivalent objects (**Fig. 1a**).

The senescence axis transfers almost without loss. The four contrasts that define it load at a median ρ of +0.727 (range 0.619–0.761); four fully independent senescence contrasts, from a different tissue (MRC5 lung fibroblasts) and a different inducer (20 Gy X-irradiation), load at a median of +0.656 (0.526–0.675).

The UVA axis does not. Its two defining contrasts load at a median of +0.807, but six independent UV-injury contrasts load at a median of **−0.017** (−0.104 to +0.276). What does load on it are photoprotective interventions (median +0.281, 0.229–0.305; Wilcoxon versus injury p = 0.0087). The single injury contrast that loads is the earliest timepoint available, UVA at 6 h (+0.276); the same study at 24 h gives +0.017, repeated exposure −0.030, chronic exposure −0.005, and UVB at 1 h and 4 h −0.104 and −0.089 (**Fig. 1b**).

We therefore treat this axis as an **early photo-adaptive response axis** rather than a UV-damage axis, and state throughout that accumulated or chronic UV injury does not load on it. This reading was generated post hoc and was subsequently put at risk of refutation in a preregistered test (below).

### A conservative effect-size anchor from an unreplicated intervention

The local dataset comprises four paired-end libraries: UVA-irradiated human dermal fibroblasts at 0 h, and the same cells 24 h after exposure to conditioned medium from parental fibroblasts (HDF-CM), from cells in the transient reprogramming phase (Repro-CM), or from established iPSCs (iPSC-CM). All four were realigned here under identical parameters (GENCODE v41, STAR 2.7.3a); input was 14.9–16.8 M read pairs, unique mapping 93.6–94.7%, and per-base mismatch 0.24–0.26% (**Fig. 2a**).

**There is one biological replicate per condition.** No differential-expression statistic is assigned to any gene. Instead we define a conservative, comparator-consistent effect-size anchor: for each expressed gene, the signed minimum of its log₂ fold change against HDF-CM and against iPSC-CM, set to zero wherever the two disagree in direction (**Fig. 2b**). Of 12,319 expressed genes, 8,921 (72.4%) agree in direction and 2,204 reach |score| ≥ log₂(1.25).

The anchor is stable to analytical choice. Across 124 recomputations — comparator choice, leave-one-study-out, removal of nuisance gene programmes, proliferation adjustment, and 100 binomial redraws at 50% read depth — D ranges from 0.241 to 0.546 and never crosses zero (**Fig. 2c**). The weakest value arises when iPSC-CM is the only comparator.

Critically, the decoupling is a property of the Repro-CM condition rather than of the construction. Applying the identical signed-minimum rule around each condition in turn gives D = **+0.461** for Repro-CM, **−0.368** for HDF-CM and **−0.090** for iPSC-CM, with the same ordering after proliferation adjustment (**Fig. 2d**).

### The decoupling index and its independence from proliferative state

Over the 2,237 comparator-consistent genes shared with both axes, the anchor correlates at ρ = +0.289 with the UVA axis and ρ = −0.172 with the senescence axis, giving D = **0.460** (expression-matched permutation, 10,000 draws, p < 1 × 10⁻⁴; null 95% interval −0.040 to +0.071). Leave-one-study-out gives 0.401–0.478, and the senescence relationship survives regression of the UVA axis out of it (ρ = −0.195).

To test whether D is a proliferation artefact we estimated a per-gene proliferation loading — the correlation of each gene with a cell-cycle meta-gene — across 345 independent fibroblast samples spanning nine treatments, two oxygen tensions and a mitochondrial mutation. The anchor itself correlates only weakly with this loading (ρ = +0.108). Residualising all three axes on it leaves D essentially unchanged (0.461 → **0.468**); restricting to proliferation-neutral genes (|loading| < 0.2, n = 715) *increases* it to **0.546**, while proliferation-linked genes (|loading| ≥ 0.4, n = 861) give 0.400 (**Fig. 3a**). Binning genes into loading deciles, D is positive in all ten bins (0.226–0.641), with every bootstrap interval above zero (**Fig. 3b**).

For context on why this test matters: repeating the age-association analysis across 1,169 Reactome pathways in 104 healthy donors, 650 pathways track donor age but only 125 (19.2%) survive proliferation adjustment, and survival is set by collinearity — 77% for |r| < 0.5, 13% for 0.7–0.8, and 5–9% above 0.8. Splicing, pre-mRNA processing and mRNA export modules sit at r = 0.89–0.94 and lose their age association entirely. We report this as a limit on what cultured-fibroblast expression data can resolve, not as evidence against an independent splicing axis of ageing; in this regime proliferation may be a mediator rather than a confounder, and adjusting for it would then be over-adjustment.

### Two module families carry the decoupling

At pathway level, the composite test that the decoupling hypothesis actually predicts — the local enrichment score against the difference between the two axes — gives ρ = **+0.476** across 553 shared Reactome pathways (p = 1.3 × 10⁻³²). The two components taken separately behave very differently: against the senescence axis alone, ρ = −0.528 (p = 4.3 × 10⁻⁴¹); against the UVA axis alone, ρ = +0.079 (p = 0.063, **not significant**). We report the non-significant component in the same figure (**Fig. 4a**).

Decomposing by module explains the asymmetry (**Fig. 4b**). Most modules that are up in the local anchor are *negative* on the UVA axis — translation (−4.1), NMD (−4.1), rRNA processing (−3.0), cell cycle (−7.8) — and so cannot contribute to decoupling. Only two families are positive on both the local and UVA axes while negative on senescence:

- **nuclear pre-mRNA processing** (splicing, 3′-end processing, export): local +2.1 to +3.5, senescence −4.6 to −7.4, UVA +0.5 to +1.4 (not significant);
- **HSF1 / UPR-PERK proteostatic stress response**: local +2.7, UVA +2.0 to +3.2, senescence −0.1 to −1.0.

The two halves of the decoupling are therefore carried by different modules. ROS detoxification, the module a photoprotection narrative would reach for, is the most UVA-aligned of all (+3.8) but is also *up* in senescence (+2.2) and is not a decoupling module.

### Senescence alters splicing outcome, and part of that is separable from arrest

Because pre-mRNA processing emerged as one carrier, we asked whether senescence changes splicing *outcome* rather than only the expression of splicing machinery. Using uniformly reprocessed junction counts, and after downsampling every library to 19.5 M junction reads, deep replicative senescence differs from proliferating fibroblasts at 651 alternative junctions (FDR < 0.05; 168 with |ΔPSI| ≥ 0.10) and early senescence at 772 (231). Two senescent states differ from each other at **zero** junctions, and in a label-permutation null on the depth-matched early-versus-proliferating comparison the true split gives 772 while all nine alternative 3-versus-3 splits give zero (**Fig. 5a**).

Whether this is senescence or merely arrest requires a quiescence control at matched population doubling. The only such dataset available provides quiescent and stress-induced senescent fibroblasts at population doubling 15; using the full 2 × 2 state-by-compound design to obtain a six-versus-six main effect, 21–25 junctions reach FDR < 0.05, and in 298 permutations of the state label within compound strata the observed count exceeds every permutation (empirical p = 0.0033, the permutation floor) (**Fig. 5b**). Proportionally the senescence-specific component is about five-fold smaller than the total: 0.14% of testable junctions versus 0.67% against proliferating cells (**Fig. 5c**).

The local intervention cannot be evaluated at this level. Among the three 24 h conditions, Repro-CM has the *highest* intronic read fraction (5.56% versus 5.19% and 4.93%) — the opposite of improved splicing efficiency — and per-gene intron retention tracks the expression anchor (ρ = −0.25), so it is not independent of expression. Junction-level ΔPSI candidates from these libraries carry an empirical false-discovery rate of 0.71–0.74 against a binomial null (**Fig. 5d**). Answering this question for a secretome requires n ≥ 3 per condition.

### Preregistered held-out generalisation

Fourteen tests were preregistered across four families (hashes in Methods), each corrected within itself. Eight succeeded and **none refuted the model** (**Fig. 6**).

**Family 1** (six tests) returned PARTIAL. A repair-aligned component replicated in an independent photoprotection study at both timepoints (ρ = +0.122 and +0.098, BH-FDR 3 × 10⁻⁴). Reversal of a chronological-age signature did **not** replicate (ρ = −0.032, BH-FDR 0.197) and disappears under proliferation adjustment (+0.008); we exclude chronological-age reversal from the claim. Two reprogramming boundary tests, which predicted null, held.

**Family 2** (two tests) returned CONFIRMED. In an independent senescence dataset — lung rather than dermal fibroblasts, X-irradiation rather than replicative or oxidative senescence — the anchor opposes senescence at ρ = −0.127 (BH-FDR 1 × 10⁻⁴), and the co-primary proliferation-residualised version, specified in advance as the harder test, also succeeds at ρ = −0.101. The effect size sits inside the range of the four discovery references (−0.044 to −0.163) with no shrinkage, and is unchanged across four arms including RIG-I and MDA5 knockouts.

**Family 3** (five tests) returned PARTIAL. Two of four photoprotection arms succeeded (ρ = +0.107 and +0.072) but both came from one study, and the locked criterion required one from each; we report this as not supported and did not relax the criterion. The failing study's own UV injury contrast was itself inert (ρ = −0.030, p = 0.59). The fifth test put the post-hoc reading of the UVA axis at risk: a three-day UVA exposure in dermal fibroblasts was predicted to load at |ρ| < 0.15, with ρ > 0.40 set as refuting. It returned **−0.017**, and two further independent injury contrasts from the same download pointed the same way, bringing the independent UV-injury set to eight, none of which loads on the axis (**Fig. 6e**).

**Family 4** (one test) returned SUCCESS with important caveats. A third independent photoprotection study gave ρ = +0.373, three-fold larger than any other. That study has no untreated arm, so the contrast could not be residualised on injury — a limitation fixed in advance. The compound is not a mild protectant: its median |log₂FC| is 0.363, two- to four-fold larger than other agents and about half the size of a UVB injury, with 50% of genes differentially expressed. Its signature does not resemble a reversal of UVB injury (ρ = +0.009 against an independent UVB injury signature). We record it as an unexplained outlier.

### A compendium of 60 perturbations, and what alignment requires

Scoring 60 fibroblast perturbations assembled from data already in hand on the same two axes (**Fig. 7**) shows the axis behaving as intended. Nine senescence contrasts occupy the negative pole (D = −0.50 to −0.73), including four fully independent arms. Donor-age signatures follow at −0.18 to −0.32. All eleven photoprotection contrasts from studies that did not build the axis have D ≥ 0 (+0.003 to +0.428); the two negative rescue points come from the study that defines the UVA axis, where a negative value is structural.

The internal control that matters most is contact inhibition: growth arrest without senescence sits at ρ_senescence = **+0.018** and D = **+0.168**, against a senescence mean of +0.67 and −0.59. Arrest alone does not move the senescence axis. This is a third, independent route to the same conclusion as the per-gene residualisation and the preregistered co-primary test.

Finally we asked what makes a perturbation resemble a given signature. Perturbation magnitude does not: across 59 signatures, ρ between median |log₂FC| and alignment with the anchor is **−0.181** (p = 0.20), and within photoprotection alone −0.204; the twelve largest perturbations on the map align at −0.157 to +0.181 (**Fig. 8a**). What does predict alignment is module content, and specifically the two families identified in Fig. 4b: HSF1 heat-shock response (ρ = +0.660, p = 4 × 10⁻⁸) and capped pre-mRNA processing (+0.604, p = 7 × 10⁻⁷), followed by splicing (+0.560) and UPR/PERK (+0.435) (**Fig. 8b**). Controlling for magnitude leaves HSF1 at +0.577 and UPR/PERK at +0.364, and within the fifteen agent-alone arms — perturbations applied without any injury — UPR/PERK still predicts alignment (ρ = +0.532, p = 0.044), with bioenergetic stressors at the top (oligomycin +0.308, 2-deoxyglucose +0.127, SURF1 deficiency +0.094) and quiescence at the bottom (−0.109). Part of this relationship is expected, since a perturbation sharing a signature's dominant modules will correlate with it; what is new is that the modules are named and that magnitude, ECM content and cell-cycle content are excluded as the explanation.

### Position of the reprogramming-phase secretome

On this map, the response of UVA-damaged fibroblasts to Repro-CM sits at D = **+0.460**, third of 60 and **first of the 39 perturbations independent of both reference axes** (**Fig. 7a**). The two contrasts above it are the UVA studies that define the positive axis and are therefore circular. Its ρ_UVA of +0.289 places it among the photoprotection group (+0.229 to +0.305) rather than the injury group (median −0.017).

We state the limits of this placement plainly. The intervention is unreplicated; the primary transcriptome has been published¹; no causal or mediation experiment was performed; and the axis does not measure chronological ageing, which these data explicitly fail to reverse.

---

## Discussion

Three claims survive the tests reported here. First, a transcriptional axis separating the senescence programme from acute adaptive stress can be constructed from public fibroblast data and is not a restatement of proliferative state: it survives per-gene residualisation on an externally estimated proliferation loading, is stronger in proliferation-neutral genes, passes a preregistered co-primary test in held-out data, and places growth arrest without senescence on the opposite side from senescence. Second, the axis is carried by two module families rather than one — proteostatic stress signalling and nuclear pre-mRNA processing — and the same two modules determine which of 60 independent perturbations resembles a given signature, while perturbation magnitude does not. Third, the response of UVA-damaged dermal fibroblasts to a reprogramming-phase secretome occupies the highest position on this axis among perturbations independent of the references.

Equally important is what the tests removed. The axis built from two UVA studies does not generalise to independent UV injury; eight independent injury contrasts, spanning 1 h to chronic exposure and both wavebands, load at or below zero. Calling it a UVA-damage axis would have been wrong, and the alternative explanation it invited — that a positive loading reflects residual damage — dissolves once the axis is seen not to carry damage at all. Reversal of a chronological-age signature failed in held-out data and vanishes under proliferation adjustment. Transcriptional mimicry of cell-intrinsic partial reprogramming is absent across three independent studies, with one naive-state arm as a reported exception. Longevity-intervention signatures give no support and disagree in sign between two references from the same source. Splicing outcome is established for the senescence reference but cannot be measured for the secretome at n = 1.

The proliferation problem deserves emphasis beyond this dataset. In cultured dermal fibroblasts from 104 donors, four fifths of age-associated pathways lose that association once proliferation is adjusted for, and whether a pathway survives is set by its collinearity with proliferation rather than by biology. Any claim that an intervention "reverses a fibroblast ageing signature" should state whether it survives this adjustment; conversely, in a system where proliferative capacity is itself part of the ageing phenotype, such adjustment may remove the causal path, and the direction of that argument should be declared rather than assumed.

The main limitation is unavoidable with the data in hand. The intervention that motivates the work is represented by one library per condition. Every downstream result is a property of a single fold-change vector, and no amount of public validation changes that; what the validation establishes is that this vector behaves consistently against many independent references, and that the same construction applied to the two comparator conditions does not reproduce it. A replicated secretome experiment — ideally ribodepleted and paired-end, at n ≥ 3, which would also make splicing outcome measurable — is the natural next step. A second limitation is that the positive pole of the axis rests on two defining studies; the photoprotection alignment now has three independent replications, but the axis itself would benefit from additional independent UVA data.

Finally, one observation remains unexplained. A single photoprotective compound aligns at ρ = +0.373, three-fold above every other, and neither its perturbation magnitude nor its module content accounts for it. We report it rather than omit it.

---

## Methods

### Local dataset and its provenance
Four paired-end RNA-seq libraries of UVA-irradiated (15 J/cm²) primary human dermal fibroblasts: 0 h after irradiation, and 24 h after exposure to HDF-CM, Repro-CM or iPSC-CM. These libraries and the accompanying phenotypic and *C. elegans* lifespan experiments are reported separately¹; the present study is a reanalysis and re-derives none of those results. Reads were adapter- and quality-trimmed (Trim Galore) and aligned with STAR 2.7.3a to GRCh38 with GENCODE v41 annotation (`--outFilterType BySJout`, sjdbOverhang 100), identically for all four samples. Reverse-stranded gene counts were used throughout.

### Conservative comparator-consistent anchor
Counts were TMM-normalised; genes with CPM ≥ 1 in at least two samples were retained (12,319 genes); log₂(CPM + 0.5) was used. For each gene, d_H = log₂CPM(Repro) − log₂CPM(HDF) and d_I = log₂CPM(Repro) − log₂CPM(iPSC); the anchor score is sign(d_H) × min(|d_H|, |d_I|) when the two agree in direction and non-zero, and 0 otherwise. **No p-value is assigned to any gene**, as there is one biological replicate per condition.

### Reference axes
Senescence axis: rank-normalised average of four independent contrasts (longitudinal replicative senescence; deep replicative senescence; late-versus-early passage; oxidative stress-induced premature senescence). UVA axis: rank-normalised average of two UVA-versus-control contrasts. Public processed data were used throughout; contrast definitions and sample selections were verified against series metadata by hand.

### Decoupling index and null model
D = ρ_Spearman(anchor, UVA axis) − ρ_Spearman(anchor, senescence axis) over genes present in all three. The null shuffles the anchor score **within deciles of mean expression** (10,000 draws), so significance is not driven by expression-level structure. Genes are not biological replicates: all reported p and FDR values are the significance of a *signature association*, never of a treatment effect.

### Proliferation loading
A cell-cycle meta-gene (Reactome "Cell Cycle, Mitotic" + "DNA Replication", 511 genes) was scored per sample across 345 independent fibroblast samples, and each gene's Spearman correlation with that meta-gene taken as its proliferation loading. Residualisation regresses each axis on this loading before recomputing D.

### Splicing outcome
Exon–exon junction counts from a uniformly reprocessed public resource were used. PSI was computed within splice-site anchored groups with ≥ 2 alternatives and adequate coverage in all samples of a contrast; logit-PSI was tested with moderated t-statistics (limma) and BH correction. Depth matching used binomial downsampling to the smallest library. Label-permutation nulls permuted the group label (within strata where a second factor existed). Local intron retention used unambiguous intronic and exonic blocks (GENCODE v41; blocks overlapping any other gene's exons or any overlapping gene span were removed) with split-aware read counting.

### Perturbation compendium
Sixty fibroblast contrasts were assembled from data already held, spanning senescence, donor age, UV injury, photoprotection, reprogramming media, metabolic and culture perturbations, and the local secretome response. Signatures contributing to a reference axis are flagged and shown as positive controls, never as independent evidence.

### Preregistration
Four preregistration documents were written and SHA-256-locked before the corresponding data were downloaded, and were not revised afterwards; clarifications made before results were seen are appended to each document's change log.

| Family | Tests | Document SHA-256 (first 8) | Outcome |
|---|---:|---|---|
| 1 — photoprotection, age, reprogramming boundary | 6 | `41e4cf4b` | PARTIAL |
| 2 — held-out senescence | 2 | `6e5a7703` | CONFIRMED |
| 3 — photoprotection replication, UVA axis | 5 | `315e68ea` | PARTIAL |
| 4 — third photoprotection study | 1 | `c542d32f` | SUCCESS (limited) |

Each family was BH-corrected within itself. No closed family was re-corrected when a later family was added; this is stated because pooling all fourteen would retrospectively change already-fixed FDRs.

### Statistics and reporting
Spearman correlations throughout. Forest-plot bars are expression-matched permutation 95% **null** intervals, not confidence intervals, and are labelled as such in every figure. Exclusion rules (gene-symbol mapping rate among expressed genes ≥ 70%; ≥ 1,000 shared genes; constructible contrast) were fixed in advance.

### Data and code availability
All public accessions are listed in Supplementary Table 1. Derived signatures, scored compendium, permutation nulls and per-figure source tables are provided as Supplementary Data. Analysis scripts are deposited at [repository to be created].

---

## Figure legends

**Fig. 1 | The two reference axes transfer to independent data very differently.**
**a**, Every compendium contrast scored on the axis to which it belongs, grouped by whether it contributed to that axis. **b**, Independent UV-injury contrasts ordered by exposure duration and waveband.

**Fig. 2 | The local anchor: quality, construction, robustness and specificity.**
**a**, Alignment metrics for all four local libraries. **b**, Per-gene log₂ fold change against each comparator; the anchor retains only direction-consistent genes. **c**, D under 124 sensitivity recomputations. **d**, The same construction applied around each of the three conditions in turn.

**Fig. 3 | The decoupling is independent of proliferative state.**
**a**, D under four prespecified proliferation-adjustment conditions. **b**, D within deciles of per-gene proliferation loading, with 500-fold bootstrap intervals.

**Fig. 4 | Pathway architecture.**
**a**, Local enrichment score against the difference between the two axes, 553 pathways. **b**, Module-by-axis map of expression-matched permutation z-scores.

**Fig. 5 | Splicing outcome.**
**a**, Senescence versus proliferating fibroblasts, depth-matched. **b**, Senescence versus quiescence at matched population doubling, with the label-permutation null. **c**, Proportion of testable junctions reaching significance against each reference. **d**, Local intron retention; descriptive only, n = 1 per condition.

**Fig. 6 | Preregistered held-out generalisation.**
**a**, Family 1, six primary tests. **b**, The repair-aligned component across three independent studies. **c**, Reprogramming boundary tests, primary and sensitivity arms. **d**, Family 2, held-out senescence, unadjusted and proliferation-adjusted. **e**, Family 3 test HO7 and the eight independent UV-injury contrasts.

**Fig. 7 | Sixty fibroblast perturbations on one axis.**
**a**, Ranked by D; open symbols contributed to a reference axis. **b**, The same perturbations in the two-axis plane.

**Fig. 8 | What alignment requires.**
**a**, Perturbation magnitude against alignment with the anchor. **b**, Module content against alignment.

---

## Reference

1. Shim, V., Ahn, K. S., Kang, J., Kim, H. K., Yoon, D. S., Jo, S., Park, C., Lee, S. G., Kang, J. H., Heo, S. Y., Park, J., Hong, K., Kang, K. *et al.* Secretomes from the transient iPSC reprogramming phase reverse cellular photoaging and extend *Caenorhabditis elegans* longevity. *J. Tissue Eng.* (submitted; manuscript JTE-Aug-26-0256). **[complete on acceptance]**

*Remaining references to be added: senescence hallmarks; SASP; proliferation confounding in senescence transcriptomics; splicing and ageing; partial reprogramming; photobiomodulation; the source publication of each compendium dataset.*
