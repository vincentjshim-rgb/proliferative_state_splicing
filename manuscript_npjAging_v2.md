---
title: "Decoupling senescence from acute stress: a proliferation-independent transcriptional axis across 60 human fibroblast perturbations"
subtitle: "Draft v2 · target journal npj Aging (Article)"
date: "12 September 2026"
---

# Abstract

Transcriptional signatures of cellular senescence are dominated by proliferative arrest, and interventions that restore growth are therefore readily mistaken for interventions that oppose senescence. Across 1,169 pathways in 104 fibroblast donors, we find that whether an age-associated pathway retains its association after adjustment for proliferation is determined by its collinearity with proliferation (77% of pathways at |r| < 0.5 versus 5–9% at |r| > 0.8) rather than by its biology. To address this, we assembled two reference axes from public human fibroblast transcriptomes by random-effects meta-analysis and defined a decoupling index, *D*, as the difference between the loading of a signature on each (pooled ρ +0.284, 95% CI 0.227 to 0.340 and −0.109, −0.155 to −0.062; difference in Fisher-*z* 0.402, 95% CI 0.324 to 0.479). *D* is robust to the confound: residualisation on a proliferation loading estimated in 345 independent samples leaves it unchanged (0.461 to 0.468), and it is greater among proliferation-neutral genes (0.546). Scoring 60 perturbations on this axis positions nine senescence contrasts at the negative pole (*D* = −0.50 to −0.73), whereas contact inhibition, which produces arrest without senescence, is positive (*D* = +0.168). Fourteen preregistered held-out tests across four correction families yielded eight successes and no refutations. Two module families account for the axis, proteostatic stress signalling and nuclear pre-mRNA processing, and the same modules predict which perturbations resemble a given signature; perturbation magnitude does not. Applied to a published transcriptome of UVA-damaged fibroblasts exposed to a reprogramming-phase secretome, the axis ranks that response first among 39 perturbations independent of both references. This intervention is unreplicated, and these data do not reverse a chronological ageing signature.

# Introduction

Cellular senescence, one of the hallmarks of ageing^1^, is defined by a stable proliferative arrest together with a characteristic secretory and chromatin state^2,3^. In transcriptomic practice, the arrest dominates: cell-cycle, ribosome-biogenesis and general biosynthetic genes covary and account for much of the variance separating senescent from proliferating cells^4^. This creates a recurrent inferential hazard, because any intervention that restores proliferative capacity will anti-correlate with a senescence signature irrespective of whether it engages senescence biology. Claims that a compound, a secretome or a partial reprogramming protocol opposes senescence are consequently difficult to distinguish from claims that it promotes growth.

The magnitude of this confound is substantial. Among 1,169 Reactome pathways whose expression tracks donor age in cultured dermal fibroblasts, only 19% retain that association after adjustment for a proliferation module score, and retention is governed almost entirely by collinearity with proliferation rather than by pathway identity. Splicing-factor expression, whose decline with age is widely reported^5,6^, correlates with proliferation at *r* ≈ 0.89 in this setting and cannot be separated from it. We regard this as a limitation of cultured-fibroblast expression data rather than as evidence against an independent splicing axis of ageing, and note that in a system where proliferative capacity is itself a component of the ageing phenotype, proliferation may act as a mediator rather than a confounder.

A second, less frequently examined hazard concerns reference signatures. Meta-signatures derived from a small number of studies are routinely reused as though they generalised^7^, yet their transfer to independent data is seldom assessed. We show below that the two axes used here behave very differently in this respect, and that the difference materially changes what each can be said to measure.

We therefore set out to construct a quantitative axis separating the senescence programme from acute stress programmes, to determine whether that separation is independent of proliferative state, and to establish the boundaries of what the axis captures. We then applied the axis to a specific question in skin ageing: the effect of a secretome collected during the transient reprogramming phase^8–12^ on the transcriptome of UVA-damaged human dermal fibroblasts. That transcriptome, together with the accompanying phenotypic and *Caenorhabditis elegans* lifespan data, is reported separately^13^; the present study is a reanalysis and does not re-derive those findings.

Two design commitments constrain what follows. First, every generalisation test was preregistered: contrasts, directional predictions, success and refutation criteria, exclusion rules and multiplicity plans were specified and hash-locked before the corresponding data were retrieved, and were not revised thereafter. Four such families were executed and each was corrected within itself; no closed family was re-corrected when a later family was added. Second, negative and unsupported results are reported in full, because they define the boundary of the claim.

# Results

## The two reference axes differ markedly in their transferability

We constructed a senescence axis by rank-averaging four independent contrasts — longitudinal replicative senescence, deep replicative senescence, late-versus-early passage, and oxidative stress-induced premature senescence — and a UVA axis from two UVA-versus-control studies (Methods). Scoring every contrast in our perturbation compendium on the axis to which it belongs reveals that the two axes are not equivalent constructs (Fig. 2D).

Across the six reference contrasts, the anchor is positively associated with the UVA studies and negatively associated with the senescence studies. A random-effects meta-analysis of the Fisher-transformed correlations gives a pooled ρ of +0.284 (95% CI 0.227 to 0.340) for the UVA axis and −0.109 (95% CI −0.155 to −0.062) for the senescence axis (Fig. 2A). Between-study heterogeneity is substantial in both (*I*² = 94% and 93%; *Q* = 16.1, df = 1 and *Q* = 44.8, df = 3, both *P* < 0.0001), which is expected at these gene counts and is the reason a random-effects rather than a fixed-effect model is used throughout. The difference between the two pooled estimates is 0.402 in Fisher-*z* units (95% CI 0.324 to 0.479; *Z* = 10.2, *P* = 2 × 10⁻²⁴).

The same comparison by gene-set enrichment gives a directionally identical answer with a different null. The 236 genes with an anchor score of at least log~2~1.25 are enriched at the top of both UVA rankings (NES = +2.77 and +3.13; BH-adjusted *P* = 1 × 10⁻²⁷ and 3 × 10⁻¹⁸) and at the bottom of two of the three usable senescence rankings (NES = −1.47, *P* = 0.010 and NES = −1.83, *P* = 0.008), while the 287 anchor-down genes are enriched at the top of the deep-senescence ranking (NES = +2.23, *P* = 5 × 10⁻¹⁰) (Fig. 2B, C). Late-passage senescence (GSE191055) is null in both directions and is reported as such.

The senescence axis transfers with minimal attenuation. The four contrasts that define it load at a median ρ of +0.727 (range 0.619–0.761), whereas four fully independent senescence contrasts, derived from a different tissue (MRC5 lung fibroblasts) and a different inducer (20 Gy X-irradiation), load at a median of +0.656 (0.526–0.675).

By contrast, the UVA axis does not transfer. Its two defining contrasts load at a median of +0.807, yet six independent UV-injury contrasts load at a median of **−0.017** (−0.104 to +0.276). The contrasts that do load are photoprotective interventions (median +0.281, range 0.229–0.305; Wilcoxon rank-sum versus injury, *P* = 0.0087). The only injury contrast that loads is the earliest timepoint available, UVA at 6 h (+0.276); the same study at 24 h gives +0.017, repeated exposure −0.030, chronic exposure −0.005, and UVB at 1 h and 4 h gives −0.104 and −0.089, respectively (Fig. 6D).

We therefore interpret this axis as an early photo-adaptive response axis rather than a UV-damage axis, and state throughout that accumulated or chronic UV injury does not load on it. This interpretation was derived post hoc and was subsequently exposed to refutation in a preregistered test (see below).

![](<public_data_tierA/derived/figures_journal/Fig1.png>){width=6.3in}

**Fig. 1 | Study design, data quality and construction of the effect-size anchor.** **A**, Design. Adult primary human dermal fibroblasts were irradiated with UVA and then exposed for 24 h to conditioned medium from unirradiated fibroblasts (HDF-CM), from the transient reprogramming phase (Repro-CM) or from iPSC (iPSC-CM); one library per condition. **B**, Principal-component analysis of the 2,000 most variable genes. **C**, Spearman correlation between libraries, average-linkage dendrogram. **D**, MA plot of the Repro-CM versus HDF-CM comparison. Dashed lines mark |log~2~FC| = log~2~1.25. **E**, The two comparator fold changes plotted against each other; shaded quadrants are direction-consistent. Counts are genes per quadrant. No *P* value is assigned to any gene because there is a single library per condition.


## A conservative effect-size anchor derived from an unreplicated intervention

The local dataset comprises four paired-end libraries: UVA-irradiated human dermal fibroblasts at 0 h, and the same cells 24 h after exposure to conditioned medium from parental fibroblasts (HDF-CM), from cells in the transient reprogramming phase (Repro-CM), or from established iPSCs (iPSC-CM). All four libraries were realigned here under identical parameters (GENCODE v41, STAR 2.7.3a). Input ranged from 14.9 to 16.8 M read pairs, unique mapping from 93.6% to 94.7%, and per-base mismatch rate from 0.24% to 0.26%. The four libraries separate along PC1, which carries 87% of the variance, and their pairwise Spearman correlations range from 0.877 to 0.989 (Fig. 1B, C).

**There is one biological replicate per condition**, and no differential-expression statistic is therefore assigned to any gene. We instead defined a conservative, comparator-consistent effect-size anchor: for each expressed gene, the signed minimum of its log~2~ fold change against HDF-CM and against iPSC-CM, set to zero wherever the two comparators disagree in direction (Fig. 1D, E). Of 12,319 expressed genes, 8,921 (72.4%) agree in direction and 2,204 reach |score| ≥ log~2~(1.25).

The anchor proved stable to analytical choice. Across 120 recomputations spanning comparator selection, leave-one-study-out, removal of nuisance gene programmes, proliferation adjustment, and 100 binomial redraws at 50% read depth, *D* ranged from 0.241 to 0.546 and did not cross zero (Fig. 3F). The lowest value arose when iPSC-CM served as the sole comparator.

Importantly, the decoupling is attributable to the Repro-CM condition rather than to the construction itself. Applying the identical signed-minimum rule around each condition in turn gives *D* = **+0.461** for Repro-CM, **−0.368** for HDF-CM and **−0.090** for iPSC-CM, with the same ordering after proliferation adjustment (Fig. 3G).

![](<public_data_tierA/derived/figures_journal/Fig2.png>){width=6.3in}

**Fig. 2 | The two reference axes are not equivalent constructs.** **A**, Random-effects meta-analysis (DerSimonian–Laird on Fisher-transformed correlations) of the association between the anchor and each reference contrast. Squares are study estimates with area proportional to inverse-variance weight, horizontal lines are 95% confidence intervals, diamonds are pooled estimates, and *Q* and *I*² quantify between-study heterogeneity. **B**, Gene-set enrichment of the anchor-derived gene sets against each rank-averaged reference axis. Tick rows mark set members. **C**, Normalised enrichment scores per reference study. Asterisks denote BH-adjusted GSEA *P*: **\****, < 10⁻⁴; ***, < 10⁻³; **, < 0.01; *, < 0.05. GSE93535 is omitted because 39% of its fold changes are saturated at the cap of an FPKM reconstruction, so a rank-based enrichment score is undefined for it. **D**, Loading of every compendium contrast on the axis to which it belongs, grouped by relation to that axis. Boxes show median and interquartile range; each point is one contrast.


## The decoupling index is independent of proliferative state

Across the 2,237 comparator-consistent genes shared with both axes, the anchor correlates at ρ = +0.289 (95% CI 0.250 to 0.326) with the UVA axis and ρ = −0.172 (95% CI −0.211 to −0.131) with the senescence axis (Fig. 3A, B). The two axes are themselves nearly orthogonal (ρ = +0.043), so the difference between these dependent correlations is interpretable: *D* = **0.460** (95% CI 0.401 to 0.519), and Williams' test for two dependent correlations sharing one variable gives *t*(2,234) = 16.7, *P* < 10⁻¹⁶ (Fig. 3C). Because genes within a signature are correlated, this *P* value is anti-conservative as a statement about independent units; an expression-decile-matched permutation of the anchor, which preserves that dependence, places *D* far outside its null (null 95% interval −0.040 to +0.071, *P* < 1 × 10⁻⁴). Leave-one-study-out yields 0.401–0.478, and the senescence relationship persists after regressing out the UVA axis (ρ = −0.195).

To determine whether *D* reflects proliferative state, we estimated a per-gene proliferation loading — the correlation of each gene with a cell-cycle meta-gene — across 345 independent fibroblast samples spanning nine treatments, two oxygen tensions and a mitochondrial mutation. The anchor correlates only weakly with this loading (ρ = +0.108). The partial correlations holding it constant are +0.315 (95% CI 0.277 to 0.352) with the UVA axis and −0.144 (95% CI −0.185 to −0.103) with the senescence axis, a difference of 0.459 — indistinguishable from the unadjusted estimate (Fig. 3D). Residualising all three axes on the loading likewise left *D* unchanged (0.461 to **0.468**); restricting the analysis to proliferation-neutral genes (|loading| < 0.2, *n* = 715) increased it to **0.546**, whereas proliferation-linked genes (|loading| ≥ 0.4, *n* = 861) gave 0.400. When genes were binned into loading deciles, *D* remained positive in all ten bins (0.226–0.641), with every bootstrap interval above zero (Fig. 3E).

![](<public_data_tierA/derived/figures_journal/Fig3.png>){width=6.3in}

**Fig. 3 | The decoupling index, its independence from proliferative state, and its specificity.** **A**, **B**, Gene-level association between the anchor and each reference axis over the 2,237 comparator-consistent genes measured in all studies. Shading is a two-dimensional count density; the line is a least-squares fit with its 95% confidence band. **C**, The two correlations and their difference, *D*, with 95% confidence intervals; the difference is tested with Williams' test for two dependent correlations sharing one variable. **D**, Partial Spearman correlation after removing a per-gene proliferation loading estimated in 345 independent fibroblast samples. Error bars are 95% confidence intervals. **E**, *D* within deciles of that loading; error bars are 500-fold bootstrap intervals. **F**, *D* under 120 recomputations spanning comparator choice, leave-one-study-out, nuisance-programme removal, proliferation handling and 100 binomial redraws at half depth. **G**, The identical signed-minimum construction applied around each of the three conditions in turn, before and after proliferation adjustment; only the Repro-CM centring is positive.


## Two module families account for the decoupling

At the pathway level, the composite test that the decoupling hypothesis predicts — the local enrichment score against the difference between the two axes — gives ρ = **+0.476** across 553 shared Reactome pathways (*P* = 1.3 × 10^−32^). The two components behave very differently when examined separately: against the senescence axis alone, ρ = −0.528 (*P* = 4.3 × 10^−41^); against the UVA axis alone, ρ = +0.079 (*P* = 0.063, not significant). We present the non-significant component in the same figure (Fig. 4A).

Decomposition by module accounts for this asymmetry (Fig. 4B). Most modules that are elevated in the local anchor are negative on the UVA axis — translation (−4.1), nonsense-mediated decay (−4.1), rRNA processing (−3.0) and cell cycle (−7.8) — and cannot therefore contribute to decoupling. Only two families are positive on both the local and UVA axes while negative on senescence: nuclear pre-mRNA processing (splicing, 3′-end processing and export; local +2.1 to +3.5, senescence −4.6 to −7.4, UVA +0.5 to +1.4, not significant), and the HSF1/UPR-PERK proteostatic stress response^14^ (local +2.7, UVA +2.0 to +3.2, senescence −0.1 to −1.0). The two halves of the decoupling are thus carried by distinct modules. Of note, reactive oxygen species detoxification, the module a photoprotection account would predict, is the most UVA-aligned of all (+3.8) but is also elevated in senescence (+2.2) and does not qualify as a decoupling module.

![](<public_data_tierA/derived/figures_journal/Fig4.png>){width=6.3in}

**Fig. 4 | Two module families account for the decoupling.** **A**, Each of 553 Reactome pathways shared by the local data and both axes, plotted as its local enrichment score against the difference between the two axes. Line and band are a least-squares fit with its 95% confidence interval; the five most extreme pathways are labelled. **B**, Expression-matched permutation *z* for each module on each signature, with BH-adjusted significance. Only nuclear pre-mRNA processing and the HSF1/UPR-PERK proteostatic stress response are positive on both the local and UVA axes while negative on senescence. Donor age is shown to demonstrate that it does not behave as the UVA axis does. **C**, Leading Reactome pathways in the local anchor by normalised enrichment score; point area is −log~10~ BH-adjusted *P*.


## Senescence alters splicing outcome, and part of this is separable from arrest

Because nuclear pre-mRNA processing emerged as one of the two carriers, we next asked whether senescence alters splicing *outcome* rather than only the expression of the splicing machinery. Using uniformly reprocessed junction counts, and after downsampling every library to 19.5 M junction reads, deep replicative senescence differed from proliferating fibroblasts at 651 alternative junctions (FDR < 0.05; 168 with |ΔPSI| ≥ 0.10) and early senescence at 772 (231 with |ΔPSI| ≥ 0.10). Two senescent states differed from one another at zero junctions. In a label-permutation null applied to the depth-matched early-versus-proliferating comparison, the true split yielded 772 whereas all nine alternative three-versus-three splits yielded zero (Fig. 5A).

Distinguishing senescence from arrest requires a quiescence control at matched population doubling. The only available dataset of this design provides quiescent and stress-induced senescent fibroblasts at population doubling 15. Using the full 2 × 2 state-by-compound design to obtain a six-versus-six main effect, 21–25 junctions reached FDR < 0.05, and across 298 permutations of the state label within compound strata the observed count exceeded every permutation (empirical *P* = 0.0033, the permutation floor; Fig. 5b). In proportional terms the senescence-specific component is approximately five-fold smaller than the total: 0.14% of testable junctions versus 0.67% against proliferating cells (Fig. 5C).

The local intervention cannot be evaluated at this level. Among the three 24 h conditions, Repro-CM shows the highest intronic read fraction (5.56% versus 5.19% and 4.93%), the opposite of improved splicing efficiency, and per-gene intron retention tracks the expression anchor (ρ = −0.25) and is therefore not independent of expression. Junction-level ΔPSI candidates from these libraries carry an empirical false-discovery rate of 0.71–0.74 against a binomial null (Fig. 5D). Addressing this question for a secretome will require *n* ≥ 3 per condition.

![](<public_data_tierA/derived/figures_journal/Fig5.png>){width=6.3in}

**Fig. 5 | Senescence alters splicing outcome, and part of the change is separable from arrest.** **A**, Alternative junctions in GSE109700 after downsampling every library to 19.5 M junction reads; logit-PSI tested with moderated *t*-statistics and Benjamini–Hochberg correction. Two senescent states differ from one another at zero junctions. **B**, Label-permutation null for the senescence main effect at matched population doubling in GSE93535 (2 × 2 state-by-compound design, six versus six). **C**, The same counts as a proportion of testable junctions. **D**, Intronic read fraction in the four local libraries. Descriptive only: *n* = 1 per condition.


## Preregistered held-out generalisation

Fourteen tests were preregistered across four families (hashes in Methods), each corrected within itself. Eight succeeded and none refuted the model.

Family 1 (six tests) returned a partial result. A repair-aligned component replicated in an independent photoprotection study at both timepoints (ρ = +0.122 and +0.098; BH-FDR 3 × 10^−4^). Reversal of a chronological ageing signature did not replicate (ρ = −0.032; BH-FDR 0.197) and was abolished by proliferation adjustment (+0.008); we accordingly exclude chronological-age reversal from the claim. Two reprogramming boundary tests, both of which predicted a null result, held (Fig. 6A).

The repair-aligned component has now been observed in six interventions across three independent studies, with mechanisms as distinct as a herbal extract, a metabolite, photobiomodulation^15^ and a peptide-fusion protein (ρ = +0.072 to +0.127; Fig. 6B). In one of these studies the anchor is negatively correlated with the UVB injury axis itself (ρ = −0.094 and −0.169), which excludes residual injury as an explanation for the positive rescue association.

Family 2 (two tests) returned a confirmed result. In an independent senescence dataset — lung rather than dermal fibroblasts, and X-irradiation rather than replicative or oxidative senescence — the anchor opposes senescence at ρ = −0.127 (BH-FDR 1 × 10^−4^), and the co-primary proliferation-residualised test, specified in advance as the more demanding of the two, also succeeded at ρ = −0.101. The effect size falls within the range of the four discovery references (−0.044 to −0.163) without attenuation, and is unchanged across four arms including RIG-I and MDA5 knockouts (Fig. 6C).

Family 3 (five tests) returned a partial result. Two of four photoprotection arms succeeded (ρ = +0.107 and +0.072) but both originated from a single study, whereas the locked criterion required one from each; we therefore report this test as not supported and did not relax the criterion. Notably, the failing study's own UV injury contrast was itself inert (ρ = −0.030, *P* = 0.59). The fifth test exposed the post-hoc reading of the UVA axis to refutation: a three-day UVA exposure in dermal fibroblasts was predicted to load at |ρ| < 0.15, with ρ > 0.40 specified as refuting. It returned **−0.017**, and two further independent injury contrasts obtained in the same retrieval pointed in the same direction, bringing the independent UV-injury set to eight, none of which loads on the axis (Fig. 6D).

Cell-intrinsic partial reprogramming was not recapitulated. Both primary contrasts fell within the prespecified |ρ| = 0.15 boundary, and a mutant control behaved identically; a single naive-state sensitivity arm exceeded the boundary and is reported rather than omitted (Fig. 6A).

Family 4 (one test) returned a success with important qualifications. A third independent photoprotection study gave ρ = +0.373, three-fold larger than any other. That study lacks an untreated arm, so the contrast could not be residualised on injury, a limitation specified in advance. The compound is not a mild protectant: its median |log~2~FC| is 0.363, two- to four-fold greater than other agents and approximately half the magnitude of a UVB injury, with 50% of genes differentially expressed. Its signature does not resemble a reversal of UVB injury (ρ = +0.009 against an independent UVB injury signature). We report it as an unexplained outlier.

![](<public_data_tierA/derived/figures_journal/Fig6.png>){width=6.3in}

**Fig. 6 | Preregistered held-out generalisation.** **A**, All thirteen preregistered primary tests with 95% confidence intervals, grouped by correction family; shading marks the prespecified null zone. Each family was Benjamini–Hochberg-corrected within itself and no closed family was re-corrected. **B**, The repair-aligned component in six interventions across three independent studies, after removing each study's own injury component. **C**, Held-out senescence (family 2) in four independent arms, before and after proliferation adjustment; shading is the range of the four discovery references. **D**, Loading of every UV contrast on the UVA reference axis; the preregistered test HO7 is shown in amber.


## A compendium of 60 perturbations

Scoring 60 fibroblast perturbations on the same two axes shows the axis behaving as intended (Fig. 7A). *D* differs across perturbation classes (Kruskal–Wallis χ²(5) = 29.8, *P* = 1.6 × 10⁻⁵), and senescence contrasts are separated from every other class (Wilcoxon rank-sum *P* = 1.5 × 10⁻⁵). Nine senescence contrasts occupy the negative pole (*D* = −0.50 to −0.73), including four fully independent arms. Donor-age signatures follow at −0.18 to −0.32. All eleven photoprotection contrasts from studies that did not contribute to the axis have *D* ≥ 0 (+0.003 to +0.428); the two negative rescue points derive from the study that defines the UVA axis, where a negative value is structural.

The most informative internal control is contact inhibition^16^. Growth arrest without senescence sits at ρ~senescence~ = **+0.018** and *D* = **+0.168**, against a senescence mean of +0.67 and −0.59, respectively. Arrest alone therefore does not displace a signature along the senescence axis. This constitutes a third, independent route to the same conclusion reached by per-gene residualisation and by the preregistered co-primary test (Fig. 7B).

![](<public_data_tierA/derived/figures_journal/Fig7.png>){width=6.3in}

**Fig. 7 | Sixty fibroblast perturbations placed on the same two axes.** **A**, *D* by perturbation class. Open symbols contributed to a reference axis and are positive controls rather than independent evidence; the omnibus test excludes the singleton Repro-CM group. **B**, The same perturbations in the plane of the two axes; *D* is the vertical distance below the dashed identity line. Contact inhibition, which produces arrest without senescence, lies at ρ~senescence~ = +0.018. **C**, The twelve highest and eight lowest of the 39 perturbations independent of both axes, with 95% confidence intervals.


## Alignment is determined by module content, not by perturbation magnitude

Finally, we asked what determines whether a perturbation resembles a given signature. Perturbation magnitude does not: across 59 signatures, the correlation between median |log~2~FC| and alignment with the anchor is **−0.181** (*P* = 0.20), and within photoprotection alone −0.204. The twelve largest perturbations on the map align at −0.157 to +0.181 (Fig. 8A).

What does predict alignment is module content, and specifically the two families identified in Fig. 4B: HSF1 heat-shock response (ρ = +0.660, *P* = 4 × 10^−8^) and capped pre-mRNA processing (+0.604, *P* = 7 × 10^−7^), followed by splicing (+0.560) and UPR/PERK (+0.435) (Fig. 8B). Controlling for magnitude leaves HSF1 at +0.577 and UPR/PERK at +0.364, and within the fifteen agent-alone arms — perturbations applied in the absence of any injury — UPR/PERK continues to predict alignment (ρ = +0.532, *P* = 0.044), with bioenergetic stressors at the top (oligomycin +0.308, 2-deoxyglucose +0.127, SURF1 deficiency +0.094) and quiescence at the bottom (−0.109). We note that part of this relationship is expected, since a perturbation sharing a signature's dominant modules will correlate with it; what is established here is that the modules are identified and that magnitude, extracellular matrix content and cell-cycle content are excluded as explanations.

![](<public_data_tierA/derived/figures_journal/Fig8.png>){width=6.3in}

**Fig. 8 | Alignment is set by module content, not by perturbation magnitude.** **A**, Each compendium signature plotted as the size of its own transcriptional response against its correlation with the anchor; the anchor is excluded. Line and band are a least-squares fit with its 95% confidence interval. **B**, Correlation between module content and alignment across the same signatures, with 95% confidence intervals and BH-adjusted significance. The two families identified in Fig. 4B are the strongest predictors; perturbation magnitude is not.


## Position of the reprogramming-phase secretome

On this map, the response of UVA-damaged fibroblasts to Repro-CM occupies *D* = **+0.460**, third of 60 and first among the 39 perturbations independent of both reference axes. The two contrasts above it are the UVA studies that define the positive axis and are therefore circular. Its ρ~UVA~ of +0.289 places it within the photoprotection group (+0.229 to +0.305) rather than the injury group (median −0.017).

We state the limits of this placement explicitly. The intervention is unreplicated; the primary transcriptome has been published^13^; no causal or mediation experiment was performed; and the axis does not measure chronological ageing, which these data explicitly fail to reverse.

# Discussion

Three conclusions are supported by the analyses reported here. First, a transcriptional axis separating the senescence programme from acute adaptive stress can be constructed from public fibroblast data and does not simply restate proliferative state: it survives per-gene residualisation on an externally estimated proliferation loading, is stronger among proliferation-neutral genes, passes a preregistered co-primary test in held-out data, and places growth arrest without senescence on the opposite side of the axis from senescence itself. Second, the axis is carried by two module families rather than one — proteostatic stress signalling and nuclear pre-mRNA processing — and these same modules determine which of 60 independent perturbations resembles a given signature, whereas perturbation magnitude does not. Third, the response of UVA-damaged dermal fibroblasts to a reprogramming-phase secretome occupies the highest position on this axis among perturbations independent of the references.

The results that were removed by these tests are equally informative. The axis constructed from two UVA studies does not generalise to independent UV injury: eight independent injury contrasts, spanning 1 h to chronic exposure and both wavebands, load at or below zero. Describing it as a UVA-damage axis would therefore have been incorrect, and the alternative explanation that this invited — that a positive loading reflects residual damage — is precluded once the axis is shown not to carry damage. Reversal of a chronological ageing signature^17,18^ failed in held-out data and was abolished by proliferation adjustment. Transcriptional mimicry of cell-intrinsic partial reprogramming was absent across three independent studies, with one naive-state arm reported as an exception. Longevity-intervention signatures provided no support and disagreed in sign between two references derived from the same source^19,20^. Splicing outcome is established for the senescence reference but cannot be measured for the secretome at *n* = 1.

The proliferation confound warrants emphasis beyond the present dataset. In cultured dermal fibroblasts from 104 donors, four-fifths of age-associated pathways lose that association once proliferation is adjusted for, and retention is determined by collinearity with proliferation rather than by pathway identity. Any claim that an intervention reverses a fibroblast ageing signature should therefore state whether it survives this adjustment. Conversely, in a system where proliferative capacity is itself part of the ageing phenotype, such adjustment may remove the causal path, and the direction of that argument should be declared rather than assumed.

The principal limitation is intrinsic to the available data. The intervention that motivates this work is represented by a single library per condition. Every downstream result is therefore a property of one fold-change vector, and no amount of public validation alters that fact. What the validation does establish is that this vector behaves consistently against many independent references, and that the identical construction applied to the two comparator conditions does not reproduce it. A replicated secretome experiment — ideally ribodepleted and paired-end at *n* ≥ 3, which would additionally render splicing outcome measurable — is the natural next step. A second limitation is that the positive pole of the axis rests on two defining studies; although the photoprotection alignment now has three independent replications, the axis itself would benefit from additional independent UVA data.

One observation remains unexplained. A single photoprotective compound aligns at ρ = +0.373, three-fold above every other tested intervention, and neither its perturbation magnitude nor its module content accounts for this. We report it rather than omit it, and note that resolving it will require a study of that compound that includes an untreated arm.

# Methods

## Local dataset and provenance

Four paired-end RNA-seq libraries of UVA-irradiated (15 J cm^−2^) primary human dermal fibroblasts were analysed: 0 h after irradiation, and 24 h after exposure to HDF-CM, Repro-CM or iPSC-CM. These libraries, together with the accompanying phenotypic and *C. elegans* lifespan experiments, are reported separately^13^; the present study is a reanalysis and re-derives none of those results. Reads were adapter- and quality-trimmed (Trim Galore 0.6.4, Cutadapt 2.8) and aligned with STAR 2.7.3a^37^ to GRCh38 with GENCODE v41 annotation^38^ (`--outFilterType BySJout`, sjdbOverhang 100), identically for all four samples. Reverse-stranded gene counts were used throughout.

## Conservative comparator-consistent anchor

Counts were TMM-normalised^39,40^ and genes with CPM ≥ 1 in at least two samples were retained (12,319 genes); log~2~(CPM + 0.5) was used. For each gene, *d*~H~ = log~2~CPM(Repro) − log~2~CPM(HDF) and *d*~I~ = log~2~CPM(Repro) − log~2~CPM(iPSC). The anchor score is sign(*d*~H~) × min(|*d*~H~|, |*d*~I~|) where the two agree in direction and are non-zero, and 0 otherwise. No *P* value is assigned to any gene, as there is a single biological replicate per condition.

## Reference axes

The senescence axis is the rank-normalised average of four independent contrasts (longitudinal replicative senescence; deep replicative senescence; late-versus-early passage; oxidative stress-induced premature senescence). The UVA axis is the rank-normalised average of two UVA-versus-control contrasts. Public processed data were used throughout, and contrast definitions and sample selections were verified against series metadata manually. Every series used in this study, together with its source publication^10,12,17,18,21–36^, is listed in Supplementary Table 1.

## Decoupling index and null model

*D* = ρ~Spearman~(anchor, UVA axis) − ρ~Spearman~(anchor, senescence axis), computed over genes present in all three. The null shuffles the anchor score within deciles of mean expression (10,000 draws), so that significance is not driven by expression-level structure. Genes are not biological replicates; all reported *P* and FDR values denote the significance of a signature association and never of a treatment effect.

## Proliferation loading

A cell-cycle meta-gene (Reactome^41^ "Cell Cycle, Mitotic" and "DNA Replication", 511 genes) was scored per sample across 345 independent fibroblast samples, and each gene's Spearman correlation with that meta-gene was taken as its proliferation loading. Residualisation regresses each axis on this loading before *D* is recomputed.

## Splicing outcome

Exon–exon junction counts from recount3^42^, a uniformly reprocessed public resource, were used. PSI was computed within splice-site-anchored groups containing at least two alternatives and with adequate coverage in all samples of a contrast. Logit-PSI was tested using moderated *t*-statistics (limma^43^) with Benjamini–Hochberg correction^44^. Depth matching used binomial downsampling to the smallest library. Label-permutation nulls permuted the group label, within strata where a second factor was present. Local intron retention was quantified over unambiguous intronic and exonic blocks (GENCODE v41), from which blocks overlapping any other gene's exons or any overlapping gene span were removed, using split-aware read counting^45^.

## Perturbation compendium

Sixty fibroblast contrasts were assembled from data already held, spanning senescence, donor age, UV injury, photoprotection, reprogramming media, metabolic and culture perturbations, and the local secretome response. Signatures contributing to a reference axis are flagged and presented as positive controls rather than as independent evidence.

## Preregistration

Four preregistration documents were written and SHA-256-locked before the corresponding data were retrieved, and were not revised thereafter. Clarifications made before results were examined are appended to each document's change log.

| Family | Tests | SHA-256 (first 8) | Outcome |
|---|---:|---|---|
| 1 — photoprotection, age, reprogramming boundary | 6 | `41e4cf4b` | Partial |
| 2 — held-out senescence | 2 | `6e5a7703` | Confirmed |
| 3 — photoprotection replication, UVA axis | 5 | `315e68ea` | Partial |
| 4 — third photoprotection study | 1 | `c542d32f` | Success (limited) |

Each family was Benjamini–Hochberg-corrected^44^ within itself. No closed family was re-corrected when a later family was added; this is stated explicitly because pooling all fourteen tests would retrospectively alter FDR values that had already been fixed.

## Statistics and reporting

Signature association was measured as a Spearman correlation and reported with a 95% confidence interval obtained by the Fisher *z* transformation, where the interval width reflects the number of shared genes and not the number of biological replicates. Correlations were combined across studies by random-effects meta-analysis (DerSimonian–Laird on Fisher-transformed estimates), with Cochran's *Q* and *I*² reported for heterogeneity and the contrast between axes taken as the difference of the two pooled *z* values with its standard error. The difference between the two correlations that defines *D* was tested with Williams' test for two dependent correlations sharing one variable^46^; because genes within a signature are not independent, this *P* value is anti-conservative, and an expression-decile-matched permutation of the anchor score (10,000 draws) is reported alongside it as a dependence-preserving calibration. Gene-set enrichment used the adaptive multilevel procedure implemented in fgsea^47,48^, with normalised enrichment scores and Benjamini–Hochberg-adjusted *P* values; references whose ranking statistic contained more than 20% tied values were excluded from enrichment analysis and this is stated where it applies. Proliferation adjustment is reported both as a partial Spearman correlation and as residualisation of each axis on the per-gene loading. Group comparisons used Kruskal–Wallis and Wilcoxon rank-sum tests, and differential splicing used moderated *t*-statistics on logit-PSI with Benjamini–Hochberg correction^43,44^. Asterisks in figures denote **\****, *P* < 10⁻⁴; ***, *P* < 10⁻³; **, *P* < 0.01; *, *P* < 0.05; n.s., not significant. Exclusion rules — gene-symbol mapping rate among expressed genes ≥ 70%, at least 1,000 shared genes, and a constructible contrast — were fixed in advance.

All *P* and FDR values in this study quantify the significance of an association between two signatures, computed over genes. None of them is a test of a treatment effect between biological replicates, and the local intervention, at one library per condition, supports no such test.

## Data and code availability

All public accessions are listed in Supplementary Table 1. Derived signatures, the scored compendium, permutation nulls and per-figure source tables are provided as Supplementary Data. Analysis scripts will be deposited in a public repository on acceptance.

# References

1. López-Otín, C., Blasco, M. A., Partridge, L., Serrano, M. & Kroemer, G. Hallmarks of aging: an expanding universe. *Cell* **186**, 243–278 (2023).
2. Gorgoulis, V. *et al.* Cellular senescence: defining a path forward. *Cell* **179**, 813–827 (2019).
3. Coppé, J.-P., Desprez, P.-Y., Krtolica, A. & Campisi, J. The senescence-associated secretory phenotype: the dark side of tumor suppression. *Annu. Rev. Pathol.* **5**, 99–118 (2010).
4. Hernandez-Segura, A. *et al.* Unmasking transcriptional heterogeneity in senescent cells. *Curr. Biol.* **27**, 2652–2660.e4 (2017).
5. Harries, L. W. *et al.* Human aging is characterized by focused changes in gene expression and deregulation of alternative splicing. *Aging Cell* **10**, 868–878 (2011).
6. Deschênes, M. & Chabot, B. The emerging role of alternative splicing in senescence and aging. *Aging Cell* **16**, 918–933 (2017).
7. Saul, D. *et al.* A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nat. Commun.* **13**, 4827 (2022).
8. Ocampo, A. *et al.* In vivo amelioration of age-associated hallmarks by partial reprogramming. *Cell* **167**, 1719–1733.e12 (2016).
9. Lu, Y. *et al.* Reprogramming to recover youthful epigenetic information and restore vision. *Nature* **588**, 124–129 (2020).
10. Gill, D. *et al.* Multi-omic rejuvenation of human cells by maturation phase transient reprogramming. *eLife* **11**, e71624 (2022).
11. Browder, K. C. *et al.* In vivo partial reprogramming alters age-associated molecular changes during physiological aging in mice. *Nat. Aging* **2**, 243–253 (2022).
12. Lu, J. Y. *et al.* Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming. *Cell* **188**, 5895–5911.e17 (2025).
13. Shim, V. *et al.* Secretomes from the transient iPSC reprogramming phase reverse cellular photoaging and extend *Caenorhabditis elegans* longevity. *J. Tissue Eng.* (submitted; manuscript JTE-Aug-26-0256).
14. Labbadia, J. & Morimoto, R. I. The biology of proteostasis in aging and disease. *Annu. Rev. Biochem.* **84**, 435–464 (2015).
15. Hamblin, M. R. Mechanisms and mitochondrial redox signaling in photobiomodulation. *Photochem. Photobiol.* **94**, 199–212 (2018).
16. Coller, H. A., Sang, L. & Roberts, J. M. A new description of cellular quiescence. *PLoS Biol.* **4**, e83 (2006).
17. Fleischer, J. G. *et al.* Predicting age from the transcriptome of human dermal fibroblasts. *Genome Biol.* **19**, 221 (2018).
18. Tsitsipatis, D. *et al.* Transcriptomes of human primary skin fibroblasts of healthy individuals reveal age-associated mRNAs and long noncoding RNAs. *Aging Cell* **22**, e13915 (2023).
19. Tyshkovskiy, A. *et al.* Identification and application of gene expression signatures associated with lifespan extension. *Cell Metab.* **30**, 573–593.e8 (2019).
20. Tyshkovskiy, A. *et al.* Universal transcriptomic hallmarks of mammalian ageing and mortality. *Nature* **654**, 173–188 (2026).
21. De Cecco, M. *et al.* L1 drives IFN in senescent cells and promotes age-associated inflammation. *Nature* **566**, 73–78 (2019).
22. Kim, H. S. *et al.* Transcriptomic analysis of human dermal fibroblast cells reveals potential mechanisms underlying the protective effects of visible red light. *J. Dermatol. Sci.* **94**, 276–283 (2019).
23. Montoni, A. *et al.* Chronic UVA1 irradiation of human dermal fibroblasts: persistence of DNA damage and validation of a cell culture-based model of photoaging. *J. Invest. Dermatol.* **139**, 1821–1824.e3 (2019).
24. Liu, X. *et al.* Reprogramming roadmap reveals route to human induced trophoblast stem cells. *Nature* **586**, 101–107 (2020).
25. Sturm, G. *et al.* OxPhos defects cause hypermetabolism and reduce lifespan in cells and in patients with mitochondrial diseases. *Commun. Biol.* **6**, 22 (2023).
26. Hasegawa, T. *et al.* Cytotoxic CD4^+^ T cells eliminate senescent cells by targeting cytomegalovirus antigen. *Cell* **186**, 1417–1431.e20 (2023).
27. Yang, T. *et al.* UV radiation-induced peptides in frog skin confer protection against cutaneous photodamage through suppressing MAPK signaling. *MedComm* **5**, e625 (2024).
28. Li, Y. *et al.* Effects of rice fermentation and its bioactive components on UVA-induced oxidative stress and senescence in dermal fibroblasts. *Photochem. Photobiol.* **101**, 392–403 (2025).
29. Yan, X. *et al.* Anti-aging and rejuvenating effects and mechanism of Dead Sea water in skin. *Int. J. Cosmet. Sci.* **46**, 307–317 (2024).
30. Fontana, G. A. *et al.* UVA irradiation promotes ROS-mediated formation of the common deletion in mitochondrial DNA. *Life* **16**, 577 (2026).
31. Victorelli, S. *et al.* Mitochondrial RNA cytosolic leakage drives the SASP. *Nat. Commun.* **16**, 10992 (2025).
32. Tanaka, H. *et al.* Nucleosome stability safeguards cell identity, stress resilience and healthy aging. *bioRxiv* https://doi.org/10.1101/2025.09.17.676776 (2025). Preprint.
33. Zhang, X. *et al.* Sauchinone attenuates UVB-induced photoaging by suppressing oxidative stress and ferroptosis through activation of the Keap1/Nrf2 pathway. *Pharm. Biol.* **64**, 764–782 (2026).
34. Zhong, J. *et al.* ST6GAL1 is a functional regulator of UVA-induced photoaging in human dermal fibroblasts. *Cells* **15**, 1497 (2026).
35. Yo, K. & Rünger, T. M. UVA and UVB induce different sets of long noncoding RNAs. *J. Invest. Dermatol.* **137**, 769–772 (2017).
36. Lämmermann, I. *et al.* Blocking negative effects of senescence in human skin fibroblasts with a plant extract. *npj Aging Mech. Dis.* **4**, 4 (2018).
37. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013).
38. Frankish, A. *et al.* GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Res.* **51**, D942–D949 (2023).
39. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* **11**, R25 (2010).
40. Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* **26**, 139–140 (2010).
41. Ragueneau, E. *et al.* The Reactome Knowledgebase 2026. *Nucleic Acids Res.* **54**, D673–D681 (2026).
42. Wilks, C. *et al.* recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol.* **22**, 323 (2021).
43. Ritchie, M. E. *et al.* limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* **43**, e47 (2015).
44. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).
45. Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* **26**, 841–842 (2010).
46. Steiger, J. H. Tests for comparing elements of a correlation matrix. *Psychol. Bull.* **87**, 245–251 (1980).
47. Subramanian, A. *et al.* Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad. Sci. USA* **102**, 15545–15550 (2005).
48. Korotkevich, G., Sukhov, V. & Sergushichev, A. Fast gene set enrichment analysis. *bioRxiv* https://doi.org/10.1101/060012 (2019). Preprint.

# Supplementary Table 1 | Public datasets

Twenty GEO series and three non-GEO resources were used. "Defines axis" marks the four contrasts that constitute a reference axis; those contrasts are positive controls and are excluded wherever independence is claimed. Family numbers refer to the four preregistered test families (Methods).

| Accession | Biological system | Assay | Design | Role | Contrasts | Ref. |
|---|---|---|---|---|---:|---:|
| GSE240226 | Human dermal fibroblast | RNA-seq | 4 arms, *n* = 3 | UVA axis (defines); rescue discovery | 3 | 28 |
| GSE302943 | hTERT BJ5ta fibroblast | RNA-seq (NextSeq 500) | 2 arms, *n* = 5 | UVA axis (defines) | 1 | 30 |
| GSE89005 | Human dermal fibroblast | Array (GPL16956) | Sham *n* = 3, UVA *n* = 2 per arm | Independent UV injury | 3 | 35 |
| GSE125429 | Human dermal fibroblast, 4 donor strains | Array (GPL13607) | 4 complete pairs | Independent UV injury, chronic | 1 | 23 |
| GSE109700 | LF1 lung fibroblast | RNA-seq | 3 states, *n* = 3 | Senescence axis (defines); splicing analysis | 2 | 21 |
| GSE179848 | HC1–HC4 dermal fibroblast, longitudinal | RNA-seq | 345 samples, 4 donors | Senescence axis (defines); proliferation loading; 12 perturbations | 12 | 25 |
| GSE191055 | Human dermal fibroblast | RNA-seq (NovaSeq 6000) | P27 vs P4, *n* = 4 | Senescence axis (defines) | 1 | 26 |
| GSE93535 | Human dermal fibroblast | RNA-seq | 4 arms, *n* = 3 | Senescence axis; quiescence contrast | 2 | 36 |
| GSE113957 | Primary dermal fibroblast, 97 donors | RNA-seq (FPKM) | Age 22–89 y | Chronological age | 1 | 17 |
| GSE226189 | Primary skin fibroblast, 82 donors | RNA-seq | Age 22–89 y | Chronological age | 1 | 18 |
| GSE165177 | Fibroblast, transient MPTR reprogramming | RNA-seq | Day 13, donor-paired | Reprogramming boundary | 1 | 10 |
| GSE297233 | Fibroblast from a 96-year-old donor, OSK | RNA-seq | *n* = 4 | Held out, family 1 (HO1) | 2 | 12 |
| GSE307377 | Primary dermal fibroblast, old vs young | RNA-seq | *n* = 9 | Held out, family 1 (HO2) | 1 | 32 |
| GSE116968 | NHDF, red-light photoprotection | RNA-seq | *n* = 6 | Held out, family 1 (HO3a, HO3b) | 8 | 22 |
| GSE149694 | Fibroblast reprogramming time course | RNA-seq | *n* = 4–5 | Held out, family 1 (HO4a, HO4b) | 6 | 24 |
| GSE306957 | MRC5 lung fibroblast, X-irradiation senescence, RIG-I and MDA5 knockouts | RNA-seq | *n* = 5–7 per arm | Held out, family 2 (HO5a, HO5b) | 4 | 31 |
| GSE240486 | NHDF, osmolyte photoprotection | RNA-seq | 6 arms, *n* = 3 | Held out, family 3 (HO6a, HO6b) | 5 | 29 |
| GSE222414 | WS1 fetal skin fibroblast, peptide-fusion photoprotection | RNA-seq | 4 arms, *n* = 4 | Held out, family 3 (HO6c, HO6d) | 3 | 27 |
| GSE329475 | Human dermal fibroblast, three-day UVA | RNA-seq | 2 arms, *n* = 3 | Held out, family 3 (HO7) | 1 | 34 |
| GSE326951 | HFF foreskin fibroblast, sauchinone photoprotection | RNA-seq | 2 arms, *n* = 3 | Held out, family 4 (HO8) | 1 | 33 |
| This study | Adult human dermal fibroblast, UVA, conditioned medium | RNA-seq | *n* = 1 per condition | Local anchor | 1 | 13 |

| Non-GEO resource | Version | Use | Ref. |
|---|---|---|---:|
| recount3 | Release of 2026-08 | Uniformly processed junction counts for GSE109700 and GSE93535 | 42 |
| Reactome | v97 | Pathway gene sets, module definitions | 41 |
| Tyshkovskiy *et al.* | 2019 and 2026 supplementary tables | External lifespan and age signatures | 19, 20 |
