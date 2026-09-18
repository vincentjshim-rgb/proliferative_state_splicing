---
title: "Transcriptional reproducibility is a property of the perturbation, not of its effect size: a cross-study benchmark of 76 human fibroblast perturbations"
subtitle: "Draft v1 · target journal npj Aging (Article)"
date: "12 September 2026"
---

# Abstract

Reference gene-expression signatures are routinely derived from one or two studies and then reused as though they generalised. We asked, for a single cell type, which kinds of perturbation actually produce a transcriptional state that survives a change of laboratory. We assembled 76 perturbation contrasts in human fibroblasts from 29 GEO series and one local dataset, spanning senescence, donor age, ultraviolet injury, photoprotective interventions, reprogramming media, metabolic and culture perturbations, and cell-free secretome preparations, and scored all 2,768 evaluable pairs by Spearman correlation over shared genes. Within a single experiment the median correlation between two arms of the same class is +0.201; between laboratories it is +0.004. The size of that collapse differs 20-fold between classes. Senescence retains 46% of its within-study coherence between studies (median ρ = +0.387, 95% CI 0.294 to 0.547), whereas ultraviolet injury retains 5% (+0.009, −0.005 to 0.029) and reprogramming media retain none (−0.017, −0.036 to 0.008). This is not a consequence of effect size. Across classes the correlation between effect magnitude and reproducibility is negative (ρ = −0.286); reprogramming media produce the largest transcriptional responses of any class and the least reproducible ones. At the level of individual pairs, magnitude explains 6.6% of the variance in similarity and adding perturbation class raises this to 30.1% (F(6,282) = 15.8, P = 1.2 × 10⁻¹⁵), and the ordering is unchanged when pairs are restricted to a matched magnitude range. Secretome interventions are intermediate: strongly reproducible within an experiment (ρ = +0.33 to +0.68) and weakly so between (+0.143, 0.104 to 0.168). We give the operational consequence — only senescence and donor age support a transferable reference axis in this cell type — and apply it to a reprogramming-phase secretome, which opposes the senescence axis at ρ = −0.172 (95% CI −0.211 to −0.131), third of 15 secretome contrasts.

# Introduction

Much of molecular ageing research is organised around reference signatures. A set of genes, or a ranked fold-change vector, is derived from one experiment and then used to score other experiments: to ask whether a compound reverses photoageing, whether an intervention opposes senescence, or whether a cell has been rejuvenated. The practice is productive and it is also fragile, because it assumes that the original experiment captured something that recurs elsewhere.

That assumption is tested unevenly. For cellular senescence it has been examined repeatedly, and the answer is nuanced: consensus gene lists derived by intersecting several models are short¹, overlap at the level of individual genes is limited across inducers², and single-cell work shows pronounced cell-type specificity³. The same literature nonetheless reports that global transcriptional variation across inducers aligns along a shared axis², which is what makes senescence signatures usable in practice. For ultraviolet injury the assumption is rarely tested at all. A conserved set of UV-responsive genes has been proposed from keratinocyte data⁴ and used as a biomarker panel⁵, and a large applied literature scores candidate photoprotective agents against a UV signature obtained in the same study. Whether a UV signature obtained in one laboratory predicts a UV response in another has not, to our knowledge, been measured.

The analogous question is now becoming urgent for secretomes. Circulating and secreted proteins are the basis of a rapidly growing family of ageing biomarkers⁶⁻⁸, and cell-free conditioned medium and extracellular-vesicle preparations are being developed as interventions. Those efforts rest on the premise that a secretome produces a characteristic and repeatable state in the cells that receive it. That premise has not been checked across laboratories.

A general difficulty underlies all three cases. When a signature fails to transfer, the natural explanation is that the effect was small and the measurement noisy. Attenuation is real, and a weak perturbation will correlate poorly with anything. If reproducibility simply tracked effect size, a benchmark of reproducibility would be a benchmark of statistical power and would carry no biological information.

We therefore fixed the cell type — human fibroblasts — and varied only the kind of perturbation, assembling every fibroblast contrast we could obtain on comparable terms, and measured cross-study reproducibility on one scale. We then asked directly whether the differences between classes are explained by how large their transcriptional effects are. They are not, and the classes separate in a way that has immediate consequences for which reference signatures can be built and which cannot.

Two features of the study should be stated at the outset. First, every hypothesis test that could have been tuned after seeing data was preregistered and hash-locked; five such documents were written and none was revised. Second, an earlier version of this analysis built a "UVA reference axis" from two UV studies, and that axis is among the objects the benchmark rejects. We report that episode because it illustrates exactly the failure the benchmark is designed to detect.

# Results

## A cross-study benchmark of fibroblast perturbations

We assembled 76 perturbation contrasts profiled in human fibroblasts (Fig. 1A). Twenty-nine independent GEO series contributed 75 of them and one unreplicated local experiment contributed the last. Each contrast was reduced to a genome-wide log₂ fold-change vector, and every pair was scored by Spearman correlation over genes measured in both, with at least 1,000 shared genes required. Of 2,850 possible pairs, 2,768 (97%) were evaluable, with a median of 11,972 shared genes.

Contrasts were assigned to seven classes on the basis of what was done to the cells, not on the basis of their transcriptional content: senescence, donor age, ultraviolet injury, photoprotection or rescue, reprogramming media, metabolic and culture perturbations, and cell-free secretome preparations. The local reprogramming-phase secretome is kept separate and used only in the final section. Class sizes range from 4 to 16 contrasts drawn from 2 to 8 independent studies (Fig. 1C). The clustered similarity matrix shows one dominant block and little else (Fig. 1B).

![](<public_data_tierA/derived/figures_benchmark/Fig1.png>){width=6.3in}

**Fig. 1 | Design of the benchmark and the similarity landscape.** **A**, Analysis design. A single cell type is held fixed and only the kind of perturbation varies; reproducibility is defined as the median correlation between contrasts of the same class obtained in different laboratories. **B**, All 76 contrasts ordered by average-linkage clustering of 1 − ρ; the bar beneath shows class membership, coloured as in **C**. Grey cells are pairs with fewer than 1,000 shared genes. **C**, Number of contrasts and of independent studies per class.

## Reproducibility differs twenty-fold between perturbation classes

Pooling all classes, the median correlation between two contrasts of the same class is +0.201 when they come from the same experiment and +0.004 when they come from different laboratories. Reproducibility is therefore the exception rather than the rule, and the interesting quantity is how much of the within-experiment agreement survives the move between experiments.

That quantity differs sharply (Fig. 2A, B). Senescence is the most reproducible class: nine contrasts from four studies, spanning replicative exhaustion, late passage, oxidative stress-induced premature senescence and X-irradiation in both dermal and lung fibroblasts, agree at a median between-study ρ of **+0.387** (95% CI 0.294 to 0.547), which is 46% of their within-study agreement of +0.843. Donor age follows at +0.220 (0.122 to 0.383). Secretome preparations reach +0.143 (0.104 to 0.168).

The remaining classes do not transfer. Photoprotective interventions agree at +0.031 (0.008 to 0.087) despite agreeing at +0.566 within a study, a retention of 5%. Ultraviolet injury agrees at **+0.009** (−0.005 to 0.029) across 42 between-study pairs. Metabolic and culture perturbations reach +0.004. Reprogramming media are the extreme case: nine contrasts agree at +0.742 within a study and at **−0.017** (−0.036 to 0.008) between studies, retaining nothing at all.

The failure of the ultraviolet class is not a platform artefact. Restricting to the RNA-seq contrasts alone, and so excluding the two microarray studies, gives a median of −0.007. Inspecting the 42 individual pairs shows that exactly one is substantially positive, at +0.348, and that pair consists of the two studies from which we had previously constructed a "UVA reference axis" (Fig. 2C). The axis was built on the single coherent pair in the class and, as reported below, did not transfer to any of the remaining eight independent UV contrasts.

![](<public_data_tierA/derived/figures_benchmark/Fig2.png>){width=6.3in}

**Fig. 2 | Cross-study reproducibility differs twenty-fold between classes.** **A**, Median Spearman ρ between contrasts of the same class obtained in different studies, with 2,000-fold bootstrap 95% confidence intervals; the number of contributing pairs is given at the left. **B**, The same classes compared with their within-study agreement; "retained" is the between-study median expressed as a percentage of the within-study median. Donor-age contrasts have no within-study pairs. **C**, All individual between-study pairs for the two extreme classes. Boxes show median and interquartile range.

## Reproducibility is not explained by effect size

If the classes differed only in how strongly they perturb the transcriptome, attenuation alone would generate the ordering above. They do not (Fig. 3).

At the class level the relationship runs the wrong way: the correlation between the median effect magnitude of a class, measured as the median |log₂FC| of its contrasts, and its between-study reproducibility is **ρ = −0.286** (P = 0.53) (Fig. 3A). Reprogramming media perturb the transcriptome more than any other class, with a median |log₂FC| of 0.638, and are the least reproducible. Ultraviolet injury perturbs it more than secretomes do (0.246 versus 0.200) and is sixteen times less reproducible.

At the level of individual pairs, magnitude does carry some signal — larger perturbations do correlate slightly better (β = 0.348, t = 4.50, P = 1 × 10⁻⁵) — but it is a minor term. A model with magnitude alone explains 6.6% of the variance in pairwise similarity; adding perturbation class raises this to 30.1%, and the improvement is decisive (F(6,282) = 15.8, P = 1.2 × 10⁻¹⁵) (Fig. 3B).

The most direct control restricts the comparison to pairs whose magnitudes lie inside the range spanned by the secretome class, 0.11 to 0.32 in |log₂FC|, leaving 204 pairs matched for perturbation size. The ordering is preserved: secretome +0.143, photoprotection +0.042, reprogramming +0.011, ultraviolet injury +0.007, metabolic +0.006 (Fig. 3C).

Reproducibility is therefore a property of what was done to the cells, not of how much it moved them.

![](<public_data_tierA/derived/figures_benchmark/Fig3.png>){width=6.3in}

**Fig. 3 | Reproducibility is not a function of effect size.** **A**, Class-level relationship between median effect magnitude and between-study reproducibility. **B**, Pair-level relationship across 290 between-study pairs; the dashed line and band are a least-squares fit with its 95% confidence interval, and the inset reports the nested-model comparison. **C**, Pairs restricted to the magnitude range spanned by the secretome class, so that the classes being compared are matched for perturbation size.

## Secretome effects are reproducible within an experiment and largely not between

Fourteen secretome contrasts were assembled from eight studies in which a cell-free preparation — conditioned medium, an EV-depleted medium fraction, small or non-small extracellular vesicles, or extruded nanovesicles — was applied to recipient fibroblasts with a matched control arm. Sources included senescent fibroblasts, bone-marrow and adipose mesenchymal stromal cells, trophoblast stem cells, endothelial cells and platelet-rich fibrin.

Within a single experiment these preparations produce clearly reproducible states. Two arms of the same study agree at a median ρ of +0.382, and every study tested individually is coherent: immortalised versus primary MSC secretome +0.681, cytokine-primed versus resting MSC-conditioned medium +0.565, four fractions of the same bone-marrow MSC preparation +0.364 (range 0.327 to 0.811), and conditioned medium from senescent cells at two exposure times +0.325 (Fig. 4A, B). Between studies the same class agrees at +0.143 (Wilcoxon P = 1 × 10⁻⁵).

The pipeline therefore detects secretome effects; they simply do not survive a change of laboratory and source. We cannot separate two readings of this result with the present design, and state both: it may be that "secretome" is not a coherent perturbation category, or that secretomes from different source cells genuinely do different things. All eight studies used a different source, so the design confounds the two.

Multidimensional scaling of the full matrix shows the same structure from another direction (Fig. 4C). The eigenvalue spectrum is flat — 18.4%, 14.7%, 9.6% — so no single direction dominates the space. The leading axis nevertheless tracks the senescence direction (|r| = 0.74) and not the ultraviolet direction, consistent with senescence being the only strongly transferable class.

![](<public_data_tierA/derived/figures_benchmark/Fig4.png>){width=6.3in}

**Fig. 4 | Secretome effects are reproducible within an experiment and weakly so between experiments.** **A**, All secretome pairs, split by whether the two contrasts come from the same experiment. **B**, Agreement between arms within each contributing study; bars span the range of pairs, the point is the median, and the dashed line marks the between-study median. **C**, Classical multidimensional scaling of all 76 contrasts.

## What the benchmark permits, and a worked application

The operational reading is a rule about which reference axes can be built. Taking a lower confidence bound of ρ = 0.10 as the minimum for a usable reference, only senescence and donor age qualify in this cell type; photoprotection, ultraviolet injury, metabolic perturbation and reprogramming media do not (Fig. 5B). A UV reference signature built in one laboratory should not be expected to score a UV experiment performed in another, and the same holds for reprogramming media, whose large and internally consistent effects are laboratory-specific.

We applied the surviving axis to the question that motivated this work. The local contrast is a transcriptome of UVA-damaged human dermal fibroblasts exposed for 24 h to conditioned medium collected during the transient reprogramming phase — the window in which epigenetic age declines before somatic identity is lost⁹⁻¹¹ — scored against conditioned medium from parental fibroblasts and from established iPSC. There is one library per condition, so the contrast is an effect-size anchor and no gene is assigned a differential-expression statistic.

Placed on the senescence axis together with every secretome contrast in the benchmark, the reprogramming-phase secretome loads at **ρ = −0.172** (95% CI −0.211 to −0.131), third most negative of 15 (Fig. 5A). Ten of the fifteen are negative, so opposing the senescence axis is common among secretomes rather than distinctive; the resting MSC-conditioned medium contrast is positive at +0.203. We state the limits plainly: the local intervention is unreplicated, its primary transcriptome is published elsewhere¹², no causal experiment was performed, and the axis measures the senescence direction and not chronological age.

![](<public_data_tierA/derived/figures_benchmark/Fig5.png>){width=6.3in}

**Fig. 5 | What the benchmark permits, and a worked application.** **A**, Every secretome contrast in the benchmark, and the local reprogramming-phase secretome, scored on the senescence reference axis; error bars are Fisher-*z* 95% confidence intervals. Shading marks the half-plane opposing senescence. **B**, Classes ranked by between-study reproducibility, with a lower confidence bound of ρ = 0.10 taken as the minimum for a usable reference axis.

## Preregistration

Five preregistration documents were written and SHA-256-locked before the corresponding data were retrieved, covering seventeen tests in five correction families; none was revised after the fact. Families 1 to 4 concerned the generalisation of an earlier framework and are reported in full in Supplementary Note 1, including the test that retired the UVA axis: a three-day UVA exposure in dermal fibroblasts was predicted to load on that axis at |ρ| < 0.15, with ρ > 0.40 specified as refuting, and returned −0.017.

Family 5 concerned the secretome class. Its primary test (HS1) predicted that secretome contrasts would not form a coherent class, with |median between-study ρ| < 0.15 as the prediction and > 0.35 as refuting. The observed value was +0.144, satisfying the prediction but lying 0.006 below the threshold; we report it as a boundary value rather than as a clean confirmation. **Its preregistered positive control (HS2) failed.** Conditioned medium from oncogene-induced senescent fibroblasts, applied to recipient fibroblasts, loaded on the senescence axis at ρ = +0.026 at day 4 and −0.067 at day 10, against a prediction of ρ > +0.20. Under the rule fixed in the same document, a failed positive control suspends interpretation of HS1, and we suspend it. The within-study control reported above serves the same purpose and passes decisively, but it was defined after the data were seen and is labelled post hoc throughout.

The failure of HS2 is itself informative: paracrine senescence induced by a secretome does not resemble replicative senescence at the whole-transcriptome level, which is consistent with the inducer-specificity documented within single laboratories².

# Discussion

The central result is a dissociation. Within one experiment, every class of perturbation we examined produces a coherent transcriptional state, with median within-study correlations from +0.07 to +0.84. Between experiments, most of that coherence disappears, and how much disappears is set by the kind of perturbation rather than by its magnitude. Senescence keeps nearly half of it; reprogramming media, which move the transcriptome hardest, keep none.

This has a direct methodological consequence. The practice of deriving a reference signature from one or two studies is safe for senescence and for donor age in fibroblasts, and unsafe for ultraviolet injury, for photoprotective interventions, for reprogramming media, and for metabolic perturbations. In those classes a signature will score contrasts from the same laboratory well and contrasts from elsewhere at approximately zero. We encountered this directly: an earlier version of this work built a UVA reference axis from the only coherent pair among 42 UV comparisons, interpreted its failure to load independent UV injury as a biological finding about early photo-adaptation, and exposed that interpretation to a preregistered test which it survived only because the prediction had been set to a null. The benchmark supplies the simpler explanation.

Why senescence should be transferable and ultraviolet injury not is worth stating as a hypothesis rather than a conclusion. Senescence is a terminal, self-reinforcing state: cells converge on it and remain there, so the measured contrast is a difference between two stable attractors. Ultraviolet injury, photoprotection and reprogramming media are transients sampled at an arbitrary time after an arbitrary dose, and the studies in our compendium differ in waveband, fluence, timepoint and medium composition. On that account reproducibility would be predicted by whether the perturbation has an endpoint, not by how hard it pushes — which is what we observe, though we have not tested it directly.

For secretome-based ageing biomarkers the result is a caution with a number attached. Cell-free preparations do produce reproducible states in recipient fibroblasts, at +0.33 to +0.68 within an experiment, and those states are largely specific to the preparation and the laboratory, agreeing at +0.14 elsewhere. Any claim that a secretome induces a characteristic recipient-cell programme should therefore be demonstrated across sources and laboratories rather than within one. Our design cannot distinguish whether the category is incoherent or whether different source cells simply produce different signals, because every study used a different source; resolving that would require several laboratories applying secretomes from the same source.

The limitations are substantial. Class assignment is ours, and a different grouping would give different numbers; we assigned classes by experimental manipulation rather than by transcriptional content to avoid circularity, but the choice is not unique. Class sizes are unequal, and donor age rests on six pairs. Reproducibility is measured as a rank correlation over the whole transcriptome, which is the most forgiving of the usual metrics; gene-list overlap would give lower values everywhere and would not change the ordering. We measured one cell type, and it does not follow that ultraviolet responses are irreproducible in keratinocytes or in skin. Finally, the local intervention that motivated the work remains a single library per condition; the benchmark tells us which axis it may be placed on, not how large its effect is.

# Methods

## Signature assembly

Seventy-six contrasts were assembled in human fibroblasts from 29 GEO series and one local experiment (Supplementary Table 1). Where processed count matrices were available, counts were TMM-normalised¹³,¹⁴ and a contrast was computed as the difference in mean log₂(CPM + 0.5) between arms, retaining genes with CPM ≥ 1 in at least two samples. Where a study deposited only a differential-expression table, the published log₂ fold change was used. Gene identifiers were mapped to HGNC symbols; where a study deposited Ensembl identifiers only, symbols were obtained from the other deposited matrices in the same batch. Contrast definitions and sample-to-arm assignments were verified manually against series metadata for every study.

The local contrast comprises four paired-end libraries of UVA-irradiated (15 J cm⁻²) primary human dermal fibroblasts: 0 h after irradiation, and 24 h after exposure to conditioned medium from parental fibroblasts, from cells in the transient reprogramming phase, or from established iPSC. Reads were trimmed (Trim Galore 0.6.4, Cutadapt 2.8) and aligned with STAR 2.7.3a¹⁵ to GRCh38 with GENCODE v41 annotation¹⁶. Because there is one biological replicate per condition, no gene is assigned a differential-expression statistic; the contrast is the signed minimum of the two comparator log₂ fold changes, set to zero where they disagree in direction.

## Similarity and reproducibility

Every pair of contrasts was scored by Spearman correlation over genes measured in both, requiring at least 1,000 shared genes; pairs below that threshold are treated as missing. Reproducibility of a class is the median correlation over its pairs drawn from different GEO series; within-study agreement is the median over its pairs drawn from the same series. Confidence intervals are 2,000- to 3,000-fold bootstrap percentile intervals of the median. Correlations of individual contrasts with a reference axis are reported with Fisher-*z* 95% confidence intervals, whose width reflects the number of shared genes and not the number of biological replicates.

## Effect-size control

The magnitude of a contrast is the median |log₂FC| over its measured genes. At the class level, magnitude and reproducibility were compared by Spearman correlation. At the pair level, the similarity of a pair was regressed on the geometric mean of the two magnitudes, and this model was compared with one additionally containing perturbation class by nested *F* test. The magnitude-matched analysis restricts all classes to the |log₂FC| range spanned by the secretome class.

## Reference axes

The senescence reference axis is the rank-averaged mean of four independent senescence contrasts. An earlier UVA axis, built by the same procedure from two UVA studies, is reported here only as an example of a reference that the benchmark rejects and is not used to support any claim.

## Preregistration and statistics

Five preregistration documents were SHA-256-locked before the corresponding data were retrieved (`41e4cf4b`, `6e5a7703`, `315e68ea`, `c542d32f`, `55bf016d`). Each family was Benjamini–Hochberg-corrected within itself¹⁷ and no closed family was re-corrected when a later family was added. Group comparisons used Wilcoxon rank-sum and Kruskal–Wallis tests. Every *P* value in this study quantifies the significance of an association between two signatures, computed over genes; none is a test of a treatment effect between biological replicates.

## Data and code availability

All accessions are listed in Supplementary Table 1. The 76 × 76 similarity matrix, all contrast vectors, the class assignments and the preregistration documents are provided as Supplementary Data. Analysis code will be deposited in a public repository on acceptance.

# References

1. Casella, G. *et al.* Transcriptome signature of cellular senescence. *Nucleic Acids Res.* **47**, 7294–7305 (2019).
2. Bridge, J. E. *et al.* A transcriptomic analysis reveals shared and inducer-specific expression patterns of cellular senescence. *bioRxiv* https://doi.org/10.64898/2026.04.20.719721 (2026). Preprint.
3. Sanborn, M. A. *et al.* Unveiling the cell-type-specific landscape of cellular senescence through single-cell transcriptomics using SenePy. *Nat. Commun.* **16**, 1884 (2025).
4. Shen, Y. *et al.* Epigenetic and genetic dissections of UV-induced global gene dysregulation in skin cells through multi-omics analyses. *Sci. Rep.* **7**, 42646 (2017).
5. Saul, D. *et al.* A new gene set identifies senescent cells and predicts senescence-associated pathways across tissues. *Nat. Commun.* **13**, 4827 (2022).
6. Wang, Y. *et al.* Organ-specific proteomic aging clocks predict disease and longevity across diverse populations. *Nat. Aging* **6**, 162–180 (2026).
7. Coenen, L. *et al.* An extracellular matrix aging clock based on circulating matrisome proteins predicts biological aging and disease. *Aging Cell* **25**, e70474 (2026).
8. Xiao, H. *et al.* Proteomic aging clocks in epidemiological studies: advances, applications and prospects. *Nat. Aging* **6**, 970–986 (2026).
9. Olova, N., Simpson, D. J., Marioni, R. E. & Chandra, T. Partial reprogramming induces a steady decline in epigenetic age before loss of somatic identity. *Aging Cell* **18**, e12877 (2019).
10. Ocampo, A. *et al.* In vivo amelioration of age-associated hallmarks by partial reprogramming. *Cell* **167**, 1719–1733.e12 (2016).
11. Gill, D. *et al.* Multi-omic rejuvenation of human cells by maturation phase transient reprogramming. *eLife* **11**, e71624 (2022).
12. Shim, V. *et al.* Secretomes from the transient iPSC reprogramming phase reverse cellular photoaging and extend *Caenorhabditis elegans* longevity. *J. Tissue Eng.* (submitted; manuscript JTE-Aug-26-0256).
13. Robinson, M. D. & Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* **11**, R25 (2010).
14. Robinson, M. D., McCarthy, D. J. & Smyth, G. K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* **26**, 139–140 (2010).
15. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* **29**, 15–21 (2013).
16. Frankish, A. *et al.* GENCODE: reference annotation for the human and mouse genomes in 2023. *Nucleic Acids Res.* **51**, D942–D949 (2023).
17. Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995).
18. Ein-Dor, L., Kela, I., Getz, G., Givol, D. & Domany, E. Outcome signature genes in breast cancer: is there a unique set? *Bioinformatics* **21**, 171–178 (2005).
19. Venet, D., Dumont, J. E. & Detours, V. Most random gene expression signatures are significantly associated with breast cancer outcome. *PLoS Comput. Biol.* **7**, e1002240 (2011).
20. Hernandez-Segura, A. *et al.* Unmasking transcriptional heterogeneity in senescent cells. *Curr. Biol.* **27**, 2652–2660.e4 (2017).

# Supplementary Table 1 | Datasets in the benchmark

Twenty-nine GEO series and one local experiment contributed 76 contrasts. Class is assigned by experimental manipulation, not by transcriptional content.

| Accession | Class | Contrasts | Contrast(s) |
|---|---|---:|---|
| GSE113957 | donor age | 1 | donor age, per decade |
| GSE179848 | donor age | 1 | late vs early passage |
| GSE226189 | donor age | 1 | donor age, per decade |
| GSE307377 | donor age | 1 | old vs young dermal fibroblast |
| GSE116968 | metabolic / culture | 2 | uv_vs_control 1h; uv_vs_control 4h |
| GSE179848 | metabolic / culture | 11 | 2-Deoxyglucose vs Control (21% O2); betahydroxybutyrate vs Control (21% O2); Contact_Inhibition vs Control (21% O2); DEX vs Control (21% O2); Galactos |
| GSE116968 | photoprotection / rescue | 6 | prered_vs_uv 1h; prered_vs_uv 4h; Pre-Red vs UV, injury-residualised, 1h; Pre-Red vs UV, injury-residualised, 4h; red_vs_control 1h; red_vs_control 4h |
| GSE222414 | photoprotection / rescue | 2 | UVB + TAT-UIFSP5 vs UVB (rescue); UVB + TAT-UIFSP6 vs UVB (rescue) |
| GSE240226 | photoprotection / rescue | 2 | Maifuyin rescue vs UVA; succinate rescue vs UVA |
| GSE240486 | photoprotection / rescue | 4 | 0.06% Osmoter alone, no UV; 0.02% Osmoter alone, no UV; UV + 0.06% Osmoter vs UV (rescue); UV + 0.02% Osmoter vs UV (rescue) |
| GSE326951 | photoprotection / rescue | 1 | UVB + sauchinone vs UVB (rescue) |
| GSE93535 | photoprotection / rescue | 1 | compound 1201 in quiescent |
| this study | Repro-CM | 1 | reprogramming-phase secretome vs HDF-CM and iPSC-CM |
| GSE149694 | reprogramming | 6 | 5iLAF-D13 vs Fibroblast-D3; Fibroblast-D7 vs Fibroblast-D3; NHSM-D13 vs Fibroblast-D3; Primed-D13 vs Fibroblast-D3; RSeT-D13 vs Fibroblast-D3; t2iLGoY |
| GSE165177 | reprogramming | 1 | transient reprogramming (MPTR) |
| GSE297233 | reprogramming | 2 | O4YRSK +dox vs -dox, day 4; OSK +dox vs -dox, day 4 |
| GSE139563 | secretome | 2 | senescent-cell CM, day 4 vs day 0; senescent-cell CM, day 10 vs day 0 |
| GSE251807 | secretome | 4 | BM-MSC conditioned medium vs DMEM; BM-MSC EV-depleted medium vs DMEM; BM-MSC small EVs vs DMEM; BM-MSC non-small EVs vs DMEM |
| GSE266052 | secretome | 2 | cytokine-primed hMSC-CM vs untreated (both TGF-beta); resting hMSC-CM vs untreated (both TGF-beta) |
| GSE268248 | secretome | 1 | PRF serum vs untreated (gingival, boundary) |
| GSE279804 | secretome | 1 | ADSC nanovesicles vs control |
| GSE282054 | secretome | 1 | hTSC secretome vs ES-CM |
| GSE293186 | secretome | 1 | endothelial-cell EVs vs control, 72 h |
| GSE306748 | secretome | 2 | primary MSC secretome vs control; immortalised MSC secretome vs control |
| GSE109700 | senescence | 2 | deep replicative senescence; early replicative senescence |
| GSE191055 | senescence | 1 | P27 vs P4 |
| GSE306957 | senescence | 4 | DDX58_KO senescent vs proliferating; IFIH1_KO senescent vs proliferating; pLenti senescent vs proliferating; siScramble senescent vs proliferating |
| GSE93535 | senescence | 2 | SIPS vs quiescent; NA |
| GSE125429 | UV injury | 1 | chronic UVA vs control |
| GSE134533 | UV injury | 1 | HDF UVB vs no UV (2 donors) |
| GSE222414 | UV injury | 1 | UVB vs control (injury) |
| GSE240226 | UV injury | 1 | UVA vs control |
| GSE240486 | UV injury | 1 | UV vs control (injury) |
| GSE302943 | UV injury | 1 | cumulative UVA vs control |
| GSE329475 | UV injury | 1 | UVA, 3 consecutive days vs NC (injury) |
| GSE89005 | UV injury | 3 | single UVA 6 h; single UVA 24 h; repeated UVA 24 h |
