---
title: "Cover letter · GeroScience"
date: "18 September 2026"
---

[Date]

Editor-in-Chief
*GeroScience*
American Aging Association

Dear Editor,

We submit for your consideration the manuscript **"Splicing-factor expression in cultured human fibroblasts reports proliferative state and the coupling weakens in skin"** as an Original Article.

**What the study asks.** Splicing-factor transcripts fall with age in cultured fibroblasts, and the field reads that fall as a signal of ageing. Those genes are also transcribed in proportion to growth. We asked how much of the age signal measured in these cells is a readout of how fast the cells divide, where that coupling stops, and what survives it.

**What we found.** Using a resource in which the division rate of every culture was measured by cell counting rather than inferred from the transcriptome, expression of 177 pre-mRNA processing genes follows measured division rate across 328 libraries (ρ = 0.71; mixed-model β = 0.69 with the cell line as a random effect). In 107 adult fibroblast donors, proliferation accounts for 70% of the age effect on these genes (63 to 73% across four cohort definitions), and after adjustment a decline of about 0.04 standard deviations per decade remains in those donors. We do not claim that this residue replicates. It halves again when a second gene-expression programme enters the model, and in the GTEx cultures — where no age effect was visible before adjustment — the adjusted decline is a suppression effect that 15 of 25 ageing programmes acquire, the mitotic cell-cycle positive control among them. We therefore report no age effect independent of proliferation. Across 63 contrasts from 23 studies the change in splicing genes tracks the cell-cycle change (ρ = 0.86), and eight of twenty-three Reactome ageing programmes follow measured division rate at |ρ| > 0.5, four of them as closely as the cell-cycle genes themselves. A preregistered comparison in GTEx locates the boundary: the coupling is strong in cultured fibroblasts (partial ρ = 0.74) and weak at both skin sites (0.23 and 0.06), while in those same donors the same proliferation score still recovers the mitotic cell-cycle programme at 0.67 and 0.60 — so the score still registers proliferation in that tissue. In skeletal muscle that score does not recover the mitotic cell-cycle programme (ρ = 0.03), so muscle cannot test the coupling and we do not use it as evidence. Post hoc, range restriction does not account for the fall, and cell composition accounts for much of what is left in skin: adjusting for keratinocyte, fibroblast and immune markers takes splicing from 0.23 to 0.10 in sun-exposed skin and from 0.06 to −0.06 in non-exposed skin while leaving the cell-cycle programme where it was (0.670 to 0.666, 0.597 to 0.586).

**Why GeroScience.** The practical consequences fall squarely in this journal's readership. Interventions are screened in cultured fibroblasts and scored by the restoration of these programmes. No class of intervention we examined shows a detectable class-level shift in splicing-factor expression beyond what its effect on the cell cycle predicts, although for reprogramming media that test has little power and individual conditions depart from the prediction in both directions. Senescence scores computed in culture depend on which panel is used, and a curated senescence pathway scores *higher* in faster-dividing cultures. We give the numbers a laboratory can act on: report the division rate beside the signature, state the panel and its cell-cycle content, and report the direction in which proliferation moves with age in the cohort used.

**How the work was done.** The study is computational and uses public data throughout; no data were generated for it. Analyses added after the main results were known are labelled post hoc wherever they appear, and those that are not on the main line of argument are reported in full as Supplementary Results. The GTEx comparison is one of the study's main results and it was preregistered: the hypotheses, definitions, exclusions and success criteria were written down and the document's SHA-256 digest recorded **before** the data were downloaded, and the analysis script verifies that digest before it runs. We then obtained the skin and skeletal-muscle projects from recount3 and analysed 2,715 samples, one per donor per tissue — cultured fibroblasts from 511 donors, sun-exposed skin from 739, non-exposed skin from 626 and skeletal muscle from 839. All six preregistered outcomes are reported. Four were supported, one could not be evaluated because there was no unadjusted age effect to attenuate, and one — a positive control asking whether collagen formation falls with age in sun-exposed skin — failed. Under the rule fixed in advance, the failure of that collagen positive control (H6) bars us from using age effects measured **within** GTEx tissues as evidence anywhere in the paper. It does not touch the tissue comparison itself, which asks how strongly the splicing genes track proliferation in each tissue and is the result we report. A second preregistered analysis, of methylation clocks in the cultures whose division rate was counted, met one of its four criteria; it is reported in the main text and in Supplementary Note 2, and we make no claim about clocks.

**What we report against ourselves.** Five findings weaken easy readings of the claim, and each is given with its number. First, the attenuation is not specific to splicing genes: the same adjustment applied to 1,000 expression-matched random gene sets removes a median of 56% of their age effect, so the 70% we measure is ordinary, and what no random set reproduces is the size of the association on either side of the adjustment. Second, the residue that survives adjustment depends on the cohort's technical covariates — −0.040 per decade (*P* = 0.005) with cell repository, sequencing instrument and sex in the model, and −0.022 (*P* = 0.11) with sequencing depth alone. Third, the 31 donors aged 83 and over all come from one cell repository, so age and cell source cannot be separated in that stratum; the principal donor-level estimates are therefore also given for donors aged 20 to 82 (Supplementary Fig. S2d and S3b, c; Supplementary Tables 1 and 5). Fourth, when division is perturbed by the treatments of the counted resource, splicing-factor expression follows the transcriptional proliferation programme more closely than the counted rate, and two perturbations — oligomycin and growth at 3% oxygen — leave a residue that neither measure explains. Fifth, an index of splicing outcome that we had expected to behave differently from splicing-factor abundance showed an age association in the deposited sample set (*P* = 0.002) and none in the adult cohort (*P* = 0.51); we report that as a negative result, in full in the Supplementary Results, rather than removing it.

**Related work.** A manuscript with overlapping authorship, *[title]* (Shim, V. *et al.*), is under consideration at the *Journal of Tissue Engineering* (JTE-Aug-26-0256). A copy accompanies this submission. The two share no data: that work is an experimental study, and the present manuscript is a reanalysis of public data that cites it in the Discussion as the earlier work that prompted the question. No finding of it is tested here.

**Availability.** All datasets are public and listed with accessions in Supplementary Table 4. Analysis code, gene lists, contrast definitions, both preregistration documents with their locked SHA-256 digests, and the thirteen Supplementary Data tables are public at https://github.com/vincentjshim-rgb/proliferative_state_splicing. We will archive the accepted release with a Zenodo DOI.

The manuscript is not under consideration elsewhere, all authors have approved the submission, and we declare no competing interests [confirm].

Thank you for considering our work.

Yours sincerely,

[Corresponding author]
[Affiliation]
[Email]

---

**Suggested reviewers.** [To be completed. Expertise to aim for: (1) fibroblast replicative ageing and the Cellular Lifespan Study; (2) transcriptomic or epigenetic clocks and their confounders; (3) splicing regulation in ageing; (4) senescence marker panels. Avoid anyone from the groups whose datasets carry the analysis, and anyone at the authors' institutions.]

**Opposed reviewers.** [Optional.]
