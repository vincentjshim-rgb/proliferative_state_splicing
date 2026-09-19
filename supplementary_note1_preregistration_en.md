---
title: "Supplementary Note 1. Preregistration of the culture-versus-tissue comparison in GTEx"
subtitle: "English translation of `preregistration_gtex_culture_tissue_boundary_ko.md`, with outcomes and deviations"
date: "15 September 2026"
---

# What this document is

The GTEx analysis reported in Fig. 6 was planned before the data were seen. The plan below was written in Korean on 13 September 2026, and its SHA-256 digest was recorded in `public_data_tierA/logs/preregistration_sha256.txt` before any GTEx expression file was downloaded or opened; the analysis script `scripts/run_gtex_boundary.R` recomputes the digest and stops if it does not match. The digest is

`61972b991dc6de770478cd7f6df39a5319e11107ef2993aac2faeee8eb8b774e`

and it still matches the current file, which has not been amended since it was locked.

This was an internal record. It was **not** deposited with a public registry such as OSF or AsPredicted, so a reader cannot verify independently when it was written; what can be verified is that the document the analysis ran against is the document reproduced here, and that its stated criteria are the criteria against which the outcomes below are reported, including the two that the analysis failed. The translation is complete rather than summarised. Gene names, dataset identifiers, statistical terms and numbers are unchanged; where the Korean is ambiguous it is translated literally and a bracketed translator's note is added.

---

# Translation

## Preregistration: culture-versus-tissue boundary test (GTEx v8 via recount3)

- Written: 2026-09-13
- Status: **written and locked before the GTEx expression data were downloaded or opened**
- Locking method: record the SHA-256 of this file in `public_data_tierA/logs/preregistration_sha256.txt`, then run the analysis script.

### 1. Purpose

The manuscript's present claim is as follows. In cultured human fibroblasts, splicing-factor expression and eight ageing programmes (chromatin organisation, cellular senescence, pre-mRNA processing, mRNA splicing, interferon signalling, DNA double-strand break repair, telomere maintenance, base excision repair) are readouts of division rate, whereas collagen formation and the unfolded protein response are not. The limitations section stated that "it says nothing about post-mitotic tissue". This test establishes that boundary from data rather than by assumption. GTEx is the only public resource in which **cultured fibroblasts**, **skin tissue** and **skeletal muscle (post-mitotic tissue)** can be compared under one definition in the same donor population.

### 2. Data (processed only)

| File | Content | Size |
|---|---|---:|
| `gtex.gene_sums.SKIN.G026.gz` + 3 metadata files | Cells – Cultured fibroblasts, Skin – Sun Exposed (Lower leg), Skin – Not Sun Exposed (Suprapubic) | about 196 MB |
| `gtex.gene_sums.MUSCLE.G026.gz` + 3 metadata files | Muscle – Skeletal | about 81 MB |

About 277 MB in total. FASTQ and junction files will not be downloaded. Splicing **outcome** (unannotated junctions) will therefore not be computed in GTEx, and only programme-level expression analysis will be performed.

### 3. Definitions fixed before seeing results

- **Gene ID → symbol**: use the mapping already used in the project (Geneid → Gene_name from `GSE282054_raw_counts.txt.gz`), unchanged.
- **Normalisation**: separately per tissue type. edgeR `filterByExpr` → TMM → log₂ CPM → gene-wise *z*.
- **Programme scores**: the list of 24 Reactome sets in `scripts/run_hallmark_audit.R`, unchanged. Mean *z* of the genes expressed in that tissue; excluded if fewer than 15.
- **Splicing set**: the 96 genes in `conserved_core/age_down_splicing_core.txt`.
- **Proliferation score**: mean *z* of the 20 markers in `scripts/run_outcome_vs_expression.R` (MKI67, CCNB1, CCNA2, CDK1, TOP2A, BIRC5, BUB1, PLK1, AURKA, TYMS, RRM2, PCNA, MCM2–MCM7, TK1, UBE2C).
- **Prevention of circularity**: when a programme is compared with the proliferation score or adjusted for proliferation, the 20 proliferation markers contained in that programme are first removed from the programme score. For the mitotic cell cycle set (the positive control), the adjusted result is not interpreted.
- **Age**: the midpoint of GTEx's released age bands (20–29, …, 70–79), i.e. 25, 35, …, 75.
- **Technical covariates**: RIN, ischaemic time, log₁₀ library size. Fields not present in the recount3 metadata will not be used, and that fact will be recorded.
- **Biological covariates**: sex, Hardy scale of death circumstances (missing values as a separate category).
- **Sample exclusions**: samples not belonging to the four tissue types above; missing age; the lowest 1% of library size within a tissue type. Where one donor has several samples of one tissue type, the one with the largest library size.
- **Definition of the correlation**: the primary measure is the **partial Spearman ρ** after regressing the technical covariates out of both variables (GTEx is post-mortem and surgical material, so technical variation is large). The raw Spearman ρ is a sensitivity analysis.
- **Age model**: `score ~ age + sex + Hardy + RIN + ischaemic time + log10 library size`; the proliferation-adjusted model adds the proliferation score to this. Attenuation = 1 − β_adj/β.

### 4. Primary hypotheses and decision criteria

| | Hypothesis | Success | Failure |
|---|---|---|---|
| **H1** culture replication | Partial ρ between the splicing-96 score and the proliferation score in GTEx cultured fibroblasts | ρ ≥ 0.5 | ρ < 0.3 (0.3–0.5 is partial replication) |
| **H2** coupling weaker in tissue | The splicing-96–proliferation ρ in lower-leg skin tissue is smaller than in cultured fibroblasts | Fisher *z* difference, one-sided *P* < 0.05 | no difference, or the opposite |
| **H3** post-mitotic tissue boundary | Partial ρ between splicing-96 and proliferation in skeletal muscle | ρ < 0.3 | ρ ≥ 0.5 |
| **H4** attenuation of the age effect in culture | In cultured fibroblasts, the age effect on splicing-96 is attenuated after adjustment for proliferation | attenuation ≥ 50% | attenuation < 25%. If the unadjusted age *P* ≥ 0.05, report as **not evaluable** |
| **H5** collagen dissociation in culture | Partial ρ between collagen formation and proliferation in cultured fibroblasts | \|ρ\| < 0.3 | \|ρ\| ≥ 0.5 |
| **H6** positive control for tissue ageing | Collagen formation falls with age in lower-leg skin tissue | β < 0, *P* < 0.05 | β ≥ 0 or *P* ≥ 0.05 (in which case the whole tissue age analysis is interpreted as underpowered) |

Note on H2: because the epidermis of skin contains proliferating basal keratinocytes, the coupling may persist in tissue and H2 may fail. In that case the conclusion is restricted to "skin tissue is not a post-mitotic tissue", and the judgement about the boundary is left to H3.

### 5. Exploratory analyses (no success criteria)

- **S1** Compare the partial ρ with proliferation of the 24 programmes side by side in culture, the two skin sites and skeletal muscle (the tissue version of Fig. 7A). [Translator's note: "Fig. 7A" is the panel number in the draft current at the time of writing; it is Fig. 5a in the revised manuscript.]
- **S2** Compare the age slope of collagen formation between sun-exposed (lower leg) and non-exposed (suprapubic) skin.
- **S3** Donor-level concordance of programme scores between a donor's cultured fibroblasts and skin tissue.
- **S4** Skin tissue composition sensitivity: add keratinocyte marker (KRT5, KRT14, KRT1, KRT10) and dermal fibroblast marker (COL1A1, COL1A2, DCN, LUM, PDGFRA) scores as covariates.
- **S5** Full table of age effects for the 24 programmes × tissue types, before and after adjustment for proliferation.

### 6. Multiple testing

- H1–H6: success is decided by the effect-size criteria above. The *P* value of each hypothesis is reported together after Benjamini–Hochberg correction within the six primary hypotheses.
- S1, S5: within each tissue type, the proliferation family and the age family are BH-corrected separately (as in the existing audit).

### 7. Reporting principles

- Report both successes and failures. Failed hypotheses will not be removed from the manuscript.
- Do not change definitions, thresholds or sample-exclusion rules after seeing results. If a change is unavoidable, record the date and reason in the amendment history below and report **the result under the original definition as well**.
- The material analysed is reprocessed public data and no causal claim is made. GTEx age is given in ten-year bands, so age effects are interpreted only at band resolution.

### Amendment history

(none)

---

# Outcomes

Recorded in `results_gtex_culture_tissue_boundary_ko.md`. Final samples after the preregistered exclusions: cultured fibroblasts 511 donors (492 with complete technical covariates), lower-leg skin 739, suprapubic skin 626, skeletal muscle 839; 388 donors contributed both cultured fibroblasts and lower-leg skin.

| | Preregistered criterion | Estimate obtained | *P* | BH-FDR | Verdict |
|---|---|---|---:|---:|---|
| **H1** | ρ ≥ 0.5 succeeds; ρ < 0.3 fails | partial ρ = +0.739 in cultured fibroblasts (*n* = 492) | 1.2 × 10⁻⁸⁵ | 7.1 × 10⁻⁸⁵ | **supported** |
| **H2** | one-sided *P* < 0.05 for the Fisher *z* difference | lower-leg skin +0.225 (*n* = 739) versus culture +0.739; *z* = −11.95 | 3.2 × 10⁻³³ (one-sided) | 9.6 × 10⁻³³ | **supported** |
| **H3** | ρ < 0.3 succeeds; ρ ≥ 0.5 fails | partial ρ = −0.115 in skeletal muscle (*n* = 839) | 8.5 × 10⁻⁴ | 1.7 × 10⁻³ | **supported** |
| **H4** | attenuation ≥ 50% succeeds; not evaluable if the unadjusted age *P* ≥ 0.05 | unadjusted age effect −0.025 per decade | 0.33 | 0.33 | **not evaluable** |
| **H5** | \|ρ\| < 0.3 succeeds; \|ρ\| ≥ 0.5 fails | partial ρ = +0.049 for collagen formation in culture | 0.28 | 0.33 | **supported** |
| **H6** | β < 0 and *P* < 0.05 succeeds | collagen formation in lower-leg skin +0.027 per decade, i.e. rising | 0.018 | 0.027 | **not supported (opposite direction)** |

Because H6 failed, the preregistered rule was applied: age effects in the GTEx tissues are treated as not validated by a positive control and are not used to support any claim in the manuscript. This applies to the tissue columns of S5 and to the S2 comparison of sun-exposed with non-exposed skin. It was noted at the time, without changing that rule, that the failure was a significant effect in the opposite direction rather than an absence of power.

---

# Deviations and post hoc additions

In the order in which they occurred. The preregistration's own amendment history is empty; the first two items below are recorded in the outcome record, and the rest are later additions, each labelled post hoc.

1. **Before the analysis was run — the *P* value attached to H4.** The preregistration did not state which *P* value belonged to H4. The draft script attached the proliferation-adjusted age *P*; before any result was seen this was changed to the **unadjusted** age *P*, which is the quantity the preregistered rule uses to decide whether H4 is evaluable. The decision rule itself was unchanged. Recorded in `results_gtex_culture_tissue_boundary_ko.md`.

2. **Before the analysis was run — a failed first execution.** The first run stopped before R started because the output directory did not exist. The directory was created and the same script was run again. No definition changed.

3. **The splicing set is the 96-gene version, fixed by the preregistration.** After the preregistration was locked, a fourth donor cohort (GSE226189) was removed from the study on 2026-09-13 because its processing batches coincide with age strata, which changed the study-wide splicing set from 96 genes to 177. The preregistered verdicts above stand on the 96-gene file, which was not rewritten. Repeating the same analysis with the 177-gene set is a **post hoc sensitivity analysis**: cultured fibroblasts +0.766, lower-leg skin +0.203, skeletal muscle −0.130, and the culture age effect −0.024 per decade unadjusted (*P* = 0.35) and −0.052 adjusted (*P* = 0.001). Every verdict is the same and H4 remains not evaluable.

4. **Post hoc — 32 off-trend cultures.** After the outcomes were known, 32 of the 511 cultures were found to lie off the main trend at low proliferation; they are enriched in particular RNA-extraction and sequencing batches and have lower RIN. Excluding them raises the raw correlation from 0.72 to 0.89. The preregistered estimate includes them, and the exclusion is reported only as an observation.

5. **Post hoc — range restriction and variance-matched subsampling** (2026-09-15, `results_gtex_posthoc_ko.md`, `public_data_tierA/derived/gtex_posthoc/`). Added in response to the objection that proliferation varies less between tissue donors than between cultures. The regression slope, which range restriction does not affect, falls from 0.554 in culture to 0.167 in lower-leg skin, 0.072 in suprapubic skin and −0.127 in muscle; a Thorndike case-II correction predicts ρ = 0.653 for lower-leg skin against the 0.225 observed; and cultures subsampled to the skin's proliferation spread give a median partial ρ of 0.872 (95% range 0.834–0.910) with a median-centred window and 0.770 (0.701–0.857) with the tails down-weighted rather than discarded.

6. **Post hoc — adjustment for cell composition** (2026-09-15). This extends the preregistered exploratory analysis S4 in three ways, all of them after the outcomes were known: the marker panels were enlarged (keratinocyte KRT5, KRT14, KRT1, KRT10, TP63, IVL, LOR; fibroblast COL1A1, COL1A2, DCN, LUM, PDGFRA, **PDGFRB**; immune PTPRC; and for muscle myofibre ACTA1, MYH1, MYH2, MYH7, CKM and satellite-cell PAX7, MYF5); the fibroblast panel thus keeps the five genes the preregistration named, PDGFRA among them, and **adds PDGFRB**; and the adjustment was applied to the splicing–proliferation coupling itself, whereas S4 as written applied it to age effects. With composition in the covariates the partial ρ between the splicing-96 score and the proliferation score falls from +0.225 (*P* = 6.6 × 10⁻¹⁰) to +0.101 (*P* = 0.006) in lower-leg skin (*n* = 739 donors), from +0.056 (*P* = 0.17) to −0.058 (*P* = 0.15) in suprapubic skin (*n* = 626), and from −0.115 (*P* = 8.5 × 10⁻⁴) to −0.040 (*P* = 0.25) in muscle (*n* = 839), while culture is essentially unchanged (+0.739 to +0.722, *n* = 492; `public_data_tierA/derived/gtex_posthoc/composition_adjusted.tsv`). The fibroblast panel was revised once, on the same day and after the first post hoc run: that run had PDGFRB in place of the preregistered PDGFRA and gave +0.111 (*P* = 0.003) in lower-leg skin, −0.052 (*P* = 0.19) in suprapubic skin, −0.021 (*P* = 0.55) in muscle and +0.721 in culture, after which PDGFRA was restored and PDGFRB kept as an addition, with no other change to markers, covariates, samples or model; the values reported above and in the main text are from the revised panel, and the first-run table is kept as `composition_adjusted_first_definition.tsv`. The preregistered S4 result on the age effect, which was run with the preregistered panel, stands as recorded (lower-leg skin splicing-96 age effect −0.032, *P* = 0.037, becoming −0.013, *P* = 0.32, with keratinocyte and fibroblast markers), and is not used as evidence, since it is a tissue age effect.

7. **Post hoc — genotype-driven positive control for donor concordance** (2026-09-15). The preregistered exploratory analysis S3 found no donor-level concordance between cultured fibroblasts and skin for any of the 25 programme scores (ρ = −0.10 to +0.13, none at FDR < 0.05). To establish that the design can detect donor identity at all, genes whose expression is strongly genotype-dependent were tested in the same 388 donors and the same matrices: *GSTM1* ρ = +0.750 (*P* = 3 × 10⁻⁷¹) and *ERAP2* ρ = +0.749 (*P* = 5 × 10⁻⁷¹), neither explained by sex, together with seven sex-chromosome genes at ρ = 0.666 to 0.736. This control was not preregistered.

8. **Post hoc — redefinition of the comparison cohort** (2026-09-15). Fig. 6c compares the GTEx cultures with the fibroblast donor cohort GSE113957. That cohort was redefined after the outcomes were known, on grounds independent of the GTEx analysis: ten progeria donors (nine of them in the 142-sample analysed set) and 26 donors under 20 years of age were removed, and cell repository, sequencing instrument and sex were added as covariates, leaving 107 normal adults. For the preregistered 96-gene set this changes the comparison cohort's estimates from −0.093 per decade unadjusted and −0.031 adjusted, as recorded at the time, to −0.128 unadjusted (*P* = 4 × 10⁻⁶) and −0.033 adjusted (95% CI −0.064 to −0.002, *P* = 0.04). Nothing in the GTEx analysis itself was recomputed.

9. **Post hoc — combined residual age effect** (2026-09-16, `public_data_tierA/derived/residual_meta/residual_age_effect_meta.tsv`). The proliferation-adjusted age effects on the preregistered 96-gene set in the GTEx cultures (−0.053 per decade, *n* = 492 donors) and in the 107-donor cohort (−0.033) were combined by inverse variance: −0.042 per decade (95% CI −0.065 to −0.019, 599 donors). Unadjusted estimates were not combined. This was not preregistered; it uses the cultures only, not the tissues to which the rule that followed the failure of H6 applies. Withdrawn at revision (2026-09-18): the pooled value is no longer reported anywhere in the manuscript or its tables; the two estimates are shown side by side in Fig. 6c and Supplementary Table 5 and are not combined, because one is an attenuation and the other a suppression effect in cultures with no unadjusted age effect, and with two studies the between-study variance is not estimable.

10. **Post hoc — specificity of the coupling in GTEx cultures** (2026-09-16, `public_data_tierA/derived/machinery_specificity/gtex_culture_partials.tsv`; Supplementary Fig. S5b, Supplementary Table 5). The comparison of splicing definitions with translation, Metabolism of RNA and rRNA processing, made first against the counted division rate, was repeated against the 20-marker proliferation score in the 511 GTEx cultures. The asymmetry seen against the counted division rate was not reproduced: with the other score held constant, translation keeps a partial ρ of 0.51 to 0.59 and the splicing definitions 0.28 to 0.41. This was not preregistered, and it uses the cultures only.

Two further discrepancies are recorded for completeness. The preregistration estimated the download at about 277 MB; the files actually retrieved total 265 MB, with SHA-256 digests in `public_data_tierA/gtex/gtex_sha256.txt`. And the panel numbering referred to in S1 is that of the draft current on 2026-09-13, not of the revised manuscript.
