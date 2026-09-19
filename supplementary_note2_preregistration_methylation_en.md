---
title: "Supplementary Note 2. Preregistration of the methylation clock comparison"
subtitle: "English translation of `preregistration_methylation_clock_ko.md`, with outcomes and deviations"
date: "16 September 2026"
---

# What this document is

The methylation analysis reported in Supplementary Fig. S7 was planned before the data were seen. The plan below was written in Korean on 16 September 2026, and its SHA-256 digest was recorded before the methylation matrix of GSE179847 was downloaded; the analysis script `scripts/revision/run_methylation_clock.R` recomputes the digest and stops if it does not match. The digest is

`be51a7550c982e2071ebdb80f36ab7098965351c6e14a65434e8df275505d929`

and it still matches the current file, which has not been amended since it was locked.

As with Supplementary Note 1, this was an internal record and was not deposited with a public registry, so a reader cannot verify independently when it was written; what can be verified is that the document the analysis ran against is the document reproduced here, and that its criteria are the criteria against which the outcomes below are reported, including the three the analysis failed. The translation is complete rather than summarised.

---

# Translation

## Preregistration: the relationship between methylation clocks and the transcriptomic proliferation readout

This document was written **before downloading** GSE179847 (Cellular Lifespan Study DNA methylation, EPIC array, 479 samples). The SHA-256 is recorded, then the data are downloaded, and the analysis script checks the digest before it runs. After locking, additions are made only to the amendment history at the end of the document.

### 1. Background and the question

This study has shown that transcriptional programmes in cultured fibroblasts, splicing factors among them, are readouts of the **measured division rate**. That methylation clocks run fast in culture (Sturm 2019) and that cell division causes methylation loss in late-replicating domains (Endicott 2022) are already known. **"The clock follows division too" is therefore not itself a new question.**

The new question is one. **In the same cultures, are the two biomarker families — the transcriptomic proliferation readout and the methylation clock — measuring the same variable?** GSE179847 and GSE179848 come from the same cell lines and the same passages, so they can be matched directly at the sample level.

### 2. Data

- **GSE179847**: `GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz` (about 2.0 GB). Rows are CpGs, columns are per-sample beta values and detection *P* values. Sample identifiers are sentrix IDs, linked to GSM and `unique_variable_name` through the GEO sample list.
- **GSE179848**: the RNA-seq data already held. 328 libraries carry a measured division rate.
- Clock coefficients: the Horvath 2013 multi-tissue clock (353 CpGs). The skin and blood clock (391) and PhenoAge (513) are used only if their coefficient tables can be obtained from a public source; if they cannot, that fact is recorded.

### 3. Definitions fixed in advance

- **DNAm age**: computed with each clock's published coefficients. For Horvath 2013 the published inverse transformation (adult basis) is used as it stands. Missing CpGs are replaced by the cohort mean of that CpG and **the number replaced is reported**.
- **Clock acceleration**: the residual of DNAm age regressed on days in culture within a cell line.
- **Measured division rate**: `24 / population_doubling_time_uhours_per_division` (divisions per day). This name is not used for any other dataset.
- **Transcriptomic proliferation score**: the mean z of the existing 20 markers (MKI67 and the rest).
- **Splicing score**: the mean z of the existing 177-gene pre-mRNA processing set.
- **Matching**: only samples with exactly the same `unique_variable_name` are treated as a pair. If no exact match exists, matching is restricted to samples identical in cell line, treatment, oxygen and passage, and the number matched is reported.
- **Note on normalisation**: this analysis does not reproduce the original paper's pipeline (BMIQ and the rest), so **absolute age estimates are approximations.** Absolute values are not interpreted; only correlations and residuals are.

### 4. Hypotheses and success criteria

| | Hypothesis | Success criterion |
|---|---|---|
| **M1** (sanity check) | In untreated healthy cell lines, DNAm age increases with cumulative divisions | ρ ≥ 0.3, *P* < 0.05 |
| **M2** | Clock acceleration is associated with the **measured division rate** | \|ρ\| ≥ 0.3, *P* < 0.05 |
| **M3** (the main test) | In the same cultures, the transcriptomic proliferation score is associated with DNAm age | \|ρ\| ≥ 0.3, *P* < 0.05 |
| **M4** | The association between the splicing score and DNAm age shrinks when the measured division rate is adjusted for | the partial correlation falls by at least 50% |

- If M1 fails, the clock computation is taken to be wrong and **M2–M4 are not interpreted.** In that case the failure is recorded as it stands.
- If any hypothesis comes out in the opposite direction it is reported as it stands. Criteria are not changed after results are seen.
- Multiple testing is corrected across the four hypotheses with Benjamini–Hochberg.

### 5. Exclusion criteria

- Measurements with a detection *P* value above 0.01 are treated as missing.
- Samples missing 30% or more of the clock CpGs are excluded and their number reported.
- Samples with no division rate, or a rate of zero or below, are excluded (the same rule as the RNA analysis).

### 6. Reporting principles

- That absolute DNAm ages are approximations is stated in the figure and in the text.
- This analysis occupies **no more than one panel and one paragraph** in the paper. The claim is limited to one thing: whether the two biomarker families read the same variable in culture.
- Endicott 2022 and Sturm 2019 are named as prior work, and what this analysis re-confirms is distinguished from what it adds.

### 7. Amendment history

(none)

---

# Outcomes

The analysis was run on 16 September 2026 with `scripts/revision/run_methylation_clock.R`, which verified the digest above before it started. Of 479 arrays, none was missing more than 30% of clock CpGs; the detection filter set 0.018% of measurements to missing and the 99 remaining gaps were filled with the cohort mean of the CpG. Clock CpGs available: 334 of 353 (Horvath 2013), 391 of 391 (skin and blood), 513 of 513 (PhenoAge), 63 of 69 (Hannum). After the division-rate exclusion, 471 methylation samples from nine cell lines and nineteen treatments remained, and exact matching on `unique_variable_name` gave 109 cultures with both a clock and transcriptional scores, from lines HC1, HC2 and HC4.

| | Test | ρ | *P* | BH-FDR | n | Criterion met |
|---|---|---|---|---|---|---|
| M1 | DNAm age vs cumulative doublings, untreated healthy lines | +0.356 | 8 × 10⁻⁴ | 0.0015 | 86 | **yes** |
| M2 | clock acceleration vs measured division rate | −0.014 | 0.76 | 0.76 | 471 | no |
| M3 | proliferation score vs DNAm age, paired cultures | +0.162 | 0.093 | 0.124 | 109 | no |
| M4 | splicing vs DNAm age with division rate held constant | −0.338 | 3 × 10⁻⁴ | 0.0013 | 109 | no |

M1 was met, so M2–M4 are interpreted under the preregistered rule. M2 and M3 failed on the effect-size criterion fixed in advance (|ρ| ≥ 0.3), not on significance alone: the 95% confidence interval for M3 is −0.03 to +0.34, so a correlation of the preregistered size is not excluded by these 109 cultures, and a larger paired set could change the verdict. M4 asked for attenuation of at least 50%; the association did not attenuate but changed sign, from +0.15 unadjusted to −0.34 with division rate held constant (an attenuation of −131%), so the criterion is not met. That sign change is recorded and not interpreted.

Supporting numbers, none of them preregistered tests. Within the 109 paired cultures the transcriptional scores follow the measured division rate as they do throughout the study (splicing ρ = +0.77, proliferation +0.69, and the two agree with each other at +0.82), while DNAm age follows it at +0.41 before time in culture is removed and at zero after. Three further clocks, computed the same way, give acceleration–rate correlations of +0.18 for skin and blood (*P* = 7 × 10⁻⁵), +0.14 for Hannum (*P* = 0.003) and +0.06 for PhenoAge (*P* = 0.18); all are below the preregistered 0.3.

# Deviations and post hoc additions

1. **The Hannum clock was added.** The preregistration named Horvath 2013, the skin and blood clock and PhenoAge. Hannum's clock was computed as well, because its coefficient table was obtained at the same time. It is reported only as a supporting comparison in Supplementary Fig. S7d and carries none of the preregistered verdicts.

2. **How "untreated healthy cell lines" was operationalised for M1.** The preregistration did not define the phrase. It was taken as the six healthy-donor lines (HC1–HC6) under the study's own `Control` treatment label at 21% oxygen, which gives n = 86. The SURF1 patient lines and every treated or hypoxic culture are excluded from M1 only; M2 uses all 471.

3. **Clock acceleration was fitted with a common time slope.** "The residual of DNAm age regressed on days in culture within a cell line" was implemented as a single regression on days with a cell-line term, which gives every line the same slope. Fitting a separate slope per line instead gives ρ = +0.092 against the measured division rate (*P* = 0.045), still well below the preregistered 0.3, and the correlation within individual lines ranges from −0.74 to +0.40. The verdict on M2 is the same either way.

4. **Post hoc — the confidence interval reported for M3.** The preregistration fixed a threshold but did not ask for an interval. The Fisher-*z* interval is given in the manuscript and above so that the failure is not read as evidence of absence.

5. **Reporting limit** (recorded 2026-09-17). The preregistration allowed one panel and one paragraph. The analysis is reported as one four-panel supplementary figure (Supplementary Fig. S7: the sanity check M1, the tests M2 and M3, and the supporting clocks; M4 is not drawn, and its estimate is in Supplementary Table 8) so that every outcome is visible, with one Results paragraph, one sentence of Discussion, one Methods subsection and Supplementary Table 8. The claim remains limited as preregistered, and no claim about clocks is made.
