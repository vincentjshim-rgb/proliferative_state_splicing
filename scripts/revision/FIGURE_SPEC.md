# Figure specification for the revised manuscript (2026-09-15)

Direction (a) was chosen: splicing **outcome** (use of unannotated junctions) and the
*COL1A2* event leave the core claims. The manuscript's single claim becomes

> In cultured human fibroblasts the abundance of splicing-factor transcripts is a
> readout of proliferative state; the coupling does not transfer to skin or muscle.

Ten figures become **eight main figures and two supplementary figures**. Numbering
changes, so every `save_fig(..., "FigN.png")` call must be renamed accordingly.

| new | old | script | question |
|---|---|---|---|
| Fig. 1 | 1 | `make_sa_fig1.R` | do splicing-machinery transcripts fall with donor age and time in culture? |
| Fig. 2 | 2 | `make_bio_fig2.R` | do they follow the measured division rate? |
| Fig. 3 | 3 | `make_sa_fig5gene.R` | which genes keep an age association after adjustment? |
| Fig. 4 | 4 | `make_sa_fig4new.R` | what does the adjustment remove, and what does the outcome measure show? |
| Fig. 5 | 6 | `make_sa_fig3.R` | does the coupling generalise across 63 contrasts? |
| Fig. 6 | 7 | `make_sa_fighall.R` | which ageing programmes are division readouts? |
| Fig. 7 | 8 | `make_sa_figgtex.R` | is the coupling a property of culture? (preregistered GTEx) |
| Fig. 8 | 9 | `make_bio_fig6.R` | do interventions raise splicing factors as they raise division? |
| Fig. S1 | 5 | `make_splicing_fig.R` | event-level splicing: why it is not interpreted |
| Fig. S2 | new | `make_sa_figgtex.R` | post hoc tests of range restriction and cell composition (GTEx) |
| Fig. S3 | 10 | `make_sashimi_local.R` | read-level coverage at *COL1A2* |

Supplementary items are numbered in order of first citation (Supplementary Results 1-8, tables 1-9, figures S1-S10; main figures 1-6 since 2026-09-17);
`scripts/artifact_build/build_review_bundle.sh` and the reading pages follow that order.

## Cohort (affects Figs. 1c, 3, 4, and Fig. 7e)

`scripts/revision/cohort_gse113957.R`. GSE113957 is no longer the published
142-sample set. Nine Hutchinson-Gilford progeria donors and 26 normal donors under
20 are removed; **the primary cohort is 107 normal donors aged 20 to 96**, and every
donor model carries log sequencing depth, cell repository (Coriell AG vs other),
sequencing instrument and sex as covariates. Sensitivity cohorts: all 133 normal
donors, normal adults 20-82 (n = 76), and the published 142.

Recomputed results live in `public_data_tierA/derived/cohort_revised/`:
`<cohort>/sample_metrics.tsv`, `<cohort>/gene_level.tsv`, `<cohort>/heatmap_matrix.tsv`,
`<cohort>/heatmap_sample_order.tsv`, `sensitivity_metrics.tsv`, `sensitivity_genes.tsv`.

Headline numbers, primary cohort (n = 107):

- machinery (177 genes) age effect **−0.133 per decade, P = 6 × 10⁻⁷**, after
  proliferation **−0.040, P = 0.005** (70% of the effect lost)
- ρ(age, proliferation) = **−0.46**; r(machinery, proliferation) = **0.89**
- gene level: **137** of 177 genes age-associated, **19** keep a decline after
  adjustment, 0 acquire one; coupling predicts how much is kept (ρ = **−0.45**)
- use of unannotated junctions (outcome): **no age association** (β = +1.5 × 10⁻⁵,
  P = 0.51); in the published cohort it was P = 0.002; in adults 20-82 it runs the
  other way (P = 0.012). **This is now reported as a negative result.**
- no gene keeps a decline in all four cohort definitions (published 40, normal 56,
  primary 19, adults 20-82 zero)

## Other recomputed inputs

- `public_data_tierA/derived/revision_stats/` — `fig2_per_line_correlations.tsv`
  (7 cell lines, ρ 0.32-0.88), `fig2_mixed_models.tsv` (β = 0.688, P = 3 × 10⁻⁴⁶),
  `fig1d_culture_time.tsv` (exact P: HC3 0.0028, HC4 0.10), `fig6_contrasts_revised.tsv`
  (63 contrasts, ρ = 0.864 [0.78, 0.92], R² = 0.923), `fig6_class_residuals.tsv`
  (raw P 0.22-0.85, BH 0.85), `fig6_null_sets.tsv` (1000 expression-matched sets:
  median ρ 0.27, max 0.58; observed 0.865), `fig6_other_programmes.tsv`
- `public_data_tierA/derived/hallmark_revised/` — `division_rate_cellcycle_removed.tsv`,
  `donor_age_primary_cohort.tsv`, `senescence_dissection_cellcycle_union.tsv`
- `public_data_tierA/derived/gtex_posthoc/` — range restriction, composition
  adjustment, donor-identity positive control, sample counts
- `public_data_tierA/derived/interventions_revised/` — reprogramming n = 7
  (ρ = 0.89 [0.43, 0.98], P = 0.007), secretome n = 14 (ρ = 0.63 [0.14, 0.87], P = 0.017)
- `public_data_tierA/derived/psi_cohorts/` — event counts by cohort, *COL1A2* PSI by
  age band and repository, local junction read counts

## Typography rules (all figures, from the journal's figure guidelines)

`scripts/sciadv_theme.R` was rewritten; use its helpers and do not hand-format.

- panel letters: **lower-case bold a, b, c** — `lab_grid(..., labels = c("a","b"))`
  (it lower-cases anyway)
- `rp(r, p)` prints `ρ = 0.71` / `P = 1 × 10⁻⁵⁰`; `rpn(r, p, n)` adds `n = 328`;
  `pfmt(p)` alone for a P value; `num(x)` for a number with a real minus sign;
  `num_axis()` as a scale label function; constants `RHO`, `R2`, `MINUS`, `TIMES`
- never print `rho`, `R2`, `9e-04`, `P = < 0.001`, or a hyphen for a negative number
- base font 8 pt (`theme_sa()` default); **nothing below 7 pt** anywhere, including
  gene labels, leader-line labels and axis text
- no red-green contrast: use `PASS` (blue) and `FAIL` (orange) for verdicts
- every panel that reports a correlation gives its *n*, and every *P* says what unit
  it is over (donors, samples, genes, contrasts)
- the 177-gene score has ONE name everywhere: **"pre-mRNA processing score"**
  (axis label) and in prose "splicing-factor expression"; Fig. 7 uses the
  preregistered 96-gene version and says so in the legend
- render at 600 dpi, full width 183 mm, and **look at the PNG** afterwards: fix any
  overlapping text, clipped labels or leader lines pointing at nothing

## Named splicing factors

Figs. 2b, 3a and 8c label "the splicing factors named in the ageing literature".
The three lists used to differ. Read `scripts/revision/named_splicing_factors.txt`
if it exists (one gene per line, sourced list) and fall back to the list in
`scripts/revision/run_interventions.R` otherwise.
