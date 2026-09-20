# Splicing-factor expression in cultured human fibroblasts reports proliferative state

Analysis code, gene lists, contrast definitions, preregistrations and derived tables for the
manuscript of the same name (submitted to *GeroScience*).

The study is a reanalysis of public data. **No sequencing data were generated for it**, and no
raw data are redistributed here: every dataset is obtainable from its own repository under the
accessions listed in Supplementary Table 4 of the manuscript.

## What is in this repository

| Path | Contents |
|---|---|
| `scripts/` | All analysis and figure code. `scripts/revision/` holds the current analyses; `scripts/artifact_build/` builds the manuscript, the supplementary PDF and the reading artifacts. |
| `supplementary_data/` | The fourteen machine-readable Supplementary Data tables cited in the manuscript. |
| `public_data_tierA/derived/figures_sciadv/` | The rendered figures (Fig. 1–7, Supplementary Fig. S1–S8), 600 dpi PNG. |
| `manuscript_splicing_v2.md`, `manuscript_splicing_v2_ko.md` | The manuscript, English and Korean. |
| `supplementary_results_en.md`, `supplementary_results_ko.md` | Supplementary Results 1–8 in full. |
| `preregistration_*.md`, `supplementary_note*_en.md` | The two preregistrations as locked, their English translations, outcomes and recorded deviations. |
| `PACKAGES.tsv` | Every R package present in the analysis library, with its version. |
| `CLAUDE.md` | The project's standing analysis rules, its negative results and its revision history. |

## Preregistration

Two analyses were preregistered as internal records — not deposited with a public registry, which
the manuscript states plainly. Each document was locked by its SHA-256 digest before the data were
downloaded, and the analysis script recomputes the digest and stops if it does not match:

| Document | SHA-256 | Checked by |
|---|---|---|
| `preregistration_gtex_culture_tissue_boundary_ko.md` | `61972b99…` | `scripts/run_gtex_boundary.R` |
| `preregistration_methylation_clock_ko.md` | `be51a755…` | `scripts/revision/run_methylation_clock.R` |

Outcomes, including the two hypotheses that failed and the one that could not be evaluated, are in
Supplementary Notes 1 and 2 and Supplementary Table 9.

## Reproducing the analysis

Requirements: R 4.5.1 with the packages in `PACKAGES.tsv`, Python 3, pandoc and LibreOffice for the
document builds. The R scripts expect an analysis library at `analysis_r_lib/` relative to the
project root; point `.libPaths()` at your own library instead if you prefer.

```bash
export PROJECT_ROOT="$(pwd)"

# figures (each script writes one numbered figure into public_data_tierA/derived/figures_sciadv/)
Rscript scripts/make_sa_fig1.R        # Fig. 1
Rscript scripts/make_bio_fig2.R       # Fig. 2
Rscript scripts/make_sa_fig3transcript.R # Fig. 3
Rscript scripts/make_sa_fig3merged.R  # Fig. 4
Rscript scripts/make_sa_fig3.R        # Fig. 5
Rscript scripts/make_sa_fighall.R     # Fig. 6
Rscript scripts/make_sa_figgtex.R     # Fig. 7, Supplementary Fig. S8, Supplementary Table 9

# documents
bash scripts/artifact_build/build_supplementary_pdf.sh
bash scripts/artifact_build/build_reading_artifacts.sh
```

The figure scripts read derived tables under `public_data_tierA/derived/`, which are produced by the
analysis scripts in `scripts/revision/` from the public datasets. Those inputs are not in this
repository; download them under the accessions in Supplementary Table 4 first.

### Retired scripts are guarded

Twenty-eight superseded scripts write to filenames the current figures use. Each begins with a
`stop()` guard and will refuse to run. Do not remove the guards: one of them would overwrite the
GTEx boundary figure with output from an analysis that still contains a withdrawn contrast. To run
one deliberately, set `ALLOW_RETIRED_FIGURE_SCRIPT=1`.

## Data availability

Accessions for every public dataset are in Supplementary Table 4. GTEx v8 was obtained through
recount3; the DNA methylation arrays are GSE179847; the donor cohort excluded from the study is
GSE226189. Derived tables underlying each figure are in `supplementary_data/`.

## Related manuscript

A manuscript with overlapping authorship is under consideration at the *Journal of Tissue
Engineering* (JTE-Aug-26-0256). The present study shares no data with it and cites it in the
Discussion as the earlier work that prompted the question.

## Licence

Code is released under the MIT Licence (`LICENSE`). Gene lists and derived tables are released under
CC BY 4.0. The underlying public datasets remain under the terms of their own repositories.
