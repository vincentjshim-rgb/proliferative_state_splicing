"""Assemble the Supplementary Information as one document: Supplementary Results 1-8,
figures S1-S8 with their legends, Tables 1-10, Notes 1 and 2 (the preregistrations in
English) and the list of
Supplementary Data files. Written as markdown, converted to docx by pandoc and to
PDF by LibreOffice, because the journal asks for a single PDF.

Run through scripts/artifact_build/build_supplementary_pdf.sh.

Every table is read from the derived TSVs at build time, and the statements in the table
headings that depend on them (arms modelled, genes passing the filter, counts in the
negative control) are computed here, not typed in. Supplementary Table 3 (b) and (f) read
reviewer_sensitivities/subset_division_rate.tsv and covariate_attribution.tsv, which
scripts/revision/run_covariate_attribution.R writes.
"""
import os
import pathlib, re, sys, csv

WS = pathlib.Path(os.environ.get("PROJECT_ROOT", pathlib.Path(__file__).resolve().parents[2]))
D = WS / "public_data_tierA/derived"
OUT = pathlib.Path(sys.argv[1])          # working directory for the build
FIG = D / "figures_sciadv"




SHORT_METRIC = {"machinery": "splicing-factor expression",
                "machinery96": "splicing-factor expression (96-gene set)",
                "unannot_reads": "unannotated junction reads",
                "unannot_junc": "distinct unannotated junctions",
                "entropy": "entropy of junction usage",
                "spliceo": "core spliceosome (127 genes)",
                "rest": "capping, 3'-end, export (50 genes)"}

SHORT_RULE = {"comparator only (sensitivity)": "comparator only"}

SHORT_MARKERS = {"fibro+immune": "fibroblast, immune",
                 "kerat+fibro+immune": "keratinocyte, fibroblast, immune",
                 "myofibre+satellite+fibro+immune": "myofibre, satellite cell, fibroblast, immune"}

SHORT_SET = {"splicing 96": "96", "splicing 177": "177"}

SHORT_GROUP = {"sex": "sex chromosome", "eqtl": "genotype-dependent"}

SHORT_COHORT = {"107 normal adults 20+ (primary)": "adults 20+ (primary)",
                "133 normal donors": "all normal (133)",
                "76 normal adults 20-82": "adults 20-82",
                "deposited 142 samples": "deposited (142)"}


TISSUES = ["culture", "legskin", "pubskin", "muscle"]
TISSUE_NAME = {"culture": "cultured fibroblasts", "legskin": "skin, sun-exposed",
               "pubskin": "skin, not exposed", "muscle": "skeletal muscle"}


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table(head, body):
    """A pipe table. Negative numbers get a real minus sign, as in the text.

    A literal | inside a cell ("median |change|", "splicing | comparator") is escaped:
    unescaped, it split the header into extra cells, and pandoc then dropped the last
    column names and shifted the others.

    The dashes of the separator carry relative column widths, which pandoc applies to any
    table with a row longer than its line limit (equal dashes gave equal widths, so a column
    of names was as narrow as a column of numbers). A page holds about 65 characters across
    (WIDTH). Columns of short unbroken values (numbers, flags) come first: they get the width
    of their longest value, or of the longest word of their heading if the page has room for
    that, so that no number is broken across lines. Columns of text share what is left in
    proportion to their longest cell, and wrap; in a table too wide for the page they are
    squeezed to a minimum before the numbers are."""
    WIDTH, PAD, SQUEEZED = 65, 2, 6.5
    raw_head = [str(h) for h in head]
    body = [[re.sub(r"^-(?=[\d.])", "−", str(c)) for c in r] for r in body]
    n = len(raw_head)
    col = [[r[j] for r in body] for j in range(n)]
    longest = [max([len(c) for c in col[j]] + [0]) for j in range(n)]
    word = [max([len(t) for c in col[j] for t in c.split()] + [0]) for j in range(n)]
    hword = [max([len(t) for t in raw_head[j].split()] + [0]) for j in range(n)]
    fixed = [longest[j] <= 10 and all(" " not in c for c in col[j]) for j in range(n)]
    text = [j for j in range(n) if not fixed[j]]
    floor = [min(max(word[j], hword[j]), 12) + PAD for j in range(n)]
    want = [max(min(longest[j], 30) + PAD, floor[j]) for j in range(n)]
    for need in ([max(longest[j], hword[j]) + PAD for j in range(n)],
                 [max(longest[j], 3) + PAD for j in range(n)]):
        units = [need[j] if fixed[j] else floor[j] for j in range(n)]
        if sum(units) <= WIDTH:
            break
    spare = WIDTH - sum(units)
    if spare >= 0:                                   # everything fits: text takes the spare
        room = [0 if fixed[j] else want[j] - floor[j] for j in range(n)]
        if sum(room) > 0:
            units = [units[j] + min(room[j], spare * room[j] / sum(room)) for j in range(n)]
    else:                                            # too wide: squeeze the text, then the rest
        left = WIDTH - sum(need[j] for j in range(n) if fixed[j])
        if left >= SQUEEZED * len(text):
            share = sum(floor[j] for j in text)
            units = [need[j] if fixed[j] else left * floor[j] / share for j in range(n)]
        else:
            scale = (WIDTH - SQUEEZED * len(text)) / sum(need[j] for j in range(n) if fixed[j])
            units = [need[j] * scale if fixed[j] else SQUEEZED for j in range(n)]
    dash = [max(3, round(2 * u)) for u in units]
    esc = lambda c: c.replace("|", "\\|")
    out = ["| " + " | ".join(esc(h) for h in raw_head) + " |",
           "|" + "|".join("-" * d for d in dash) + "|"]
    out += ["| " + " | ".join(esc(c) for c in r) + " |" for r in body]
    return "\n".join(out)


def pivot(path):
    """one row per programme, one column per tissue"""
    r = read(path)
    progs = []
    for x in r:
        if x["programme"] not in progs:
            progs.append(x["programme"])
    body = []
    for pr in progs:
        cells = []
        for t in TISSUES:
            v = [x for x in r if x["programme"] == pr and x["tissue"] == t]
            cells.append(f"{float(v[0]['rho_prolif']):.2f}" if v else "")
        body.append([pr] + cells)
    return table(["programme"] + [TISSUE_NAME[t] for t in TISSUES], body)


def sources(path):
    seen = {}
    for x in read(path):
        seen.setdefault(x["accession"], x["source_publication"])
    return table(["accession", "source publication"], [[k, v] for k, v in seen.items()])


def tsv(path, cols=None, rename=None, rows=None, fmt=None, short=None, derive=None):
    """read a TSV and return markdown rows, optionally selecting and renaming columns;
    derive maps a column name to a function of the row, for cells built from several fields"""
    r = read(path)
    if rows:
        r = [x for x in r if rows(x)]
    cols = cols or list(r[0].keys())
    head = [(rename or {}).get(c, c.replace("_", " ")) for c in cols]
    body = []
    for x in r:
        cells = []
        for c in cols:
            if derive and c in derive:
                cells.append(str(derive[c](x)))
                continue
            v = x.get(c, "")
            if short and c in short:
                v = short[c].get(v, v)
            if c == "pct_lost" and float(x.get("p", 1)) >= 0.05:
                v = "-"
            if fmt and c in fmt:
                try:
                    v = fmt[c](float(v))
                except (ValueError, TypeError):
                    pass
            cells.append(str(v))
        body.append(cells)
    return table(head, body)


def sig(n, d=2):
    return f"{float(n):.{d}g}"


def num(n, d=3):
    return f"{float(n):.{d}f}"


def span(a, b, d=3):
    """an interval as one cell, 'lower to upper', which keeps wide tables within the page"""
    return " to ".join(num(v, d).replace("-", "−") for v in (a, b))


def signed(v, d=3):
    return ("−" if v < 0 else "+") + f"{abs(v):.{d}f}"


def quantile(v, q):
    """linear interpolation between order statistics (the default of R's quantile)"""
    v = sorted(v)
    h = (len(v) - 1) * q
    lo = int(h)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (h - lo)


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        rk, i = [0.0] * len(v), 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                rk[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return rk
    ra, rb = rank(a), rank(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    cov = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    return cov / (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5


def negative_control(path, ref_path):
    """Summary of the random-set negative control for the proliferation adjustment, computed
    from the 1,000 draws; the splicing set's own row comes from the subset table."""
    finite = lambda x: all(x[k] == x[k] and abs(x[k]) != float("inf")
                           for k in ("beta", "p", "beta_adj", "p_adj", "lost"))
    r = [{k: float(v) for k, v in x.items()} for x in read(path)]
    r = [x for x in r if finite(x)]
    ref = [x for x in read(ref_path) if x["metric"] == "machinery"][0]
    rb, rp, rba, rpa, rl = (float(ref[k]) for k in ("beta", "p", "beta_adj", "p_adj", "pct_lost"))
    s = [x for x in r if x["p"] < 0.05]
    lost = [x["lost"] for x in s]
    as_much = sum(v >= rl for v in lost)
    body = [
        ["random gene sets drawn", f"{len(r)}"],
        ["sets with an age association of their own (unadjusted P < 0.05)", f"{len(s)}"],
        ["of those, increasing with age / decreasing with age",
         f"{sum(x['beta'] > 0 for x in s)} / {sum(x['beta'] < 0 for x in s)}"],
        ["share of their age effect removed by adjustment, %: median", f"{quantile(lost, .5):.0f}"],
        ["the same, interquartile range", f"{quantile(lost, .25):.0f} to {quantile(lost, .75):.0f}"],
        ["the same, 2.5th to 97.5th percentile", f"{quantile(lost, .025):.0f} to {quantile(lost, .975):.0f}"],
        [f"of those, sets losing at least as much as the splicing set ({rl:.1f}%)",
         f"{as_much} of {len(s)} ({100 * as_much / len(s):.0f}%)"],
        ["of those, sets keeping an age effect after adjustment (adjusted P < 0.05)",
         f"{sum(x['p_adj'] < 0.05 for x in s)} of {len(s)}"],
        ["unadjusted age effect per decade, range over all sets",
         f"{signed(min(x['beta'] for x in r))} to {signed(max(x['beta'] for x in r))}"],
        ["smallest unadjusted P over all sets", sig(min(x["p"] for x in r))],
        ["smallest adjusted P over all sets", sig(min(x["p_adj"] for x in r))],
        ["splicing set, unadjusted age effect per decade (P)", f"{rb:.3f} ({sig(rp)})"],
        ["splicing set, adjusted age effect per decade (P)", f"{rba:.3f} ({sig(rpa)})"],
        ["random sets with an unadjusted effect at least as large as the splicing set's",
         f"{sum(abs(x['beta']) >= abs(rb) for x in r)}"],
        ["random sets with an unadjusted P at least as small", f"{sum(x['p'] <= rp for x in r)}"],
        ["random sets with an adjusted P at least as small", f"{sum(x['p_adj'] <= rpa for x in r)}"],
    ]
    return table(["quantity", "value"], body)


def hiseq_note(path):
    """how the two entangled technical covariates of the primary cohort overlap"""
    r = read(path)
    hi = [x for x in r if "HiSeq" in x["instr"]]
    ag = sum(x["repo"] == "AG" for x in hi)
    return (f"{len(hi)} of the {len(r)} libraries were sequenced on the HiSeq 2500, "
            f"{'all' if ag == len(hi) else str(ag)} of them from the AG repository")


def gtex_age(path, tissue):
    """age effects of every programme in one GTEx tissue, before and after adjustment"""
    r = [x for x in read(path) if x["tissue"] == tissue]
    n = sorted({x["n"] for x in r})
    body = [[x["programme"], num(x["beta_age"]), sig(x["p_age"]), sig(x["FDR_age"]),
             num(x["beta_age_adj"]), sig(x["p_age_adj"])] for x in r]
    return (f"*{TISSUE_NAME[tissue].capitalize()}* ({' to '.join(n)} donors).\n\n"
            + table(["programme", "per decade", "P", "FDR", "per decade (adj)", "P (adj)"], body))


def concordance_range(path):
    rho = [float(x["rho"]) for x in read(path)]
    n = sorted({x["n_donors"] for x in read(path)})
    return (f"rho = {signed(min(rho), 2)} to {signed(max(rho), 2)} for the {len(rho)} programme "
            f"scores in the same {n[0]} donors")


def legend(tag):
    """pull a supplementary legend out of the manuscript so the two never drift"""
    s = (WS / "manuscript_splicing_v2.md").read_text()
    m = re.search(r"\*\*Supplementary " + tag + r"\.(.*?)\*\*(.*?)(?:\n\n)", s, re.S)
    if not m:
        return ""
    return ("**Supplementary " + tag + ". " + m.group(1).strip() + "** "
            + m.group(2).strip()).replace("\n", " ")


title = re.search(r'^title: "(.*)"', (WS / "manuscript_splicing_v2.md").read_text(), re.M).group(1)


def md_file(name, demote=True):
    """a project markdown file without its front matter, headings pushed down one level"""
    x = (WS / name).read_text()
    if x.startswith("---"):
        x = x.split("---", 2)[2].lstrip("\n")
    return re.sub(r"^# ", "## ", x, flags=re.M) if demote else x


## statements in the table headers that depend on the tables are derived from them
_arms = read(D / "treatment_perturbation/treatment_counts.tsv")
_out = [x["cond"] for x in _arms if x["modelled"] != "TRUE"]
arms_note = ("no arm fell below these limits" if not _out
             else "the arms left out were " + ", ".join(_out))

_sel = [x for x in read(D / "machinery_specificity/counted_libraries_partials.tsv")
        if x["splicing_definition"] == "selected_177"]
_rate = read(D / "division_rate/sample_division_rate.tsv")
matrix_note = (f"{_sel[0]['n_genes_splicing']} of the 177 genes pass it, and there the 177-gene "
               f"score correlates with the measured rate at "
               f"{float(_sel[0]['rho_splicing_whole_set']):.2f} "
               f"({spearman([float(x['splice']) for x in _rate], [float(x['rate']) for x in _rate]):.2f} "
               f"for the score of Fig. 2a)")

n_secretome = sum(x["class"] == "secretome" for x in read(FIG / "supp_table1.tsv"))

figs = "\n\n".join(f"![]({FIG}/FigS{i}.png){{width=6.3in}}\n\n" + legend(f"Fig. S{i}") for i in range(1, 9))

parts = [f"""# Supplementary Information

**{title}**

[Authors and affiliations to be completed before submission]

This file contains Supplementary Results 1 to 8, Supplementary Figures S1 to S8,
Supplementary Tables 1 to 10, Supplementary Notes 1 and 2 and the list of Supplementary Data
files.
Reference numbers are those of the main reference list. Every number here is reproduced by
the scripts named in the Code availability statement.

# Supplementary Results

{md_file("supplementary_results_en.md")}

# Supplementary Figures

{figs}

# Supplementary Tables

## Supplementary Table 1. Cohort sensitivity in the donor series

Every donor-level result under four definitions of the GSE113957 cohort. Age effects
are per decade of donor age; "lost" is the proportion of the unadjusted effect removed
by adjustment for proliferation.

{tsv(FIG / "supp_table4.tsv",
     cols=["cohort", "n", "rho_age_proliferation", "machinery_beta", "machinery_P",
           "machinery_beta_adj", "machinery_P_adj", "pct_of_age_effect_lost"],
     rename={"cohort": "cohort", "rho_age_proliferation": "rho age-prol",
             "machinery_beta": "effect", "machinery_P": "P",
             "machinery_beta_adj": "adj effect", "machinery_P_adj": "adj P",
             "pct_of_age_effect_lost": "% lost"},
     short={"cohort": SHORT_COHORT})}

Gene-level counts and the outcome index in the same four cohorts.

{tsv(FIG / "supp_table4.tsv",
     cols=["cohort", "unannotated_junction_beta", "unannotated_junction_P",
           "genes_age_associated", "genes_keeping_decline", "psi_age_events"],
     rename={"unannotated_junction_beta": "unann. effect", "unannotated_junction_P": "unann. P",
             "genes_age_associated": "genes, age", "genes_keeping_decline": "genes kept",
             "psi_age_events": "PSI events"},
     short={"cohort": SHORT_COHORT})}

## Supplementary Table 2. Treatment perturbation in the counted resource

The splicing score is the 177-gene score of Fig. 2. All of this table is post hoc.

**(a) Treatment arms.** Libraries and cell lines per arm. Arms with fewer than five
libraries or present in only one cell line were to be left out of the models; {arms_note}.

{tsv(D / "treatment_perturbation/treatment_counts.tsv",
     cols=["cond", "n", "lines", "modelled"],
     rename={"cond": "treatment", "n": "libraries", "lines": "cell lines"})}

**(b) Residual effect of each treatment**, with cell line, oxygen tension and days
grown in every model and Benjamini–Hochberg correction within each model. A treatment
is counted as keeping a residual effect when its coefficient is significant at
BH-FDR < 0.05 within the model, and as accounted for by the adjustment named in the
first column otherwise. Counts in the text and filled symbols in Supplementary Fig. S1c
use BH-FDR < 0.05.

{tsv(D / "treatment_perturbation/treatment_coefficients.tsv",
     cols=["response", "term", "beta", "lo", "hi", "p", "FDR"],
     rename={"response": "response and adjustment", "term": "treatment",
             "beta": "effect", "lo": "95% CI lower", "hi": "95% CI upper", "p": "P"},
     fmt={"beta": lambda v: num(v, 3), "lo": lambda v: num(v, 3),
          "hi": lambda v: num(v, 3), "p": lambda v: sig(v), "FDR": lambda v: sig(v)})}

**(c) Matched treated-minus-control differences**, median per arm, against control
cultures of the same cell line and oxygen tension nearest in days grown.

{tsv(D / "treatment_perturbation/matched_summary.tsv",
     cols=["cond", "n", "d_rate", "d_splice", "d_prolif", "d_days"],
     rename={"cond": "treatment", "n": "cultures", "d_rate": "delta division rate",
             "d_splice": "delta splicing", "d_prolif": "delta proliferation",
             "d_days": "median |delta days|"},
     fmt={"d_rate": lambda v: num(v, 3), "d_splice": lambda v: num(v, 2),
          "d_prolif": lambda v: num(v, 2), "d_days": lambda v: num(v, 1)})}

**(d) Oxygen tension.** Cultures at 3% oxygen against cultures of the same cell line at
21% oxygen: medians of the differences (3% minus 21%), with the number of cultures in
which the 3% value is the lower one. The unmatched rows reproduce a pooled comparison
and are given so that its composition can be seen.

{tsv(D / "treatment_perturbation/oxygen_contrast.tsv",
     cols=["contrast", "line", "n", "n_21", "d_rate", "d_splice", "d_prolif", "d_days",
           "n_rate_lower", "n_splice_lower"],
     rename={"n": "n at 3%", "n_21": "n at 21%", "d_rate": "delta rate",
             "d_splice": "delta splicing", "d_prolif": "delta proliferation",
             "d_days": "median |delta days|", "n_rate_lower": "rate lower",
             "n_splice_lower": "splicing lower"},
     fmt={"d_rate": lambda v: num(v, 3), "d_splice": lambda v: num(v, 2),
          "d_prolif": lambda v: num(v, 2), "d_days": lambda v: num(v, 1)})}

The oxygen term in models of control cultures, with cell line and days grown as covariates.

{tsv(D / "treatment_perturbation/oxygen_models.tsv",
     cols=["sample_set", "response", "held_constant", "n", "n_3", "beta", "lo", "hi", "p"],
     rename={"sample_set": "libraries", "held_constant": "held constant", "n_3": "n at 3%",
             "beta": "effect of 3% oxygen", "lo": "95% CI lower", "hi": "95% CI upper", "p": "P"},
     fmt={"beta": lambda v: num(v, 3), "lo": lambda v: num(v, 3), "hi": lambda v: num(v, 3),
          "p": lambda v: sig(v)})}

## Supplementary Table 3. Sensitivity analyses in the donor cohort

**(a) Technical quality covariates.** Uniquely mapped percentage, multi-mapping
percentage, intronic signal and exon assignment rate added to the cohort model.

{tsv(D / "reviewer_sensitivities/qc_covariates.tsv",
     short={'metric': SHORT_METRIC},
     fmt={"beta": lambda v: f"{v:.3g}", "p": lambda v: sig(v),
          "beta_adj": lambda v: f"{v:.3g}", "p_adj": lambda v: sig(v),
          "pct_lost": lambda v: f"{v:.0f}"})}

**(b) Core spliceosome versus the rest of the set.** The 127 genes belonging to
Reactome mRNA splicing, and the 50 genes of capping, 3'-end processing and export.

{tsv(D / "reviewer_sensitivities/core_spliceosome_subset.tsv",
     short={'metric': SHORT_METRIC},
     fmt={"beta": lambda v: num(v, 3), "p": lambda v: sig(v),
          "beta_adj": lambda v: num(v, 3), "p_adj": lambda v: sig(v),
          "pct_lost": lambda v: f"{v:.0f}"})}

The same three scores against the measured division rate in the 328 counted libraries
(Spearman correlation with its 95% confidence interval), all scored in the expression
matrix of the programme audit.

{tsv(D / "reviewer_sensitivities/subset_division_rate.tsv",
     cols=["gene_set", "genes_measured", "n_libraries", "rho_rate", "lo", "hi", "p"],
     rename={"gene_set": "gene set", "genes_measured": "genes measured",
             "n_libraries": "libraries", "rho_rate": "rho with division rate",
             "lo": "95% CI lower", "hi": "95% CI upper", "p": "P"},
     fmt={"rho_rate": lambda v: num(v, 2), "lo": lambda v: num(v, 2),
          "hi": lambda v: num(v, 2), "p": lambda v: sig(v)})}

**(c) Base coverage versus read equivalents.** Gene sums divided by read length
(51 or 75 nucleotides in this series).

{tsv(D / "reviewer_sensitivities/read_length_conversion.tsv",
     short={'metric': SHORT_METRIC},
     fmt={"beta": lambda v: num(v, 3), "p": lambda v: sig(v),
          "beta_adj": lambda v: num(v, 3), "p_adj": lambda v: sig(v),
          "pct_lost": lambda v: f"{v:.0f}"})}

**(d) Transcriptome-wide age predictor (post hoc).** Ten-fold cross-validation repeated
five times in the 107 donors; the permuted-age null is twenty repetitions.

{tsv(D / "age_predictor/age_predictor_cv.tsv",
     cols=["model", "r", "lo", "hi", "MAE", "r_with_proliferation"],
     rename={"lo": "95% CI lower", "hi": "95% CI upper",
             "r_with_proliferation": "r with proliferation score"},
     fmt={"r": lambda v: num(v, 2), "lo": lambda v: num(v, 2), "hi": lambda v: num(v, 2),
          "MAE": lambda v: num(v, 1), "r_with_proliferation": lambda v: num(v, 2)})}

**(e) Negative control for the proliferation adjustment.** The unadjusted and the
proliferation-adjusted cohort models applied to 1,000 random gene sets matched to the
splicing set for expression level and containing no splicing or cell-cycle genes. The
share removed is given for the random sets that have an age association of their own,
because it is undefined without one; effects are in the units of the splicing score.

{negative_control(D / "reviewer_sensitivities/adjustment_negative_control.tsv",
                  D / "reviewer_sensitivities/core_spliceosome_subset.tsv")}

**(f) Which technical covariate the residual depends on.** The proliferation-adjusted
age effect of the 177-gene score in the 107 donors, per decade with its 95% confidence
interval, under subsets of the cohort's technical covariates; the share removed is
relative to the unadjusted model that carries the same covariates. The two technical
terms overlap: {hiseq_note(D / "cohort_revised/primary/sample_metrics.tsv")}.

{tsv(D / "reviewer_sensitivities/covariate_attribution.tsv",
     cols=["covariates", "beta_unadj", "p_unadj", "beta_adj", "ci", "p", "pct_removed"],
     rename={"beta_unadj": "unadjusted", "p_unadj": "P", "beta_adj": "adjusted",
             "ci": "95% CI (adjusted)", "p": "P (adj)", "pct_removed": "% removed"},
     derive={"ci": lambda x: span(x["lo"], x["hi"])},
     fmt={"beta_unadj": lambda v: num(v, 3), "p_unadj": lambda v: sig(v),
          "beta_adj": lambda v: num(v, 3), "p": lambda v: sig(v),
          "pct_removed": lambda v: f"{v:.0f}"})}

## Supplementary Table 4. Datasets and contrasts

All 63 contrasts, with the class assigned by experimental manipulation, the change in
the cell-cycle and splicing sets, and the source publication of each series. The
secretome class holds {n_secretome} contrasts.

{tsv(FIG / "supp_table1.tsv",
     cols=["accession", "class", "contrast", "cellcycle_change", "splicing_change"],
     rename={"cellcycle_change": "cell cycle", "splicing_change": "splicing"})}

Source publication of each series.

{sources(FIG / "supp_table1.tsv")}

**Intervention classes predicted out of sample (post hoc).** Median absolute residual
against median absolute change, from the in-sample fit of all 63 contrasts and from fits
that exclude the class being predicted.

{tsv(D / "interventions_revised/prediction_accuracy.tsv",
     cols=["scheme", "group", "n_contrasts", "n_series", "median_abs_change",
           "median_abs_residual", "share_predicted", "n_within_0.05", "n_beyond_0.10"],
     rename={"n_contrasts": "contrasts", "n_series": "series", "median_abs_change": "median |change|",
             "median_abs_residual": "median |residual|", "share_predicted": "share predicted",
             "n_within_0.05": "within 0.05", "n_beyond_0.10": "beyond 0.10"},
     fmt={"median_abs_change": lambda v: num(v, 3), "median_abs_residual": lambda v: num(v, 3),
          "share_predicted": lambda v: num(v, 2)})}

Mean class residual with its 95% confidence interval and the shift detectable with 80% power.

{tsv(D / "interventions_revised/class_residuals_in_and_out_of_sample.tsv",
     cols=["scheme", "class", "unit", "n", "mean_residual", "lo", "hi",
           "detectable_shift_80pct_power", "class_median_abs_change"],
     rename={"mean_residual": "mean residual", "lo": "95% CI lower", "hi": "95% CI upper",
             "detectable_shift_80pct_power": "detectable shift",
             "class_median_abs_change": "class median |change|"},
     fmt={"mean_residual": lambda v: num(v, 3), "lo": lambda v: num(v, 3), "hi": lambda v: num(v, 3),
          "detectable_shift_80pct_power": lambda v: num(v, 3),
          "class_median_abs_change": lambda v: num(v, 3)})}

## Supplementary Table 5. Specificity of the coupling, and the specification range of the residual age effect

All of this table is post hoc.

**(a) The counted libraries.** Correlation with the measured division rate across the
{_sel[0]["n"]} counted libraries for each splicing definition and each comparator, the
partial correlation of each when the other is held constant, their difference, and the
Pearson correlation between the two scores.
Every set, the 177-gene set included, is scored in the expression matrix of the programme
audit so that all sets pass one filter; {matrix_note}. No comparison shares a gene: shared
genes were removed from both sides for the unselected definitions and from the comparator
only for the selected set, and where a definition is nested in Metabolism of RNA its genes
were removed from that comparator.

{tsv(D / "machinery_specificity/counted_libraries_partials.tsv",
     cols=["splicing_definition", "comparator", "n_genes_splicing", "n_genes_comparator",
           "rho_splicing", "rho_comparator", "partial_splicing", "partial_comparator",
           "asymmetry", "pct_of_selected_asymmetry", "r_between"],
     rename={"splicing_definition": "splicing definition", "n_genes_splicing": "genes, splicing",
             "n_genes_comparator": "genes, comparator", "rho_splicing": "rho splicing",
             "rho_comparator": "rho comparator", "partial_splicing": "splicing | comparator",
             "partial_comparator": "comparator | splicing",
             "pct_of_selected_asymmetry": "% of selected asymmetry",
             "r_between": "r between scores"},
     fmt={"rho_splicing": lambda v: num(v, 2), "rho_comparator": lambda v: num(v, 2),
          "partial_splicing": lambda v: num(v, 2), "partial_comparator": lambda v: num(v, 2),
          "asymmetry": lambda v: num(v, 2), "pct_of_selected_asymmetry": lambda v: f"{v:.0f}",
          "r_between": lambda v: num(v, 2)})}

The comparisons affected by the choice of rule, repeated under the other rule: the
selected set with shared genes removed from both sides, and the other definitions with
shared genes removed from the comparator only.

{tsv(D / "machinery_specificity/shared_rule_sensitivity.tsv",
     cols=["splicing_definition", "comparator", "shared_rule", "n_shared_removed",
           "n_genes_splicing", "n_genes_comparator", "partial_splicing", "partial_comparator",
           "asymmetry"],
     rename={"splicing_definition": "splicing definition", "shared_rule": "shared genes removed from",
             "n_shared_removed": "shared genes", "n_genes_splicing": "genes, splicing",
             "n_genes_comparator": "genes, comparator", "partial_splicing": "splicing | comparator",
             "partial_comparator": "comparator | splicing"},
     short={"shared_rule": SHORT_RULE},
     fmt={"partial_splicing": lambda v: num(v, 3), "partial_comparator": lambda v: num(v, 3),
          "asymmetry": lambda v: num(v, 3)})}

**(b) The donor cohort.** Age effect per decade in the 107 adult donors under each model.

{tsv(D / "machinery_specificity/donor_age_models.tsv",
     cols=["splicing_definition", "comparator", "side", "metric", "model", "n_genes", "beta", "lo", "hi", "p"],
     rename={"splicing_definition": "splicing definition", "n_genes": "genes",
             "beta": "per decade", "lo": "95% CI lower", "hi": "95% CI upper", "p": "P"},
     fmt={"beta": lambda v: num(v, 3), "lo": lambda v: num(v, 3),
          "hi": lambda v: num(v, 3), "p": lambda v: sig(v)})}

**(c) GTEx cultures.** The same comparison against the 20-marker transcriptional
proliferation score in 511 cultures. The asymmetry is not reproduced.

{tsv(D / "machinery_specificity/gtex_culture_partials.tsv",
     cols=["splicing_definition", "comparator", "n", "rho_splicing", "rho_comparator",
           "partial_splicing", "partial_comparator", "asymmetry"],
     rename={"splicing_definition": "splicing definition", "rho_splicing": "rho splicing",
             "rho_comparator": "rho comparator", "partial_splicing": "splicing | comparator",
             "partial_comparator": "comparator | splicing"},
     fmt={"rho_splicing": lambda v: num(v, 2), "rho_comparator": lambda v: num(v, 2),
          "partial_splicing": lambda v: num(v, 2), "partial_comparator": lambda v: num(v, 2),
          "asymmetry": lambda v: num(v, 2)})}

**(d) The two cohorts, not pooled.** Proliferation-adjusted age effects of the two
fibroblast cohorts, reported side by side; they are not combined, because one is an
attenuation and the other a suppression effect in cultures with no unadjusted age effect.

{tsv(D / "residual_meta/residual_age_effect_meta.tsv",
     rows=lambda x: "combined" not in x["cohort"],
     cols=["gene_set", "cohort", "n", "beta", "lo", "hi", "p"],
     rename={"gene_set": "gene set", "beta": "per decade", "lo": "95% CI lower",
             "hi": "95% CI upper", "p": "P"},
     fmt={"beta": lambda v: num(v, 4), "lo": lambda v: num(v, 4), "hi": lambda v: num(v, 4),
          "p": lambda v: sig(v)})}

**(e) The four definitions of the fibroblast cohort**, reported as sensitivity and not
pooled, since they are not independent of one another.

{tsv(D / "residual_meta/cohort_definition_sensitivity.tsv",
     cols=["metric", "cohort", "n", "beta_adj", "lo_adj", "hi_adj", "p_adj"],
     rename={"beta_adj": "adjusted per decade", "lo_adj": "95% CI lower",
             "hi_adj": "95% CI upper", "p_adj": "P"},
     fmt={"beta_adj": lambda v: num(v, 3), "lo_adj": lambda v: num(v, 3),
          "hi_adj": lambda v: num(v, 3), "p_adj": lambda v: sig(v)})}

## Supplementary Table 6. Ageing programmes against measured division rate and donor age

First: correlation with the measured division rate across 328 counted libraries, with
all genes and after removing every gene belonging to the Reactome cell-cycle union.
Second: the same programmes against donor age in the 107-donor primary cohort, before
and after adjustment for proliferation (adj), each with its nominal P value and BH-FDR.

{tsv(D / "hallmark_revised/division_rate_cellcycle_removed.tsv",
     cols=["programme", "n_genes", "n_cc_removed", "rho_rate", "rho_rate_noCC", "FDR_noCC"],
     rename={"n_genes": "genes", "n_cc_removed": "cell-cycle genes removed",
             "rho_rate": "rho, all genes", "rho_rate_noCC": "rho, after removal",
             "FDR_noCC": "FDR"},
     fmt={"rho_rate": lambda v: num(v, 2), "rho_rate_noCC": lambda v: num(v, 2),
          "FDR_noCC": lambda v: sig(v)})}

{tsv(D / "hallmark_revised/donor_age_primary_cohort.tsv",
     cols=["programme", "rho_age", "beta_decade", "p_age", "FDR_age", "beta_adj", "p_adj", "FDR_adj"],
     rename={"rho_age": "rho with age", "beta_decade": "per decade", "p_age": "P",
             "FDR_age": "FDR", "beta_adj": "per decade (adj)", "p_adj": "P (adj)",
             "FDR_adj": "FDR (adj)"},
     fmt={"rho_age": lambda v: num(v, 2), "beta_decade": lambda v: num(v, 3),
          "p_age": lambda v: sig(v), "FDR_age": lambda v: sig(v),
          "beta_adj": lambda v: num(v, 3), "p_adj": lambda v: sig(v),
          "FDR_adj": lambda v: sig(v)})}

## Supplementary Table 7. Senescence panels in cultured fibroblasts

Post hoc. Cell-cycle content is the share of a panel's measured genes belonging to the
Reactome cell-cycle union. Correlations are with the measured division rate across 328
counted libraries; age effects are per decade in the 107 adult donors. The percentage of
the age effect lost is meaningful only for the five panels that had one.

{tsv(D / "senescence_panels/senescence_panel_audit.tsv",
     cols=["panel", "genes_measured", "pct_cellcycle", "rho_rate", "rho_rate_noCC",
           "median_gene_rho", "beta", "FDR", "beta_adj", "FDR_adj", "pct_lost"],
     rename={"genes_measured": "measured", "pct_cellcycle": "% cell cycle",
             "rho_rate": "rho rate", "rho_rate_noCC": "rho after removal",
             "median_gene_rho": "median gene rho", "beta": "age per decade",
             "beta_adj": "adjusted", "pct_lost": "% of age effect lost"},
     fmt={"pct_cellcycle": lambda v: f"{v:.0f}", "rho_rate": lambda v: num(v, 2),
          "rho_rate_noCC": lambda v: num(v, 2), "median_gene_rho": lambda v: num(v, 2),
          "beta": lambda v: num(v, 3), "FDR": lambda v: sig(v),
          "beta_adj": lambda v: num(v, 3), "FDR_adj": lambda v: sig(v),
          "pct_lost": lambda v: f"{v:.0f}"})}

## Supplementary Table 8. Preregistered methylation clock tests

The four hypotheses fixed in Supplementary Note 2, each with the criterion it was given.
M1 and M2 are Spearman correlations across methylation samples; M3 and M4 are across the
109 cultures carrying both a clock and transcriptional scores. "Met" is the preregistered
criterion, which for M1 to M3 is |rho| >= 0.3 with a significant P value and for M4 is an
attenuation of at least 50% when the measured division rate is held constant.

{tsv(D / "methylation_clock/preregistered_tests.tsv",
     cols=["test", "rho", "p", "FDR", "n", "supported"],
     rename={"rho": "rho", "p": "P", "supported": "met"},
     fmt={"rho": lambda v: num(v, 3), "p": lambda v: sig(v), "FDR": lambda v: sig(v)})}

## Supplementary Table 9. The preregistered GTEx hypotheses and their outcomes

The six hypotheses as they were locked before the data were downloaded, with the
criterion fixed in advance, the value observed and the verdict. Four were supported; H4
could not be evaluated because the GTEx cultures had no unadjusted age effect to
attenuate; and H6 — a positive control asking whether collagen formation falls with age
in sun-exposed skin — was not supported, which under the preregistered rule bars age
effects measured within GTEx tissues from supporting any claim in this paper. The full
document, its English translation and its ten recorded deviations are Supplementary
Note 1.

{tsv(D / "gtex_boundary/preregistered_scorecard.tsv",
     cols=["hypothesis", "test", "criterion", "observed", "verdict"])}

## Supplementary Table 10. Programmes against proliferation and donor age in GTEx

**(a) Partial correlation with proliferation.** Partial Spearman correlation of each of
the twenty-five programme scores with the proliferation score in each tissue, after
regressing out RNA integrity, ischaemic time and library size (the values of Fig. 6b).

{pivot(D / "gtex_boundary/programme_by_tissue.tsv")}

**(b) Age effects before and after adjustment for proliferation.** Change in each
programme score per decade of donor age, from linear models with sex, Hardy death
classification, RNA integrity, ischaemic time and library size as covariates, and the
same with the proliferation score added (adj). FDR is Benjamini–Hochberg within a
tissue and was computed for the unadjusted effects. The preregistered positive control
for age effects in tissue failed (H6), so under the preregistered rule the age effects
in skin and muscle are listed for completeness and support no claim.

{gtex_age(D / "gtex_boundary/programme_by_tissue.tsv", "culture")}

{gtex_age(D / "gtex_boundary/programme_by_tissue.tsv", "legskin")}

{gtex_age(D / "gtex_boundary/programme_by_tissue.tsv", "pubskin")}

{gtex_age(D / "gtex_boundary/programme_by_tissue.tsv", "muscle")}

**(c) Range restriction (post hoc).** Spread of the proliferation score, the partial
correlation of the splicing score with it, and the regression slope of the splicing score
on it with its 95% confidence interval, which a narrower spread does not affect; for the
preregistered 96-gene set and for the 177-gene set.

{tsv(D / "gtex_posthoc/range_restriction.tsv",
     cols=["tissue", "set", "sd_prolif", "partial_rho", "p", "slope", "ci"],
     rename={"tissue": "tissue (donors)", "set": "genes in set", "sd_prolif": "SD of proliferation",
             "partial_rho": "partial rho", "p": "P", "ci": "95% CI of slope"},
     short={"set": SHORT_SET},
     derive={"tissue": lambda x: TISSUE_NAME[x["tissue"]] + " (" + x["n"] + ")",
             "ci": lambda x: span(x["slope_lo"], x["slope_hi"])},
     fmt={"sd_prolif": lambda v: num(v, 3), "partial_rho": lambda v: num(v, 3),
          "p": lambda v: sig(v), "slope": lambda v: num(v, 3)})}

The partial correlation expected in each tissue if the narrower spread of proliferation
were the only difference from culture, against the one observed (96-gene set).

{tsv(D / "gtex_posthoc/range_restriction_prediction.tsv",
     cols=["tissue", "sd_ratio", "predicted_rho_if_only_range", "observed_rho"],
     rename={"sd_ratio": "SD relative to culture",
             "predicted_rho_if_only_range": "rho expected from range restriction alone",
             "observed_rho": "rho observed"},
     short={"tissue": TISSUE_NAME},
     fmt={"sd_ratio": lambda v: num(v, 3), "predicted_rho_if_only_range": lambda v: num(v, 3),
          "observed_rho": lambda v: num(v, 3)})}

Cultured fibroblasts subsampled to the spread of proliferation in sun-exposed skin: a
window centred on the median, and kernel weights that down-weight the tails without
discarding them; the median partial correlation over the draws with its 2.5th and 97.5th
percentiles.

{tsv(D / "gtex_posthoc/variance_matched_summary.tsv",
     cols=["scheme", "target_sd", "donors_available", "draws", "achieved_sd_median",
           "rho_median", "range"],
     rename={"target_sd": "target SD", "donors_available": "cultures available",
             "achieved_sd_median": "SD achieved (median)", "rho_median": "partial rho (median)",
             "range": "2.5th to 97.5th percentile"},
     derive={"range": lambda x: span(x["rho_lo"], x["rho_hi"])},
     fmt={"target_sd": lambda v: num(v, 3), "achieved_sd_median": lambda v: num(v, 3),
          "rho_median": lambda v: num(v, 3)})}

**(d) Adjustment for cell composition (post hoc).** The splicing–proliferation partial
correlation with the technical covariates only, and with the marker scores named in the
third column added to them. Block (f) repeats the adjustment across programmes, which is
what separates a composition effect on splicing from a composition effect on the whole
proliferation axis.

{tsv(D / "gtex_posthoc/composition_adjusted.tsv",
     cols=["tissue", "set", "markers", "rho_tech", "p_tech", "rho_tech_comp", "p_comp"],
     rename={"tissue": "tissue (donors)", "set": "genes in set", "markers": "marker scores added",
             "rho_tech": "partial rho, technical", "p_tech": "P",
             "rho_tech_comp": "partial rho, with composition", "p_comp": "P (composition)"},
     short={"set": SHORT_SET, "markers": SHORT_MARKERS},
     derive={"tissue": lambda x: TISSUE_NAME[x["tissue"]] + " (" + x["n_comp"] + ")"},
     fmt={"rho_tech": lambda v: num(v, 3), "p_tech": lambda v: sig(v),
          "rho_tech_comp": lambda v: num(v, 3), "p_comp": lambda v: sig(v)})}

**(e) Positive control for donor concordance (post hoc).** Spearman correlation between a
donor's cultured fibroblasts and the same donor's sun-exposed skin for genes whose
expression depends strongly on genotype or on sex. The column before the note is the
product of the gene's correlations with sex in the two sample types, that is, the
concordance that sex alone would produce. For comparison, the preregistered exploratory
analysis gave {concordance_range(D / "gtex_boundary/donor_concordance_culture_vs_legskin.tsv")}
(Supplementary Fig. S8c).

{tsv(D / "gtex_posthoc/donor_identity_positive_control.tsv",
     cols=["gene", "group", "n", "rho", "p", "ceiling_from_sex", "note"],
     rename={"n": "donors", "p": "P", "ceiling_from_sex": "rho expected from sex alone"},
     short={"group": SHORT_GROUP},
     fmt={"rho": lambda v: num(v, 3), "p": lambda v: sig(v),
          "ceiling_from_sex": lambda v: num(v, 3)})}

**(f) The same composition adjustment for six sets (post hoc).** Run on the splicing
sets alone, block (d) cannot distinguish an adjustment that removes splicing's coupling
from one that flattens the proliferation signal. Repeated for the mitotic cell-cycle
positive control, collagen formation and the two unselected Reactome splicing definitions,
the two separate: the cell-cycle programme and collagen formation do not move, while every
definition of the splicing set falls (the values of Fig. 6d). The fibroblast marker score
contains COL1A1 and COL1A2, so collagen formation is not an independent check.

{tsv(D / "gtex_posthoc/composition_adjusted_all_programmes.tsv",
     cols=["tissue", "set", "rho_tech", "rho_tech_comp", "drop_abs"],
     rename={"tissue": "tissue", "set": "programme", "rho_tech": "partial rho, technical",
             "rho_tech_comp": "partial rho, with composition", "drop_abs": "difference"},
     short={"set": SHORT_SET},
     derive={"tissue": lambda x: TISSUE_NAME[x["tissue"]]},
     fmt={"rho_tech": lambda v: num(v, 3), "rho_tech_comp": lambda v: num(v, 3),
          "drop_abs": lambda v: num(v, 3)})}

"""]

for i, (f, ttl) in enumerate((("supplementary_note1_preregistration_en.md", "Preregistration for the GTEx analysis"),
                              ("supplementary_note2_preregistration_methylation_en.md", "Preregistration for the methylation clock analysis")), 1):
    parts.append(f"# Supplementary Note {i}. {ttl}\n\n" + md_file(f) + "\n\n\n\n")

parts.append("""# Supplementary Data

Machine-readable tables accompany this file.

| File | Deposited as | Contents |
|---|---|---|
| Supplementary Data 1 | `SupplementaryData1_contrasts.tsv` | The 63 contrasts: change in the 177-gene splicing set, in the cell-cycle set and in the splicing set without the 23 shared genes, with the class of each contrast (`leave_out_defining_series.tsv` in the same folder gives the refit without the three set-defining series) |
| Supplementary Data 2 | `SupplementaryData2_gtex_programmes_by_tissue.tsv` | Programme scores against proliferation and donor age in every GTEx tissue, with age effects before and after adjustment |
| Supplementary Data 3 | `SupplementaryData3_gene_level_by_cohort.tsv` | For each of the four cohort definitions: the proliferation–age correlation, the machinery–proliferation correlation, and the number of genes with an age association before and after adjustment for proliferation |
| Supplementary Data 4 | `SupplementaryData4_programme_division_rate.tsv` | Division-rate correlations for every programme, with and without cell-cycle genes |
| Supplementary Data 5 | `SupplementaryData5_residual_specification.tsv` | The residual age effect under every model specification — splicing definition, comparator programme and covariate set (Supplementary Table 5) |
| Supplementary Data 6 | `SupplementaryData6_composition_adjusted.tsv` | Partial correlation with proliferation before and after cell-composition marker scores, for six sets in four tissues (Fig. 6d, Supplementary Fig. S8b) |
| Supplementary Data 7 | `SupplementaryData7_library_scores_vs_counted_rate.tsv` | Counted division rate, 20-marker proliferation score and the 177- and 96-gene splicing scores for the 328 counted libraries (Figs 2, 5a and 7b) |
| Supplementary Data 8 | `SupplementaryData8_intervention_prediction.tsv` | Out-of-sample prediction of the intervention classes: contrasts, series, median change and residual under each prediction scheme (Supplementary Fig. S4, Supplementary Results 5) |
| Supplementary Data 9 | `SupplementaryData9_senescence_panels.tsv` | The seven senescence panels: gene counts, cell-cycle content, correlation with the counted division rate with and without cell-cycle genes, and the age effect before and after adjustment (Supplementary Fig. S6, Supplementary Table 7) |
| Supplementary Data 10 | `SupplementaryData10_methylation_tests.tsv` | The four preregistered methylation clock tests with their correlations, *P* values, BH-FDR and verdicts (Supplementary Table 8) |
""")

(OUT / "supplementary.md").write_text("\n".join(parts))
print("supplementary.md written:", len("\n".join(parts)), "characters")
