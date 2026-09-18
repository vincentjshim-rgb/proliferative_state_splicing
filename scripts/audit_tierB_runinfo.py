#!/usr/bin/env python3
"""Join ENA run metadata to GEO sample metadata and summarize Tier-B raw size.

This script never downloads sequence reads.  It consumes only the small ENA
filereport TSV snapshots stored under public_data_tierA/metadata/tierB_runinfo.
"""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path


ROOT = Path("public_data_tierA")
RUNINFO = ROOT / "metadata" / "tierB_runinfo"
AUDIT = ROOT / "derived" / "audit"

PROJECTS = {
    "GSE191055": "PRJNA789540",
    "GSE109700": "PRJNA431791",
    "GSE302943": "PRJNA1292871",
    "GSE240226": "PRJNA1003002",
    "GSE165177": "PRJNA693500",
    "GSE179848": "PRJNA745334",
    "GSE226189": "PRJNA939148",
    "GSE93535": "PRJNA361062",
}


def parse_semicolon_int(value: str) -> int:
    return sum(int(part) for part in value.split(";") if part.strip())


with (AUDIT / "all_geo_samples.tsv").open() as handle:
    geo_rows = list(csv.DictReader(handle, delimiter="\t"))

sample_lookup: dict[tuple[str, str], dict[str, str]] = {}
experiment_lookup: dict[tuple[str, str], dict[str, str]] = {}
for row in geo_rows:
    relation = row.get("relation", "")
    biosample = re.search(r"BioSample:\s*https?://[^| ]+/([^ |]+)", relation)
    experiment = re.search(r"SRA:\s*https?://[^| ]+[?&]term=([^ |]+)", relation)
    if biosample:
        sample_lookup[(row["dataset"], biosample.group(1))] = row
    if experiment:
        experiment_lookup[(row["dataset"], experiment.group(1))] = row

joined: list[dict[str, object]] = []
summary: list[dict[str, object]] = []
for dataset, project in PROJECTS.items():
    path = RUNINFO / f"{project}.tsv"
    with path.open() as handle:
        run_rows = list(csv.DictReader(handle, delimiter="\t"))
    total_bytes = 0
    mapped = 0
    layouts = set()
    for run in run_rows:
        fastq_bytes = parse_semicolon_int(run.get("fastq_bytes", ""))
        total_bytes += fastq_bytes
        layouts.add(run.get("library_layout", ""))
        geo = sample_lookup.get((dataset, run.get("sample_accession", "")))
        if geo is None:
            geo = experiment_lookup.get((dataset, run.get("experiment_accession", "")))
        if geo is not None:
            mapped += 1
        joined.append({
            "dataset": dataset,
            "bioproject": project,
            "geo_accession": geo.get("accession", "") if geo else "",
            "geo_title": geo.get("title", "") if geo else "",
            "sample_accession": run.get("sample_accession", ""),
            "experiment_accession": run.get("experiment_accession", ""),
            "run_accession": run.get("run_accession", ""),
            "library_layout": run.get("library_layout", ""),
            "read_count": run.get("read_count", ""),
            "base_count": run.get("base_count", ""),
            "fastq_bytes": fastq_bytes,
            "fastq_GiB": round(fastq_bytes / 1024**3, 4),
        })
    summary.append({
        "dataset": dataset,
        "bioproject": project,
        "runs": len(run_rows),
        "mapped_to_geo_samples": mapped,
        "layout": ";".join(sorted(layouts)),
        "fastq_bytes": total_bytes,
        "fastq_GiB": round(total_bytes / 1024**3, 2),
    })

AUDIT.mkdir(parents=True, exist_ok=True)
for filename, rows in [
    ("tierB_run_level_inventory.tsv", joined),
    ("tierB_project_size_summary.tsv", summary),
]:
    with (AUDIT / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_subset(name: str, dataset: str, selected_rows: list[dict[str, object]], rationale: str) -> dict[str, object]:
    return {
        "priority_subset": name,
        "dataset": dataset,
        "runs": len(selected_rows),
        "unique_geo_samples": len({r["geo_accession"] for r in selected_rows}),
        "layout": ";".join(sorted({str(r["library_layout"]) for r in selected_rows})),
        "fastq_bytes": sum(int(r["fastq_bytes"]) for r in selected_rows),
        "fastq_GiB": round(sum(int(r["fastq_bytes"]) for r in selected_rows) / 1024**3, 2),
        "recommended_working_disk_GiB_3x": round(sum(int(r["fastq_bytes"]) for r in selected_rows) / 1024**3 * 3, 1),
        "rationale": rationale,
    }


by_dataset = {dataset: [r for r in joined if r["dataset"] == dataset] for dataset in PROJECTS}
subsets: list[dict[str, object]] = []
subsets.append(add_subset(
    "S1_P4_P27_only", "GSE191055",
    [r for r in by_dataset["GSE191055"] if str(r["geo_title"]).startswith(("P4_", "P27_"))],
    "Eight HDF passage samples; exclude old/young skin samples that are biological n=1 each",
))
subsets.append(add_subset(
    "S2_UVA_heldout_all", "GSE302943", by_dataset["GSE302943"],
    "Ten single-end BJ-5ta samples; moderate-size held-out UVA splicing test; lock filename/paper labels",
))
subsets.append(add_subset(
    "S3_SIPS_all", "GSE93535", by_dataset["GSE93535"],
    "All 12 GEO samples represented by 13 runs; SIPS and rescue exact-event validation",
))
subsets.append(add_subset(
    "S4_UVA_rescue_all", "GSE240226", by_dataset["GSE240226"],
    "All 12 samples; UVA injury plus two rescue contrasts",
))
subsets.append(add_subset(
    "S5_senescence_all", "GSE109700", by_dataset["GSE109700"],
    "All nine paired-end samples; deep sequencing makes this expensive despite small n",
))
subsets.append(add_subset(
    "S6_MPTR_day13", "GSE165177",
    [r for r in by_dataset["GSE165177"]
     if re.match(r"^O[23]_(transiently_reprogrammed|negative_control)_13days_exp[12]$", str(r["geo_title"]))],
    "Published principal day-13 window; four matched pairs but only O2/O3 donors",
))

longitudinal_candidates = [
    r for r in by_dataset["GSE179848"]
    if re.match(r"^HC[1-4]_sp[23]_Normal_Control_ox21_p\d+$", str(r["geo_title"]))
]
longitudinal_selected: list[dict[str, object]] = []
for donor in ("HC1", "HC2", "HC3", "HC4"):
    donor_rows = [r for r in longitudinal_candidates if str(r["geo_title"]).startswith(donor + "_")]
    donor_rows.sort(key=lambda r: int(re.search(r"_p(\d+)$", str(r["geo_title"])).group(1)))
    longitudinal_selected.extend(donor_rows[:2] + donor_rows[-2:])
subsets.append(add_subset(
    "S7_longitudinal_first2_last2", "GSE179848", longitudinal_selected,
    "Sixteen normal-control samples preserving four donor-level early/late contrasts; not the 345-sample project",
))

age_sorted = sorted(by_dataset["GSE226189"], key=lambda r: int(re.search(r"AGE(\d+)_", str(r["geo_title"])).group(1)))
age_extremes = age_sorted[:10] + age_sorted[-10:]
subsets.append(add_subset(
    "S8_age_extremes_10plus10", "GSE226189", age_extremes,
    "Twenty donor extremes as a cost-reduced exact-event sensitivity set; processed counts remain primary for age regression",
))

with (AUDIT / "tierB_priority_subset_size_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(subsets[0]))
    writer.writeheader()
    writer.writerows(subsets)

with (AUDIT / "tierB_runinfo_sha256.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["file", "bytes", "sha256"])
    writer.writeheader()
    for path in sorted(RUNINFO.glob("*.tsv")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        writer.writerow({"file": str(path), "bytes": path.stat().st_size, "sha256": digest})

print(f"Wrote {len(joined)} run records and {len(summary)} project summaries")
