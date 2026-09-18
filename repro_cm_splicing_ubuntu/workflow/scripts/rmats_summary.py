#!/usr/bin/env python3
"""Create ΔPSI candidate tables from rMATS output generated with --statoff."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable


EVENT_TYPES = ("SE", "A5SS", "A3SS", "MXE", "RI")
COORDINATE_COLUMNS = {
    "SE": (
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),
    "A5SS": (
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "flankingES",
        "flankingEE",
    ),
    "A3SS": (
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "flankingES",
        "flankingEE",
    ),
    "MXE": (
        "1stExonStart_0base",
        "1stExonEnd",
        "2ndExonStart_0base",
        "2ndExonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),
    "RI": (
        "riExonStart_0base",
        "riExonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def parse_counts(value: str) -> list[int]:
    if value is None or value == "":
        return []
    return [int(float(item)) for item in value.split(",")]


def parse_levels(value: str) -> list[float]:
    if value is None or value == "":
        return []
    levels: list[float] = []
    for item in value.split(","):
        if item not in {"NA", "NaN", "nan", ""}:
            levels.append(float(item))
    return levels


def calculate_levels(inclusion: list[int], skipping: list[int], inc_len: float, skip_len: float):
    if len(inclusion) != len(skipping):
        raise ValueError("inclusion and skipping count vectors have different lengths")
    levels: list[float] = []
    for inc_count, skip_count in zip(inclusion, skipping):
        inc_norm = inc_count / inc_len if inc_len > 0 else 0.0
        skip_norm = skip_count / skip_len if skip_len > 0 else 0.0
        denominator = inc_norm + skip_norm
        levels.append(inc_norm / denominator if denominator > 0 else math.nan)
    return levels


def finite_mean(values: Iterable[float]) -> float:
    finite = [value for value in values if not math.isnan(value)]
    return sum(finite) / len(finite) if finite else math.nan


def format_float(value: float) -> str:
    return "NA" if math.isnan(value) else f"{value:.6f}"


def event_key(event_type: str, row: dict[str, str]) -> str:
    components = [event_type, row.get("chr", ""), row.get("strand", "")]
    components.extend(row.get(column, "") for column in COORDINATE_COLUMNS[event_type])
    return ":".join(components)


def load_event_rows(rmats_dir: Path, event_type: str, count_mode: str):
    final_path = rmats_dir / f"{event_type}.MATS.{count_mode}.txt"
    if final_path.is_file():
        return read_tsv(final_path)

    definition_path = rmats_dir / f"fromGTF.{event_type}.txt"
    raw_path = rmats_dir / f"{count_mode}.raw.input.{event_type}.txt"
    if not definition_path.is_file() or not raw_path.is_file():
        raise FileNotFoundError(
            f"missing rMATS files for {event_type}/{count_mode}: "
            f"{definition_path} ; {raw_path}"
        )
    definition_fields, definitions = read_tsv(definition_path)
    raw_fields, raw_rows = read_tsv(raw_path)
    by_id = {row["ID"]: row for row in definitions}
    merged: list[dict[str, str]] = []
    for raw_row in raw_rows:
        definition = by_id.get(raw_row["ID"])
        if definition is None:
            raise ValueError(f"rMATS raw event ID not found in definitions: {raw_row['ID']}")
        merged.append({**definition, **raw_row})
    fields = definition_fields + [field for field in raw_fields if field not in definition_fields]
    return fields, merged


def summarize_mode(
    rmats_dir: Path,
    contrast: str,
    group1: str,
    group2: str,
    count_mode: str,
    min_reads: int,
    min_abs_delta: float,
    statistical_test: bool,
    max_fdr: float,
):
    summarized: list[dict[str, str]] = []
    source_fields: list[str] = []
    for event_type in EVENT_TYPES:
        fields, rows = load_event_rows(rmats_dir, event_type, count_mode)
        for field in fields:
            if field not in source_fields and field not in {"PValue", "FDR"}:
                source_fields.append(field)
        for row in rows:
            inc1 = parse_counts(row.get("IJC_SAMPLE_1", ""))
            skip1 = parse_counts(row.get("SJC_SAMPLE_1", ""))
            inc2 = parse_counts(row.get("IJC_SAMPLE_2", ""))
            skip2 = parse_counts(row.get("SJC_SAMPLE_2", ""))
            inc_len = float(row["IncFormLen"])
            skip_len = float(row["SkipFormLen"])

            levels1 = parse_levels(row.get("IncLevel1", ""))
            levels2 = parse_levels(row.get("IncLevel2", ""))
            if not levels1:
                levels1 = calculate_levels(inc1, skip1, inc_len, skip_len)
            if not levels2:
                levels2 = calculate_levels(inc2, skip2, inc_len, skip_len)
            mean1 = finite_mean(levels1)
            mean2 = finite_mean(levels2)
            delta = mean1 - mean2 if not math.isnan(mean1) and not math.isnan(mean2) else math.nan
            totals1 = [inc + skip for inc, skip in zip(inc1, skip1)]
            totals2 = [inc + skip for inc, skip in zip(inc2, skip2)]
            min_total1 = min(totals1) if totals1 else 0
            min_total2 = min(totals2) if totals2 else 0
            p_value_text = row.get("PValue", "NA") or "NA"
            fdr_text = row.get("FDR", "NA") or "NA"
            try:
                fdr_value = float(fdr_text)
            except ValueError:
                fdr_value = math.nan
            if statistical_test and math.isnan(fdr_value):
                raise ValueError(
                    f"statistical_test=true but rMATS returned no numeric FDR for {event_type}"
                )
            passes = (
                not math.isnan(delta)
                and abs(delta) >= min_abs_delta
                and min_total1 >= min_reads
                and min_total2 >= min_reads
                and (not statistical_test or fdr_value <= max_fdr)
            )

            clean_source = {
                key: value for key, value in row.items() if key not in {"PValue", "FDR"}
            }
            summarized.append(
                {
                    "contrast": contrast,
                    "group1": group1,
                    "group2": group2,
                    "event_type": event_type,
                    "event_key": event_key(event_type, row),
                    "mean_psi_group1": format_float(mean1),
                    "mean_psi_group2": format_float(mean2),
                    "delta_psi_group1_minus_group2": format_float(delta),
                    "abs_delta_psi": format_float(abs(delta)) if not math.isnan(delta) else "NA",
                    "min_total_junction_reads_group1": str(min_total1),
                    "min_total_junction_reads_group2": str(min_total2),
                    "passes_exploratory_filter": "yes" if passes else "no",
                    "p_value": p_value_text if statistical_test else "NA",
                    "fdr": fdr_text if statistical_test else "NA",
                    "statistical_test": "run" if statistical_test else "not_run_no_biological_replicates",
                    **clean_source,
                }
            )
    derived_fields = [
        "contrast",
        "group1",
        "group2",
        "event_type",
        "event_key",
        "mean_psi_group1",
        "mean_psi_group2",
        "delta_psi_group1_minus_group2",
        "abs_delta_psi",
        "min_total_junction_reads_group1",
        "min_total_junction_reads_group2",
        "passes_exploratory_filter",
        "p_value",
        "fdr",
        "statistical_test",
    ]
    output_fields = derived_fields + [field for field in source_fields if field not in derived_fields]
    summarized.sort(
        key=lambda row: (
            row["abs_delta_psi"] == "NA",
            -float(row["abs_delta_psi"]) if row["abs_delta_psi"] != "NA" else 0,
            row["event_type"],
            row.get("geneSymbol", ""),
        )
    )
    candidates = [row for row in summarized if row["passes_exploratory_filter"] == "yes"]
    return output_fields, summarized, candidates


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_rmats(
    rmats_dir: Path,
    contrast: str,
    group1: str,
    group2: str,
    summary_dir: Path,
    min_reads: int,
    min_abs_delta: float,
    statistical_test: bool = False,
    max_fdr: float = 0.05,
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for mode in ("JC", "JCEC"):
        fields, all_rows, candidates = summarize_mode(
            rmats_dir,
            contrast,
            group1,
            group2,
            mode,
            min_reads,
            min_abs_delta,
            statistical_test,
            max_fdr,
        )
        write_tsv(summary_dir / f"all_events.{mode}.tsv", fields, all_rows)
        write_tsv(summary_dir / f"candidates.{mode}.tsv", fields, candidates)
        counts[f"all_{mode}"] = len(all_rows)
        counts[f"candidates_{mode}"] = len(candidates)
    return counts
