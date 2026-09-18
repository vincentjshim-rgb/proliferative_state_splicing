#!/usr/bin/env python3
"""Snakemake script: combine per-contrast exploratory rMATS candidate tables."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


contrast_metadata = {
    row["contrast_id"]: row for row in read_tsv(Path(snakemake.input.contrasts))
}
all_rows: list[dict[str, str]] = []
for path in snakemake.input.candidates:
    rows = read_tsv(Path(path))
    for row in rows:
        metadata = contrast_metadata[row["contrast"]]
        all_rows.append(
            {
                "priority": metadata["priority"],
                "contrast_interpretation": metadata["interpretation"],
                **row,
            }
        )

long_fields: list[str] = []
for row in all_rows:
    for field in row:
        if field not in long_fields:
            long_fields.append(field)

long_path = Path(snakemake.output.long)
ranked_path = Path(snakemake.output.ranked)
report_path = Path(snakemake.output.report)
for output in (long_path, ranked_path, report_path):
    output.parent.mkdir(parents=True, exist_ok=True)

with long_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=long_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(all_rows)

grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
for row in all_rows:
    grouped[row["event_key"]].append(row)

ranked_rows: list[dict[str, str]] = []
for key, rows in grouped.items():
    primary_rows = [row for row in rows if row["priority"] == "primary"]
    deltas = [float(row["delta_psi_group1_minus_group2"]) for row in rows]
    primary_deltas = [
        float(row["delta_psi_group1_minus_group2"]) for row in primary_rows
    ]
    primary_direction_consistent = (
        "yes"
        if len(primary_deltas) >= 2
        and (all(value > 0 for value in primary_deltas) or all(value < 0 for value in primary_deltas))
        else "no"
    )
    min_coverage = min(
        min(
            int(row["min_total_junction_reads_group1"]),
            int(row["min_total_junction_reads_group2"]),
        )
        for row in rows
    )
    representative = rows[0]
    ranked_rows.append(
        {
            "event_key": key,
            "event_type": representative["event_type"],
            "GeneID": representative.get("GeneID", ""),
            "geneSymbol": representative.get("geneSymbol", ""),
            "chr": representative.get("chr", ""),
            "strand": representative.get("strand", ""),
            "primary_contrast_support_count": str(len(primary_rows)),
            "all_contrast_support_count": str(len(rows)),
            "primary_direction_consistent": primary_direction_consistent,
            "max_abs_delta_psi": f"{max(abs(value) for value in deltas):.6f}",
            "minimum_junction_coverage_seen": str(min_coverage),
            "contrasts": ";".join(row["contrast"] for row in rows),
            "delta_psi_values": ";".join(
                f"{row['contrast']}={row['delta_psi_group1_minus_group2']}" for row in rows
            ),
            "evidence_level": "exploratory_N_of_1",
        }
    )

ranked_rows.sort(
    key=lambda row: (
        -int(row["primary_contrast_support_count"]),
        row["primary_direction_consistent"] != "yes",
        -int(row["all_contrast_support_count"]),
        -float(row["max_abs_delta_psi"]),
        row["geneSymbol"],
    )
)
ranked_fields = [
    "event_key",
    "event_type",
    "GeneID",
    "geneSymbol",
    "chr",
    "strand",
    "primary_contrast_support_count",
    "all_contrast_support_count",
    "primary_direction_consistent",
    "max_abs_delta_psi",
    "minimum_junction_coverage_seen",
    "contrasts",
    "delta_psi_values",
    "evidence_level",
]
with ranked_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=ranked_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(ranked_rows)

double_primary = sum(
    1
    for row in ranked_rows
    if row["primary_contrast_support_count"] == "2"
    and row["primary_direction_consistent"] == "yes"
)
report_path.write_text(
    "\n".join(
        [
            "Repro-CM exploratory splicing candidate catalog",
            "",
            f"Candidate rows across all contrasts: {len(all_rows)}",
            f"Unique splice events: {len(ranked_rows)}",
            f"Events supported by both primary Repro-CM contrasts with consistent direction: {double_primary}",
            "",
            "Ranking is based on cross-comparison consistency, absolute delta PSI, and junction coverage.",
            "No differential-splicing P value or FDR was calculated because there is one biological sample per condition.",
            "The 24 h versus 0 h comparisons confound treatment with elapsed time.",
            "Candidates require independent junction RT-PCR/ddPCR or replicated RNA-seq validation.",
            "",
        ]
    ),
    encoding="utf-8",
)
