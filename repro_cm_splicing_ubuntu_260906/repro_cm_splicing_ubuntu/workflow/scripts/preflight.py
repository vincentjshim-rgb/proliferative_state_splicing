#!/usr/bin/env python3
"""Validate paired FASTQ inputs and record an immutable input inventory."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import statistics
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import yaml


REQUIRED_COLUMNS = {
    "sample_id",
    "condition",
    "biological_replicate",
    "r1",
    "r2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--json", required=True, dest="json_out")
    parser.add_argument("--profiles", required=True)
    parser.add_argument("--checksums", required=True)
    parser.add_argument("--skip-sha256", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"manifest is missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("manifest has no samples")
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("sample_id values must be unique")
    return rows


def normalize_read_id(header: str) -> str:
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def read_record(handle, path: Path, record_number: int):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\r\n")
    plus = handle.readline()
    quality = handle.readline().rstrip("\r\n")
    if not sequence or not plus or quality == "":
        raise ValueError(f"truncated FASTQ record {record_number} in {path}")
    if not header.startswith("@"):
        raise ValueError(f"invalid FASTQ header at record {record_number} in {path}")
    if not plus.startswith("+"):
        raise ValueError(f"invalid FASTQ plus line at record {record_number} in {path}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"sequence/quality length mismatch at record {record_number} in {path}"
        )
    return header, sequence


def profile_pair(r1: Path, r2: Path, scan_records: int) -> tuple[dict, dict]:
    lengths1: list[int] = []
    lengths2: list[int] = []
    with gzip.open(r1, "rt", encoding="ascii", errors="strict") as h1, gzip.open(
        r2, "rt", encoding="ascii", errors="strict"
    ) as h2:
        for idx in range(1, scan_records + 1):
            rec1 = read_record(h1, r1, idx)
            rec2 = read_record(h2, r2, idx)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError(
                    f"R1/R2 record-count mismatch within first {scan_records} records: "
                    f"{r1} ; {r2}"
                )
            if normalize_read_id(rec1[0]) != normalize_read_id(rec2[0]):
                raise ValueError(
                    f"R1/R2 identifiers differ at record {idx}: "
                    f"{rec1[0].strip()} ; {rec2[0].strip()}"
                )
            lengths1.append(len(rec1[1]))
            lengths2.append(len(rec2[1]))
    if not lengths1:
        raise ValueError(f"empty FASTQ pair: {r1} ; {r2}")
    return length_summary(lengths1), length_summary(lengths2)


def length_summary(lengths: list[int]) -> dict:
    counts = Counter(lengths)
    modal_length, modal_count = sorted(
        counts.items(), key=lambda item: (-item[1], -item[0])
    )[0]
    return {
        "records_scanned": len(lengths),
        "min_read_length": min(lengths),
        "max_read_length": max(lengths),
        "modal_read_length": modal_length,
        "modal_fraction": round(modal_count / len(lengths), 6),
        "mean_read_length": round(statistics.fmean(lengths), 3),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mkdir_for(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def main() -> int:
    args = parse_args()
    manifest_path = Path(args.manifest).resolve()
    config_path = Path(args.config).resolve()
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    scan_records = int(config["analysis"].get("fastq_scan_records", 10000))
    rows = read_manifest(manifest_path)

    errors: list[str] = []
    warnings: list[str] = []
    profiles: list[dict] = []
    checksums: list[dict] = []
    maximum_lengths: list[int] = []
    condition_counts = Counter(row["condition"] for row in rows)

    contrasts_path = Path(config["contrasts"]).resolve()
    if not contrasts_path.is_file():
        errors.append(f"contrasts file not found: {contrasts_path}")
    else:
        with contrasts_path.open(newline="", encoding="utf-8") as handle:
            contrast_rows = list(csv.DictReader(handle, delimiter="\t"))
        for contrast in contrast_rows:
            for field in ("group1", "group2"):
                condition = contrast[field]
                if condition not in condition_counts:
                    errors.append(
                        f"contrast {contrast['contrast_id']} references absent condition: {condition}"
                    )
            if config["analysis"]["rmats"].get("statistical_test", False):
                if any(condition_counts.get(contrast[field], 0) < 2 for field in ("group1", "group2")):
                    errors.append(
                        f"statistical_test=true requires at least two biological samples in both "
                        f"groups of {contrast['contrast_id']}"
                    )

    for condition, count in sorted(condition_counts.items()):
        if count == 1:
            warnings.append(
                f"{condition}: one biological sample; statistical differential splicing is unsupported"
            )

    for row in rows:
        sample_id = row["sample_id"]
        r1 = Path(row["r1"]).expanduser()
        r2 = Path(row["r2"]).expanduser()
        for mate, path in (("R1", r1), ("R2", r2)):
            if not path.is_absolute():
                errors.append(f"{sample_id} {mate}: path must be absolute: {path}")
            if not path.is_file():
                errors.append(f"{sample_id} {mate}: file not found: {path}")
            elif not os.access(path, os.R_OK):
                errors.append(f"{sample_id} {mate}: file is not readable: {path}")
            elif path.suffix != ".gz":
                warnings.append(f"{sample_id} {mate}: expected .gz input: {path}")
        if any(not path.is_file() for path in (r1, r2)):
            continue
        try:
            r1_profile, r2_profile = profile_pair(r1, r2, scan_records)
        except (OSError, EOFError, UnicodeError, ValueError) as exc:
            errors.append(f"{sample_id}: {exc}")
            continue

        for mate, path, profile in (
            ("R1", r1, r1_profile),
            ("R2", r2, r2_profile),
        ):
            maximum_lengths.append(profile["max_read_length"])
            profiles.append(
                {
                    "sample_id": sample_id,
                    "condition": row["condition"],
                    "biological_replicate": row["biological_replicate"],
                    "mate": mate,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    **profile,
                }
            )
            checksums.append(
                {
                    "sample_id": sample_id,
                    "mate": mate,
                    "path": str(path),
                    "sha256": "SKIPPED" if args.skip_sha256 else sha256(path),
                }
            )

    if maximum_lengths and len(set(maximum_lengths)) > 1:
        warnings.append(
            "maximum read lengths differ across files; rMATS will use the maximum "
            "with --variable-read-length"
        )
    nominal_read_length = max(maximum_lengths) if maximum_lengths else None

    json_out = Path(args.json_out)
    profiles_out = Path(args.profiles)
    checksums_out = Path(args.checksums)
    for output in (json_out, profiles_out, checksums_out):
        mkdir_for(output)

    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest_path),
        "config": str(config_path),
        "sample_count": len(rows),
        "condition_counts": dict(condition_counts),
        "biological_replicates_per_condition": dict(condition_counts),
        "scan_records_requested": scan_records,
        "nominal_read_length": nominal_read_length,
        "errors": errors,
        "warnings": warnings,
        "analysis_status": "blocked" if errors else "ready_for_exploratory_analysis",
        "statistical_warning": (
            "The present manifest has one biological sample per condition. "
            "Do not report differential-splicing P values or FDR."
        ),
    }
    json_out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    profile_fields = [
        "sample_id",
        "condition",
        "biological_replicate",
        "mate",
        "path",
        "size_bytes",
        "records_scanned",
        "min_read_length",
        "max_read_length",
        "modal_read_length",
        "modal_fraction",
        "mean_read_length",
    ]
    with profiles_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=profile_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(profiles)

    with checksums_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "mate", "sha256", "path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(checksums)

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    print(f"Preflight passed for {len(rows)} paired-end samples")
    print(f"Detected nominal read length: {nominal_read_length}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
