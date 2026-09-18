#!/usr/bin/env python3
"""Snakemake script: run rMATS without statistics and summarize ΔPSI candidates."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.params.script_dir).resolve()))
from rmats_summary import summarize_rmats


def executable() -> str:
    for name in ("rmats.py", "rmats"):
        found = shutil.which(name)
        if found:
            return found
    raise RuntimeError("rMATS executable was not found in PATH")


def absolute_paths(values) -> list[str]:
    return [str(Path(value).resolve()) for value in values]


contrast = str(snakemake.wildcards.contrast)
group1 = str(snakemake.params.group1)
group2 = str(snakemake.params.group2)
bams1 = absolute_paths(snakemake.input.b1)
bams2 = absolute_paths(snakemake.input.b2)
gtf = str(Path(snakemake.input.gtf).resolve())
profile_json = Path(snakemake.input.profile_json)
with profile_json.open(encoding="utf-8") as handle:
    profile = json.load(handle)

configured_read_length = snakemake.params.read_length
if str(configured_read_length).lower() == "auto":
    read_length = int(profile["nominal_read_length"])
else:
    read_length = int(configured_read_length)

summary_dir = Path(snakemake.output.all_jc).parent
contrast_dir = summary_dir.parent
contrast_dir.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(prefix=".rmats_work_", dir=contrast_dir) as temp_name:
    work_dir = Path(temp_name)
    od = work_dir / "output"
    tmp = work_dir / "tmp"
    od.mkdir()
    tmp.mkdir()
    b1_file = work_dir / "b1.txt"
    b2_file = work_dir / "b2.txt"
    b1_file.write_text(",".join(bams1) + "\n", encoding="utf-8")
    b2_file.write_text(",".join(bams2) + "\n", encoding="utf-8")

    command = [
        executable(),
        "--b1",
        str(b1_file),
        "--b2",
        str(b2_file),
        "--gtf",
        gtf,
        "-t",
        "paired",
        "--readLength",
        str(read_length),
        "--variable-read-length",
        "--libType",
        str(snakemake.params.library_type),
        "--nthread",
        str(snakemake.threads),
        "--od",
        str(od),
        "--tmp",
        str(tmp),
        "--task",
        "both",
        "--individual-counts",
    ]
    if not bool(snakemake.params.statistical_test):
        command.append("--statoff")
    if bool(snakemake.params.allow_clipping):
        command.append("--allow-clipping")
    if bool(snakemake.params.novel_splice_sites):
        command.append("--novelSS")

    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write("COMMAND\n" + " ".join(command) + "\n\n")
        log_handle.flush()
        subprocess.run(
            command,
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            env=os.environ.copy(),
        )

    temporary_summary = work_dir / "summary"
    counts = summarize_rmats(
        rmats_dir=od,
        contrast=contrast,
        group1=group1,
        group2=group2,
        summary_dir=temporary_summary,
        min_reads=int(snakemake.params.min_reads),
        min_abs_delta=float(snakemake.params.min_abs_delta),
        statistical_test=bool(snakemake.params.statistical_test),
        max_fdr=float(snakemake.params.max_fdr),
    )

    report = temporary_summary / "README.txt"
    report.write_text(
        "\n".join(
            [
                f"Contrast: {contrast}",
                f"Group 1: {group1}",
                f"Group 2: {group2}",
                f"Delta PSI direction: mean({group1}) - mean({group2})",
                f"Read length used by rMATS: {read_length}",
                f"Library type: {snakemake.params.library_type}",
                (
                    "Statistical test: RUN"
                    if bool(snakemake.params.statistical_test)
                    else "Statistical test: NOT RUN (--statoff)"
                ),
                (
                    f"FDR filter: <= {snakemake.params.max_fdr}"
                    if bool(snakemake.params.statistical_test)
                    else "Reason: the present design has one biological sample per condition."
                ),
                f"Exploratory filter: abs(delta PSI) >= {snakemake.params.min_abs_delta} and",
                f"  total inclusion+skipping junction reads >= {snakemake.params.min_reads}",
                "  in every sample of both groups.",
                f"JC events: {counts['all_JC']}",
                f"JC candidates: {counts['candidates_JC']}",
                f"JCEC events: {counts['all_JCEC']}",
                f"JCEC candidates: {counts['candidates_JCEC']}",
                "Do not interpret this table as statistically significant differential splicing.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    final_raw = contrast_dir / "raw"
    if final_raw.exists():
        shutil.rmtree(final_raw)
    shutil.move(str(od), str(final_raw))
    summary_dir.mkdir(parents=True, exist_ok=True)
    for source in temporary_summary.iterdir():
        destination = summary_dir / source.name
        os.replace(source, destination)
