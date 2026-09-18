#!/usr/bin/env python3
"""Snakemake script: build a STAR index with the detected nominal read length."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


with Path(snakemake.input.profile_json).open(encoding="utf-8") as handle:
    profile = json.load(handle)
configured_read_length = snakemake.params.read_length
if str(configured_read_length).lower() == "auto":
    read_length = int(profile["nominal_read_length"])
else:
    read_length = int(configured_read_length)
sjdb_overhang = max(1, read_length - 1)

index_dir = Path(snakemake.output.index)
if index_dir.exists():
    shutil.rmtree(index_dir)
index_dir.mkdir(parents=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)

command = [
    "STAR",
    "--runMode",
    "genomeGenerate",
    "--runThreadN",
    str(snakemake.threads),
    "--genomeDir",
    str(index_dir),
    "--genomeFastaFiles",
    str(Path(snakemake.input.fasta).resolve()),
    "--sjdbGTFfile",
    str(Path(snakemake.input.gtf).resolve()),
    "--sjdbOverhang",
    str(sjdb_overhang),
]
with log_path.open("w", encoding="utf-8") as log_handle:
    log_handle.write("COMMAND\n" + " ".join(command) + "\n\n")
    log_handle.flush()
    subprocess.run(command, check=True, stdout=log_handle, stderr=subprocess.STDOUT)

(index_dir / "index_metadata.json").write_text(
    json.dumps(
        {
            "read_length": read_length,
            "sjdb_overhang": sjdb_overhang,
            "reference_fasta": str(Path(snakemake.input.fasta).resolve()),
            "annotation_gtf": str(Path(snakemake.input.gtf).resolve()),
        },
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
