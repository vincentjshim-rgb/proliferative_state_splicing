#!/usr/bin/env python3
"""Record tool and host versions used for a workflow run."""

from __future__ import annotations

import datetime
import platform
import shutil
import subprocess
from pathlib import Path


COMMANDS = {
    "python": ["python", "--version"],
    "fastqc": ["fastqc", "--version"],
    "fastp": ["fastp", "--version"],
    "multiqc": ["multiqc", "--version"],
    "STAR": ["STAR", "--version"],
    "samtools": ["samtools", "--version"],
    "featureCounts": ["featureCounts", "-v"],
    "rMATS": ["rmats.py", "--version"],
}


lines = [
    f"generated_at_utc\t{datetime.datetime.now(datetime.timezone.utc).isoformat()}",
    f"platform\t{platform.platform()}",
    f"machine\t{platform.machine()}",
]
for name, command in COMMANDS.items():
    executable = shutil.which(command[0])
    if not executable:
        lines.append(f"{name}\tNOT_FOUND")
        continue
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    message = (completed.stdout + " " + completed.stderr).strip().replace("\n", " | ")
    lines.append(f"{name}\t{message}\texit_code={completed.returncode}")

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("\n".join(lines) + "\n", encoding="utf-8")
