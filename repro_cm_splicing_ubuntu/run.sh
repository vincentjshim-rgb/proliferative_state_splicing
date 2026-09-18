#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${project_root}"

if ! command -v snakemake >/dev/null 2>&1; then
    echo "ERROR: snakemake is not available in PATH." >&2
    echo "Run ./setup_runner_env.sh and then use ./run_with_env.sh" >&2
    exit 1
fi

workflow_cores="${CORES:-16}"
snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --cores "${workflow_cores}" \
    --software-deployment-method conda \
    --conda-prefix .snakemake/conda \
    --rerun-incomplete \
    --printshellcmds \
    "$@"
