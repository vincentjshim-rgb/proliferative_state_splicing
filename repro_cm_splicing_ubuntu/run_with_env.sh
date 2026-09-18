#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner_prefix="${project_root}/.runner_env"

if [[ ! -x "${runner_prefix}/bin/snakemake" ]]; then
    echo "ERROR: runner environment not found: ${runner_prefix}" >&2
    echo "Run ./setup_runner_env.sh first." >&2
    exit 1
fi

export PATH="${runner_prefix}/bin:${PATH}"
exec "${project_root}/run.sh" "$@"
