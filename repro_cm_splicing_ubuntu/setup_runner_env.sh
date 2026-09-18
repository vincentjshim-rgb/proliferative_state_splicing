#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runner_prefix="${project_root}/.runner_env"
environment_file="${project_root}/envs/snakemake.yaml"

if command -v micromamba >/dev/null 2>&1; then
    if [[ -d "${runner_prefix}" ]]; then
        micromamba env update -y -p "${runner_prefix}" -f "${environment_file}"
    else
        micromamba create -y -p "${runner_prefix}" -f "${environment_file}"
    fi
elif command -v mamba >/dev/null 2>&1; then
    if [[ -d "${runner_prefix}" ]]; then
        mamba env update -p "${runner_prefix}" -f "${environment_file}" --prune
    else
        mamba env create -p "${runner_prefix}" -f "${environment_file}"
    fi
elif command -v conda >/dev/null 2>&1; then
    if [[ -d "${runner_prefix}" ]]; then
        conda env update -p "${runner_prefix}" -f "${environment_file}" --prune
    else
        conda env create -p "${runner_prefix}" -f "${environment_file}"
    fi
else
    echo "ERROR: micromamba, mamba, or conda is required." >&2
    echo "Install one package manager, then rerun this script." >&2
    exit 1
fi

echo "Runner environment is ready: ${runner_prefix}"
echo "Next: ./setup_reference.sh"
