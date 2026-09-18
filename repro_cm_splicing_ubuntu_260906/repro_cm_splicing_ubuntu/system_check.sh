#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="${project_root}/config/samples.tsv"

echo "Project root: ${project_root}"
echo "Kernel: $(uname -srmo)"
echo "CPU threads: $(getconf _NPROCESSORS_ONLN)"
echo
free -h
echo
df -h "${project_root}"
echo
echo "FASTQ files:"
tail -n +2 "${manifest}" | while IFS=$'\t' read -r sample condition replicate r1 r2; do
    for path in "${r1}" "${r2}"; do
        if [[ -r "${path}" ]]; then
            stat --printf='%n\t%s bytes\n' "${path}"
        else
            echo "MISSING_OR_UNREADABLE: ${path}"
        fi
    done
done
