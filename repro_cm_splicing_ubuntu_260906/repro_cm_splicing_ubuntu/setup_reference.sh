#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
release="50"
reference_dir="${project_root}/resources/gencode_v${release}"
base_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${release}"
fasta_gz="${reference_dir}/GRCh38.primary_assembly.genome.fa.gz"
gtf_gz="${reference_dir}/gencode.v${release}.primary_assembly.annotation.gtf.gz"
fasta="${fasta_gz%.gz}"
gtf="${gtf_gz%.gz}"

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl is required to download GENCODE." >&2
    exit 1
fi
if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip is required to validate and decompress GENCODE." >&2
    exit 1
fi

mkdir -p "${reference_dir}"

download_if_missing() {
    local url="$1"
    local destination="$2"
    if [[ -s "${destination}" ]]; then
        echo "Using existing file: ${destination}"
        return
    fi
    local partial="${destination}.part"
    curl --fail --location --retry 5 --retry-delay 5 --output "${partial}" "${url}"
    mv "${partial}" "${destination}"
}

download_if_missing \
    "${base_url}/GRCh38.primary_assembly.genome.fa.gz" \
    "${fasta_gz}"
download_if_missing \
    "${base_url}/gencode.v${release}.primary_assembly.annotation.gtf.gz" \
    "${gtf_gz}"

gzip -t "${fasta_gz}"
gzip -t "${gtf_gz}"

if [[ ! -s "${fasta}" ]]; then
    gzip -dc "${fasta_gz}" > "${fasta}.part"
    mv "${fasta}.part" "${fasta}"
fi
if [[ ! -s "${gtf}" ]]; then
    gzip -dc "${gtf_gz}" > "${gtf}.part"
    mv "${gtf}.part" "${gtf}"
fi

printf '%s\n' \
    "GENCODE release ${release}" \
    "${base_url}/GRCh38.primary_assembly.genome.fa.gz" \
    "${base_url}/gencode.v${release}.primary_assembly.annotation.gtf.gz" \
    > "${reference_dir}/source_urls.txt"

sha256sum "${fasta_gz}" "${gtf_gz}" "${fasta}" "${gtf}" \
    > "${reference_dir}/reference_sha256.txt"

echo "Reference files are ready in ${reference_dir}"
echo "The STAR index will be built by Snakemake using the detected read length."
