#!/usr/bin/env bash
set -euo pipefail

root_dir="${1:-public_data_tierA}"
manifest="${root_dir}/download_manifest.tsv"
validation_log="${root_dir}/logs/validation.tsv"
checksum_log="${root_dir}/logs/sha256.tsv"
inventory_log="${root_dir}/logs/file_inventory.tsv"

if [[ ! -f "${manifest}" ]]; then
  printf 'Manifest not found: %s\n' "${manifest}" >&2
  exit 1
fi

mkdir -p "${root_dir}/logs"
printf 'dataset\ttarget\tbytes\tvalidator\tstatus\n' > "${validation_log}"

failures=0
while IFS=$'\t' read -r dataset category role url target; do
  candidate="${root_dir}/${target}.part"
  final="${root_dir}/${target}"

  if [[ -f "${candidate}" ]]; then
    input="${candidate}"
  elif [[ -f "${final}" ]]; then
    input="${final}"
  else
    printf '%s\t%s\t0\tmissing\tFAIL\n' "${dataset}" "${target}" >> "${validation_log}"
    failures=$((failures + 1))
    continue
  fi

  bytes="$(stat -c '%s' "${input}")"
  validator="unknown"
  status="FAIL"
  case "${target}" in
    *.gz)
      validator="gzip_test"
      if gzip -t "${input}"; then status="PASS"; fi
      ;;
    *.tar)
      validator="tar_list"
      if tar -tf "${input}" >/dev/null; then status="PASS"; fi
      ;;
    *.xlsx|*.zip|*.docx)
      validator="zip_test"
      if unzip -tqq "${input}"; then status="PASS"; fi
      ;;
    *)
      validator="nonempty"
      if [[ -s "${input}" ]]; then status="PASS"; fi
      ;;
  esac

  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${dataset}" "${target}" "${bytes}" "${validator}" "${status}" >> "${validation_log}"
  if [[ "${status}" != "PASS" ]]; then
    failures=$((failures + 1))
  fi
done < <(tail -n +2 "${manifest}")

if (( failures > 0 )); then
  printf '%s file(s) failed validation; no .part files were finalized.\n' "${failures}" >&2
  exit 2
fi

while IFS=$'\t' read -r dataset category role url target; do
  candidate="${root_dir}/${target}.part"
  final="${root_dir}/${target}"
  if [[ -f "${candidate}" ]]; then
    mv "${candidate}" "${final}"
  fi
done < <(tail -n +2 "${manifest}")

printf 'sha256\tbytes\ttarget\n' > "${checksum_log}"
printf 'dataset\tcategory\trole\tbytes\ttarget\n' > "${inventory_log}"
while IFS=$'\t' read -r dataset category role url target; do
  final="${root_dir}/${target}"
  checksum="$(sha256sum "${final}" | cut -d' ' -f1)"
  bytes="$(stat -c '%s' "${final}")"
  printf '%s\t%s\t%s\n' "${checksum}" "${bytes}" "${target}" >> "${checksum_log}"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "${dataset}" "${category}" "${role}" "${bytes}" "${target}" >> "${inventory_log}"
done < <(tail -n +2 "${manifest}")

printf 'Validated and finalized %s files.\n' "$(($(wc -l < "${manifest}") - 1))"
