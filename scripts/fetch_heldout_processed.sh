#!/usr/bin/env bash
# Figure 6 held-out datasets. 사전등록 문서(repro_cm_preregistration_phase2_ko.md) 고정 후 다운로드.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
OUT=$WS/public_data_tierA/heldout
META=$OUT/metadata
mkdir -p "$OUT" "$META"
get () { # $1=url $2=dest
  local d="$2"; [ -s "$d" ] && { echo "skip $(basename $d)"; return; }
  if curl -fsSL -m 1200 -o "$d.part" "$1"; then mv "$d.part" "$d"
    echo "ok   $(basename $d)  $(stat -c%s "$d") bytes"
  else rm -f "$d.part"; echo "FAIL $(basename $d)"; fi
}
declare -A F=(
 [GSE297233]="GSE297233_raw_counts_matrix.csv.gz"
 [GSE307377]="GSE307377_raw.txt.gz GSE307377_tpm.txt.gz"
 [GSE116968]="GSE116968_Processed_data_NHDF_RNAseq_1_post_1h_.xlsx GSE116968_Processed_data_NHDF_RNAseq_2_Post_4h_.xlsx"
 [GSE149694]="GSE149694_RAW.tar filelist.txt"
)
for G in GSE297233 GSE307377 GSE116968 GSE149694; do
  P=$(echo $G | sed 's/...$/nnn/')
  for f in ${F[$G]}; do
    [ "$f" = "filelist.txt" ] && dest="$OUT/${G}_filelist.txt" || dest="$OUT/$f"
    get "https://ftp.ncbi.nlm.nih.gov/geo/series/${P}/${G}/suppl/${f}" "$dest"
  done
  get "https://ftp.ncbi.nlm.nih.gov/geo/series/${P}/${G}/soft/${G}_family.soft.gz" "$META/${G}_family.soft.gz"
done
( cd "$OUT" && sha256sum *.gz *.tar *.xlsx 2>/dev/null > heldout_sha256.txt || true )
echo "TOTAL: $(du -sh "$OUT" | cut -f1)"
