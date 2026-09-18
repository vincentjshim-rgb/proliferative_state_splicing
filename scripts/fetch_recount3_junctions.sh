#!/usr/bin/env bash
# recount3 junction 자료 다운로드. FASTQ 재정렬을 대체한다.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
OUT=$WS/public_data_tierA/recount3
mkdir -p "$OUT"
BASE=https://duffel.rail.bio/recount3/human/data_sources/sra
# SRP144355 GSE113957 fibroblast donors age 1-96 (n=143)
# SRP096629 GSE93535  quiescent vs SIPS at matched PD15 (+1201 rescue)
# SRP131506 GSE109700 deep replicative senescence
for SRP in SRP144355 SRP096629 SRP131506; do
  T=${SRP: -2}
  for kind in "junctions/$T/$SRP/sra.junctions.$SRP.ALL.MM.gz" \
              "junctions/$T/$SRP/sra.junctions.$SRP.ALL.RR.gz" \
              "junctions/$T/$SRP/sra.junctions.$SRP.ALL.ID.gz" \
              "metadata/$T/$SRP/sra.recount_project.$SRP.MD.gz" \
              "metadata/$T/$SRP/sra.sra.$SRP.MD.gz" \
              "metadata/$T/$SRP/sra.recount_qc.$SRP.MD.gz" \
              "metadata/$T/$SRP/sra.recount_seq_qc.$SRP.MD.gz" \
              "gene_sums/$T/$SRP/sra.gene_sums.$SRP.G026.gz"; do
    f=$(basename "$kind")
    [ -s "$OUT/$f" ] && { echo "skip $f"; continue; }
    if curl -fsSL -m 900 -o "$OUT/$f.part" "$BASE/$kind"; then
      mv "$OUT/$f.part" "$OUT/$f"; echo "ok   $f  $(stat -c%s "$OUT/$f") bytes"
    else
      rm -f "$OUT/$f.part"; echo "MISS $f"
    fi
  done
done
( cd "$OUT" && sha256sum *.gz > recount3_sha256.txt )
echo "TOTAL: $(du -sh "$OUT" | cut -f1)"
