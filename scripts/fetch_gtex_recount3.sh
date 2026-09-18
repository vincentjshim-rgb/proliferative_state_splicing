#!/usr/bin/env bash
# GTEx v8 (recount3) processed gene sums + metadata for the culture-vs-tissue boundary test.
# Processed files only: no FASTQ, no junctions.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
OUT=$WS/public_data_tierA/gtex
mkdir -p "$OUT"
BASE=https://duffel.rail.bio/recount3/human/data_sources/gtex
for P in SKIN MUSCLE; do
  T=${P: -2}
  for kind in "gene_sums/$T/$P/gtex.gene_sums.$P.G026.gz" \
              "metadata/$T/$P/gtex.gtex.$P.MD.gz" \
              "metadata/$T/$P/gtex.recount_qc.$P.MD.gz" \
              "metadata/$T/$P/gtex.recount_project.$P.MD.gz"; do
    f=$(basename "$kind")
    [ -s "$OUT/$f" ] && { echo "skip $f"; continue; }
    if curl -fsSL -m 1800 -o "$OUT/$f.part" "$BASE/$kind"; then
      mv "$OUT/$f.part" "$OUT/$f"
      if gzip -t "$OUT/$f"; then echo "ok   $f  $(stat -c%s "$OUT/$f") bytes  gzip-ok"
      else echo "BAD  $f  gzip integrity failed"; mv "$OUT/$f" "$OUT/$f.corrupt"; fi
    else
      rm -f "$OUT/$f.part"; echo "MISS $f"
    fi
  done
done
( cd "$OUT" && sha256sum *.gz > gtex_sha256.txt )
for P in SKIN MUSCLE; do T=${P: -2}
  for kind in gene_sums/$T/$P/gtex.gene_sums.$P.G026.gz metadata/$T/$P/gtex.gtex.$P.MD.gz \
              metadata/$T/$P/gtex.recount_qc.$P.MD.gz metadata/$T/$P/gtex.recount_project.$P.MD.gz; do
    printf "GTEx_v8_recount3_%s\tgtex\tculture_tissue_boundary\t%s\tgtex/%s\n" "$P" "$BASE/$kind" "$(basename $kind)" \
      >> $WS/public_data_tierA/download_manifest.tsv
  done
done
echo "TOTAL: $(du -sh "$OUT" | cut -f1)"
