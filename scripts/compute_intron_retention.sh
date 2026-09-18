#!/usr/bin/env bash
# 로컬 4개 시료의 intron retention burden.
# spliced read가 intron을 가로지르는 것을 intronic coverage로 오계수하지 않도록 -split 필수.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
BAMDIR=$WS/public_data_tierA/derived/local_bam_v41
ANN=$WS/public_data_tierA/derived/ir_annotation_v41
OUT=$WS/public_data_tierA/derived/intron_retention
mkdir -p "$OUT"

# BAM 헤더 순서로 genome file을 만들고 BED를 같은 순서로 재정렬한다.
samtools view -H "$BAMDIR/HDF.Aligned.sortedByCoord.out.bam" \
  | awk -F'\t' '/^@SQ/{n="";l=""; for(i=2;i<=NF;i++){if($i~/^SN:/)n=substr($i,4); if($i~/^LN:/)l=substr($i,4)} print n"\t"l}' \
  > "$OUT/genome.txt"
for r in exonic_unique intronic_unique; do
  bedtools sort -faidx "$OUT/genome.txt" -i "$ANN/$r.bed" > "$OUT/$r.gsorted.bed"
done
echo "regions: exonic=$(wc -l < "$OUT/exonic_unique.gsorted.bed")  intronic=$(wc -l < "$OUT/intronic_unique.gsorted.bed")"

count_one () {   # $1 = sample label
  local s=$1 bam="$BAMDIR/$1.Aligned.sortedByCoord.out.bam"
  for r in exonic_unique intronic_unique; do
    [ -s "$OUT/${s}.${r}.counts" ] && continue
    bedtools coverage -counts -split -sorted -g "$OUT/genome.txt" \
      -a "$OUT/$r.gsorted.bed" -b "$bam" > "$OUT/${s}.${r}.counts.tmp"
    mv "$OUT/${s}.${r}.counts.tmp" "$OUT/${s}.${r}.counts"
  done
  echo "  counted $s"
}
for s in HDF REP IPS UVA0; do count_one "$s" & done
wait
echo "ALL COUNTS DONE"
