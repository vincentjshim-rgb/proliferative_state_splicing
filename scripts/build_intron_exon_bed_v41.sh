#!/usr/bin/env bash
# GENCODE v41 (로컬 BAM과 동일 annotation)에서 명확한 exonic/intronic 구간을 만든다.
# 다른 유전자의 exon 또는 유전자 span과 겹치는 구간은 모두 제거한다.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
GTF=/home/shim/Downloads/project/project_shim/RNA_seq/gencode.v41.chr_patch_hapl_scaff.annotation.gtf
OUT=$WS/public_data_tierA/derived/ir_annotation_v41
mkdir -p "$OUT"

echo "[1/6] exons"
awk -F'\t' 'BEGIN{OFS="\t"} $3=="exon" && $1~/^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $9~/gene_type "protein_coding"/ {
   match($9,/gene_name "[^"]+"/); g=substr($9,RSTART+11,RLENGTH-12);
   print $1,$4-1,$5,g,".",$7 }' "$GTF" | sort -k1,1 -k2,2n > "$OUT/exons.bed"

echo "[2/6] gene spans"
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $1~/^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $9~/gene_type "protein_coding"/ {
   match($9,/gene_name "[^"]+"/); g=substr($9,RSTART+11,RLENGTH-12);
   print $1,$4-1,$5,g,".",$7 }' "$GTF" | sort -k1,1 -k2,2n > "$OUT/genes.bed"

echo "[3/6] unambiguous exonic blocks (one gene only)"
bedtools merge -i "$OUT/exons.bed" -c 4 -o distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} $4 !~ /,/ && $3-$2>=50 {print $1,$2,$3,$4}' > "$OUT/exonic_unique.bed"

echo "[4/6] regions where gene spans overlap"
bedtools merge -i "$OUT/genes.bed" -c 4 -o distinct \
  | awk -F'\t' 'BEGIN{OFS="\t"} $4 ~ /,/ {print $1,$2,$3}' > "$OUT/gene_overlap.bed"

echo "[5/6] intronic blocks: gene span minus ANY exon minus overlapping gene spans"
bedtools merge -i "$OUT/exons.bed" | cut -f1-3 > "$OUT/exons_merged_all.bed"
bedtools subtract -a "$OUT/genes.bed" -b "$OUT/exons_merged_all.bed" \
  | bedtools subtract -a - -b "$OUT/gene_overlap.bed" \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3-$2>=50 {print $1,$2,$3,$4}' > "$OUT/intronic_unique.bed"

echo "[6/6] summary"
for f in exonic_unique intronic_unique; do
  printf "  %-18s blocks=%s  genes=%s  bp=%s\n" "$f" \
    "$(wc -l < "$OUT/$f.bed")" "$(cut -f4 "$OUT/$f.bed" | sort -u | wc -l)" \
    "$(awk '{s+=$3-$2} END{print s}' "$OUT/$f.bed")"
done
