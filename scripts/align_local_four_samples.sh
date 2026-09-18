#!/usr/bin/env bash
# 로컬 4개 시료(15J_0hr 포함)를 동일 파라미터로 정렬한다.
# 원본 프로젝트 디렉터리는 읽기 전용으로 취급하고 출력은 workspace에만 쓴다.
set -euo pipefail
WS=/home/shim/Downloads/develop_shim_project_260906
SRC=/home/shim/Downloads/project/project_shim/RNA_seq
IDX=$SRC                                  # GENCODE v41 GRCh38, sjdbOverhang 100
OUT=$WS/public_data_tierA/derived/local_bam_v41
TRIM=$OUT/trim
mkdir -p "$OUT" "$TRIM"
THREADS=20
# samtools 1.3.1 fails to index STAR 2.7.3a BAMs ("corrupted or unsorted") although
# quickcheck passes; use the 1.21 build in the conda env instead.
SAMTOOLS=/home/shim/anaconda3/envs/qiime2-amplicon-2024.10/bin/samtools
[ -x "$SAMTOOLS" ] || SAMTOOLS=samtools

# 1) 15J_0hr 만 trimming (나머지 3개는 2022년 Trim Galore 결과 재사용)
if [ ! -s "$TRIM/15J_0hr.R1_val_1.fq.gz" ]; then
  echo "[$(date +%T)] trim_galore 15J_0hr"
  trim_galore --paired --gzip --cores 4 --output_dir "$TRIM" \
    "$WS/00.RNA-seq_samle/15J_0hr.R1.fastq.gz" "$WS/00.RNA-seq_samle/15J_0hr.R2.fastq.gz"
fi

align () {           # $1=label  $2=R1  $3=R2
  local lab=$1 r1=$2 r3=$3
  if [ -s "$OUT/${lab}.Aligned.sortedByCoord.out.bam" ]; then
    echo "[$(date +%T)] $lab already aligned; skip"; return
  fi
  echo "[$(date +%T)] STAR $lab"
  STAR --runThreadN $THREADS \
       --genomeDir "$IDX" \
       --readFilesIn "$r1" "$r3" \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --limitBAMsortRAM 30000000000 \
       --outFileNamePrefix "$OUT/${lab}."
  "$SAMTOOLS" index -@ 8 "$OUT/${lab}.Aligned.sortedByCoord.out.bam"
  echo "[$(date +%T)] $lab done"
}

align HDF  "$SRC/01.trim/HDF_CM_24h.R1_val_1.fq.gz"  "$SRC/01.trim/HDF_CM_24h.R2_val_2.fq.gz"
align REP  "$SRC/01.trim/Rep_CM_24h.R1_val_1.fq.gz"  "$SRC/01.trim/Rep_CM_24h.R2_val_2.fq.gz"
align IPS  "$SRC/01.trim/iPSC_CM_24h.R1_val_1.fq.gz" "$SRC/01.trim/iPSC_CM_24h.R2_val_2.fq.gz"
align UVA0 "$TRIM/15J_0hr.R1_val_1.fq.gz"            "$TRIM/15J_0hr.R2_val_2.fq.gz"

echo "[$(date +%T)] ALL DONE"
ls -la "$OUT"/*.bam
