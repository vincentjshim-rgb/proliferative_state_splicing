# Tier B raw RNA-seq 다운로드 의사결정표

작성일: 2026-09-11  
상태: ENA run metadata만 조회함. 아래 FASTQ는 아직 다운로드하지 않음.

## 1. 정확한 project 규모

| dataset | BioProject | run | layout | 전체 FASTQ |
|---|---|---:|---|---:|
| GSE191055 | PRJNA789540 | 37 | paired | 47.35 GiB |
| GSE109700 | PRJNA431791 | 9 | paired | 115.70 GiB |
| GSE302943 | PRJNA1292871 | 10 | single | 11.00 GiB |
| GSE240226 | PRJNA1003002 | 12 | paired | 26.04 GiB |
| GSE165177 | PRJNA693500 | 95 | single | 54.76 GiB |
| GSE179848 | PRJNA745334 | 345 | paired | 1,648.81 GiB |
| GSE226189 | PRJNA939148 | 82 | paired | 1,055.58 GiB |
| GSE93535 | PRJNA361062 | 13 | single | 17.36 GiB |

project 전체 다운로드는 GSE191055, GSE179848, GSE226189에서 불필요하다. sample-level subset을 사용한다.

## 2. 분석에 필요한 subset과 실제 크기

| 우선 | subset | run/sample | FASTQ | 권장 작업공간(약 3배) | 판정 |
|---:|---|---:|---:|---:|---|
| 1 | GSE191055 P4 n=4 + P27 n=4 | 8/8 | 2.06 GiB | 6.2 GiB | 다운로드 1순위 |
| 2 | GSE302943 UVA n=5 + control n=5 | 10/10 | 11.00 GiB | 33.0 GiB | 다운로드 1순위 |
| 3 | GSE93535 Q/Q1201/SIPS/SIPS1201 | 13 run/12 sample | 17.36 GiB | 52.1 GiB | 1–2번 signal 후 |
| 4 | GSE240226 네 arm 전체 | 12/12 | 26.04 GiB | 78.1 GiB | UVA rescue splicing이 필요할 때 |
| 5 | GSE165177 day-13 matched set | 8/8 | 4.40 GiB | 13.2 GiB | donor n=2라 낮은 우선순위 |
| 6 | GSE109700 전체 | 9/9 | 115.70 GiB | 347.1 GiB | 비용 대비 후순위 |
| 7 | GSE179848 HC1–HC4 first2/last2 | 16/16 | 88.50 GiB | 265.5 GiB | processed 결과로 충분, raw 보류 |
| 8 | GSE226189 youngest10+oldest10 | 20/20 | 227.01 GiB | 681.0 GiB | raw 보류 |

## 3. 권장 다운로드 wave

### Wave B1

- GSE191055 P4/P27: 2.06 GiB
- GSE302943 UVA/control: 11.00 GiB

합계 FASTQ 13.06 GiB, 최소 작업공간 약 39.2 GiB다. 이 두 자료로 replicative senescence와 UVA의 exact splicing event가 모두 검증되는지 먼저 판정한다.

### Wave B2 — B1에 재현 가능한 event/module이 있을 때만

- GSE93535: 17.36 GiB
- 또는 GSE240226: 26.04 GiB

SIPS rescue가 중심이면 GSE93535, UVA rescue가 중심이면 GSE240226를 선택한다. 둘을 동시에 받을 필요는 없다.

### 보류

GSE109700, GSE179848, GSE226189 raw는 processed gene-level 결과로 현재 질문에 충분하고 저장공간 비용이 크다. Figure 3 후보가 명확해진 뒤 junction read로 꼭 확인해야 할 경우에만 받는다.

## 4. raw보다 먼저 필요한 내부 자료

1. 정상 `iPSC_CM_24h.R1.fastq.gz`와 제공자 checksum
2. 로컬 library strandedness/sample sheet
3. extracellular miRNA raw count 또는 FASTQ와 replicate 정보
4. antibody-array raw spot intensity와 blank/duplicate 정보

특히 1번이 없으면 Tier B raw를 받아도 Repro-CM-specific splicing이라는 중심 질문을 완성할 수 없다.

## 5. 감사 파일

- project summary: `public_data_tierA/derived/audit/tierB_project_size_summary.tsv`
- run-level accession/byte: `public_data_tierA/derived/audit/tierB_run_level_inventory.tsv`
- subset size: `public_data_tierA/derived/audit/tierB_priority_subset_size_summary.tsv`
- 조회한 ENA metadata snapshot: `public_data_tierA/metadata/tierB_runinfo/`

