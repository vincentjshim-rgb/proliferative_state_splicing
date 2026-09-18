# Repro-CM 후속 연구: Tier A 실행 결과와 논문 스토리 재판정

작성일: 2026-09-11  
범위: AlphaGenome·공간 오믹스 제외, 신규 wet-lab 없이 가능한 범위  
상태: Tier A processed data 분석 완료, Tier B raw FASTQ 미다운로드

## 1. 최종 판단

현재 결과는 다음 두 주장을 지지하지 않는다.

1. 서로 다른 UVA 조건에 공통인 하나의 보편적 photoaging 전사 시그니처가 존재한다.
2. partial reprogramming의 전사 반응은 일반적으로 senescence의 반대 방향이다.

대신 자료가 지지하는 더 정확한 결론은 다음이다.

> 인간 섬유아세포의 UVA 반응은 조사량·반복 횟수·회복 시간에 따라 급성 손상, 누적 손상, 만성 적응 상태로 갈라지며, chronological aging, replicative/SIPS senescence, partial reprogramming도 하나의 “젊음–노화” 축으로 환원되지 않는다.

따라서 후속 논문의 가장 방어력 있는 중심 질문은 아래와 같다.

> Repro-CM 연관 반응은 어떤 UVA/aging 상태축을 실제로 되돌리며, 어떤 축에는 영향을 주지 않는가?

이 질문은 단순한 “회춘 입증”보다 약하지만 과학적으로 더 새롭고, 현재 공개자료가 드러낸 이질성을 정면으로 다룬다. 다만 로컬 Repro-CM RNA-seq의 iPSC-CM R1 파일이 손상되어 있어, **현재 단계에서는 Repro-CM 자체에 이 결론을 연결할 수 없다.**

## 2. Tier A 실행 및 QC 결과

| 항목 | 결과 |
|---|---:|
| manifest에 고정된 원본 파일 | 39개 |
| 검증된 압축 용량 | 399,545,781 byte, 381.04 MiB |
| 형식·압축 무결성 통과 | 39/39 |
| SHA-256 기록 | 39/39 |
| GEO series | 11개 |
| GEO sample metadata | 744 sample |
| 설치·확인한 R package | limma, edgeR, fgsea, AnnotationDbi, org.Hs.eg.db, metafor, readxl, data.table |

검증 기록은 `public_data_tierA/logs/validation.tsv`, checksum은 `public_data_tierA/logs/sha256.tsv`, 분석 역할과 caveat는 `public_data_tierA/analysis_registry.tsv`에 고정했다.

`public_data_tierA/logs/rejected/GPL16956_family.soft.gz.incomplete`는 GPL 전체 SOFT 다운로드 재시도 중 끊긴 472.6 MiB 파일이다. 분석 입력이 아니며 rejected 폴더에 격리했다. 대신 필요한 platform sequence table은 정상 파일로 확보했다. 이 파일은 사용자 확인 없이 삭제하지 않았다.

### 로컬 선행연구 자료의 중단 조건

`iPSC_CM_24h.R1.fastq.gz`는 전체 gzip CRC 검사에 실패한다. R2와 보이는 read 수가 같아도 원본 동일성을 보장하지 못하므로 사용 금지다. 이 파일을 정상 원본으로 교체하고 제공자 checksum을 확인하기 전에는 다음을 확정할 수 없다.

- Repro-CM 대 iPSC-CM 비교
- Repro-CM-specific gene rank
- Repro-CM-specific splicing event
- Figure 4–6에서 사용할 Repro-CM 수용세포 response

## 3. Figure 1 결과: “재현성”의 정확한 의미

### 직접 답변

Figure 1에 다른 public data를 쓰는 것은 맞다. 그러나 공개자료에는 Repro-CM 처치군이 없으므로 **Repro-CM 효과의 직접 재현**이 아니라 **독립 UVA biology에 대한 외적 일치성**을 검증한다.

### 공개 UVA 코호트 결과

| 비교 | 공유 유전자 | Spearman ρ | 해석 |
|---|---:|---:|---|
| GSE125429 만성 5주 UVA vs GSE240226 급성 3 h UVA | 11,658 | -0.0049 | 공통 전유전체 방향 없음 |
| GSE302943 누적 5일 UVA vs GSE125429 만성 UVA | 10,538 내 overlap | 0.0003 | 일치하지 않음 |
| GSE302943 누적 UVA vs GSE240226 급성 UVA | 10,538 내 overlap | 0.3469 | 뚜렷한 일치 |
| Reactome: GSE240226 vs GSE302943 | 공통 pathway | 0.4062 | pathway 수준에서도 일치 |
| Reactome: GSE125429 vs GSE240226 | 공통 pathway | -0.0358 | 일치하지 않음 |

GSE89005 sensitivity 분석도 같은 구조를 보였다.

- single UVA 6 h vs GSE240226 acute: gene ρ=0.184, Reactome ρ=0.457
- single UVA 6 h vs GSE302943 cumulative: gene ρ=0.237, Reactome ρ=0.186
- single UVA 6 h vs GSE125429 chronic: gene ρ=0.0075, Reactome ρ=0.026
- repeated UVA 24 h는 위 세 연구와 거의 일치하지 않거나 약한 반대 방향

따라서 Figure 1에서 모든 UVA를 합쳐 meta-signature를 만들면 논리 오류다. 주분석은 다음처럼 층화해야 한다.

1. acute/early injury: GSE240226와 GSE89005 single 6 h
2. cumulative injury: GSE302943
3. chronic/adapted state: GSE125429와 GSE89005 repeated arm

Repro-CM rank가 복구되면 이 세 축에 각각 투영해야 한다. 일부 축만 반대 방향이면 “photoaging 전반을 역전”이 아니라 해당 축만 완화한다고 써야 한다.

### GSE302943 label 오류

GEO sample title/treatment field와 제출 count filename의 UVA/control 배치가 서로 반대였다. 논문 Figure S6의 열 라벨과 보고된 MMP1 증가·COL1A1 감소·mitochondrial gene 감소를 재현하는 쪽은 filename 라벨이었다. 따라서 filename+논문 Figure S6를 primary로 고정했고 GEO field는 submission metadata 오류로 처리했다. 근거는 [원 논문](https://pmc.ncbi.nlm.nih.gov/articles/PMC13117788/)과 [GEO sample record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9114511)에 남아 있다.

### GSE89005 probe annotation

GPL16956에는 공식 gene annotation이 없어 58,944개 60-mer probe를 GRCh38/GENCODE v50에 재매핑했다. MAPQ≥30, aligned length≥55 nt, NM≤2, 단일 gene 조건을 통과한 probe는 13,903개였고 protein-coding gene용 probe는 5,957개였다. 회수율이 낮으므로 GSE89005는 sensitivity 자료로만 사용한다.

## 4. Figure 2 결과: aging, senescence, MPTR은 동일 축이 아니다

### 재현성이 높은 부분

senescence model 사이에는 비교적 안정적인 일치가 있었다.

| 비교 | gene-level ρ | Reactome ρ |
|---|---:|---:|
| GSE109700 early vs deep senescence | 0.828 | — |
| GSE109700 deep vs GSE191055 P27 | 0.279 | 0.429 |
| GSE109700 deep vs GSE93535 SIPS | 0.278 | 0.688 |
| GSE191055 P27 vs GSE93535 SIPS | 0.365 | 0.718 |

즉 replicative senescence와 stress-induced premature senescence를 하나의 넓은 축으로 모델링할 근거는 있다.

### chronological aging의 불일치

- GSE113957 성인 22–89세(n=97)와 GSE226189 22–89세(n=82): gene ρ=0.188이지만 Reactome ρ=-0.167
- GSE113957 adult age와 deep senescence: gene ρ=0.450, Reactome ρ=0.815
- GSE226189 age와 deep senescence: gene ρ=0.044, Reactome ρ=-0.245
- GSE179848 donor-level late passage(n=4 donor)는 GSE113957 및 senescence와 양의 방향이나, GSE226189와는 약함

GSE113957의 1–96세 전체를 쓰면 아동 발달이 adult aging에 섞이므로 primary 분석은 공통 범위 22–89세로 제한했다. GSE179848의 345 sample은 독립 표본이 아니므로 HC1–HC4 네 donor의 early–late 차이만 분석 단위로 사용했다.

이 결과는 chronological age score를 하나만 만들면 cohort-specific proliferation/culture effect를 biological aging으로 오인할 수 있음을 보여준다.

### MPTR day-specific 재분석

처음에는 GSE165177의 성공한 MPTR 13쌍을 day 10/13/15/17에 걸쳐 paired model로 합쳤다. 그러나 Gill 논문의 핵심 window는 13일이므로 시간점별 donor 평균 효과를 다시 계산했다. 원 논문은 [eLife 논문](https://doi.org/10.7554/eLife.71624)과 [PMC 전문](https://pmc.ncbi.nlm.nih.gov/articles/PMC9023058/)에서 확인했다.

| day | matched pair | 고유 donor | 추론 수준 |
|---:|---:|---:|---|
| 10 | 3 | 3 | 소표본 descriptive |
| 13 | 4 | 2(O2, O3) | descriptive only |
| 15 | 2 | 1(O3) | descriptive only |
| 17 | 4 | 3 | 소표본 descriptive |

13일 MPTR의 logFC는 GSE113957 adult age와 거의 무관했다(ρ=0.0076). 그러나 late replicative state(ρ=0.163), early senescence(ρ=0.243), deep senescence(ρ=0.214), P27(ρ=0.374), SIPS(ρ=0.213)와는 같은 방향이었다.

13일 marker도 단순한 anti-senescence를 지지하지 않는다.

- SOX2 +8.18, POU5F1 +4.94
- MMP1 -2.35, COL1A1 +1.02: matrix/photoaging marker는 복구 방향
- MKI67 -1.63, LMNB1 -1.73, CDKN1A +0.99: 일반적인 anti-senescence 방향과 불일치
- NANOG/LIN28A 약 +0.11: 강한 pluripotency 전환은 아님

따라서 MPTR은 collagen/MMP 축을 회복시키면서 cell-cycle arrest/stress 축은 senescence와 비슷할 수 있다. “부분 리프로그래밍=senescence reversal”로 놓는 것은 오류다.

## 5. Figure 5 사전 점검: longevity와의 단순 수렴도 성립하지 않음

Gladyshev/Tyshkovskiy 2026 signature와 공개 UVA·aging·senescence effect를 screening했다. 이는 아직 Ensembl ortholog를 엄밀히 적용하기 전의 gene-symbol 기반 go/no-go 분석이다.

- 누적 UVA GSE302943은 age-adjusted rodent mortality와 양의 상관(ρ=0.194), rodent maximum-lifespan slope와 약한 음의 상관(ρ≈-0.041)을 보였다.
- acute UVA GSE240226은 age-adjusted mortality와 약한 양의 상관(ρ=0.091)이었다.
- chronic GSE125429는 aging/mortality/lifespan signature와 사실상 무관했다.
- deep senescence는 rodent age(ρ=0.186), mortality(ρ=0.177), human multi-tissue age(ρ=0.217)와 양의 상관이었다.
- 그러나 maximum-lifespan slope는 일부 senescence 자료와도 양의 상관을 보여, “aging의 정확한 반대=longevity”로 해석할 수 없었다.
- Tyshkovskiy 2019의 소수 curated gene set에 대한 GSEA는 현재 자료에서 FDR을 통과하지 못했다.

따라서 Figure 5는 Repro-CM이 longevity를 예측한다는 주 Figure로 쓰면 안 된다. 엄밀한 ortholog mapping과 local Repro-CM rank가 준비된 뒤 **분자적 수렴 여부를 묻는 보조 Figure**로만 유지한다. C. elegans 수명 연장을 사람 섬유아세포 전사 방향의 인과적 확인으로 쓰는 것도 금지한다.

## 6. 수정된 Figure 1–6 스토리

| Figure | 핵심 질문 | 현재 판정 | 허용되는 결론 |
|---|---|---|---|
| 1 | UVA 반응은 재현되는가? | 부분 통과 | UVA response는 regimen-dependent이며 acute/cumulative 축은 일부 재현됨 |
| 2 | age·senescence·MPTR이 한 축인가? | 분석 완료, 단일축 가설 기각 | 이 상태들은 분리해야 하며 ECM 회복과 senescence reversal도 다름 |
| 3 | Repro-CM-specific splicing/RBP가 있는가? | 중단 | iPSC R1 복구 및 replicate가 있는 raw cohort 검증 후 candidate 가능 |
| 4 | secretome/miRNA가 response를 설명하는가? | 중단 | 원 intensity/count가 있을 때 network-consistent candidate만 가능 |
| 5 | longevity biology와 수렴하는가? | 약한/혼합 | molecular convergence의 제한적 sensitivity 분석만 가능 |
| 6 | 어떤 결론이 held-out 자료에 남는가? | 선행 Figure 의존 | leave-dataset-out과 negative control을 통과한 축만 최종 score로 제시 |

권장 제목의 방향은 다음과 같다.

> Regimen-resolved transcriptomic triangulation separates photoaging, cellular senescence and partial reprogramming states in human fibroblasts

Repro-CM 파일이 복구되고 외부축에 선택적으로 일치하면 부제 또는 마지막 결과에 Repro-CM case study를 넣는다. 복구되지 않으면 Repro-CM을 제목·중심 결론에 넣지 않는 편이 안전하다.

## 7. 논리 오류 점검

| 위험한 추론 | 왜 오류인가 | 교정 |
|---|---|---|
| public UVA data가 Repro-CM 효과를 재현 | public cohort에 Repro-CM 처치가 없음 | external concordance라고 표현 |
| 여러 UVA 연구를 바로 합침 | dose·시간·반복·cell state가 서로 다르고 실제 방향도 불일치 | regimen별 study effect 후 계층화 |
| p가 매우 작으므로 효과가 큼 | 유전자 수가 많으면 작은 ρ도 매우 작은 p를 만듦 | ρ, 방향 일치율, dataset-level 재현성 중심 |
| 345 sample을 n=345로 처리 | 같은 donor의 longitudinal 반복 | donor difference n=4가 생물학적 단위 |
| MPTR 13일 n=4 donor | 실제로 O2/O3 두 donor의 두 experiment | n donor=2로 표시, descriptive only |
| aging과 senescence를 같은 score로 합침 | 두 age cohort와 pathway 결과가 충돌 | age, senescence, proliferation, ECM을 별도 축으로 유지 |
| lifespan signature와 반대면 수명 연장 | 조직·종·개입 맥락이 다르고 상관은 인과가 아님 | molecular convergence만 서술 |
| network path가 mediator를 입증 | 입력 edge와 seed가 문헌 편향을 가짐 | prioritized candidate, not mediator |
| 로컬 n=1 DEG의 FDR 사용 | biological variance를 추정할 수 없음 | effect rank만 생성, 통계추론은 public replicate에서 |
| 0 h vs 24 h를 치료 효과로 해석 | 자연 회복·시간·배지 효과가 혼재 | 24 h CM comparator끼리만 primary 비교 |

## 8. in silico-only 논문의 실제 선례와 차이

in silico 중심 논문은 Aging Cell에 실릴 수 있다. 예를 들어 다음 연구들은 공개자료 통합, 대규모 표본, 독립 검증 또는 새로운 모델을 주된 기여로 삼았다.

- 2,539명 혈액 transcriptome meta-analysis 후 별도 3,535명에서 module을 검증한 [de Magalhães group의 Aging Cell 연구](https://onlinelibrary.wiley.com/doi/10.1111/acel.12160)
- GTEx·TCGA와 senescence meta-signature로 조직별 aging–cancer–senescence 관계를 분리한 [Aging Cell 연구](https://onlinelibrary.wiley.com/doi/10.1111/acel.13041)
- 6,471 raw RNA-seq에서 QC 후 3,060 sample로 인간 aging transcriptome 모델을 만든 [Aging Cell 연구](https://onlinelibrary.wiley.com/doi/10.1111/acel.13280)
- 3,176 human skeletal-muscle methylome/transcriptome을 통합한 [Aging Cell meta-analysis](https://onlinelibrary.wiley.com/doi/full/10.1111/ACEL.13859)

반면 자주 인용되는 Gladyshev 사례는 순수 in silico 선례가 아니다.

- [Tyshkovskiy et al., Cell Metabolism 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6907080/)은 공개자료뿐 아니라 신규 RNA-seq, metabolomics 및 동물실험을 포함했다.
- [Tyshkovskiy et al., Nature 2026](https://www.nature.com/articles/s41586-026-10542-3)은 11,000개 이상의 transcriptome을 통합했지만 ITP mouse intervention RNA-seq를 새로 생산했다.

[Aging Cell 공식 범위](https://onlinelibrary.wiley.com/page/journal/14749726/homepage/forauthors.html)는 mechanistic aging biology와 중요성을 요구하며, 큰 새 진전이 없는 purely descriptive study는 범위 밖이라고 명시한다. 그러므로 “public DEG 여러 개를 합친 논문”으로는 부족하고, regimen-resolved state decomposition, 사전 고정한 통계, 독립 검증, 공개 코드라는 명확한 방법론적·생물학적 기여가 필요하다.

## 9. Aging Cell 게재 가능성

### 현재 그대로

낮다. 핵심 이유는 손상된 comparator, local n=1, 직접 Repro-CM replication 부재, Figure 3–4의 원자료 부재, 단순 longevity convergence 실패다.

### 다음 조건을 충족한 완전 in silico 원고

중간 수준까지 올릴 수 있다.

1. iPSC-CM R1 정상 파일과 checksum 확보
2. local analysis에서 p/FDR을 쓰지 않고 두 time-matched comparator에 일관된 rank만 정의
3. Figure 1을 regimen-resolved UVA state model로 전환
4. Figure 2에서 age/senescence/proliferation/ECM/pluripotency를 분리
5. 최소 두 raw cohort에서 exact splicing 또는 RNA-processing module을 재현
6. 모든 tuning과 held-out validation을 분리
7. secretome raw 값이 없으면 Figure 4의 mediator 언어를 포기
8. Figure 5를 수명 예측이 아닌 ortholog-aware molecular convergence로 제한
9. analysis container, code, manifest, sample inclusion table, full numerical source data 공개

최종적으로 Aging Cell 가능성을 만드는 것은 Figure 수가 아니라, 공개 데이터의 이질성을 이용해 기존의 단일 회춘축 해석을 수정하고 Repro-CM이 어느 상태축에 선택적으로 작용하는지를 정량화하는 새 개념이다.

## 10. 결과 파일

- Figure 1: `public_data_tierA/derived/figure1_public_uva/`
- Figure 2 aging/senescence: `public_data_tierA/derived/figure2_public_aging/`
- Figure 2 MPTR: `public_data_tierA/derived/figure2_mptr/`
- Figure 5 screening: `public_data_tierA/derived/figure5_longevity/`
- GEO/sample audit: `public_data_tierA/derived/audit/`
- 실행 코드: `scripts/`

