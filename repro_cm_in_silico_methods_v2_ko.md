# Repro-CM 후속 연구: Figure 1–6 상세 Methods v2

작성일: 2026-09-11  
원칙: AlphaGenome·공간 오믹스 제외, 연구별 effect-first 분석, biological unit 보존

## 1. 공통 분석 원칙

1. 서로 다른 연구의 expression matrix를 직접 합치지 않는다. 각 연구 안에서 effect와 uncertainty를 먼저 추정한다.
2. GEO sample 수가 아니라 donor/culture replicate/paired unit을 생물학적 n으로 사용한다.
3. 로컬 Repro-CM 자료는 조건당 n=1이므로 p-value·FDR을 만들지 않는다.
4. discovery에서 선택한 threshold와 parameter는 held-out data를 보기 전에 고정한다.
5. gene-level, pathway-level, direction consistency를 함께 보고하고 p-value만으로 결론을 내리지 않는다.
6. “rejuvenation”, “reversal”, “mediator”, “lifespan prediction”은 아래의 사전 조건을 통과할 때만 사용한다.

공통 gene statistic은 연구 내 logFC와 양측 p-value로부터 다음 signed z를 만든다.

\[
z_g = \mathrm{sign}(\mathrm{logFC}_g)\Phi^{-1}(1-p_g/2)
\]

연구 간 평가는 Spearman correlation, 같은 방향 비율, preranked GSEA/fgsea를 사용한다. 다중 비교는 Figure family별 Benjamini–Hochberg FDR로 통제한다.

## 2. Figure 1 — regimen-resolved UVA validation

### 2.1 자료와 contrast

| 자료 | biological unit | primary contrast |
|---|---|---|
| GSE125429 | HDF strain 4개 | paired chronic UVA − paired control |
| GSE240226 | replicate n=3/arm | acute UVA − control; 각 rescue − UVA |
| GSE302943 | replicate n=5/arm | cumulative UVA − control |
| GSE89005 | sham n=3, UVA n=2/arm | single 6 h, single 24 h, repeated 24 h를 각각 sham과 비교 |

### 2.2 전처리

- GSE125429: 제출된 log2 block만 primary로 사용하고 paired limma model에 donor strain을 포함한다.
- GSE240226/GSE302943: raw gene count를 edgeR filtering·TMM normalization 후 quasi-likelihood 또는 limma-voom으로 분석한다.
- GSE302943: group은 count filename과 논문 Figure S6를 기준으로 고정한다. GEO title을 사용한 반대 라벨 결과는 metadata-error sensitivity로 남긴다.
- GSE89005: 60-mer probe를 GRCh38/GENCODE v50에 minimap2로 재매핑한다. MAPQ≥30, aligned length≥55, NM≤2, 단일 gene probe만 포함한다. limma로 각 regimen을 별도 분석한다.

### 2.3 가설 검정 순서

1. GSE125429와 GSE240226를 discovery로 사용하되 먼저 일치성을 검사한다.
2. gene 및 Reactome NES 상관이 사전 기준을 통과할 때만 공통 UVA meta-rank를 생성한다.
3. 통과하지 않으면 acute/cumulative/chronic 세 strata를 독립 축으로 유지한다. 현재 결과에서는 이 경로가 확정됐다.
4. GSE302943은 held-out cumulative validation, GSE89005는 platform sensitivity로 사용한다.
5. GSE240226 rescue effect가 injury effect와 반대인지 평가한다. 같은 UVA comparator를 공유하므로 두 effect의 상관 p-value는 보조 descriptive 결과로 표시한다.

### 2.4 Repro-CM anchor 연결

정상 iPSC R1을 확보한 뒤 24 h 세 조건의 log2(CPM+0.5)를 계산한다.

\[
d_H(g)=x_{Repro}(g)-x_{HDF}(g)
\]

\[
d_I(g)=x_{Repro}(g)-x_{iPSC}(g)
\]

두 비교의 부호가 같을 때만 다음 score를 부여한다.

\[
S(g)=\mathrm{sign}(d_H)\min(|d_H|,|d_I|)
\]

부호가 다르면 S(g)=0이다. 이 rank를 UVA 세 축에 각각 투영한다. direct replication이라는 표현은 사용하지 않는다.

## 3. Figure 2 — state-axis decomposition

### 3.1 chronological age

- GSE113957: healthy donor만 포함하고 primary는 22–89세 n=97이다. `expression ~ age/10 + sex + platform`으로 decade effect를 추정한다. 1–96세 전체는 developmental sensitivity다.
- GSE226189: 82 donor에서 `count ~ age/10 + sex`를 사용한다. passage/culture covariate가 충분히 보고되지 않았다는 한계를 명시한다.

### 3.2 senescence

- GSE109700: proliferating, early, deep 각 n=3. early−proliferating와 deep−proliferating를 별도 산출한다.
- GSE191055: HDF P4 n=4 대 P27 n=4의 publisher DE table을 사용한다. young/old skin n=1은 제외한다.
- GSE93535: Q, Q1201, SIPS, SIPS1201 각 n=3. SIPS−Q와 SIPS1201−SIPS를 분석한다. Tier A에서는 논문 Cuffdiff table의 비교 가능 gene universe만 사용하므로 결과는 held-out sensitivity다.
- GSE179848: Normal, Control, 21% oxygen, HC1–HC4만 포함한다. donor마다 마지막 두 timepoint 평균에서 처음 두 timepoint 평균을 빼고 네 donor difference를 one-sample moderated model에 넣는다. sample 345개를 n으로 사용하지 않는다.

### 3.3 partial reprogramming

- all-day model: successful transiently reprogrammed sample과 같은 donor·day·experiment의 negative control을 짝지은 13 pair에 `expression ~ pair + condition` limma model을 적용한다. 이는 generic MPTR response다.
- day-specific model: 각 pair의 transient−negative log2-RPM 차이를 만들고 같은 donor의 experiment replicate를 먼저 평균한다.
- day 13은 four pair이나 고유 donor가 O2/O3 두 명이므로 p/FDR 없이 mean logFC, donor SD, 방향만 보고한다.
- day 10/17은 donor 3, day 15는 donor 1이다. 어느 경우도 “몇 년 젊어짐”을 계산하지 않는다.

### 3.4 다축 score

최소 다음 축을 따로 유지한다.

- chronological age: GSE113957, GSE226189
- replicative/SIPS senescence: GSE109700, GSE191055, GSE93535, GSE179848
- proliferation/cell cycle
- ECM/photoaging: COL1A1, COL1A2, MMP family 및 curated ECM pathway
- SASP/inflammation
- apoptosis/stress
- fibroblast identity
- pluripotency

Rejuvenation-like라는 단어는 age와 senescence가 모두 반대이고, identity가 보존되며, pluripotency와 proliferation만으로 설명되지 않을 때만 허용한다. 현재 MPTR 결과에는 적용하지 않는다.

## 4. Figure 3 — alternative splicing과 RBP

### 4.1 시작 gate

- 로컬 iPSC-CM R1 gzip integrity pass
- R1/R2 전체 record 및 ID pair 일치
- strandedness 확인 또는 RSeQC 추정
- 동일 reference를 쓸 public raw cohort 최소 두 개 확보

gate 이전에는 로컬 splicing p-value를 계산하지 않는다.

### 4.2 alignment 및 event quantification

- reference: GRCh38 primary assembly, GENCODE v50
- QC: FastQC/MultiQC, fastp, RSeQC infer_experiment, junction saturation, gene-body coverage
- alignment: STAR two-pass, 동일 parameter로 모든 raw cohort 처리
- event discovery: LeafCutter intron clusters와 MAJIQ local splice variations를 primary/secondary로 병행
- annotation sensitivity: rMATS SE/A5SS/A3SS/MXE/RI event
- 최소 junction support와 PSI missingness cutoff는 결과를 보기 전 고정

### 4.3 로컬 candidate score

로컬 조건당 n=1이므로 Repro−HDF와 Repro−iPSC의 ΔPSI가 같은 방향이고 두 비교 모두 |ΔPSI|≥0.10인 event를 candidate로 한다. 0.05와 0.15는 sensitivity cutoff다. 이것은 통계적 DEG가 아니라 effect-size candidate다.

### 4.4 public validation

- first wave: GSE191055 P4/P27와 GSE302943 UVA/control
- second wave: GSE93535 SIPS/rescue 또는 GSE240226 UVA/rescue
- event coordinate는 동일 GRCh38 annotation으로 재처리해 직접 맞춘다.
- 정확히 같은 event가 없을 때만 gene-level splicing module을 보조로 쓴다.
- replicate cohort에서는 donor-aware beta-binomial 또는 tool-native model로 ΔPSI, CI, FDR을 산출한다.
- 두 독립 cohort에서 같은 방향이거나, 하나의 discovery와 하나의 held-out에서 재현될 때만 robust event로 부른다.

### 4.5 RBP triangulation

1. event 인접 intron/exon에서 sequence motif enrichment
2. ENCODE eCLIP target overlap
3. 해당 RBP expression 또는 regulon activity 방향
4. RBP degree와 expression을 맞춘 permutation

세 근거가 일치해도 candidate regulator이며 mediator로 부르지 않는다.

## 5. Figure 4 — secretome/miRNA network

### 5.1 필수 내부 입력

- antibody array background-corrected spot intensity, blank, duplicate, exposure 정보
- extracellular miRNA raw FASTQ 또는 sample-level raw count
- biological replicate, batch, CM collection day와 pooling 정보

PDF heatmap 또는 선택된 factor 이름만으로 abundance weight를 만들지 않는다. 원수치가 없으면 unweighted exploratory network로 강등한다.

### 5.2 network 구성

- protein ligand/receptor와 intracellular signaling: version 고정 OmniPath human directed/signed edges
- pathway annotation: Reactome v97
- miRNA target: experimentally supported high-confidence edge만 별도 signed layer로 구성
- 수용 HDF에서 발현되지 않는 receptor/target 제거
- GSE226189와 GSE179848에서 보존되는 fibroblast coexpression/context edge를 보조 가중치로 사용

### 5.3 scoring과 null

- seed-to-response signed personalized PageRank 또는 random walk with restart
- response target은 discovery Repro rank, 평가 target은 held-out UVA/senescence axis
- AUROC, AUPRC, preranked enrichment를 모두 보고
- degree-matched seed permutation, edge rewiring, leave-one-ligand-out, leave-one-miRNA-family-out 수행
- tuning cohort와 평가 cohort를 분리

TIMP2·miR-302 계열은 상위에 남을 경우 network-consistent candidate로만 제시한다.

## 6. Figure 5 — longevity/mortality convergence

### 6.1 reference

- Tyshkovskiy 2019: common lifespan-extension intervention gene/pathway 및 effect-associated signature
- Tyshkovskiy 2026: rodent chronological age, mortality, age-adjusted mortality, maximum lifespan, human multi-tissue age, module table

### 6.2 최종 분석

- Ensembl Compara에서 one-to-one ortholog만 primary로 매핑한다.
- mouse symbol 대 human symbol의 단순 대문자 일치는 screening에만 허용한다.
- aging/mortality와 maximum-lifespan slope를 각각 독립 outcome으로 다룬다. 한 축의 부호를 뒤집어 다른 축으로 정의하지 않는다.
- Repro-CM rank의 signed correlation, up/down gene-set fgsea, module score를 계산한다.
- tissue, species, intervention을 하나씩 제외하는 leave-one-group-out 결과가 유지돼야 molecular convergence로 부른다.
- lifespan extension 또는 mortality reduction을 예측한다고 쓰지 않는다.

현재 screening에서 단순 lifespan convergence가 강하지 않았으므로 Figure 5는 secondary/exploratory로 사전 지정한다.

## 7. Figure 6 — validation과 최종 candidate lock

### 7.1 후보 선택은 discovery에서만

각 후보는 다음 evidence vector를 가진다.

- Repro 두 comparator의 같은 방향
- 적어도 한 UVA stratum에서 injury-opposing 또는 rescue-concordant
- senescence/age 축별 방향
- fibroblast identity 보존
- splicing/RBP 근거
- network stability
- longevity/mortality 보조 근거

가중치는 held-out data를 보기 전에 고정한다.

### 7.2 robustness

- leave-one-dataset-out
- leave-one-donor-out가 가능한 연구에서 donor influence
- effect cutoff와 gene-set size 변화
- platform 제외: array-only, RNA-seq-only
- immortalized BJ-5ta 제외
- chronic UVA 제외/포함
- degree-matched network null
- random gene-set 및 label permutation

### 7.3 held-out 성공 기준

- 사전 고정된 방향이 held-out cohort에서 일치
- effect size가 discovery CI 또는 prediction interval과 양립
- 한 dataset 제거로 결론 부호가 바뀌지 않음
- negative-control pathway/cell state에서 같은 성능이 나오지 않음

## 8. 재현성 산출물

- 원본 URL, byte, SHA-256, acquisition date manifest
- sample inclusion/exclusion table과 이유
- 모든 model formula 및 contrast table
- R/Python environment lock 또는 container
- Figure source-data TSV
- 분석 protocol version과 deviation log
- code repository와 archival DOI

Aging Cell은 제출 시 Data Availability Statement를 요구하고, 가능한 경우 script와 분석 artefact 공개를 권고한다. 공식 요구사항은 [Aging Cell author guidelines](https://onlinelibrary.wiley.com/page/journal/14749726/homepage/forauthors.html)을 따른다.

