# Repro-CM stress–senescence decoupling 연구: 상세 분석 방법 v3

작성 기준일: 2026-09-11  
이 문서는 실제 수행된 분석과 향후 다운로드 승인 후 수행할 held-out 분석을 구분한다.

## 1. 전체 설계

분석 질문은 다음 세 단계로 고정한다.

1. UVA 손상 HDF에 대한 Repro-CM 24 h 반응이 HDF-CM 및 iPSC-CM 반응과 다른가?
2. 이 Repro-CM-specific response가 독립 public UVA 및 senescence signature와 어떤 방향으로 연결되는가?
3. 그 연결은 단순 transcriptional reversal, broad rejuvenation, longevity, 또는 stress–senescence decoupling 중 어느 모델과 가장 잘 맞는가?

Local 자료는 조건당 n=1이므로 discovery anchor로만 사용한다. 통계적 일반화는 독립 public study와 expression-matched permutation에서 얻는다.

## 2. Local mRNA-seq 입력과 QC

### 2.1 입력

- HDF-CM: paired-end trimmed FASTQ
- Repro-CM: paired-end trimmed FASTQ
- iPSC-CM: 사용자가 지정한 유효 raw R1에서 생성된 paired trimmed FASTQ
- 원본 iPSC R1: `/home/shim/Downloads/project/project_shim/RNA_seq/00.data/iPSC_CM_24h.R1.fastq.gz`

세 시료는 선행연구 설계상 UVA 조사 HDF를 각 CM으로 24 h 처리한 recipient-cell RNA-seq이다. CM 생산세포 자체의 transcriptome이 아니다.

### 2.2 무결성 및 FastQC

- raw/trimmed gzip stream 전체 검사
- paired read 파일 존재 및 STAR 입력 성공 확인
- FastQC에서 per-base sequence quality, GC content, adapter content, duplication을 기록
- duplication WARN/FAIL은 library complexity 경고로 기록하되, mapping/assignment와 함께 판단

### 2.3 STAR 정렬 및 gene count

기존 GENCODE v41 GRCh38 STAR index를 사용했다.

```text
STAR --runThreadN 12 \
  --genomeDir <GENCODE-v41-GRCh38-STAR-index> \
  --readFilesIn <R1_val_1.fq.gz> <R2_val_2.fq.gz> \
  --readFilesCommand zcat \
  --outSAMtype None \
  --quantMode GeneCounts
```

STAR `ReadsPerGene.out.tab`의 세 count column을 비교했다. reverse-stranded assignment가 84.62–86.21%, forward가 5.74–5.89%여서 column 4를 primary count로 선택했다. unstranded column 2는 sensitivity에만 사용했다.

기존 RSEM count는 STAR 입력 대비 비정상적으로 적은 count만 포함해 제외했다. 향후 원고 Methods에는 “RSEM을 사용했다”고 쓰면 안 된다.

## 3. Gene annotation, normalization, filtering

- Ensembl version suffix 제거
- `org.Hs.eg.db`의 ENSEMBL→SYMBOL mapping
- 하나의 symbol로 매핑되는 여러 Ensembl row는 count 합산
- symbol은 대문자로 통일
- `edgeR::calcNormFactors`로 TMM normalization
- non-log CPM 계산
- `CPM >= 1`인 시료가 3개 중 2개 이상인 유전자 유지

각 비교의 log2 효과는 pseudocount 0.5를 사용했다.

```text
dH(g) = log2(CPM_Repro(g)+0.5) - log2(CPM_HDF(g)+0.5)
dI(g) = log2(CPM_Repro(g)+0.5) - log2(CPM_iPSC(g)+0.5)
```

## 4. Conservative Repro-CM anchor

두 비교자의 방향이 같을 때만 유전자를 Repro-specific로 인정한다.

```text
if sign(dH) = sign(dI) and sign(dH) != 0:
    S(g) = sign(dH) × min(|dH|, |dI|)
else:
    S(g) = 0
```

signed minimum을 쓰는 이유는 어느 한 비교자에서만 큰 효과가 나온 유전자가 anchor를 지배하지 못하게 하기 위해서다.

Primary threshold:

- rank analysis: `S != 0`인 모든 comparator-consistent 유전자
- 표 또는 label용 high-effect set: `|S| >= log2(1.25)`
- stronger subset sensitivity: `|S| >= log2(1.5)` 또는 `|S| >= 1`

Local n=1이므로 edgeR/DESeq2의 group dispersion 및 DEG p-value를 계산하지 않는다.

## 5. Public effect signature 표준화

각 public contrast에서 gene별 logFC와 p-value를 얻고 다음 signed z를 계산했다.

```text
Z(g) = sign(logFC(g)) × Phi^-1(1 - p(g)/2)
```

극단 p-value는 machine minimum과 1 사이로 제한한다. 중복 symbol은 사전에 제거한다.

Primary UVA studies:

- GSE240226: acute UVA vs control
- GSE302943: cumulative UVA vs control

Primary senescence studies:

- GSE179848: donor-level late vs early replicative trajectory
- GSE109700: deep replicative senescence vs proliferating
- GSE191055: P27 vs P4
- GSE93535: SIPS vs quiescent

GSE109700 early senescence는 deep contrast와 같은 study이므로 primary meta-axis에는 중복 포함하지 않는다. chronic UVA 및 GSE89005 timepoint는 context/boundary sensitivity로만 둔다.

## 6. Dataset별 signature association

Local `S(g)`와 각 public `Z(g)`의 Spearman rho를 계산한다. 효과크기 자체의 확인을 위해 local S와 public logFC의 Spearman rho도 별도 기록한다.

### 6.1 Expression-matched permutation

Local 세 시료의 평균 log2 expression을 decile로 나눈다. 각 decile 안에서 local score rank를 5,000회 섞어 null rho를 만든다.

```text
p_empirical = (1 + #{|rho_null| >= |rho_observed|}) / (5000 + 1)
```

모든 dataset contrast의 empirical p에 BH correction을 적용한다. 이 검정은 평균발현량에 따른 검출편향을 보존한다.

## 7. Meta-axis와 decoupling index

### 7.1 Dataset 내부 rank normalization

서로 다른 platform과 statistic scale을 직접 평균하지 않는다. 각 study statistic을 inverse-normal rank로 바꾼다.

```text
R_k(g) = Phi^-1((rank_k(g)-0.5)/n_k)
```

### 7.2 Axis 정의

```text
UVA_meta(g) = mean[R_GSE240226(g), R_GSE302943(g)]

SEN_meta(g) = mean[R_GSE179848(g), R_GSE109700(g),
                   R_GSE191055(g), R_GSE93535(g)]
```

공통 유전자 중 local comparator-consistent `S != 0`만 primary test에 사용한다.

### 7.3 중심 statistic

```text
rho_U = Spearman(S, UVA_meta)
rho_S = Spearman(S, SEN_meta)
D = rho_U - rho_S
```

`D > 0`은 local response가 UVA 쪽과 더 같은 방향이고 senescence 쪽과 더 반대 방향임을 뜻한다.

Local 평균발현량 decile 안에서 S를 10,000회 섞어 D의 one-sided empirical p를 계산한다.

```text
p_D = (1 + #{D_null >= D_observed}) / (10000 + 1)
```

이 가설은 현재 자료를 본 뒤 정교화되었으므로 discovery p로 표기한다. 새 held-out 자료부터 confirmatory test로 전환한다.

### 7.4 UVA-adjusted senescence

```text
SEN_residual = residuals[lm(SEN_meta ~ UVA_meta)]
rho_residual = Spearman(S, SEN_residual)
```

`rho_residual < 0`이면 UVA–senescence 공통 성분을 제거해도 local anti-senescence 관계가 남는 것이다.

## 8. Rescue residual analysis

GSE240226의 두 rescue contrast를 사용한다.

- Maifuyin vs UVA
- succinate vs UVA

각 rescue signed-z rank에서 같은 연구의 UVA-injury rank로 설명되는 성분을 제거한다.

```text
RESCUE_residual = residuals[lm(rankZ(RESCUE) ~ rankZ(UVA_injury))]
rho = Spearman(S, RESCUE_residual)
```

Local expression decile 안에서 5,000회 permutation을 시행하고 두 rescue test에 BH correction한다. 두 contrast는 같은 study에서 나온 것이므로 완전히 독립이라고 표현하지 않는다.

## 9. Robustness 분석

### 9.1 Leave-one-study-out

6개 primary study를 하나씩 제외하고 UVA/SEN meta-axis 및 D를 다시 계산한다. 모든 경우 `D > 0`이어야 방향 안정성을 통과한다.

### 9.2 Comparator sensitivity

- `dH`만 이용한 Repro-vs-HDF score
- `dI`만 이용한 Repro-vs-iPSC score
- signed-min AND anchor

세 경우 각각 D를 계산한다.

### 9.3 Gene-program removal

Reactome GMT를 사용해 아래 gene union을 정의하고 제거 후 D를 다시 계산한다.

- ribosome/translation/rRNA/NMD/EIF2/interferon/RNA-processing
- cell-cycle/mitosis/DNA-replication/G1-S/G2-M
- 두 군의 합집합
- gene-symbol prefix `RPL`, `RPS`, `MRPL`, `MRPS`, `MT-`
- local |score| 상위 250 및 상위 500

### 9.4 Technical read thinning

각 local gene count를 독립적으로 `Binomial(count, 0.5)`에서 다시 표본추출하고 TMM, filtering, anchor, D를 100회 재계산한다.

기록 지표:

- full anchor와 Spearman rho
- absolute-score top-250 overlap
- rho_U, rho_S, D의 2.5/50/97.5 percentile

이것은 sequencing-depth sensitivity이지 biological bootstrap이 아니다.

## 10. Reactome pathway 분석

- Reactome GMT 사용
- local score 또는 public signed z를 decreasing order로 정렬
- `fgseaMultilevel`
- `minSize=15`, `maxSize=500`, local analysis `eps=1e-10`
- BH FDR 사용

Dataset마다 얻은 pathway NES를 pathway 이름으로 교집합 결합한다. 2개 UVA NES 평균과 4개 senescence NES 평균을 만든 뒤 local NES와 Spearman 상관한다.

Primary pathway statistic은 553개 공통 pathway 전체의 상관이다. figure의 8개 selected pathway는 해석 예시이며 독립 confirmatory test가 아니다.

중복 Reactome pathway 때문에 pathway-level p-value는 보수적으로 해석하고, 최종본에서는 leading-edge gene Jaccard overlap으로 중복 module을 축약한다.

## 11. miRNA QC와 통합

### 11.1 Count 정리

miRDeep2 output은 하나의 mature miRNA가 여러 precursor row에 반복될 수 있다. 같은 mature miRNA의 `read_count` 최대값을 취해 중복합산을 피한다. 시료별 total mature count로 CPM을 계산한다.

Family 정의:

- miR-302/367: `miR-302a/b/c/d/e/f` 및 `miR-367`
- miR-371/372/373
- let-7 family

miRNA score도 mRNA와 같은 signed-min 규칙으로 기술한다. n=1이므로 differential-expression p-value를 계산하지 않는다.

### 11.2 Validated target test

OmniPath archive에서 다음 조건을 만족하는 edge만 사용한다.

- `mirnatarget=True`
- source가 대상 reprogramming miRNA family
- source database에 `miRTarBase` 포함
- PubMed reference 존재

family별 target gene set을 만든 뒤 local mRNA rank, UVA meta, senescence meta에서 fgsea한다. Repro-enriched miRNA의 recipient target은 negative NES가 예상되지만, 현재 miR-302/367 결과는 FDR 0.399로 실패했다. 따라서 causal network를 만들지 않는다.

## 12. Age, longevity, partial-reprogramming specificity

아래 reference와 local S의 Spearman association 및 expression-matched permutation을 계산했다.

- Tyshkovskiy 2026 rodent chronological age, mortality, age-adjusted mortality, maximum lifespan, human multi-tissue age
- Tyshkovskiy 2019 longevity-intervention/lifespan gene sets
- GSE113957 및 GSE226189 fibroblast chronological age
- GSE165177 maturation-phase transient reprogramming

가설방향과 맞지 않거나 서로 재현되지 않는 결과는 음성결과로 보존한다. multiple reference를 시험한 뒤 일부만 선택해 “rejuvenation”으로 주장하지 않는다.

## 13. Figure별 통계 표기 원칙

- Forest point: observed Spearman rho
- Forest horizontal bar: confidence interval이 아니라 expression-matched **null 95% interval**이라고 명시
- `p`와 `FDR`: gene-signature association의 permutation statistic
- local anchor panel: descriptive, `n=1 per condition`을 figure 안에 표시
- pathway NES: fgsea competitive enrichment; biological replicate p-value가 아님
- study-level 결과: accession, contrast, sample 수, primary/sensitivity 역할 표시
- 모든 panel에서 exact n gene/pathway와 seed 기록

## 14. Held-out 데이터 분석 전 고정할 사항

새 파일을 열기 전에 다음을 문서와 코드에 고정한다.

1. GSE297233 primary contrast: OSK dox vs no-dox day 4; 변형 OSK는 sensitivity
2. GSE307377 primary contrast: old vs young, sex를 covariate로 포함
3. GSE116968 primary rescue: Pre-Red vs UV를 1 h와 4 h에서 각각 계산; timepoint meta는 random-effects 또는 equal-weight rank meta
4. GSE149694: day 3 fibroblast를 reference로, day 7/day 13의 trajectory를 분리
5. GSE306957: senescent control vs proliferating control로 IFN/SASP axis를 정의하고, pathway perturbation은 별도 contrast
6. primary outcome: 방향과 D 변화
7. secondary outcome: pathway NES 및 leading-edge overlap
8. FDR family를 accession별이 아니라 전체 held-out primary tests에 대해 정의

이미 결과를 본 dataset은 held-out이라고 다시 부르지 않는다.

## 15. 재현성 파일

실행 스크립트:

- `scripts/run_local_repro_anchor.R`
- `scripts/run_local_longevity_validation.R`
- `scripts/run_decoupling_validation.R`
- `scripts/audit_local_mirna.R`
- `scripts/run_mirna_target_integration.R`
- `scripts/make_manuscript_diagnostic_figures.R`

주요 결과:

- `public_data_tierA/derived/local_repro_anchor/`
- `public_data_tierA/derived/decoupling_validation/`
- `public_data_tierA/derived/local_mirna_audit/`
- `public_data_tierA/derived/manuscript_diagnostic_figures/`
- `public_data_tierA/derived/evidence_scorecard.tsv`

모든 최종 run에서 R `sessionInfo()`와 random seed를 보존한다. 최종 원고 전에는 하나의 master script/Makefile로 실행순서를 고정하고, input SHA256 manifest와 모든 figure source TSV를 함께 공개한다.
