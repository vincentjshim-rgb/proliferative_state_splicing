# Repro-CM 후속 연구: 완전 in silico 논문 스토리와 Figure 1–6 설계

작성일: 2026-09-10  
검토 대상: JTE proof, “Secretomes from the transient iPSC reprogramming phase reverse cellular photoaging and extend Caenorhabditis elegans longevity”  
범위: AlphaGenome 및 공간 오믹스 제외. 신규 실험 없이 수행 가능한 후속 연구를 설계한다.

> 업데이트: 이 문서는 분석 전 가설 설계본이다. 2026-09-11 Tier A 결과에서 universal UVA signature와 단일 rejuvenation axis가 지지되지 않았으므로, 실제 해석은 `repro_cm_tierA_execution_and_story_revision_ko.md`를 우선한다.

## 1. 결론부터

완전한 in silico 후속 논문은 가능하다. 다만 현재 네 개 RNA-seq 라이브러리를 다시 DEG·GO 분석하는 수준으로는 Aging Cell의 기준에 미치기 어렵고, 선행논문의 Figure 5–6과 중복될 위험도 크다.

가장 타당한 새 논문의 중심 질문은 다음과 같다.

> Repro-CM과 연관된 수용세포 전사·RNA-processing 상태가 독립적인 인간 섬유아세포의 UVA 손상, 세포노화, 연령 증가, 부분 리프로그래밍 자료에서 일관된 방향성을 보이며, 그 신호가 알려진 장수 관련 전사 프로그램과도 수렴하는가?

권장 중심 주장은 다음처럼 제한한다.

> Cross-study transcriptomic triangulation identifies a Repro-CM-associated RNA-processing program that is anti-concordant with fibroblast photoaging and senescence while preserving fibroblast identity.

한국어로는 다음 의미다.

> 여러 독립 연구를 교차 검증한 결과, Repro-CM 연관 RNA-processing 프로그램은 섬유아세포 정체성을 잃지 않으면서 광노화·세포노화 방향과 반대되는 분자 상태를 보인다.

이 문장은 “Repro-CM이 노화를 역전한다”, “특정 miRNA/TIMP2가 효과를 매개한다”, “수명을 연장한다”보다 좁지만, 현재 자료가 실제로 지지할 수 있는 가장 강한 주장이다.

## 2. Figure 1에 대한 직접 답변

### 공개 데이터로 효과의 재현성을 보이는 것이 맞는가?

방향은 맞지만, 정확한 표현은 “효과 재현성”이 아니라 “독립적 외부 일치성 또는 일반화 가능성”이어야 한다.

Repro-CM을 동일한 방법으로 제조하고, 동일 UVA 조건에서 HDF-CM·iPSC-CM·vehicle과 비교한 독립 공개 데이터가 있어야 Repro-CM 처치 효과의 직접 재현이라고 부를 수 있다. 현재 확인한 공개자료에는 그런 독립 코호트가 없다.

따라서 공개 데이터로 검증할 수 있는 것은 다음이다.

- 독립 UVA 자료에서 나타나는 손상 방향과 Repro-CM 연관 방향이 반대인가?
- 독립적인 UVA 보호·회복 처치가 만드는 방향과 Repro-CM 연관 방향이 같은가?
- 독립적인 연령·세포노화 신호와 반대이며, 단순 증식이나 탈분화 신호로 설명되지 않는가?
- 한 데이터셋이나 한 플랫폼에만 의존하지 않고 여러 연구에서 방향이 유지되는가?

논문과 Figure 제목에는 “replication of the Repro-CM effect”를 쓰지 않고 다음 중 하나를 쓴다.

- Independent external concordance of the Repro-CM-associated response
- Cross-cohort validation of photoaging-opposing transcriptional programs
- Generalizability across independent UVA and fibroblast-aging cohorts

## 3. 선행논문에서 무엇을 가져오고 무엇을 새롭게 해야 하는가

### 선행논문이 이미 제시한 것

JTE proof는 다음을 이미 보고한다.

- UVA 15 J/cm²를 받은 HDF에서 Repro-CM이 생존, scratch closure, apoptosis, senescence 관련 지표를 개선함
- 마우스 상처 회복과 C. elegans 평균 수명 증가
- 수용 HDF RNA-seq의 stress/proteostasis 및 cell-cycle 관련 변화
- extracellular miR-302/367 및 miR-371/372/373 enrichment, let-7 depletion
- 항체 array에서 TIMP2를 포함한 후보 단백질

또한 proof 자체가 omics 연관성은 개별 miRNA·표적·단백질의 필요성이나 충분성을 입증하지 못한다고 명시한다.

### 후속 논문에서 반복하면 안 되는 것

- 같은 네 RNA-seq sample의 PCA, DEG Venn diagram, GO plot을 새 발견처럼 다시 제시
- 같은 miRNA–mRNA 상관 네트워크를 약간 다른 database로 재작성
- 선행논문에 나온 C. elegans 수명 결과를 새 외부 검증으로 다시 계산
- N=1 비교에서 edgeR p-value 또는 FDR을 새로운 통계 근거로 사용
- 후보 TIMP2, miR-302 또는 splicing event를 in silico 연결만으로 mediator라고 선언

### 새 논문의 독립적 기여

새 논문의 주된 증거는 공개 다중 코호트에서 나와야 한다. 로컬 네 sample은 “case-study anchor”이자 순위 시그니처 생성 자료이고, 다음이 진짜 신규 기여가 된다.

1. UVA–senescence–chronological aging–partial reprogramming을 분리하는 사전 정의된 다축 점수
2. 연구별 effect size를 먼저 산출한 뒤 meta-analysis하는 재현 가능한 프레임워크
3. RNA splicing/RBP 층의 독립 raw-data 검증
4. secretome/miRNA에서 수용세포 반응으로 이어지는 fibroblast-context network model
5. 장수·사망률 시그니처와의 종간 수렴성
6. discovery와 held-out validation을 분리한 후보 우선순위화

## 4. 논문 전체 스토리

논문의 논리 흐름은 아래와 같다.

| 단계 | 앞 Figure가 남기는 질문 | 다음 Figure의 답 |
|---|---|---|
| Figure 1 | Repro-CM 연관 상태는 UVA 손상과 실제로 관련 있는가? | 독립 UVA 자료에서 손상과 반대이고 보호 처치와 같은 방향인지 검증 |
| Figure 2 | 그것이 단순 stress 감소인가, aging/senescence와도 연결되는가? | 여러 연령·세포노화·부분 리프로그래밍 코호트에서 상태를 분해 |
| Figure 3 | 유전자 발현 외에 더 구체적인 RNA 조절층이 있는가? | 재현 가능한 splicing event와 RBP 프로그램 탐색 |
| Figure 4 | extracellular cargo가 그 RNA 상태에 어떻게 연결될 수 있는가? | ligand/receptor 및 miRNA/RBP를 문맥 특이 네트워크로 우선순위화 |
| Figure 5 | 이 상태가 피부 손상을 넘어 장수 생물학과 수렴하는가? | 독립 포유류 aging/mortality/longevity 시그니처와 종간 검증 |
| Figure 6 | 어떤 결과가 안정적이고 이후 실험에서 먼저 검증할 후보인가? | leave-dataset-out, negative control, held-out cohort로 최종 후보와 한계 제시 |

이 흐름에서 “연관”이 “매개”로 바뀌는 순간은 없다. 실험이 없는 논문이므로 마지막까지 candidate regulator, network-consistent predictor, molecular convergence라는 용어를 유지한다.

## 5. Figure 1–6 상세 설계

## Figure 1. 독립 UVA 자료에서 Repro-CM 연관 회복 방향을 검증

### 목적

선행논문에서 얻은 Repro-CM 연관 수용세포 변화가 특정 한 sample 또는 분석 pipeline의 산물이 아니라 UVA biology와 방향적으로 일치함을 보인다.

### 권장 패널

- A: 선행논문과 후속 분석의 관계를 보여주는 설계도. 로컬 N=1 자료는 anchor, 공개자료는 독립 평가 자료로 명확히 분리
- B: Repro-CM 대 HDF-CM 및 Repro-CM 대 iPSC-CM의 sign-consistent effect-rank 시그니처
- C: 독립 UVA 코호트별 injury signature와의 directional GSEA
- D: 독립 UVA-protection 또는 rescue signature와의 directional GSEA
- E: rank–rank hypergeometric overlap 또는 effect-size 상관
- F: 연구를 하나씩 제외한 leave-one-study-out 결과 및 negative-control 결과

### 핵심 자료

- [GSE125429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125429): 네 개 독립 HDF strain의 만성 UVA1 대 비조사 paired array
- [GSE89005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89005): primary skin fibroblast의 single/repeated UVA, 6 h/24 h array
- [GSE240226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240226): HDF control, UVA, UVA+두 보호 처치의 RNA-seq, 각 군 반복
- [GSE302943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302943): UVA 대 control 각 n=5 RNA-seq. 단, hTERT-immortalized BJ-5ta라는 차이를 명시

### 결과를 지지한다고 판단할 조건

- Repro-CM 연관 rank가 최소 두 개의 독립 UVA injury signature와 유의하게 반대 방향
- GSE240226의 독립 UVA 보호 처치 방향과 같은 방향
- 특정 DEG cutoff 하나가 아니라 여러 rank/threshold에서 결론 유지
- Repro-CM 시그니처가 HDF-CM 또는 iPSC-CM comparator 시그니처보다 더 일관된 성능
- expression level과 gene length를 맞춘 permutation에서도 empirical FDR 통과

### 여기서 허용되는 결론

“Repro-CM-associated transcriptional response is externally concordant with attenuation of UVA injury.”

### 금지되는 결론

“Public data reproduce the Repro-CM treatment effect.”

공개자료에는 Repro-CM 처치가 없으므로 직접 재현이 아니다.

## Figure 2. rejuvenation-like state를 aging, senescence, proliferation, identity 축으로 분해

### 목적

Figure 1에서 보인 반-UVA 방향이 진짜 aging/senescence 관련 변화인지, 아니면 단순한 cell-cycle 변화, 스트레스 억제, 선택적 세포사멸 또는 fibroblast identity 상실인지를 구분한다.

### 권장 패널

- A: 분석할 상태축의 사전 정의: chronological aging, replicative senescence, UVA injury, fibroblast identity, pluripotency, proliferation, apoptosis, SASP, interferon
- B: 독립 연령 코호트에서 얻은 age effect와 Repro-CM 연관 rank의 meta-analytic concordance
- C: replicative 및 stress-induced senescence 코호트에서의 concordance
- D: MPTR time-course에서 Repro-CM 연관 상태가 어느 시점과 닮았는지
- E: fibroblast identity–pluripotency–proliferation 삼각 좌표
- F: 결과를 한눈에 보여주는 다축 radar 또는 coefficient plot. “몇 년 젊어짐” 수치는 표시하지 않음

### 핵심 자료

- [GSE113957](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113957): 건강한 133명, 1–94세 및 HGPS 10명의 dermal fibroblast RNA-seq
- [GSE226189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189): 22–89세 건강인 82명의 primary skin fibroblast RNA-seq
- [GSE179848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179848): 345 sample의 longitudinal fibroblast replicative lifespan 자료. 분석 단위는 sample 수가 아니라 donor와 longitudinal trajectory
- [GSE191055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191055): HDF passage 4 대 passage 27, 각 n=4
- [GSE109700](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109700): proliferating, early senescence, deep senescence LF1, 각 n=3
- [GSE93535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93535): HDF quiescence 및 H2O2-induced premature senescence, rescue 조건
- [GSE165177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165177): 세 donor와 여러 시점·대조를 포함한 human MPTR RNA-seq

### rejuvenation-like라는 표현을 허용할 최소 조건

다음이 동시에 만족되어야 한다.

1. chronological age 및 senescence 방향과 반대
2. fibroblast identity score 보존
3. pluripotency score의 비정상 상승 없음
4. proliferation만으로 전체 효과가 설명되지 않음
5. apoptosis 또는 low-quality cell selection만으로 설명되지 않음
6. 적어도 두 개 이상의 독립 코호트에서 방향 유지

이 조건을 만족하지 못하면 “rejuvenation-like”도 제거하고 “UVA-recovery-associated”로 범위를 좁힌다.

### aging clock 사용 원칙

서로 다른 platform과 library protocol 사이에서 clock의 절대 출력값을 “30년 젊어짐”처럼 해석하지 않는다. 훈련 분포 밖의 N=1 sample에 적용한 절대 나이는 scale artifact일 가능성이 있다. 같은 연구 안의 상대적 방향과 여러 clock/signature의 일치만 보조 결과로 쓴다. 이 원칙은 aging clock 해석에 관한 [Aging Cell 논평](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.14377)과도 맞는다.

## Figure 3. Repro-CM-specific alternative splicing과 RBP 프로그램

### 목적

유전자 수준의 broad signature를 보다 구체적인 RNA-processing 층으로 좁히되, N=1에서 통계적 유의성을 발명하지 않는다.

### 권장 패널

- A: 로컬 FASTQ QC, junction saturation, insert size, strandedness, mapping summary
- B: Repro-CM 대 HDF-CM 및 Repro-CM 대 iPSC-CM에서 방향이 같은 candidate splicing events
- C: public senescence/aging/MPTR raw RNA-seq에서 동일 event 또는 동일 splicing module의 방향 검증
- D: event 주변 RBP motif enrichment
- E: ENCODE eCLIP target enrichment와 해당 RBP expression/activity의 삼각 검증
- F: isoform의 coding domain, reading frame 또는 NMD 가능성에 대한 제한적 annotation

### 로컬 자료에서의 통계 원칙

각 조건이 biological replicate n=1이므로 p-value와 FDR을 계산하지 않는다. Repro-CM 대 두 time-matched CM comparator에서 방향이 같고, 최소 효과 크기가 큰 event만 후보로 만든다.

- 기본 후보 기준: 두 비교 모두 같은 방향
- 최소 효과: absolute delta PSI 0.15, 엄격 민감도 분석은 0.20
- 최소 junction support: 총 informative junction count 20 이상
- read/bootstrap 변화에 따른 PSI rank stability 보고
- 15J_0hr 비교는 시간과 처치가 섞이므로 보조 그림에만 둠

### 독립 검증

가능한 경우 GSE191055, GSE109700, GSE179848, GSE165177의 raw RNA-seq를 동일 reference와 pipeline으로 재처리한다. exact junction을 검증하려면 read length, paired/single-end, depth가 충분해야 한다. 짧은 single-end 자료는 gene-level 또는 exon-module 수준 보조자료로만 사용한다.

### Figure 3의 중단 규칙

독립 자료 두 개 이상에서 event 방향을 확인하지 못하면 splicing을 논문 제목과 중심 기전에서 뺀다. 이때 Figure 3은 transcript/RBP module의 탐색적 분석으로 축소하거나 supplement로 이동한다.

### 허용·금지 결론

- 허용: “RBP X and isoform Y are reproducible correlates of the Repro-CM-associated state.”
- 금지: “RBP X–isoform Y mediates Repro-CM rejuvenation.”

## Figure 4. secretome/miRNA에서 receptor–RBP–isoform으로 이어지는 문맥 특이 네트워크

### 목적

선행논문의 extracellular cargo 후보가 수용 HDF의 관찰된 RNA 상태와 연결될 수 있는지를 투명하고 검증 가능한 네트워크 모델로 평가한다.

### 왜 black-box virtual cell보다 이 방법이 적합한가

현재 자료는 CM perturbation별 biological replicate와 dose/time series가 없어 대규모 virtual cell model을 훈련하거나 공정하게 검증할 수 없다. 반면 ligand–receptor seed와 fibroblast-context network를 이용한 random walk with restart는 입력, 경로, 성능 기준을 모두 공개할 수 있다.

따라서 이를 “virtual cell이 기전을 증명했다”고 표현하지 않고, “context-specific in silico perturbation model”이라고 부른다.

이 접근의 직접적인 방법론 참고 사례인 [McNeill et al., Aging Cell 2026](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.70683)은 extracellular ligand를 receptor seed로 연결하고 문맥 특이 network에서 random walk를 수행했다. 다만 그 연구는 독립적인 proteomics, RNA-seq와 기능 검증을 갖고 있었으므로, 본 연구가 in silico만으로 같은 수준의 인과성을 주장할 수는 없다.

### 이 Figure에 필요한 내부 원자료

- antibody array의 background-corrected spot intensity, technical replicate와 blank
- CM 종류별 extracellular miRNA raw FASTQ 또는 최소한 raw count matrix와 biological replicate 정보
- 각 sample의 batch, collection day, pooling 방식

PDF의 heatmap 이미지와 선별된 factor 이름만으로 weighted main analysis를 하면 안 된다. 수치 원자료가 없으면 Figure 4는 unweighted sensitivity analysis로 낮추고 TIMP2 중심 주장을 포기한다.

### 권장 패널

- A: secreted protein/miRNA → receptor/target → fibroblast network 분석도
- B: GSE226189 또는 GSE179848로 구축한 fibroblast-context network의 보존성
- C: ligand/receptor seed별 random-walk ranking과 관찰된 수용세포 response의 enrichment
- D: high-confidence miRNA target의 signed negative edge를 추가했을 때 성능 변화
- E: degree-matched random seed, network rewiring, leave-one-ligand-out null 결과
- F: TIMP2를 포함한 후보들의 순위와 불확실성

### 성능 기준

- 관찰된 Repro-CM 연관 rank를 reference로 AUROC, AUPRC, GSEA를 모두 보고
- network degree가 높은 gene을 자동으로 뽑는 현상을 degree-matched seed permutation으로 통제
- parameter를 고른 자료와 최종 성능을 평가하는 자료를 분리
- 한 ligand를 뺐을 때 순위가 완전히 무너지면 안정 후보로 간주하지 않음
- 서로 다른 network source 또는 edge threshold에서 rank stability 확인

### 허용·금지 결론

- 허용: “TIMP2 and selected miRNA families are prioritized as network-consistent candidates.”
- 금지: “TIMP2 is necessary or sufficient for the Repro-CM effect.”

## Figure 5. 포유류 aging, mortality, longevity intervention 시그니처와의 종간 수렴

### 목적

선행논문의 C. elegans 수명 결과를 재사용해 결론을 강화하는 대신, 독립적으로 구축된 포유류 수명·사망률 시그니처와 Repro-CM 연관 상태가 분자 수준에서 수렴하는지를 검증한다.

### 권장 자료

- [Tyshkovskiy et al., Cell Metabolism 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6907080/): 17개 수명연장 intervention, 77개 control–intervention comparison의 signature tables
- [Tyshkovskiy et al., Nature 2026](https://www.nature.com/articles/s41586-026-10542-3): 11,000개 이상 transcriptome, 25개 이상 조직, 4개 포유류에서 구축한 aging, mortality, lifespan 관련 계수와 module

### 권장 패널

- A: human–mouse–worm one-to-one ortholog mapping과 분석에 남는 gene 수
- B: aging 및 mortality coefficient와의 signed enrichment
- C: maximum lifespan 및 lifespan-extending intervention signature와의 enrichment
- D: tissue/study별 forest plot
- E: ortholog 수, 발현량, gene length를 맞춘 random null
- F: Repro-CM 후보 network module과 universal aging/mortality module의 overlap

### 중요한 해석 제한

- 사람 세포의 24 h CM 반응과 mouse lifespan intervention은 동일 현상이 아니다.
- 유전자 방향의 수렴은 수명 연장 기전을 입증하지 않는다.
- 선행논문의 C. elegans +43.45% 결과는 배경적 biological context일 뿐 독립 검증 자료로 다시 세지 않는다.
- “Repro-CM extends mammalian lifespan”은 절대 쓰지 않는다.

허용되는 표현은 “molecular convergence with longevity-associated programs”이다.

## Figure 6. 안정성 선택, held-out 검증 및 다음 실험 후보

### 목적

Figure 1–5의 서로 다른 증거를 순환논리 없이 통합하고, 어떤 후보가 임계값이나 특정 데이터셋 선택에 덜 민감한지를 보여준다.

### 후보 점수의 구성

후보별로 다음의 독립 증거를 정규화해 결합한다.

1. 로컬 time-matched Repro-CM specificity
2. 독립 UVA injury와의 anti-concordance
3. 연령·senescence와의 anti-concordance
4. fibroblast identity 보존 및 pluripotency artifact 부재
5. splicing/RBP 독립 재현
6. extracellular network propagation rank
7. longevity/mortality signature 수렴

한 dataset에서 얻은 정보가 ranking과 검증에 동시에 사용되지 않도록 discovery와 validation cohort를 사전에 나눈다.

### 권장 held-out 구성

- UVA hold-out: GSE302943 또는 GSE240226의 rescue arm 중 하나
- aging hold-out: GSE307377 raw 자료. 네 young donor와 다섯 old donor이지만 현재 GEO series에서 processed matrix 제공 여부를 다운로드 전에 재확인
- senescence hold-out: GSE93535 rescue comparison
- reprogramming hold-out: GSE237269 또는 GSE297233은 sample design과 replicate를 metadata 단계에서 확인한 뒤 채택

### 권장 패널

- A: discovery–validation 분리도
- B: 각 증거축의 candidate coefficient/rank matrix
- C: dataset bootstrap과 threshold multiverse에서의 stability-selection frequency
- D: held-out cohort 성능
- E: candidate-level LINCS/CMap genetic perturbation 또는 compound signature의 보조 검증
- F: 최종 5–10개 후보, 예측 방향, 향후 loss/rescue 실험 명세

### LINCS/CMap 사용 원칙

전체 L1000 파일을 처음부터 다운로드하지 않는다. 최종 후보가 정해진 뒤 candidate-level genetic overexpression/knockdown 및 제한된 compound signature를 조회한다. replicate는 MODZ 등으로 통합하고, 세포주별 방향 일치와 무관한 gene-set negative control을 함께 보고한다.

Figure 6의 결론은 “검증할 가치가 높은 반증 가능한 후보를 만들었다”이며 “in silico로 기전을 증명했다”가 아니다.

## 6. Results 문단의 연결 문장

### Results 1에서 2로

“Repro-CM-associated changes opposed independent UVA injury signatures and aligned with unrelated UVA-protective responses. We therefore asked whether this concordance reflected a broader aging-related state or merely acute stress attenuation.”

논리: UVA 일치성에서 곧바로 rejuvenation으로 넘어가지 않고, Figure 2에서 대안을 검증한다.

### Results 2에서 3으로

“The response was anti-concordant with multiple aging and senescence programs while retaining fibroblast identity and showed only partial overlap with proliferation and pluripotency. We next examined whether a reproducible RNA-processing layer distinguished this state.”

논리: identity와 proliferation을 통제한 뒤에만 rejuvenation-like라는 제한적 표현을 허용한다.

### Results 3에서 4로

“A subset of splice events and RBP programs showed directional support across independent datasets. We therefore tested whether extracellular cargo reported in Repro-CM could prioritize these downstream programs in a fibroblast-context network.”

논리: splicing 연관성을 발견했을 뿐 extracellular cargo가 이를 유발했다고 전제하지 않는다.

### Results 4에서 5로

“Network propagation nominated candidate extracellular regulators that recapitulated parts of the observed recipient-cell response. We next tested whether the resulting modules converged with independently derived mammalian aging, mortality, and longevity-intervention signatures.”

논리: network prediction을 causal regulator로 바꾸지 않는다.

### Results 5에서 6으로

“Because cross-species convergence alone cannot establish mechanism, we integrated only nonredundant evidence streams and evaluated candidate stability in datasets withheld from prioritization.”

논리: cross-species overlap의 과잉해석을 held-out validation으로 제어한다.

## 7. 반드시 피해야 할 논리적 오류

| 오류 | 왜 문제인가 | 교정 방법 |
|---|---|---|
| 공개 UVA 일치성을 Repro-CM 효과 재현이라고 부름 | 공개자료에 Repro-CM 처치가 없음 | external concordance/generalizability로 표현 |
| 로컬 N=1에서 DEG 또는 splicing FDR 제시 | biological variance를 추정할 수 없음 | effect-rank와 descriptive candidate만 사용 |
| 15J_0hr 대 24 h CM을 처치 효과로 해석 | 시간과 처치가 완전히 confounded | Repro 대 HDF-CM/iPSC-CM의 24 h 비교를 중심으로 사용 |
| young-like signature를 rejuvenation으로 등치 | 증식, 사멸, stress 감소, 탈분화 가능 | identity·pluripotency·proliferation·apoptosis 축 동시 평가 |
| senescence score 감소를 senomorphic 효과로 단정 | senescent cell 선택적 사멸일 수 있음 | apoptosis/viability 관련 signature와 독립 평가 |
| OSKM/MPTR 유사성을 안전한 회춘으로 해석 | fibroblast identity 소실과 pluripotency 가능 | identity 보존 조건을 필수 gate로 설정 |
| miRNA 표적 역상관을 mediation으로 해석 | database prediction과 상관은 인과가 아님 | high-confidence target, null test, candidate 표현만 사용 |
| network rank를 mechanism이라고 표현 | topology와 seed bias가 큼 | degree-matched null, rewiring, held-out benchmark |
| 종간 signature overlap을 lifespan mechanism으로 해석 | 조직·종·시간·개입이 다름 | molecular convergence로 제한 |
| 같은 네 sample을 새 독립 증거로 재사용 | 선행논문과 중복·salami slicing 위험 | 재분석임을 밝히고 공개 코호트를 주 증거로 설정 |
| 여러 연구의 normalized matrix를 합쳐 DE 수행 | platform/batch가 group과 섞임 | 각 연구 내부 effect를 구한 후 random-effects meta-analysis |
| single-cell의 cell 수를 n으로 사용 | donor pseudoreplication | donor를 biological unit로 사용 |
| 후보를 고른 자료에서 같은 후보를 검증 | 순환논리와 과적합 | discovery와 held-out cohort 분리 |
| cross-platform clock 값을 “몇 년”으로 해석 | out-of-distribution 및 scale 차이 | 상대 방향과 module score만 사용 |

## 8. 특히 제외하거나 낮은 가중치로 둘 자료

- GSE274955: young, middle, old 각각 한 donor이고 두 부위가 같은 donor에서 옴. cell 수가 많아도 연령별 biological n=1이므로 핵심 replication에 부적합
- GSE203034: 한 남성 donor를 두 연령에 측정한 자료. 흥미로운 longitudinal 사례지만 population-level 일반화 불가
- GTEx skin: 피부 조직의 cell-composition 변화가 섞이므로 cultured HDF의 핵심 검증이 아니라 보조 결과
- GSE125429와 GSE89005: microarray이므로 gene-level 방향에는 유용하지만 exact splicing 검증에는 부적합
- 짧은 single-end RNA-seq: junction-level splicing의 negative result를 evidence of absence로 해석하지 않음

## 9. in silico만으로 게재된 선례에서 배워야 할 점

### Aging Cell의 computational 선례

- [Shokhirev & Johnson, Aging Cell 2021](https://onlinelibrary.wiley.com/doi/10.1111/acel.13280)은 6,000개 이상의 public raw RNA-seq를 다시 처리하고 3,060개 high-quality sample로 조직·성별·건강상태를 아우르는 aging model과 공개 도구를 만들었다. 핵심은 sample 수가 아니라 동일한 preprocessing, 외부 일반화, 재현성 한계 공개다.
- [Deng et al., Aging Cell 2023, SenoRanger](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13809)은 6개 senescence 연구 13조건을 GTEx, ARCHS4, Tabula Sapiens와 통합하고 consensus candidate와 software resource를 만들었다.
- [Santos et al., Aging Cell 2026](https://onlinelibrary.wiley.com/doi/10.1111/acel.70696)은 25,000개 이상의 공개 인간·마우스 transcriptome에 적용 가능한 entropy metric으로 광범위한 조직·질환 일반화를 보였다.

### Gladyshev 계열 연구에서 배울 점

- [Gross et al., Nature Aging 2026](https://www.nature.com/articles/s43587-026-01161-8)은 2,358개 longevity gene, human interactome, 6,442개 compound와 pAGE를 결합하고 알려진 성공·실패 intervention으로 계산 예측을 benchmark했다. 최종 산출물을 “falsifiable framework”로 표현했지 계산만으로 수명 연장을 증명했다고 하지 않았다.
- [Tyshkovskiy et al., Cell Metabolism 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6907080/)은 17개 lifespan-extending intervention과 77 comparison을 통합하고 leave-one-intervention-out 식의 일반화를 사용했다. 이 논문은 일부 신규 실험도 포함하므로 순수 public-data 논문의 선례로 과장해서는 안 된다.
- [Tyshkovskiy et al., Nature 2026](https://www.nature.com/articles/s41586-026-10542-3)은 11,000개 이상 transcriptome, 25개 이상 조직, 네 포유류를 통합하고 within-dataset centering, mixed model, leave-one-dataset/tissue-out 검증을 사용했다.

### 본 연구에 주는 현실적 교훈

in silico-only 자체가 탈락 이유는 아니다. 그러나 상위 저널의 계산 논문은 대체로 다음 중 하나 이상을 갖는다.

- 대규모 다중 코호트
- 다른 연구에도 적용 가능한 새로운 방법 또는 resource
- 명확한 discovery–validation 분리
- 성공뿐 아니라 실패 또는 neutral intervention까지 포함한 benchmark
- 코드·데이터·모형의 완전 공개

Repro-CM 하나의 N=1 dataset에 GO와 network를 붙이는 것은 이 기준과 거리가 멀다. 반대로 “aging reversal + senescence reversal + fibroblast identity preservation”의 세 축을 일반화 가능한 평가 프레임워크로 만들고 여러 독립 intervention에 benchmark한다면 계산 논문으로서의 독창성이 생긴다.

## 10. Aging Cell 게재 가능성 평가

[Aging Cell 저자 지침](https://onlinelibrary.wiley.com/page/journal/14749726/homepage/forauthors.html)은 fundamental/mechanistic aging을 중시하고, major advance가 아닌 순수 기술적·서술적 연구는 범위 밖이라고 명시한다. 또한 supporting data의 공공 repository 등록과 accession link를 요구하고, 가능하면 분석 script도 공개하도록 요구한다.

### 현재 네 sample 재분석만으로 투고

부적합하다.

- 모든 조건 biological n=1
- 0 h 대 24 h가 time-confounded
- 선행논문과 데이터·질문이 중복
- 인과 기전 검증 불가
- public accession과 원자료 공개가 proof에서 확인되지 않음

### Figure 1–6 전체 계획을 엄격히 수행한 경우

Aging Cell에 도전할 논리적 근거는 생긴다. 단, 여전히 경쟁적인 borderline 기획이다. 다음이 필수다.

- 로컬 원자료의 repository accession 및 명시적 secondary-analysis disclosure
- 손상된 iPSC-CM FASTQ 교체
- 독립 UVA, aging, senescence, MPTR 코호트의 사전 정의 및 study-level meta-analysis
- 다른 연구에도 쓸 수 있는 다축 metric 또는 network method
- held-out validation과 topology/threshold negative controls
- causal language 제거
- splicing이 재현되지 않으면 과감히 중심 주장 축소

### Nature Aging 수준

현재 Repro-CM N=1 anchor를 중심으로는 현실성이 낮다. Nature Aging을 목표로 하려면 Repro-CM 사례를 넘어 다수의 paracrine/CM intervention을 포괄하는 일반적 resource 또는 예측 모델을 만들고, 독립 성공·실패 intervention에서 강한 prospective benchmark를 보여야 한다. 실험이 없는 경우 방법론과 규모의 기준은 더 높다.

### 결과가 약할 때의 합리적 분기

- UVA 및 senescence 일반화는 강하나 splicing/network가 약함: transcriptomic meta-analysis 중심으로 축소
- splicing은 강하나 longevity 수렴이 약함: RNA-processing/photoaging 논문으로 축소
- aging보다 UVA에만 특이적: “rejuvenation”을 버리고 photo-damage recovery로 제목 변경
- 공개자료에서 방향이 불안정: 후속 논문을 중단하거나 hypothesis/resource note로 전환

## 11. 투고 전 비가역적 gate

다음 세 조건이 충족되기 전에는 본 분석을 본격 실행하거나 논문 제목을 확정하지 않는다.

1. iPSC_CM_24h.R1.fastq.gz의 정상 원본 재확보와 full gzip integrity 확인
2. 선행 연구의 RNA-seq 및 miRNA-seq public accession 확보
3. antibody-array 수치표와 extracellular miRNA count/raw 자료의 존재 여부 확인

조건 2가 불가능하면 Aging Cell의 repository 요건과 재현성 측면에서 치명적이다. proof correction이 아직 가능하다면 JTE 원고에 Data Availability와 accession을 추가하는 것이 우선이다.

## 12. 권장 제목 후보

가장 안전한 제목:

> Cross-study transcriptomic triangulation identifies a photoaging-opposing RNA-processing state associated with transient-reprogramming secretome

aging 결과가 두 개 이상 독립 코호트에서 강할 때:

> A fibroblast-context computational framework links transient-reprogramming secretome responses to photoaging, senescence, and conserved aging programs

splicing이 독립 검증될 때만:

> Conserved RNA-binding and splicing programs mark the photoaging-opposing response associated with transient-reprogramming secretome

피해야 할 제목:

- Repro-CM reverses human aging
- Virtual-cell discovery of the causal Repro-CM rejuvenation mechanism
- TIMP2/miR-302 mediates lifespan extension

## 13. 최종 판단

Figure 1은 공개 데이터로 수행하는 것이 맞다. 그러나 동일 Repro-CM 처치의 재현이 아니라 독립 UVA biology에 대한 외부 방향 검증이다.

Figure 2–6도 전산 분석으로 구성할 수 있지만, 논문의 강점은 “실험을 컴퓨터로 대체했다”가 아니라 “선행 연구의 제한된 N=1 omics를 독립적인 UVA·aging·senescence·partial-reprogramming·longevity 자료에 걸쳐 반증 가능하게 시험했다”는 데 있어야 한다.

가장 큰 리스크는 computational sophistication 부족이 아니라 다음 세 가지다.

- 같은 데이터의 중복 출판
- N=1에서의 거짓 통계적 확신
- 연관성에서 인과성·회춘·수명연장으로의 과도한 도약

이 세 가지를 명시적으로 통제하고, 외부 다중 코호트와 일반화 가능한 분석 프레임워크를 주된 신규 기여로 만들 때에만 Aging Cell 투고를 진지하게 고려할 수 있다.
