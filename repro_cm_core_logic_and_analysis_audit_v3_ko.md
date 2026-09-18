# Repro-CM 후속 in-silico 연구: 핵심 논리와 분석 점검 v3

작성 기준일: 2026-09-11  
범위: AlphaGenome·공간오믹스 제외, 신규 wet-lab 없이 가능한 후속 연구

## 1. 최종 판단

현재 데이터가 가장 강하게 지지하는 한 문장은 다음과 같다.

> **부분 리프로그래밍기 secretome에 노출된 UVA 손상 섬유아세포는 광스트레스 전사반응을 단순히 되돌리는 대신, 급성 stress/repair 반응을 유지하면서 안정적 senescence 전사궤적을 선택적으로 반대하는 상태로 이동한다.**

영문 working claim:

> **The response to a reprogramming-phase secretome is aligned with acute photostress and repair programs while consistently opposing stable senescence trajectories, indicating stress–senescence decoupling rather than simple transcriptional reversal.**

이 문장은 현재 결과를 가장 적게 과장하면서도 선행논문의 기능적 결과와 연결된다. 핵심은 `photoaging reversal`, `general rejuvenation`, `longevity`가 아니라 **stress–senescence decoupling**이다.

현재 상태에서 주장하면 안 되는 문장은 다음과 같다.

- “Repro-CM이 UVA 전사손상을 전반적으로 역전한다.” 실제로 acute/cumulative UVA signature와 양의 상관이다.
- “Repro-CM이 보편적 biological age를 낮춘다.” fibroblast age 자료 일부에서만 약한 반대 방향이고, human multi-tissue 및 rodent age signature에서는 거의 0이다.
- “Repro-CM이 보편적 longevity program을 활성화한다.” maximum-lifespan signature와 반대 방향이었다.
- “miR-302/367이 recipient cell의 변화를 매개한다.” extracellular miRNA는 n=1이고, 검증 표적 GSEA도 유의하지 않았다.
- “부분 리프로그래밍 세포의 intracellular program을 recipient가 그대로 재현한다.” GSE165177 MPTR과의 상관은 사실상 0이다.

## 2. iPSC-CM 원본 교체와 정량 결과

사용자가 지정한 `/home/shim/Downloads/project/project_shim/RNA_seq/00.data/iPSC_CM_24h.R1.fastq.gz`는 전체 gzip 검사를 통과했다. 이 파일에서 생성되어 있던 paired Trim Galore 결과를 HDF-CM 및 Repro-CM과 함께 STAR로 다시 정량했다.

기존 RSEM 산출물은 STAR 입력 read pair의 약 2–3%에 해당하는 count만 남아 있어 분석 입력에서 제외했다. 새 분석은 GENCODE v41 기반 STAR `ReadsPerGene.out.tab`의 reverse-stranded 열을 사용했다.

| 시료 | 입력 read pairs | unique mapping | reverse-strand gene assignment | mismatch |
|---|---:|---:|---:|---:|
| HDF-CM | 16,676,510 | 94.48% | 85.69% | 0.26% |
| Repro-CM | 15,342,652 | 93.67% | 84.62% | 0.26% |
| iPSC-CM | 16,843,629 | 94.66% | 86.21% | 0.26% |

세 library 모두 FastQC per-base quality와 adapter content는 PASS였다. 다만 생물학적 반복은 조건당 1개이므로, local RNA-seq에는 DEG p-value를 부여하지 않았다.

## 3. 보수적 Repro-CM anchor

TMM-CPM 후 `CPM >= 1`인 시료가 2개 이상인 유전자만 사용했다. 각 유전자에서 다음을 계산했다.

- `dH = log2CPM(Repro-CM) - log2CPM(HDF-CM)`
- `dI = log2CPM(Repro-CM) - log2CPM(iPSC-CM)`
- 두 방향이 같을 때만 `score = sign(dH) × min(|dH|, |dI|)`
- 방향이 다르면 score를 0으로 설정

그 결과 12,319개 발현 유전자 중 8,921개(72.4%)가 두 비교자에 대해 같은 방향이었다. 2,204개는 `|score| >= log2(1.25)`, 580개는 `|score| >= log2(1.5)`, 60개는 `|score| >= 1`이었다.

이 score는 **효과크기 기반의 기술적 anchor**이지, local sample 수준의 통계적 DEG가 아니다.

## 4. 핵심 수치

### 4.1 독립 데이터셋별 방향

| 축 | 데이터셋 | Repro-CM anchor와 Spearman rho | 발현량-매칭 FDR | 판정 |
|---|---|---:|---:|---|
| UVA | GSE240226 acute UVA | +0.255 | 0.00045 | 같은 방향 |
| UVA | GSE302943 cumulative UVA | +0.313 | 0.00045 | 같은 방향 |
| senescence | GSE179848 longitudinal RS | -0.100 | 0.00045 | 반대 방향 |
| senescence | GSE109700 deep RS | -0.134 | 0.00045 | 반대 방향 |
| senescence | GSE191055 passage RS | -0.044 | 0.00180 | 반대 방향 |
| senescence | GSE93535 SIPS | -0.163 | 0.00045 | 반대 방향 |

`p`가 작은 이유를 유전자 수가 많기 때문이라고만 볼 수 없도록, local 평균발현량 decile 안에서 score를 섞는 5,000회 순열을 사용했다. 그래도 유전자는 독립적인 생물학적 반복이 아니므로, 이 FDR은 **signature association**에 대한 유의성이지 치료 효과의 biological replicate p-value가 아니다.

### 4.2 중심 decoupling index

두 UVA 연구를 rank-normalize한 뒤 평균하여 UVA meta-axis를 만들고, 네 독립 senescence 연구를 같은 방식으로 합쳐 senescence meta-axis를 만들었다.

- 공통 comparator-consistent 유전자: 2,237개
- `rho(Repro, UVA meta) = +0.289`
- `rho(Repro, senescence meta) = -0.172`
- `D = rho(UVA) - rho(senescence) = 0.460`
- 발현량-층화 10,000회 순열: `empirical p < 1 × 10^-4`
- UVA 축을 회귀 제거한 senescence residual과도 `rho = -0.195`

즉 단순히 UVA와 senescence가 원래 반대이기 때문에 생긴 결과만은 아니다.

### 4.3 대안설명에 대한 민감도

| 검정 | 결과 |
|---|---|
| 6개 primary 연구를 하나씩 제외 | D = 0.401–0.478 |
| Repro vs HDF만 사용 | D = 0.419 |
| Repro vs iPSC만 사용 | D = 0.241 |
| translation/rRNA/IFN 유전자 제거 | D = 0.464 |
| cell-cycle/DNA-replication 유전자 제거 | D = 0.475 |
| 위 두 군 모두 제거 | D = 0.476 |
| local top 500 유전자 제거 | D = 0.442 |
| read count 50% thinning 100회 | median D = 0.438; 95% 범위 0.421–0.457 |

50% thinning은 sequencing/counting 안정성만 평가한다. 생물학적 반복을 대신하지 않는다.

### 4.4 “양의 UVA 상관은 잔존 손상일 뿐”이라는 반론

GSE240226에서 Maifuyin 및 succinate rescue effect는 UVA injury effect와 각각 음의 상관이었다. 각 rescue signature에서 UVA component를 회귀 제거한 뒤에도 Repro-CM anchor는 다음과 같이 양의 상관을 보였다.

- Maifuyin residual rescue: `rho = +0.127`, expression-matched `FDR < 0.0002`
- Succinate residual rescue: `rho = +0.094`, expression-matched `FDR < 0.0002`

따라서 Repro-CM이 단지 손상된 채 남아 있다는 설명보다는, **acute stress component와 repair-aligned component가 공존한다**는 설명이 더 잘 맞는다. 다만 두 rescue arm이 한 연구에 속하므로 독립 rescue 연구가 하나 더 필요하다.

### 4.5 pathway 수준

553개 공통 Reactome pathway 전체에서:

- local Repro-CM NES 대 senescence meta-NES: `rho = -0.528`, `p = 4.29 × 10^-41`
- local Repro-CM NES 대 UVA meta-NES: `rho = +0.079`, `p = 0.063`

따라서 **pathway 전체 수준의 anti-senescence 관계는 매우 강하지만, UVA와의 global pathway concordance는 유의하지 않다.** UVA와 같은 방향인 것은 rRNA modification, mature-mRNA transport, splicing, PERK, heat-stress, NGF transcription 등 선택된 stress-adaptation module이다. 선택 pathway heatmap은 설명용이고, 553-pathway 전체 상관이 검정용이다.

## 5. miRNA 및 longevity 판정

### miRNA

- mature miRNA mapping: HDF 0.728%, iPSC 3.200%, Repro 10.739%
- biological replicate: 조건당 1개
- Repro-CM에서 miR-302/367 및 miR-371/372/373 family composition이 높아 transient-reprogramming identity에는 부합
- OmniPath 내 miRTarBase 문헌근거 표적을 이용한 miR-302/367 target GSEA: `NES = -1.05`, `FDR = 0.399`

결론: 본문 기전 figure가 아니라 Supplementary의 identity/candidate 자료로 제한한다.

### longevity와 broad rejuvenation

- Tyshkovskiy 2019 lifespan-positive gene set: 유의하지 않음
- Tyshkovskiy 2026 rodent maximum-lifespan signature: `rho ≈ -0.055`, 가설과 반대
- human multi-tissue age 및 rodent chronological-age: 거의 0
- GSE165177 MPTR: 전체 `rho = 0.0004`, day 13 `rho = 0.011`, 모두 비유의

결론: 선행논문의 C. elegans lifespan 결과는 중요한 생물학적 맥락이지만, 새 transcriptomic 논문의 중심을 `longevity signature`로 잡으면 현재 데이터와 충돌한다.

## 6. Figure 1–6의 새 흐름

현재 생성된 PNG는 최종 조판 전의 diagnostic panel이다. 최종 논문은 아래 순서가 가장 논리적이다.

### Figure 1. 공개자료가 규정하는 UVA-to-senescence 문제와 재현성 경계

- A: 선행연구의 질문을 새 분석 질문으로 변환하는 schematic
- B: acute/cumulative/chronic UVA 연구 간 gene/pathway concordance
- C: 네 독립 senescence 모델 간 concordance
- D: UVA와 senescence가 동일 축이 아님을 보여주는 비교

증명하는 것: public data가 재현하는 것은 `Repro-CM 치료 효과`가 아니라 UVA 및 senescence reference state이다. Repro-CM 효과를 직접 재현하려면 독립 Repro-CM 처리자료가 필요하다.

예상 증거등급: B. UVA는 context-dependent이고 senescence 쪽 재현성이 더 강하다.

### Figure 2. 유효한 iPSC-CM을 포함한 local Repro-CM anchor

- A: FASTQ–trim–STAR QC와 strandedness
- B: Repro-vs-HDF와 Repro-vs-iPSC 효과 산점도 및 signed-min anchor
- C: 상위 유전자와 전체 Reactome enrichment
- D: parameter/strand/read-thinning stability

증명하는 것: local signature가 입력오류나 한 비교자에 의해 만들어진 것은 아니다.

예상 증거등급: 기술적 A, 생물학적 B-/F. n=1 때문에 local differential inference는 불가능하다.

### Figure 3. 중심 결과—stress–senescence decoupling

- A: 2 UVA + 4 senescence 연구 forest plot
- B: UVA meta와 senescence meta에 대한 local rho
- C: D=0.460과 10,000회 matched-null 분포
- D: UVA를 회귀 제거한 residual senescence 관계

증명하는 것: 동일 local response가 acute/cumulative UVA와 같은 방향이면서 안정적 senescence와 반복적으로 반대 방향이다.

예상 증거등급: A-. 현재 논문의 가장 강한 결과다.

### Figure 4. RNA handling/proteostasis 중심의 상태 전환

- A: 전체 553 Reactome pathway의 local-vs-senescence scatter
- B: 선택된 stress/RNA-processing pathway heatmap
- C: 중복 pathway를 축약한 leading-edge network
- D: translation/rRNA/IFN 및 cell-cycle 제거 민감도

증명하는 것: 결과가 일부 marker gene이 아니라 pathway architecture에서도 유지된다.

예상 증거등급: anti-senescence A-, 구체적 기전 B. pathway 상관만으로 분자 인과를 주장할 수 없다.

### Figure 5. repair alignment와 robustness

- A: UVA injury와 두 rescue signature의 관계
- B: UVA component 제거 후 Repro-rescue residual association
- C: leave-one-study-out
- D: comparator, gene-removal, read-thinning sensitivity

증명하는 것: positive UVA association이 단순 잔존 손상만으로 설명되지 않으며, 중심 대조가 분석 선택에 견고하다.

예상 증거등급: B+–A-. 독립 rescue 연구가 추가되면 A-에 가까워진다.

### Figure 6. held-out state-map과 반증 가능한 모델

다운로드 승인 후 새 자료를 사용한다.

- A: GSE297233/GSE149694에서 intracellular partial-reprogramming 축과 비교
- B: GSE307377에서 최근 human dermal fibroblast age 축 검증
- C: GSE306957에서 IFN/mtRNA–SASP dependency의 위치 확인
- D: GSE116968에서 독립 UVB photoprotection signature 검증
- E: `damage sensing retained → repair/RNA-quality-control engaged → stable senescence commitment opposed` 모델

증명하는 것: 모델이 이미 사용한 자료에만 맞는 사후 설명인지, 새로운 state와 intervention에도 일반화되는지 확인한다.

예상 증거등급: 아직 미평가. 이 figure가 성공해야 상위권 저널 도전이 현실적이다.

## 7. 논리 오류 점검

| 잠재 오류 | 현재 조치 | 잔여 위험 |
|---|---|---|
| public UVA 자료로 Repro-CM 효과를 “재현”했다고 표현 | reference-state 재현으로 명칭 수정 | 직접 Repro-CM 외부자료 부재 |
| 양의 UVA 상관을 protective response로 임의 해석 | 독립 rescue component를 회귀분리 | rescue 두 개가 같은 연구 |
| n=1 local RNA-seq에 DEG p-value 사용 | 효과크기 anchor만 사용 | biological variance 추정 불가 |
| ribosome/IFN/cell cycle이 상관을 독점 | 각 program 및 조합 제거 | 미측정 세포상태 confounding 가능 |
| gene 수가 많아 p만 작아짐 | expression-matched permutation 및 study-level LOO | 유전자는 완전 독립이 아님 |
| 선택 pathway만 보여 cherry-picking | 553개 전체 pathway 상관을 primary test로 설정 | Reactome pathway끼리 중복 |
| senescence를 age/longevity와 동일시 | 세 축을 분리하고 음성결과 보고 | 제목에서 rejuvenation 사용 시 과장 위험 |
| miRNA correlation을 전달·인과로 해석 | validated-target GSEA 실패를 공개하고 보조자료로 이동 | 직접 uptake/perturbation 없음 |
| 선행논문과 중복 출판 | 새 질문·새 metric·독립자료 일반화를 중심으로 구성 | local source data와 생물학적 맥락 중복은 명시 필요 |
| 사후가설을 확증적으로 표현 | 현 단계는 exploratory로 명시 | 새 held-out 자료를 분석 전에 contrast/threshold 고정해야 함 |

## 8. in-silico 논문의 선례가 주는 기준

- Gladyshev 연구진의 2026 Nature 연구는 4종, 25개 이상 조직, 11,000개 이상 transcriptome을 통합해 age/mortality/lifespan biomarker를 만들고 여러 외부 endpoint로 검증했다. 핵심은 단순 enrichment가 아니라 **새로운 정량 framework와 대규모 독립검증**이다: <https://doi.org/10.1038/s41586-026-10542-3>
- Aging Cell의 RNAge 연구는 6개 age dataset으로 score를 만들고 별도 fibroblast/brain/reprogramming 자료로 검증했다: <https://doi.org/10.1111/acel.70075>
- npj Aging의 knowledge-primed transcriptomic clock은 887개 skin transcriptome으로 예측성과 pathway 해석력을 함께 검증했다: <https://www.nature.com/articles/s41514-021-00068-5>
- PASTA는 다조직·bulk·single-cell 및 수백만 perturbation transcriptome에 적용되는 도구와 공개 package를 제공했다: <https://doi.org/10.1002/advs.76740>

따라서 “in silico만으로 가능”은 맞지만, 본 연구도 최소한 **새 metric, 사전에 고정된 held-out test, 독립 연구 간 재현성, 코드 공개**가 있어야 한다. 현재 decoupling index가 새 metric의 출발점이다.

## 9. 저널 현실성

### 현재 결과만으로

- **Aging Cell:** scope에는 정확히 맞지만 아직 미달. local n=1, 선행논문과의 자료 중복, 사후가설, 직접 causal validation 부재가 큰 약점이다. “highly significant fundamental aging biology” 수준으로 보이기 어렵다. 공식 scope: <https://www.anatsoc.org.uk/journals/aging-cell>
- **npj Aging:** 가장 적절한 상향 목표다. senescence, aging signaling, intervention을 명시적으로 다루며 기술적 건전성·결론의 강한 근거·분야 중요성을 심사기준으로 둔다. 2026-12-03 마감의 Cellular Senescence collection도 현재 열려 있다. <https://www.nature.com/npjamd/aims>, <https://www.nature.com/npjamd/for-authors-and-referees/guideforreviewers>, <https://www.nature.com/npjamd/calls-for-papers>
- **GeroScience:** 주제는 맞지만 새 기전 또는 강한 translational link가 필요해 현재는 도전적이다.
- **Aging (Albany NY) / Frontiers in Aging:** held-out validation까지 완료하면 현실적인 후보군이다. Frontiers는 public-data computational study에 적절한 validation을 명시적으로 요구한다: <https://www.frontiersin.org/journals/aging/about>
- **Mechanisms of Ageing and Development:** 현재 전략과는 맞지 않는다. 공식 설명상 기능분석과 연결되지 않은 descriptive/correlative omics는 고려하지 않는다고 명시한다: <https://shop.elsevier.com/journals/mechanisms-of-ageing-and-development/0047-6374>
- **Stem Cell Reports:** secretome source가 reprogramming cell이라는 점만으로는 부족하다. 새 stem-cell biology 또는 broad conceptual advance가 요구되므로 현 상태의 우선순위는 낮다: <https://www.isscr.org/stem-cell-reports/>

### Figure 6 held-out 검증이 성공할 경우

1차 목표는 `npj Aging`, 2차는 `GeroScience` 또는 `Aging`, 보수적 목표는 `Frontiers in Aging` 또는 genomics/systems-biology 계열 저널이 합리적이다. Aging Cell 투고는 가능하지만, 편집단계 통과 가능성을 높이려면 독립 Repro-CM 자료 또는 최소 하나의 새 기능/인과 실험이 사실상 필요하다.

## 10. 다음 다운로드 전 데이터 점검

아래 파일은 아직 다운로드하지 않았다.

| 우선순위 | 자료 | 용도 | 표본/파일 | 판단 |
|---:|---|---|---|---|
| 1 | GSE297233 | 96세 donor fibroblast의 4일 OSK/변형 OSK partial reprogramming; mesenchymal drift | 8 samples, raw-count matrix 471.5 KB | 가장 먼저 |
| 2 | GSE307377 | 2026년 공개된 독립 young/old primary dermal fibroblast age test | 4 young + 5 old; counts 1.3 MB + TPM 1.4 MB | 매우 높음 |
| 3 | GSE116968 | 독립 photoprotection: UVB와 red-light pretreatment, 1 h/4 h | 24 samples; processed XLSX 11.7 MB ×2 | rescue 반론 검증에 중요 |
| 4 | GSE149694 | human fibroblast reprogramming time course | 32 samples; processed TAR 6.5 MB | state trajectory 검증 |
| 5 | GSE306957 | MRC5 senescence에서 RIG-I/MDA5/MAVS/mtRNA–SASP 기전 | 34 samples; processed TAR 140.5 MB | IFN 해석 경계 검증 |

공식 GEO 근거:

- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297233>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307377>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116968>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149694>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306957>

권장 다운로드는 우선 1–4의 processed data 총 약 33 MB이다. FASTQ는 이 단계에서 필요 없다. 결과가 중심모델을 지지하고 원시처리 일관성이 최종 심사에 필요하다고 판단될 때만 raw FASTQ로 확장한다.

## 11. 생성 결과 위치

- 전체 증거 판정표: `public_data_tierA/derived/evidence_scorecard.tsv`
- local mRNA anchor: `public_data_tierA/derived/local_repro_anchor/`
- 중심 decoupling 결과: `public_data_tierA/derived/decoupling_validation/`
- local miRNA 품질 및 target 검정: `public_data_tierA/derived/local_mirna_audit/`
- 진단용 figure: `public_data_tierA/derived/manuscript_diagnostic_figures/`

현재 결론은 **논문이 될 만한 통계적 패턴은 확보했지만, 그 패턴의 정확한 이름은 longevity/rejuvenation이 아니라 stress–senescence decoupling이며, 상위 저널 가능성은 새 held-out Figure 6의 성공 여부에 달려 있다**는 것이다.
