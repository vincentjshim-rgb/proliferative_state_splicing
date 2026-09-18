# Repro-CM–노화성 HDF 연구의 논문 전략 및 AI 적용 타당성 검토

> 검토 대상: `aging_virtual_cell_cis_regulation_research_note.md`, 현재 FASTQ 구성, Repro-CM splicing workflow 문서  
> 기준일: 2026-09-09  
> 목적: AlphaGenome/virtual cell의 타당한 위치를 정하고, Aging Cell 수준의 논문이 되기 위한 주장–근거–실험 연결을 설계한다.

## 1. 결론부터

### 핵심 판단

1. **AlphaGenome은 현재 논문의 주 분석축으로 넣지 않는 것이 타당하다.** 현재 개입은 DNA 염기서열 변이가 아니라 서로 다른 conditioned medium(CM)이다. 따라서 관찰하려는 변화는 주로 수용 세포에서 일어나는 `trans-acting extracellular signal -> signaling/RBP state -> RNA/splicing -> phenotype`이며, AlphaGenome의 `DNA sequence/variant -> cis-regulatory output` 질문과 다르다.
2. **virtual cell은 제한적으로 타당하다.** 다만 현재 4개의 bulk RNA-seq로 모델을 학습하거나 CM 효과를 예측했다고 주장할 수는 없다. 복합 CM을 단백질·EV·후보 ligand 또는 유전자 perturbation으로 분해하고, matched single-cell perturbation data를 만든 뒤, 후보 개입을 순위화하고 **미관측 개입을 전향적으로 검증**하는 용도로 써야 한다.
3. **Stereo-seq은 현재 단층 HDF 배양 논문에는 타당성이 낮다.** 조직 niche가 없는 단일세포 배양에서 spatial omics를 붙이면 기술 전시가 된다. 3D skin equivalent, explant 또는 in vivo 피부 모델로 확장할 때만 고려한다.
4. 현재 자료는 중요한 pilot이지만 **조건당 biological sample 1개**, time-matched vehicle control 부재, healthy control 부재, CM donor/batch 효과의 미분리 때문에 Aging Cell 원저의 중심 근거로는 부족하다. 현재 자료로 허용되는 결론은 “후속 검증할 splice-junction 후보를 발굴했다”까지이다.
5. 가장 강한 논문 축은 AI 자체가 아니라 다음의 기계론이다.

> **Reprogramming-transition cell의 secretome이 UVA로 손상된 인간 피부 fibroblast의 특정 receptor–signaling–splicing regulator–isoform 축을 조절하여, fibroblast identity를 잃지 않으면서 지속 가능한 기능 회복을 유도한다.**

AI는 이 축에서 `후보 유전자/조합 우선순위화 -> prospective test`를 효율화할 때만 기여한다.

### 투고 가능성의 현실적 등급

| 연구 상태 | Aging Cell 판단 | 이유 |
|---|---|---|
| 현재 4개 bulk RNA-seq만 | 투고 비권장 | N=1, 필수 대조군 부재, 통계적 추론 불가, 기전과 기능 검증 없음 |
| 독립 반복 + senescence/ECM 기능 검증 | 여전히 약함 | 현상 재현은 되지만 active CM component와 인과기전이 없음 |
| active factor/fraction + receptor/RBP + causal isoform rescue + 3D/ex vivo 검증 | **현실적인 Aging Cell 후보** | descriptive study를 넘어 aging mechanism과 intervention causality를 제공 |
| 위 내용 + 여러 donor/노화 유도법/조직에서 일반화 + in vivo 효능·안전성 | Nature Aging 등 더 높은 범위 검토 가능 | 특정 배양모델을 넘어 일반 geroscience 원리와 조직 기능을 제시 |

정확한 게재 가능성은 결과의 효과크기·재현성·새로움에 달려 있어 지금 수치로 말할 수 없다. 다만 현재 상태는 Aging Cell 공식 scope가 배제하는 “purely descriptive” 연구에 가깝고, 완성 설계는 해당 저널이 요구하는 mechanistic aging biology에 맞출 수 있다.

---

## 2. 현재 자료가 실제로 답할 수 있는 질문

현재 등록된 네 조건은 다음과 같다.

- `15J_0hr`
- `HDF_CM_24h`
- `Rep_CM_24h`
- `iPSC_CM_24h`
- 각 조건당 biological sample 1개

현재 가장 해석 가능한 비교는 다음 두 가지의 **탐색적 후보 비교**이다.

1. `Rep_CM_24h vs HDF_CM_24h`
2. `Rep_CM_24h vs iPSC_CM_24h`

이 비교도 CM class뿐 아니라 source cell line, donor, CM production batch가 함께 바뀌었다면 “Repro-CM 고유 효과”로 확정할 수 없다. 반면 모든 `24h` 조건과 `15J_0hr`의 비교에는 treatment와 elapsed time이 동시에 달라지는 결정적 confounding이 있다.

### 현재 자료에서 가능한 표현

> 단일 pilot experiment에서 Rep-CM 처리 HDF는 HDF-CM 또는 iPSC-CM 처리 HDF와 비교하여 splice-junction usage가 다른 후보 event를 보였다.

### 현재 자료에서 불가능한 표현

- Rep-CM이 UVA-induced aging을 유의하게 되돌렸다.
- Rep-CM이 biological age를 낮췄다.
- 확인된 alternative splicing event가 기능 회복의 원인이다.
- iPSC-CM보다 Rep-CM이 우월하다.
- 특정 secreted factor 또는 EV가 효과의 원인이다.
- virtual cell이 새로운 rejuvenation intervention을 정확히 예측했다.

---

## 3. 현재 연구 프레이밍의 논리 점검

### 3.1 너무 넓은 원래 수식

원 노트의 식은 장기 연구 프로그램에는 맞지만, 한 논문의 identifiable question으로는 너무 넓다.

\[
P(\text{future state}\mid \text{genotype},\text{initial state},\text{niche},\text{perturbation},\text{time})
\]

문제는 genotype, cell state, niche, perturbation, time을 한 연구에서 모두 독립적으로 변화시키고 충분히 반복 측정하지 않으면 각 효과를 분리할 수 없다는 점이다. 또한 단순 조건부 확률은 intervention causality를 뜻하지 않는다.

현재 논문에 맞는 질문은 다음처럼 좁히는 것이 좋다.

\[
\Delta Y_{t+\tau}=f(X_t,\operatorname{do}(\text{defined CM component or gene perturbation}),\tau)
\]

- \(X_t\): UVA 손상 후의 초기 HDF 상태
- 개입: 표준화된 CM, 분획, 단일 ligand, neutralization, receptor/RBP perturbation
- \(Y\): senescence + fibroblast function + isoform outcome의 사전 정의된 묶음
- genotype과 spatial niche는 현재 core claim에서 제외

### 3.2 `variant -> cis effect -> cell-state transition -> phenotype`

이 연결은 개념적으로 옳지만 현재 실험에는 variant가 없다. 따라서 이 사슬을 현재 CM 논문의 main story로 사용하면 첫 화살표가 데이터 없이 시작된다. 장기 precision-aging 프로그램 또는 별도의 pharmacogenomic response study에는 적합하다.

### 3.3 `young-like transcriptome -> rejuvenation`

이 연결은 성립하지 않는다. young-like expression은 다음 현상으로도 생길 수 있다.

- cell-cycle 재진입 또는 증식성 세포의 선택
- 일시적 stress suppression
- 세포 사멸로 senescent subpopulation이 사라짐
- fibroblast identity의 상실 또는 부분 탈분화
- global RNA content 변화와 normalization artifact

따라서 rejuvenation은 transcriptomic direction 하나가 아니라 다음의 conjunction으로 정의해야 한다.

\[
\text{rejuvenation evidence}=
\text{aging-feature attenuation}
+\text{function recovery}
+\text{identity preservation}
+\text{durability}
+\text{safety}
\]

초기 원고에서는 “rejuvenation”보다 **attenuation of photoaging-associated phenotypes** 또는 **functional restoration**이 안전하다.

### 3.4 `alternative splicing change -> phenotype recovery`

ΔPSI는 인과성을 뜻하지 않는다. splicing 변화가 gene-expression 변화, cell cycle, RNA quality 또는 splicing-factor abundance의 결과일 수 있다. 다음 두 방향의 실험이 모두 필요하다.

- **Necessity:** CM이 효과가 있어도 isoform-switching ASO/CRISPR 또는 candidate RBP 억제로 기능 회복이 사라지는가?
- **Sufficiency/rescue:** CM 없이도 해당 isoform을 유도하면 phenotype이 개선되는가? 반대 isoform을 되돌리면 CM 효과가 소실되는가?

### 3.5 `CM treatment -> secreted factor mechanism`

recipient RNA-seq만으로 secreted factor를 지정할 수 없다. CM에는 soluble proteins, metabolites, lipids, extracellular particles/vesicles, 배지 성분, source-cell death products가 함께 있다. active fraction과 factor를 직접 분리해야 한다.

### 3.6 `24 h response -> aging reversal`

24 h는 early response를 잡는 데 좋지만 stable aging reversal을 주장하기에는 너무 짧다. 최소한 treatment washout 뒤 수일에서 1–2주간 phenotype이 유지되는지 확인해야 한다. 반대로 early time point는 receptor signaling과 primary transcription/splicing response를 찾기 위해 필요하다. 즉 `0.5–2 h`, `6 h`, `24 h`, `72 h`, washout follow-up을 역할별로 나눠야 한다.

### 3.7 `UVA photoaging -> general aging`

UVA-induced damage는 extrinsic photoaging 모델이지 intrinsic chronological aging 전체와 동일하지 않다. 결론은 처음에는 skin photoaging/HDF senescence에 한정한다. replicative senescence, irradiation-induced senescence 또는 naturally aged donor HDF에서 재현될 때 일반화 범위를 넓힌다.

### 3.8 `scRNA-seq -> splicing`

표준 3'-tag scRNA-seq는 상세 isoform/alternative-splicing 결론에 부적합하다. scRNA-seq는 responder heterogeneity와 state-transition에 사용하고, splicing은 deep paired-end bulk RNA-seq, junction RT-PCR/ddPCR, long-read RNA-seq 또는 targeted isoform sequencing으로 검증한다.

---

## 4. AlphaGenome 적용 여부

### 현재 core paper: 사용하지 않음

AlphaGenome은 1 Mb DNA sequence로 RNA, accessibility, TF binding, histone modification, contact map, splicing 관련 track과 **variant effect**를 예측한다. 그러나 다음을 입력으로 받는 모델은 아니다.

- Rep-CM/HDF-CM/iPSC-CM 처리
- UVA dose와 treatment time
- 수용 HDF의 phospho-signaling 또는 RBP abundance
- CM의 단백질·metabolite·EV 조성
- 동일 DNA 서열에서 treatment에 의해 발생한 trans-splicing 변화

AlphaGenome 논문 자체도 condition-specific variant effect, cell-context specificity, personal-genome prediction, complex phenotype 연결이 여전히 제한적이라고 명시한다. 특히 AlphaGenome이 splice junction을 예측할 수 있다는 사실과, **CM이 유도한 splicing 변화**를 예측할 수 있다는 주장은 전혀 다르다.

### AlphaGenome이 타당해지는 별도 질문

다음 조건이 생기면 companion aim으로 사용할 수 있다.

> “왜 일부 HDF donor만 Rep-CM에 반응하는가? 반응 차이가 cis-regulatory 또는 splice-regulatory variant로 설명되는가?”

필요한 설계:

1. 충분히 많은 genotyped/WGS recipient donors에서 Rep-CM response를 반복 측정한다.
2. `genotype x treatment` response-QTL 또는 heterozygote 내 allele-specific expression/splicing을 찾는다.
3. AlphaGenome으로 candidate variant의 expression, accessibility, TF-binding, splicing delta를 예측한다.
4. response cell state에서 열린 cCRE/scATAC evidence와 겹치는 후보만 남긴다.
5. base editing 또는 MPRA/minigene으로 해당 염기를 직접 검증한다.
6. 편집된 isogenic cell에서 CM response 차이를 확인한다.

이 경우에도 AlphaGenome은 causal evidence가 아니라 variant prioritization 도구다. 표본수는 effect size와 allele frequency에 대한 사전 power calculation로 정해야 하며, 현재 4개 sample과는 규모가 전혀 다르다.

### 의사결정

- **현재 논문:** AlphaGenome figure를 넣지 않는다.
- **향후 precision-response 연구:** donor heterogeneity와 genotype 자료가 확보되면 별도 논문/aim으로 전개한다.
- 후보 SNP 없이 임의 aging GWAS variant를 AlphaGenome에 넣어 CM 기전과 연결하는 것은 피한다. 이는 post hoc story가 된다.

---

## 5. Virtual cell 적용안: “virtual fibroblast”로 좁혀 사용

### 타당한 이유

이 연구에는 많은 후보 ligand, receptor, signaling node, RBP와 조합이 생길 수 있어 모두 wet-lab으로 검사하기 어렵다. 따라서 virtual cell은 다음 질문에 실제 가치를 줄 수 있다.

> “UVA-damaged HDF의 state를 기능적으로 더 건강한 방향으로 이동시키면서 fibroblast identity와 종양 억제 안전성을 유지할 가능성이 높은 정의된 perturbation은 무엇인가?”

중요한 것은 whole CM을 하나의 자연어 label로 넣는 것이 아니라 개입을 다음처럼 정의하는 것이다.

- recombinant ligand의 identity와 dose
- neutralizing antibody 또는 receptor inhibitor
- CRISPRi/a target
- splice-switching oligonucleotide
- EV-rich/EV-depleted/soluble fraction
- 처리 시간

### 현재 자료로 당장 할 수 없는 것

- 4개 bulk sample로 foundation/transition model 학습
- CM label만으로 unseen CM mixture 예측
- gene-level virtual cell output으로 isoform mechanism 증명
- public cancer cell perturbation data에서 HDF photoaging으로 zero-shot 이전한 값을 효능 근거로 사용

### 권장 데이터 생성

1. **Discovery panel:** secretome proteomics와 현재 RNA/splicing pilot로 30–100개의 ligand–receptor–TF/RBP–isoform 후보를 만든다.
2. **Matched perturbation screen:** UVA-damaged primary HDF에서 focused CRISPRi/a 또는 arrayed perturbation을 실시한다. 가능하면 non-targeting control, positive control, vehicle, Rep-CM을 같은 batch에서 처리한다.
3. **Single-cell readout:** scRNA-seq는 cell-state heterogeneity와 responder fraction을 읽는다. 핵심 splice event는 별도 targeted junction assay로 같은 조건에서 읽는다.
4. **Donor structure:** donor와 CM batch가 train/test 양쪽에 섞이지 않게 group-wise split한다.
5. **Prospective test:** 모델 학습 때 보지 않은 단일/조합 perturbation을 사전 등록된 기준으로 골라 새 donor와 새 CM batch에서 검증한다.

### 예측 목표: 단일 “young score” 금지

다목적 score가 적절하다.

\[
S(a)=w_1F_{ECM}+w_2F_{repair}-w_3SASP-w_4DNA\ damage
-w_5identity\ loss-w_6unsafe\ proliferation
\]

- $F_{ECM}$: collagen/ECM synthesis와 deposition
- $F_{repair}$: migration 또는 3D tissue repair
- SASP: 단백질 수준의 inflammatory output
- DNA damage: γH2AX/53BP1 foci, persistent DNA-SCARS 등
- identity loss: fibroblast identity/ECM program 붕괴
- unsafe proliferation: 비정상적 cell-cycle activation, anchorage-independent growth 등

가중치는 결과를 본 뒤 정하지 말고 사전 정의하거나, Pareto ranking으로 제시한다.

### 모델 평가의 필수 조건

- mean/no-change, linear, additive baseline과 비교
- held-out **donor**, **perturbation**, **CM batch**, 가능하면 **senescence trigger** 평가
- all-gene correlation만 쓰지 않고 perturbation-specific DE, pathway, responder fraction, functional endpoint 예측 평가
- batch 또는 평균 treatment shift만 학습한 모델이 아닌지 Systema/Cell-Eval 유형의 평가 사용
- top-k 모델 추천군과 non-model/random 또는 baseline 추천군의 prospective hit rate 비교
- 모델이 baseline을 이기지 못하면 결과를 숨기지 않고, 모델을 제거하거나 negative benchmarking result로 제한

### 어떤 모델을 쓸 것인가

처음부터 거대한 새 foundation model을 만들 이유는 없다.

1. pretrained embedding/State 계열을 사용할 수 있는지 기술적으로 점검한다.
2. 동일 split에서 regularized linear model, nearest-neighbor/perturbation mean, additive model을 먼저 구축한다.
3. matched HDF perturbation data로 fine-tuning 또는 작은 conditional transition model을 비교한다.
4. 성능과 prospective utility가 확인될 때만 “virtual fibroblast”라는 표현을 쓴다.

현재 field에서도 unseen perturbation/context 일반화는 어려운 문제이며, 최신 모델의 성공 주장과 독립 benchmark의 경고가 공존한다. 따라서 virtual cell은 wet-lab의 대체물이 아니라 **검증할 실험을 더 잘 고르는 시스템**으로 정의하는 것이 가장 방어 가능하다.

### Go/No-go gate

| Gate | 통과 조건 | 실패 시 조치 |
|---|---|---|
| G1: biological effect | 독립 donor/batch에서 Rep-CM phenotype 재현 | AI 단계 중단, CM protocol부터 교정 |
| G2: defined perturbations | active fraction과 후보 ligand/receptor/RBP panel 확보 | whole-CM 비교 논문으로 축소 |
| G3: predictive value | held-out donor/perturbation에서 simple baseline 초과 | virtual cell claim 제거 |
| G4: prospective utility | 새 실험에서 모델 추천 hit rate가 대조 추천보다 높음 | hypothesis generator로만 보고 |
| G5: functional causality | 추천 perturbation이 RNA뿐 아니라 function 개선 | transcriptomic state prediction으로만 제한 |

---

## 6. 권장 실험 설계

### 6.1 biological unit를 명확히 한다

이 연구에는 최소 두 종류의 독립 biological unit가 있다.

1. **recipient HDF donor**
2. **CM source line/production batch**

같은 donor의 세포를 여러 well에 나눈 것은 technical replicate이지 독립 biological replicate가 아니다. 하나의 CM pool을 여러 well에 처리해도 CM biological replicate가 늘어나지 않는다.

권장 시작점은 여러 recipient donor와 최소 3개의 독립 CM production batch이며, 최종 수는 pilot variance와 primary endpoint에 대한 power analysis로 정한다. 성별·donor age·passage는 기록하고 무작위화/블라인딩 가능한 단계는 적용한다.

### 6.2 최소 조건 행렬

| Recipient state | Treatment | 목적 |
|---|---|---|
| Healthy/young-like HDF | vehicle/base medium | 정상 기준 |
| UVA-damaged HDF | vehicle/base medium, time-matched | 자연 회복과 시간효과 |
| UVA-damaged HDF | HDF-CM | 일반 HDF paracrine control |
| UVA-damaged HDF | Rep-CM | 주 개입 |
| UVA-damaged HDF | iPSC-CM | pluripotent-state secretome comparator |
| Healthy HDF | Rep-CM | 비손상 세포의 비정상 증식/identity effect |

모든 CM은 동일한 basal medium, conditioning duration, source-cell viable count/confluence, 수집·원심분리·filtration·freeze–thaw, 총 단백질 또는 정의된 normalization 조건을 맞춘다. source-cell death와 reprogramming reagent/vector carry-over도 측정한다.

특히 Rep-CM은 “reprogramming-transition”이라는 source state가 batch마다 같은지 marker panel과 transcript/state score로 release criterion을 정해야 한다. source state가 확인되지 않으면 CM의 차이를 reprogramming stage의 차이라고 해석할 수 없다. 논문 첫 Methods와 Figure 1에서 HDF-CM, Rep-CM, iPSC-CM의 source cell identity와 collection protocol을 모호하지 않게 정의한다.

### 6.3 시간 설계

- `0 h`: 손상 후 treatment 전 baseline
- `0.5–2 h`: receptor/phospho-signaling
- `6 h`: primary transcription/RBP response
- `24 h`: 현재 pilot과 연결되는 RNA/splicing response
- `72 h`: phenotype conversion 여부
- `washout + 7–14 d`: 지속성, 재발, 이상 증식

모든 시점을 RNA-seq할 필요는 없다. early signaling은 phosphoprotein assay, 핵심 time point는 omics, 장기 시점은 phenotype으로 역할을 나누면 된다.

### 6.4 senescence/photoaging 판정

단일 marker로 정의하지 않는다. 최소 세 축을 사용한다.

1. **stable cell-cycle arrest:** EdU/Ki-67, 장기 증식 회복 여부
2. **structural/damage markers:** SA-β-gal, LMNB1, p16/p21, γH2AX/53BP1 또는 DNA-SCARS
3. **secretory/functional phenotype:** SASP protein panel, ROS/mitochondrial function, collagen I/III, MMP1, migration/ECM deposition

세포 수 감소로 marker-positive fraction이 낮아진 것을 rejuvenation으로 오해하지 않도록 viability, apoptosis와 absolute cell counts를 함께 측정한다.

### 6.5 splicing discovery와 검증

현재 rMATS workflow는 pilot 후보 생성에 적절하다. replicate 실험에서는 다음을 추가한다.

- donor-paired design과 CM-batch structure를 보존
- FDR를 포함한 replicate-aware differential splicing
- gene-expression 변화와 ΔPSI를 분리해 보고
- 가능하면 두 알고리즘 또는 junction-level orthogonal check로 robustness 확인
- IGV/Sashimi plot의 수동 확인
- independent samples에서 junction RT-PCR/ddPCR
- coding consequence가 중요하면 long-read 또는 targeted full-length sequencing
- 단백질 isoform/기능 결과 확인

상위 후보는 ΔPSI 크기만으로 고르지 말고 `coverage + donor consistency + aging concordance + RBP motif/evidence + predicted protein consequence + manipulability`로 순위화한다.

### 6.6 CM의 active component 규명

권장 순서:

1. Rep-CM, HDF-CM, iPSC-CM의 unbiased proteomics/cytokine/metabolite profiling
2. molecular-weight fractionation, heat/protease/RNase sensitivity, EV-rich와 EV-depleted fraction 비교
3. unconditioned medium을 동일 처리한 process control
4. candidate depletion/neutralization으로 necessity 확인
5. recombinant factor 또는 defined combination add-back으로 sufficiency 확인
6. recipient receptor KO/knockdown 또는 inhibitor로 receptor dependency 확인

EV를 주장한다면 단순 precipitation product를 “exosome”이라고 부르지 말고 MISEV2023의 separation, characterization, process-control 원칙을 따른다.

### 6.7 통계 모형

가능하면 CM source와 recipient donor를 교차시킨 설계를 사용한다. 일반적인 형태는 다음과 같다.

\[
y \sim \text{recipient state}*\text{treatment}*\text{time}
+(1|\text{recipient donor})+(1|\text{CM source line/batch})
\]

- 주 효과보다 `damaged state x treatment` interaction이 핵심이다.
- technical well을 biological N으로 세지 않는다.
- omics의 discovery와 confirmatory endpoint를 구분한다.
- primary endpoint와 contrast를 사전 지정한다.
- 제외 기준, batch correction, multiple testing, effect size와 confidence interval을 보고한다.

---

## 7. 권장 논문 스토리와 Figure 흐름

### 중심 가설

> Rep-CM에는 reprogramming transition에서만 풍부한 active component 조합이 있고, 이 조합은 UVA-damaged HDF의 정의된 receptor–signaling–RBP 축을 바꾸어 기능적으로 중요한 isoform program을 회복한다. 이 효과는 단순한 증식, selection 또는 탈분화가 아니며 3D 피부 환경에서도 ECM 기능을 개선한다.

후보가 확인되기 전에는 `[factor]`, `[receptor]`, `[RBP]`, `[isoform]` 자리를 비워 두고 결과에 따라 채운다. 알려진 pathway에 맞추어 후보를 억지로 선택하지 않는다.

### Figure 1. Rep-CM이 time-matched UVA-HDF phenotype을 선택적으로 개선한다

**질문:** Rep-CM 효과가 존재하고 재현되는가?

필요 자료:

- CM 생성과 처리 schematic
- healthy + UVA/vehicle + 세 CM의 완전한 대조군
- 여러 recipient donor와 여러 CM batch
- viability/apoptosis, EdU/Ki-67, SA-β-gal, LMNB1/p16/p21, SASP
- collagen/MMP/ECM 또는 migration 같은 fibroblast function

**논리 연결:** Rep-CM 처리군 개선 -> 후속 기전 탐색의 근거.

**오류 방지:** vehicle 24 h가 없으면 시간효과를 Rep-CM 효과로 오인한다. viability가 없으면 senescent-cell killing을 rejuvenation으로 오인한다. HDF-CM과 iPSC-CM은 “음성 대조”가 아니라 각각 생물학적으로 active comparator다.

**통과 기준:** 독립 donor와 CM batch에서 effect direction이 재현되고, 적어도 하나의 cell-state 축과 하나의 fibroblast-function 축이 함께 개선된다.

### Figure 2. Rep-CM이 age/photoaging-concordant RNA 및 splicing program을 조절한다

**질문:** phenotype 변화와 연결되는 분자 프로그램은 무엇인가?

필요 자료:

- replicate bulk RNA-seq 또는 total RNA-seq
- 현재 4개 pilot의 candidate catalog는 이 단계의 후보 prior로 사용
- external naturally aged HDF dataset 및 독립 senescence signatures와 방향성 비교
- ΔPSI, junction coverage, gene expression을 분리한 분석
- 상위 event의 RT-PCR/ddPCR 및 가능하면 long-read 확인

**논리 연결:** replicated phenotype -> coherent molecular program.

**오류 방지:** 젊은 방향의 PCA 이동 자체는 rejuvenation이 아니다. cell-cycle genes를 제외한 sensitivity analysis를 하고, 외부 age model은 현재 결과를 보지 않은 상태로 고정한다. acute UVA response와 chronological aging이 일치하지 않는 부분도 보고한다.

**통과 기준:** donor-consistent splice events가 있고, 최소 일부가 기능 축과 dose-response/시간순서상 연결된다.

### Figure 3. Rep-CM의 active fraction과 factor를 규명한다

**질문:** 효과를 전달하는 source-side 물질은 무엇인가?

필요 자료:

- CM proteomics/cytokine/metabolite/particle profiling
- fractionation과 process controls
- candidate neutralization/depletion
- recombinant add-back 및 조합 실험

**논리 연결:** Rep-CM 고유 phenotype -> active biochemical cause.

**오류 방지:** enrichment는 causality가 아니다. neutralization만으로는 nonspecific toxicity 가능성이 있으므로 add-back과 dose-response를 함께 본다. 하나의 factor가 전체 CM을 완전히 재현하지 못하면 조합 또는 permissive factor로 결론을 제한한다.

**통과 기준:** 최소한 necessity 또는 sufficiency가 강하게 성립하고, 이상적으로 둘 다 성립한다.

### Figure 4. Focused perturbation screen과 virtual fibroblast가 causal node를 전향적으로 선별한다

**질문:** ligand에서 isoform/phenotype으로 이어지는 핵심 recipient-side node와 유효 조합은 무엇인가?

필요 자료:

- candidate receptor, signaling node, TF/RBP의 focused CRISPRi/a 또는 arrayed perturbation
- matched scRNA-seq와 핵심 junction assay
- simple baselines 대 virtual-cell model의 held-out 성능
- 새 donor/CM batch에서 unseen perturbation의 prospective validation

**논리 연결:** active factor -> recipient network -> 검증할 mechanism.

**오류 방지:** 모델 attention/feature importance를 기전으로 간주하지 않는다. 모델이 모든 처리군의 평균 변화를 복제해도 높은 correlation이 나올 수 있으므로 perturbation-specific metric을 쓴다. 예측 후 같은 data에서 후보를 재선정하면 prospective test가 아니다.

**통과 기준:** 모델 추천이 baseline 추천보다 실제 hit enrichment를 보이고, hit가 RNA뿐 아니라 functional score를 개선한다.

**예산이 제한되면:** Figure 4의 AI 부분을 제거하고 Figure 3 후보를 직접 arrayed perturbation으로 검증한다. Aging Cell에는 AI label보다 인과기전이 중요하다.

### Figure 5. `[receptor]–[RBP]–[isoform]` 축의 necessity와 sufficiency를 증명한다

**질문:** 특정 splicing 변화가 기능 회복을 매개하는가?

필요 자료:

- early phospho-signaling과 시간 순서
- receptor inhibition/KO
- RBP protein, localization 또는 phosphorylation
- RBP knockdown/overexpression; 가능하면 CLIP/eCLIP 또는 minigene motif evidence
- splice-switching ASO/CRISPR와 isoform-specific rescue
- downstream protein/localization/activity

**논리 연결:** network hit -> direct molecular mechanism -> phenotype.

**오류 방지:** RBP 발현과 ΔPSI의 상관만으로 직접 조절을 주장하지 않는다. isoform manipulation이 total gene abundance도 크게 바꾸면 splice-specific effect와 dosage effect를 분리한다.

**통과 기준:** pathway 차단이 CM-induced isoform과 기능 효과를 모두 없애고, isoform-specific rescue가 phenotype을 복구한다.

### Figure 6. 효과의 지속성, identity 보존 및 tissue relevance를 보인다

**질문:** 이 현상이 배양 dish의 일시적 transcriptional effect를 넘는가?

필요 자료:

- washout 후 장기 추적
- fibroblast identity, karyotype/genomic stability 또는 적절한 safety surrogate
- 3D human skin equivalent 또는 human skin explant에서 ECM/조직 구조/상처회복
- 가능하면 naturally aged donor HDF 또는 다른 senescence trigger에서 일반화

**논리 연결:** 세포 기전 -> aging-relevant tissue function.

**오류 방지:** 3D model에서 bulk collagen 증가가 recipient HDF의 회복인지 단순 세포수 증가인지 구분한다. iPSC-like marker의 유도나 과도한 proliferation을 rejuvenation으로 해석하지 않는다.

**통과 기준:** 조직 수준 기능, 지속성, identity/safety가 같은 방향을 지지한다.

### 원고의 문단 흐름

**Introduction**

1. 피부 노화와 UV 손상은 fibroblast senescence, ECM 저하, altered intercellular communication을 동반하지만, 손상된 세포를 기능적으로 회복시키는 paracrine mechanism은 충분히 규명되지 않았다.
2. cellular reprogramming은 여러 age-associated feature를 낮출 수 있지만 identity loss와 안전성 문제가 있다. 따라서 reprogramming 과정의 유익한 extracellular signal만 분리하는 cell-free 접근이 매력적이다.
3. alternative splicing은 노화와 함께 변하지만, 대부분의 연구는 원인인지 결과인지 밝히지 못했다.
4. 이 연구는 reprogramming-transition secretome이 photoaged HDF의 기능을 회복시키고, 그 효과가 정의된 receptor–RBP–isoform 축에 의해 매개된다는 가설을 검증한다.

Introduction에는 실제 variant 분석이 없다면 AlphaGenome을 넣지 않는다. virtual cell도 논문의 생물학적 gap이 아니라 후보선별 방법으로 Results/Methods에서 등장시키는 편이 자연스럽다.

**Results**

`재현 가능한 phenotype -> replicated molecular program -> active source factor -> model/perturbation-guided recipient node -> causal isoform -> tissue relevance`

각 Results section 첫 문장은 바로 앞 section이 남긴 질문으로 시작한다. 예를 들어 “Rep-CM이 phenotype을 개선했으므로, 다음으로 어떤 molecular program이 이 효과와 동반되는지 물었다”와 같이 이어 간다. 다만 “동반”에서 바로 “매개”로 건너뛰지 않고, 매개성은 Figure 5의 loss/rescue 실험 뒤에만 사용한다.

**Discussion**

1. 가장 먼저 새롭게 증명한 한 가지 causal chain을 요약한다.
2. 왜 transition-state secretome이 HDF-CM 및 fully reprogrammed iPSC-CM과 다른지 source biology를 논의한다.
3. splicing이 단순 biomarker가 아니라 mediator였다는 necessity/sufficiency 근거를 설명한다.
4. virtual fibroblast는 discovery efficiency와 prospective hit enrichment만큼만 주장한다.
5. UVA model, in vitro system, CM 표준화, donor diversity와 장기 안전성의 한계를 명시한다.
6. 마지막에 defined factor/combination으로의 전환 가능성을 제시하되 organismal longevity나 임상효능으로 확장하지 않는다.

---

## 8. 결과에 따른 스토리 분기

### Branch A: 강한 splice mechanism이 나온다

가장 좋은 Aging Cell 스토리다.

`Rep-CM active factor -> receptor/signaling -> RBP activity -> causal isoform switch -> stable HDF/ECM recovery`

이 경우 alternative splicing을 title과 abstract 전면에 둔다.

### Branch B: RNA state와 function은 좋아지지만 causal splice event가 없다

splicing을 억지로 중심에 두지 않는다.

`Rep-CM active factor -> receptor/pathway -> functional recovery`

splicing은 exploratory supplement로 낮춘다. 기전의 깊이가 충분하면 Aging Cell 가능성은 남지만 active factor와 downstream pathway의 gain/loss validation이 필요하다.

### Branch C: Rep-CM과 iPSC-CM 차이는 있으나 healthy/vehicle 대비 회복이 없다

“rejuvenation” 논문이 아니다. reprogramming-stage-specific paracrine biology 또는 secretome comparison 연구로 재프레이밍한다.

### Branch D: CM 효과가 batch에 따라 불안정하다

biomarker-standardized CM manufacturing 또는 active-component identification이 먼저다. 평균 효과를 보고 논문 스토리를 진행하면 재현성 문제가 심각하다.

### Branch E: CM이 senescent cell을 선택적으로 죽인다

이는 rejuvenation이 아니라 senolytic/senomorphic 가능성이다. absolute cell counts, apoptosis, surviving-cell function으로 분리하고 스토리를 다시 세운다.

---

## 9. Aging Cell 및 더 높은 저널을 위한 체크리스트

### Aging Cell에 필요한 최소 패키지

- [ ] age/photoaging 질문이 명확하고 일반 세포보호와 구분됨
- [ ] time-matched healthy 및 damaged vehicle control
- [ ] 독립 recipient donors와 CM source/batches
- [ ] 다중 senescence marker와 fibroblast-specific functional endpoint
- [ ] active CM fraction/factor에 대한 perturbational evidence
- [ ] receptor/downstream mechanism의 gain/loss evidence
- [ ] 핵심 isoform의 orthogonal validation과 causal rescue
- [ ] cell identity, viability, unsafe proliferation 통제
- [ ] 3D skin/ex vivo 또는 강한 독립 human relevance
- [ ] 모든 omics에서 effect size, uncertainty, multiple testing과 independent validation
- [ ] “rejuvenation”의 범위를 데이터가 지지하는 수준으로 제한

### Nature Aging급으로 올라가려면 추가로 필요한 것

- UVA-HDF 한 모델을 넘어 naturally aged donor 또는 여러 senescence triggers에서 보편성
- human tissue/ex vivo와 in vivo에서 기능적 효과
- 개별 factor가 아니라 aging의 일반 원리로 확장되는 새로운 기전
- 장기 효능과 종양/재프로그래밍 안전성
- 독립 cohort 또는 외부 validation
- virtual cell을 넣는다면 실험비용을 실제로 줄이거나 unseen intervention을 성공적으로 발견한 정량적 증거

### 현재 프로젝트의 강점

- paired-end mRNA-seq여서 junction-level exploratory discovery에 적합
- HDF-CM, Rep-CM, iPSC-CM comparator가 있어 source-state-specific hypothesis를 만들 수 있음
- 분석 workflow가 N=1 한계를 명시하고 통계 검정을 끈 점은 연구 진실성 측면에서 매우 좋음
- aging–splicing–intercellular communication을 연결하는 주제는 Aging Cell scope와 잘 맞을 잠재력이 있음

### 가장 큰 위험

1. AI 용어가 biology보다 앞서는 것
2. acute stress recovery를 aging reversal로 부르는 것
3. CM batch/pseudoreplication 문제
4. source factor를 recipient transcriptome만으로 추정하는 것
5. splice event가 phenotype의 원인이라는 실험이 없는 것
6. in vitro UVA model 결과를 organismal aging/longevity로 확장하는 것

---

## 10. 실행 우선순위

### Phase 0. 현재 pilot 정리

1. 기존 rMATS exploratory workflow를 완료한다.
2. 두 primary contrast에서 같은 방향으로 보이는 event와 coverage가 좋은 후보를 만든다.
3. 후보를 “결론”이 아니라 replicated experiment의 assay-development list로 사용한다.

### Phase 1. AI보다 먼저 해야 할 실험

1. time-matched vehicle/healthy control을 포함한 독립 biological repeat
2. Rep-CM phenotype의 재현성과 CM batch effect 측정
3. senescence + ECM function + viability 동시 측정
4. top junction RT-PCR/ddPCR assay 확립

**Go 조건:** Rep-CM 효과가 donor/batch를 넘어 재현되고 phenotype과 splice candidate가 함께 움직임.

### Phase 2. 기전 압축

1. CM proteomics/fractionation
2. active factor neutralization/add-back
3. receptor와 RBP perturbation
4. causal isoform manipulation

**Go 조건:** 최소 하나의 끊김 없는 causal chain이 성립.

### Phase 3. Virtual fibroblast

1. focused perturbation panel과 matched scRNA/targeted splicing data 생성
2. simple baseline과 State/conditional model 비교
3. held-out donor/perturbation 평가
4. top predicted novel interventions의 prospective validation

**원칙:** Phase 1과 2가 실패하면 virtual-cell 적용으로 스토리를 구하지 않는다.

### Phase 4. 저널 상승 실험

1. 3D human skin/explant
2. naturally aged HDF 또는 다른 senescence trigger
3. durability와 safety
4. 필요 시 in vivo 피부 기능

---

## 11. 권장 제목 초안

기전이 완성되기 전:

> **A reprogramming-transition secretome attenuates photoaging-associated phenotypes in human dermal fibroblasts through remodeling of RNA splicing**

factor–RBP–isoform 축이 확인된 후:

> **A reprogramming-transition secretome restores dermal fibroblast function through a `[factor]–[receptor]–[RBP]–[isoform]` axis**

virtual cell이 전향적 발견에 성공한 경우에만:

> **Perturbation-guided virtual fibroblast modeling identifies a secretome-responsive splicing circuit that restores aged cell function**

“longevity”, “precision aging”, “full virtual cell”, “reversal of aging”은 현재 근거 범위를 넘으므로 title에서 제외한다.

---

## 12. 주요 근거 자료

1. AlphaGenome 원논문: Avsec et al. *Nature* (2026), “Advancing regulatory variant effect prediction with AlphaGenome.”  
   https://www.nature.com/articles/s41586-025-10014-0
2. AI Virtual Cell의 정의와 검증 원칙: Bunne et al. *Cell* (2024), “How to build the virtual cell with artificial intelligence.”  
   https://pubmed.ncbi.nlm.nih.gov/39672099/
3. State: Adduri et al. *Cell* (online 2026-08-31), perturbation response prediction across contexts.  
   https://doi.org/10.1016/j.cell.2026.07.052
4. 현재 perturbation model의 baseline 문제: Ahlmann-Eltze et al. *Nature Methods* (2025).  
   https://www.nature.com/articles/s41592-025-02772-6
5. perturbation-specific evaluation: Viñas Torné et al. *Nature Biotechnology* (2026 issue; online 2025), Systema.  
   https://www.nature.com/articles/s41587-025-02777-8
6. Aging Cell 공식 aims and scope.  
   https://onlinelibrary.wiley.com/page/journal/14749726/homepage/productinformation.html
7. Cellular senescence marker consensus: Gorgoulis et al. *Cell* (2019).  
   https://pubmed.ncbi.nlm.nih.gov/31675495/
8. Aging과 alternative splicing의 연결 및 인과성 과제: Deschênes & Chabot, *Aging Cell* (2017).  
   https://onlinelibrary.wiley.com/doi/full/10.1111/acel.12646
9. Aging Cell의 UV-HDF 기계론 논문 사례: Chen et al. *Aging Cell* (2024), METTL14–miR-100-3p–ERRFI1 axis.  
   https://onlinelibrary.wiley.com/doi/full/10.1111/acel.14123
10. human fibroblast age transcriptome reference: Fleischer et al. *Genome Biology* (2018).  
    https://link.springer.com/article/10.1186/s13059-018-1599-6
11. transient reprogramming에서 transcriptome과 fibroblast function을 함께 본 사례: Gill et al. *eLife* (2022).  
    https://doi.org/10.7554/eLife.71624
12. EV를 주장할 때의 최소 원칙: MISEV2023.  
    https://pmc.ncbi.nlm.nih.gov/articles/PMC10850029/

---

## 최종 권고

이 연구는 **“AlphaGenome + virtual cell + Stereo-seq + aging”을 모두 담는 플랫폼 논문**보다, **“reprogramming-transition secretome이 photoaged fibroblast의 causal splicing circuit와 기능을 회복한다”는 한 문장짜리 기계론 논문**으로 설계해야 성공 가능성이 높다.

- AlphaGenome: 현재 core에서 제외
- virtual cell: matched perturbation data와 prospective validation이 가능할 때만 제한적으로 포함
- Stereo-seq: 3D/ex vivo 조직 단계 전에는 제외
- 현재 RNA-seq: 후보 생성과 assay 설계에 사용
- 가장 먼저 투자할 것: 대조군, 독립 donor/CM batch, 기능 phenotype, causal perturbation

AI의 역할은 논문의 주장을 넓히는 것이 아니라, 검증 가능한 causal chain을 더 효율적으로 찾는 데 두어야 한다.
