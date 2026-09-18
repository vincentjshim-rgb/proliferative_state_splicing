---
title: "npj Aging 제출 전 모의 심사"
subtitle: "원고 v5 · Splicing-factor expression tracks proliferation in cultured human fibroblasts but not in skin or muscle"
date: "2026년 9월 15일"
---

# 한눈에 보기

**판단: 현재 원고는 투고 준비가 끝나지 않았다.** 첫 번째 기둥(배양에서 splicing factor의 양은 측정된 분열 속도를 따른다)과 세 번째 기둥(그 결합은 조직으로 옮겨지지 않는다)은 이번 재분석에서도 유지된다. 반면 두 번째 기둥(같은 세포에서 splicing outcome은 나이에 따라 늘고 증식 보정에 흔들리지 않는다)과 Fig. 5의 *COL1A2* 사례는, 공개 메타데이터만 봐도 드러나는 코호트 구성 문제를 통제하면 유지되지 않는다. 이 문제는 리뷰어가 GSE113957의 원 논문(Fleischer 2018)을 열어보기만 해도 찾을 수 있다.

아래 여섯 가지는 투고 전에 반드시 해결해야 한다. 근거 수치는 모두 이번에 원자료로 직접 다시 계산한 값이며, 원고에는 아직 반영하지 않았다.

| | 문제 | 확인한 근거 | 영향받는 곳 |
|---|---|---|---|
| 1 | 142명 코호트에 조로증(HGPS) 환자 9명과 18세 미만 31명이 포함되어 있고, 원고에는 언급이 없다 | SRA 시료 속성 `disease: HGPS` 9건(2–8세). 원고에서 "progeria/HGPS" 검색 0건 | Fig. 1A, 3, 4, 5, 8E, 초록 |
| 2 | 83세 이상 공여자 31명이 모두 Coriell AG 계열이고, splicing outcome의 연령 효과는 이 군집이 만든다 | 83세 미만만 보면 ρ = −0.05(*P* = 0.47). 세포주 출처·시퀀서를 공변량으로 넣으면 *P* = 0.68 | 두 번째 기둥, Fig. 4C–E, 5, 초록·고찰 |
| 3 | *COL1A2*의 원위 acceptor junction은 독립 정렬에서 다중 정렬 read로만 지지된다 | 로컬 STAR 네 라이브러리에서 고유 read 0개, 다중 정렬 read 0–3개, GENCODE v41 주석 없음 | Fig. 5C–D, 10, 본문 |
| 4 | "감소가 남는 39개 유전자"는 조건마다 크게 달라진다 | 31 / 56 / 64 / 1개. 다섯 조건 모두에 남는 유전자 0개 | Fig. 3D, 초록, 고찰의 권고 |
| 5 | Fig. 7 본문에 사실 오류가 있다 | "여덟 개 모두 세포주기 집합(0.62) 이상"이라고 썼지만, 넷은 0.56–0.60 | Fig. 7 본문, 초록, 고찰 |
| 6 | 투고 규정을 지키지 못한 항목이 있다 | 초록 305단어(규정 150). 참고문헌 19는 제출 중인 원고. 윤리·저자 기여·이해상충 없음. 로컬 RNA-seq 미기탁 | 전체 |

반대로, 리뷰어가 가장 먼저 물을 두 가지는 미리 확인해 두었고 **결론이 유지된다.** 328개 라이브러리의 반복 측정을 고려한 혼합 모형에서도 분열 속도와의 결합이 남는다. 64개 대비의 결합은 발현량을 맞춘 무작위 유전자 집합 1,000개 중 어느 것보다도 강하다.

# 편집자 판단 (모의)

**저자께**

주제는 이 저널의 범위에 잘 맞는다. 노화 전사체 연구에서 splicing factor의 양을 노화 지표로 읽는 관행을, 세포 계수로 잰 분열 속도와 사전등록한 조직 비교로 검증했다는 점이 편집자로서 흥미롭다. 결합 자체가 알려져 있다는 것은 서론이 이미 인정하고 있다. 기여의 크기는 "노화 주장이 나오는 계에서 이 결합이 얼마나 크고 어디서 멈추는가"에 달려 있다.

심사로 보내기 전에 두 가지를 요구할 것이다. 첫째, 초록을 150단어 이하로 줄이고 필수 항목(윤리, 데이터·코드 공개, 저자 기여, 이해상충)을 채울 것. 둘째, 관련 원고(Journal of Tissue Engineering 심사 중)의 사본을 커버레터와 함께 낼 것. 겹치는 자료는 Repro-CM 대비 하나와 read 수준 그림 하나뿐이라 중복 게재 문제는 크지 않다고 보지만, 편집자는 직접 확인한다.

심사 결과는 현재 상태라면 **major revision 이상, 거절 위험이 높음**으로 예상한다. 아래 리뷰어 의견 중 1–3번 문제는 리뷰어 한 명이라도 GSE113957 메타데이터를 열면 드러나고, 그 경우 원고의 두 번째 기둥이 무너진다.

# 직접 확인한 재분석

아래 계산은 원고 산출물(`public_data_tierA/derived/`)과 recount3 메타데이터로 다시 한 것이다. 스크립트는 `scripts/presubmission_checks/`에 있다. 원고와 그림에는 아직 반영하지 않았다.

## 142명 코호트에 조로증 환자 9명과 아동 31명이 있다

GSE113957(SRP144355) 시료 143개 중 원고는 142개를 썼다. SRA 시료 속성을 보면 정상 133명과 **HGPS 9명**(나이 2, 3, 3, 4, 5, 8, 8, 8, 8세)이다. 18세 미만은 31명이고 나이 범위는 1–96세다. HGPS 세포의 증식 점수 평균은 0.02로, 10세 미만 정상 공여자(0.47)보다 훨씬 낮다. HGPS는 *LMNA*의 cryptic splice donor가 활성화되어 생기는 병이므로 splicing outcome 분석에 넣을 수 없다.

| 분석 | 원고 (142명) | HGPS 제외 (133명) | 20세 이상 정상 (107명) |
|---|---|---|---|
| 기계 발현의 연령 효과, 보정 전 → 후 | *P* 4 × 10⁻⁶ → 0.035 (73% 감소) | 2 × 10⁻⁶ → 0.15 (82%) | 1 × 10⁻⁴ → 0.21 (82%) |
| 주석 없는 junction 비율, 보정 전 → 후 | *P* 0.002 → 0.004 | 0.002 → 0.004 | 0.02 → 0.04 |
| 증식 보정 후 감소가 남는 유전자 | 39 | 31 | 9 |

첫 번째 기둥은 HGPS를 빼면 오히려 강해진다. 원고의 공개 수치(73%, *P* = 0.035)는 정확히 재현된다.

## 83세 이상은 모두 한 세포주 은행 출신이다

세포주 ID를 보면 83세 이상 31명은 **모두 Coriell NIA Aging Cell Repository의 AG 계열**이다. 20세 미만 정상 공여자 26명 중 24명은 GM 계열이다. 나이와 세포주 출처가 거의 완전히 겹친다. 게다가 AG 계열은 모든 연령대에서 주석 없는 junction 비율이 더 높다.

| 연령대 (정상) | AG 계열 (n) | 그 외 (n) |
|---|---|---|
| 0–20세 | 80.8 (2) | 76.7 (24) |
| 20–40세 | 81.4 (12) | 76.7 (19) |
| 40–60세 | 79.3 (16) | 75.8 (4) |
| 60–83세 | 78.6 (15) | 75.2 (10) |
| 83–96세 | 84.5 (31) | — |

<small>주석 없는 junction read 비율의 중앙값 × 10⁴</small>

같은 모형을 조건별로 다시 적합하면 두 기둥의 운명이 갈린다.

| 조건 | n | 주석 없는 junction: 연령 효과 *P* (보정 전 → 후) | 기계 발현: *P* (보정 전 → 후, 감소율) |
|---|---:|---|---|
| 원고 | 142 | 0.002 → 0.004 | 4 × 10⁻⁶ → 0.035 (73%) |
| 정상 + 세포주 출처·시퀀서 공변량 | 133 | **0.68 → 0.85** | 4 × 10⁻⁹ → 5 × 10⁻⁴ (65%) |
| 정상, 83세 미만 | 102 | **0.47 → 0.95** (ρ = −0.05) | 2 × 10⁻⁴ → 0.054 (69%) |
| 정상, AG 계열, 83세 미만 | 45 | 0.019 → 0.056 (**ρ = −0.35, 반대 방향**) | 0.041 → 0.12 (72%) |
| 정상, 20–82세 | 76 | 0.061 → 0.19 (ρ = −0.20) | 0.0016 → 0.033 (64%) |

PSI event(Fig. 5)도 같다.
- 연령 연관 event는 원고 135개에서, 83세 미만만 보면 18개로 줄어든다.
- 20–82세에서는 19개이고, 증식 보정 후에는 0개다.
- *COL1A2* event 9개의 *P*는 원고 3 × 10⁻⁸–1 × 10⁻³에서, 83세 미만 3 × 10⁻⁵–0.2, 20–82세 0.002–0.6이 된다.
- *P*가 가장 작은 *COL1A2* event(PSI 약 0.04의 원위 acceptor)의 PSI 중앙값은 0–83세 구간에서 0.035–0.049로 평평하다가, 83세 이상에서 0.072로 뛴다.

본문의 "나이에 따라 일관된 몫을 잃는 작고 점진적인 재분배"는 자료와 맞지 않는다. 계단형 변화이고, 그 계단이 세포주 출처와 겹친다. 초고령의 실제 생물학일 수도 있지만 이 자료로는 구분할 수 없다. 세포주 출처를 공변량으로 넣어도 PSI 효과는 남는다(83세 이상이 모두 AG라 공변량이 이 군집을 분리하지 못하기 때문). 그러므로 결정적인 검사는 83세 이상을 뺀 분석이다.

## *COL1A2* 원위 junction은 다중 정렬 read로만 보인다

이 연구에서 독립적으로 정렬한 로컬 라이브러리 네 개(STAR 2.7.3a, GENCODE v41)의 `SJ.out.tab`에서 두 junction을 비교했다.

- 근위 junction(chr7:94,406,304–94,407,846, 주석 있음)은 고유 정렬 read 3,965–6,579개로 지지된다.
- Fig. 5에서 나이와 함께 늘어난다고 한 원위 junction(chr7:94,406,304–94,418,498)은 **주석이 없고, 고유 read가 네 라이브러리 모두 0개**다. HDF-CM과 Repro-CM에서 다중 정렬 read 2개와 3개로만 나타난다.

collagen 유전자는 Gly-X-Y 반복을 암호화하는 짧고 서로 비슷한 exon이 많아서, 정렬기가 비슷한 exon 사이에 가짜 junction을 만들기 쉽다. 리뷰어는 이 junction이 정렬 인공물인지 물을 것이다. 연결된 문제도 두 가지 있다. Fig. 5B 가운데 막대의 이름은 "annotated alternative splicing events"인데 이 event의 짝은 주석이 없다. Fig. 10의 표시 창(94,403,500–94,412,500)은 원위 acceptor 위치를 포함하지 않는다.

## 39개 유전자는 조건마다 달라진다

| 조건 | 연령 연관 | 보정 후 연관 | 감소가 남음 | 보정 후에만 연관 |
|---|---:|---:|---:|---:|
| 원고 (142명) | 129 | 46 | 39 | 7 |
| HGPS 제외 (133명) | 134 | 41 | 31 | 7 |
| + 세포주 출처·시퀀서 공변량 | 152 | 60 | 56 | 2 |
| 83세 미만 (102명) | 124 | 81 | 64 | 12 |
| 20–82세 (76명) | 88 | 1 | 1 | 0 |

다섯 조건에 모두 남는 유전자는 없다. "보정 후에도 감소가 남는 39개 유전자는 대안 지표가 될 수 있다"는 고찰의 권고는 삭제해야 한다. 한편 20–82세 성인에서는 연령 연관 88개 중 1개만 보정을 견딘다. "연령 연관의 대부분이 분열 속도로 설명된다"는 핵심 주장은 성인에서 더 강하다.

## 분열 속도 결합은 반복 측정을 고려해도 유지된다

Fig. 2A의 328개 라이브러리는 세포주 7개(HC1 110, HC2 121, HC3 29, HC4 33, SURF1 세 세포주 12·13·10)에서 반복 채취한 것이다. 원고의 ρ = 0.71, *P* = 1 × 10⁻⁵⁰은 이를 독립 표본처럼 계산했다. 다시 계산해도 결합은 남는다.

| 분석 | 추정치 | *P* |
|---|---|---|
| 세포주별 Spearman ρ | 0.71, 0.71, 0.88, 0.62, 0.76, 0.32, 0.66 | SURF1_2(n = 13)만 0.3, 나머지 < 0.05 |
| 혼합 모형, 세포주 무작위 절편 | β = 0.69 (표준화) | 3 × 10⁻⁴⁶ |
| 무작위 절편 + 기울기 | β = 0.71 | 4 × 10⁻³⁰ |
| 무작위 절편 + 배양 일수 | β = 0.54 | 3 × 10⁻¹⁹ |
| 세포주 × 처리 × 산소 평균 (41개 계열) | ρ = 0.54 | 2 × 10⁻⁴ |
| 무처리 건강 세포주, 21% O₂ (70개, 4개 세포주) | ρ = 0.87, 혼합 모형 β = 0.92 | 8 × 10⁻²⁴ |

원고에는 합친 ρ 대신 혼합 모형 추정치와 세포주별 값을 보고하면 된다.

## 64개 대비 결합은 무작위 유전자 집합보다 강하다

암 signature 연구(Venet 2011)가 보인 것처럼, 발현량이 비슷한 아무 유전자 집합이나 증식을 따라갈 수 있다. 그래서 GSE113957 발현량 구간을 맞춘 무작위 177개 유전자 집합(splicing·세포주기 유전자 제외)을 1,000번 뽑아 같은 점수를 계산했다.

- 무작위 집합: ρ 중앙값 0.27(95% 범위 −0.02–0.49, 최대 0.59), *R*² 중앙값 0.04(최대 0.40)
- 관찰값: ρ = 0.86, *R*² = 0.92. 1,000개 모두보다 크다(경험적 *P* < 0.001).
- senescence 대비 9개를 빼면 ρ = 0.80, *R*² = 0.68(n = 56). *R*² 0.92는 극단의 senescence 대비에 일부 기대고 있다.
- 연구당 한 값으로 평균하면 ρ = 0.79, *R*² = 0.94.
- 같은 방식으로 본 Reactome Metabolism of RNA(세포주기 유전자 제거, 710개)는 ρ = 0.77이다. 결합은 splicing에 한정되지 않고 RNA 처리 전반에 걸친다. extracellular matrix organisation은 ρ = −0.59다.

# 리뷰어 1 — 노화·splicing 생물학

**총평.** 측정된 분열 속도와 사전등록 조직 비교는 이 분야에 드문 설계이고, 실패를 숨기지 않은 점이 좋다. 그러나 splicing outcome과 *COL1A2* 결과는 코호트 구성에 취약하고, 몇몇 해석은 자료보다 앞서 나간다.

## Major

1. <span class="sev c">치명</span> **공여자 코호트의 구성.** 142명 코호트에 HGPS 섬유아세포 9개와 18세 미만 31명이 들어 있는데, 원고에는 이 사실이 없다. HGPS는 *LMNA*의 비정상 splicing으로 생기는 병이라 splicing outcome 분석과 섞을 수 없다. 노화 해석을 하려면 성인 정상 공여자로 분석하고, 아동을 넣는다면 따로 보여야 한다. Fig. 1C와 7B는 22–89세 성인 97명을 쓰고 Fig. 3–5와 8E는 1–96세 142명을 써서 코호트 정의도 서로 다르다.
2. <span class="sev c">치명</span> **세포주 출처와 초고령의 혼동.** 83세 이상이 모두 Coriell AG 계열이다. splicing outcome과 *COL1A2* PSI의 연령 효과는 이 군집을 빼면 사라지거나 반대 방향이 된다. "splicing outcome은 나이에 따라 늘고 증식 보정에 흔들리지 않는다"는 두 번째 핵심 주장을 이 자료로는 지지할 수 없다.
3. <span class="sev m">주요</span> **"splicing precision"이라는 이름.** 주석 없는 junction read 비율은 splicing 정확도 말고도 pre-mRNA와 intron read의 비율, RNA 분해, read 길이, 주석의 불완전성에 좌우된다. 원고는 이를 "a direct index of splicing precision"이라고 부르고, "precision"을 9번, "dysfunction"을 3번 쓴다. recount3 QC 지표(고유 매핑률, intron read 비율 등)를 공변량으로 넣고 "use of unannotated junctions"로 이름을 바꿔야 한다.
4. <span class="sev c">치명</span> **대표 사례 *COL1A2*.** 원위 junction이 독립 정렬에서 다중 정렬 read로만 지지되고, 변화는 83세 이상의 계단형이다. 고유 read만으로 PSI를 다시 계산하지 못한다면 대표 사례에서 내리는 것이 안전하다.
5. <span class="sev m">주요</span> **반대 방향 해리를 보고하지 않았다.** 원고는 oligomycin과 dexamethasone이 분열을 늦추면서 splicing factor를 올린다고만 쓴다. 같은 표(`perturbation_rate_vs_splicing.tsv`)에는 반대 방향 해리도 있다.
   - 3% 산소는 분열 속도 변화 −0.02인데 splicing 점수는 −0.40이다.
   - β-hydroxybutyrate는 −0.07과 −0.32다.
   - mitoNUITs는 −0.20과 −0.31이다.

   본문의 "glucose restriction"은 실제 처리가 2-deoxyglucose(해당과정 억제제)다.
6. <span class="sev m">주요</span> **"splicing 유전자" 집합의 구성.** 177개에는 nucleoporin(*NDC1*, *NUP43*, *NUP50*, *NUP155*), RNA polymerase II 소단위(*POLR2C*, *POLR2F*, *POLR2G*), *UBA52*, *PPP1CA*, *SEC13*처럼 splicing factor가 아닌 유전자가 들어 있다. 보정 후 증가로 뒤집히는 6개 중 4개가 nucleoporin이다. 이름을 "pre-mRNA processing genes"로 통일하고, spliceosome 핵심 유전자로 민감도를 보여야 한다.
7. <span class="sev m">주요</span> **노화 프로그램 감사(Fig. 7)의 유전자 중복.**
   - chromatin organisation, DNA 복구, telomere 유지에는 복제 의존 histone처럼 S기에 전사되는 유전자가 많다. 결합이 이들에서 오는지 senescence 집합에만 분리했으므로, 모든 프로그램에서 Reactome Cell Cycle 유전자를 뺀 결과를 보여야 한다.
   - Reactome Cellular Senescence는 연구자들이 가장 많이 쓰는 senescence 패널이 아니다. SenMayo(Saul 2022) 같은 패널로도 확인해야 "senescence 점수가 반대를 보고한다"는 주장이 공정해진다.
   - 본문의 "All eight equal or exceed the cell-cycle gene set itself (0.62)"는 틀렸다(interferon 0.60, DSB repair 0.59, telomere 0.58, BER 0.56).
8. <span class="sev m">주요</span> **"문헌이 거명하는 splicing factor"의 정의.** Fig. 2B는 28개, Fig. 3A와 9C는 25개이고 구성원이 서로 다르다. Fig. 2B에서 빨간색으로 표시한 "reported in ageing" 8개도 기준과 출처가 없다. 목록 하나를 출처와 함께 정해 모든 그림에 써야 한다.
9. <span class="sev m">주요</span> **개입 분류.** GSE149694의 "Fibroblast-D7 vs Fibroblast-D3"는 배양 대조인데 reprogramming 8개에 들어 있다(Fig. 9C의 "fibroblast d7"). GSE297233의 "O4YRSK"를 "mutant control"로 부른 것도 원 논문(Lu 2025 *Cell*)과 대조해야 한다.

## Minor

- 참고문헌 6(Gorgoulis 2019)과 7(Ocampo 2016)은 "개입을 splicing factor 회복 능력으로 선별하고 회복을 회춘으로 읽는다"는 문장을 직접 뒷받침하지 않는다.
- "almost deterministically"(2회)는 과장이다.
- 조직 splicing 문헌이 빠졌다. GTEx 조직의 연령 연관 splicing 지도(Wang 2018 *Sci. Rep.*)와 노화 근육 단백체의 spliceosome 변화(Ubaida-Mohien 2019 *eLife*)는 "근육의 splicing factor 양을 기각할 수 없다"는 고찰 문장과 직접 연결된다.
- hallmark 프로그램을 고르면서 López-Otín의 hallmark 논문(2013, 2023)을 인용하지 않았다.

# 리뷰어 2 — 통계·계산

**총평.** 전반적으로 전통적 통계를 성실하게 썼다. 다만 헤드라인 *P* 값 몇 개가 독립성 가정을 어기고, 효과 크기의 불확실성이 거의 보고되지 않았다.

## Major

1. <span class="sev m">주요</span> **반복 측정.** 328개 라이브러리는 7개 세포주의 반복 측정이므로 ρ = 0.71, *P* = 1 × 10⁻⁵⁰은 부풀려져 있다. 세포주를 무작위 효과로 둔 모형을 보고해야 한다. (사전 점검: 결론 유지, 위 표 참조)
2. <span class="sev m">주요</span> **대비 수준 결합의 특이성.** 발현량이 비슷한 무작위 집합과의 비교, senescence 대비를 뺀 결과, 같은 연구 안의 대비를 묶은 분석을 보고해야 한다. 같은 대조군을 공유하는 대비는 독립이 아니다. (사전 점검: 무작위 집합보다 강함, senescence 제외 시 *R*² 0.68)
3. <span class="sev m">주요</span> **증식 보정의 해석.** 증식 점수는 연령과 splicing 발현 사이의 매개변수다. 보정 후 남는 효과는 대리 변수의 측정오차로 과소·과대 추정될 수 있다. 같은 라이브러리의 기술적 공변(품질, 조성)이 보정에 흡수될 수도 있다. 발현량을 맞춘 무작위 유전자 집합에 같은 보정을 적용해 연령 효과가 얼마나 줄어드는지 음성 대조로 보여야 한다.
4. <span class="sev m">주요</span> **표본이 작은 상관의 *P*.**
   - HC3(n = 7)의 ρ = −0.96은 정확 검정 양측 *P* ≈ 0.003인데, 본문은 "all *P* < 0.001"이라고 쓴다.
   - HC4(n = 6)의 *P* = 0.072는 점근 근사값이다.
   - n < 10에서는 정확 검정을 써야 한다.
5. <span class="sev m">주요</span> **단측과 양측.** H2의 *P* = 3 × 10⁻³³은 사전등록대로 단측인데, Methods는 "*P* values are two-sided"라고 쓴다.
6. <span class="sev m">주요</span> **사전등록의 검증 가능성.**
   - SHA-256 고정은 내부 기록이라 리뷰어가 시점을 확인할 수 없다.
   - 문서가 한국어라 리뷰어가 읽을 수도 없다.
   - 영어 번역본을 보충 자료로 내고, 원문 해시와 함께 공개 저장소에 올린 뒤 "공개 등록 없이 내부적으로 고정했다"고 정확히 써야 한다.
   - 결과 기록에 있는 실행 전 변경(H4에 붙일 *P*)도 밝혀야 한다.
7. <span class="sev n">사소</span> **recount3 gene sums는 base coverage 합이다.** read 수로 쓴 근거를 적거나 read 길이로 환산한 민감도를 보여야 한다(GSE113957, GTEx 모두 해당).
8. <span class="sev n">사소</span> **효과 크기의 불확실성.** 주요 ρ에 95% 신뢰구간(Fisher *z*)이 없다.

## Minor

- *P* 표기가 그림(9e-04, 8e-01)과 본문(9 × 10⁻⁴)에서 다르고, "*P* = < 0.001"(Fig. 1D)과 "*P* = < 2e-16"(Fig. 7E)은 잘못된 표기다.
- "significant"를 *P* 없이 쓴 문장이 있다: "Entropy of junction usage … significant before adjustment and not after."
- Fig. 6C는 BH 보정 뒤 여섯 종류가 모두 *P* = 0.94로 같아서 이상해 보인다. 원래 *P*를 함께 적어야 한다.
- Fig. 5B는 서로 다른 척도를 한 막대그래프에 섞었다. 27%와 102%는 계수 기반이고 15%는 event 수(20/135) 기반이다.
- Fig. 9B의 secretome은 n = 14인데 본문은 "fifteen … preparations"(Repro-CM 포함)라고 쓴다.
- 같은 177개 유전자 점수의 이름이 그림마다 여섯 가지다(아래 그림 점검).

# 리뷰어 3 — 조직·GTEx·배양 생물학

**총평.** 사전등록한 배양 대 조직 비교는 원고에서 가장 설득력 있는 부분이다. 다만 "결합이 조직에서 약해진다"는 결과에는 생물학적 해석 말고 통계적·구성적 설명도 가능하며, 원고가 이를 배제하지 않았다.

## Major

1. <span class="sev m">주요</span> **범위 제한.** 조직에서는 증식 점수의 분산이 배양보다 훨씬 작으므로, 증식 1단위당 관계가 같아도 상관은 약해진다. 조직별 증식 점수 분산, 회귀 기울기, 분산을 맞춘 부분표본의 상관을 보여야 한다.
2. <span class="sev m">주요</span> **세포 구성.** 피부의 증식 신호는 표피 기저 keratinocyte에서 오고 splicing 발현은 모든 세포에서 온다. 결과 기록(S4)에서는 연령 효과에만 구성 보정을 했다. 결합 자체도 keratinocyte·섬유아세포 marker로 보정해야 한다.
3. <span class="sev m">주요</span> **Fig. 8D에 양성 대조가 없다.** "공여자의 배양 세포와 피부가 25개 프로그램 중 하나도 일치하지 않는다"는 결과는 측정이 원래 잡음이 많아서 생길 수도 있다. 유전형이 발현을 크게 좌우하는 유전자(강한 cis-eQTL·결실 다형성 유전자)나 성염색체 유전자에서는 같은 공여자의 일치가 검출된다는 것을 먼저 보여야, "배양이 공여자를 지운다"는 해석이 성립한다.
4. <span class="sev m">주요</span> **제목과 결과의 충돌.** 일광 노출 피부의 partial ρ = 0.23은 유의하다(*P* ≈ 7 × 10⁻¹⁰). "but not in skin"은 결과와 맞지 않는다.
5. <span class="sev n">사소</span> **"independent cohorts"(Fig. 8E).** GSE113957은 splicing 집합을 정의한 자료이기도 해서, 집합 정의에 대해서는 독립이 아니다. 그 쪽 추정치에는 HGPS와 아동도 포함되어 있다.
6. <span class="sev n">사소</span> **배양 이탈 표본 32개(6%)** 때문에 Fig. 8B 배양 패널이 비선형이다. 포함·제외 두 결과를 이미 legend에 적었으니 본문에도 짧게 적는 것이 좋다.
7. <span class="sev n">사소</span> **표본 수.** Fig. 8A는 배양 511명, 본문과 8E는 492명이다(공변량이 완비된 수). legend에 둘 다 적어야 한다. 비노출 피부(626명)의 n은 legend에 없다.
8. <span class="sev n">사소</span> 조직의 splicing outcome은 GTEx junction 자료로 측정할 수 있다. 할 수 없는 이유가 자료 접근이 아니라 범위라면 그렇게 적는 것이 좋다.

# 그림 점검

## 모든 그림에 공통

- **패널 표기.** 규정은 소문자 굵은 a, b, c다. 현재는 대문자 A–E다.
- **기호.** "rho", "R2"를 ρ, *R*²로 바꾸고, 음수에는 하이픈 대신 수학 minus(−)를 쓴다. 규정은 그리스 문자를 symbol 서체로 쓰라고 한다.
- ***P* 표기.** "9e-04", "8e-01", "1e-03"을 본문처럼 9 × 10⁻⁴, 0.8, 0.001로 쓴다.
- **같은 점수의 이름.** 177개 유전자 점수를 Fig. 1은 "pre-mRNA processing score", Fig. 2는 "splicing-factor expression", Fig. 4B는 "pre-mRNA processing machinery expression", Fig. 6은 "pre-mRNA processing change", Fig. 9B는 "change in splicing factors"라고 부른다. Fig. 8은 다른 판(96개)이다.
- **글자 크기.** 가장 작은 글자가 5.2–5.6 pt다(Fig. 3A 유전자명, Fig. 3D 축, 여러 지시선 라벨). 규정은 인쇄 크기에서 8 pt를 권장한다.
- **색.** Fig. 8F는 판정을 초록과 빨강으로 구분한다. 규정은 적록 대비를 피하라고 한다.
- **그림 수.** 10개다. 규정에 상한은 없지만 "불필요한 그림은 피하라"고 한다. Fig. 10은 보충 자료로 옮기는 것이 자연스럽다.

## 그림별

| 그림 | 문제 | 고칠 것 |
|---|---|---|
| Fig. 1 | A: "Donor age 1 to 96 years"(C는 22–89세). D: "*P* = < 0.001" 오표기, HC3 n = 7의 *P*. B: FDR e-표기, x축이 108까지 가는데 100 눈금이 없음 | 코호트 정의 후 A 수정. 정확 검정. 표기 통일 |
| Fig. 2 | A: 반복 측정을 합친 ρ와 *P* = 1e-50. B: 검정 방법·n이 legend에 없고 빨강 기준의 출처가 없음, 28개 목록이 Fig. 3·9와 다름. C: *HNRNPD*와 *HNRNPA1* 라벨이 겹치고 ρ에 *P*·n이 없음. D: 공여자를 무시한 직선 적합 | 혼합 모형. 목록 통일과 출처. 라벨 겹침 해소. 세포주별 선 |
| Fig. 3 | 142명(HGPS·아동 포함) 기반이라 전체 재계산 필요. D: 보정 후에만 연관이 생긴 7개가 그림에서 구분되지 않음. 가장 작은 글자 5.2 pt | 코호트 확정 후 재작성. 7개를 기호로 구분 |
| Fig. 4 | A: "actually cuts"와 "annotated", "is transcribed"와 유전자명이 겹침, "1-96 years". C: "r = −0.23" 라벨이 점과 겹침. D: "nothing lost"와 "*P* = 0.004"가 겹치고, 범례의 회색 견본이 빨강·파랑 막대와 맞지 않음. E: 83세 이상 군집 | 두 번째 기둥 결정 후 재작성 |
| Fig. 5 | B: 계수 기반(27%, 102%)과 event 수 기반(15%)을 한 그래프에 섞음, "annotated" 표기가 사실과 다름. C: 유전자명이 이탤릭이 아니고 좌표계(GRCh38)가 없음. D: 83세 이상 계단형, *P*·n 없음. E: 그림 항목(post-translational phosphorylation, IGFBP, 혈소판)과 본문 항목(collagen degradation·biosynthesis, elastic fibre, TGF-β, MET)이 다름 | 두 번째 기둥 결정 후 재작성. 본문과 그림의 항목 일치 |
| Fig. 6 | 축 "cell-cycle module change"("module" 표현, 단위 없음). A: "UVA, 3 days"와 "stress-induced senescence" 라벨이 점과 겹치고 지시선이 교차, "R2". B: 패널 제목("metabolic")이 종류 이름("metabolic / culture")과 다름. C: 모든 종류 *P* = 0.94, "+0.0" 눈금, "no class of perturbation retains an effect"는 "차이가 검출되지 않음"이 정확함 | 축 이름·단위, 라벨 재배치, 원래 *P* 병기 |
| Fig. 7 | **B: FDR ≥ 0.05인 11개 프로그램(cellular senescence 포함)의 점이 보이지 않음.** alpha = 0이 채움과 테두리를 모두 투명하게 만들어, 이름표 7개의 지시선이 빈 곳을 가리킴. B: 설명 글자가 해당 음영 상자에서 멀리 있음. A: 범례가 값 라벨(+0.21 등)과 겹침. C: "150-gene average" 라벨이 violin과 겹침. E: "*P* = < 2e-16", "*P* = 8e-01" | 채움만 투명한 빈 원으로. 범례 위치. 본문 오류와 함께 수정 |
| Fig. 8 | A: 배양 n = 511 대 본문 492. B: 피부 패널 라벨이 점과 겹침. E: *P* 라벨이 신뢰구간 선과 겹침, "142 donors, 1-96 years". F: 초록·빨강 판정색, ">=", H2가 단측임이 드러나지 않음 | n 병기, 라벨 위치, 색각이상 친화 색과 기호 |
| Fig. 9 | A: "t2iLGoY-D13"(Fig. 6은 "naive medium, day 13", C는 "naive (t2iLGoY)"). B: secretome n = 14 대 본문 fifteen. C: 배양 대조 "fibroblast d7"을 reprogramming 열에 넣음, "repro-phase CM" 열은 대부분 미측정 | 명칭 통일, n 명시, 분류 재검토 |
| Fig. 10 | 창이 원위 acceptor를 포함하지 않음. 본문의 "20 read 기준에 도달하는 대체 acceptor가 없다"는 사실이지만, 원위 junction이 다중 정렬 read로만 나타난다는 점이 빠짐. 유전자명 이탤릭 아님 | 보충 자료로 옮기고 문장 보강 |

# 참고문헌 점검

PubMed E-utilities와 Crossref로 30개의 제목, 저자, 학술지, 권, 쪽, 연도를 대조했다. **29개는 서지사항이 맞다.**

| # | 문헌 | 결과 |
|---:|---|---|
| 1 | Harries 2011 *Aging Cell* | 일치 (PMID 21668623) |
| 2 | Deschênes & Chabot 2017 *Aging Cell* | 일치 (28703423) |
| 3 | Casella 2019 *Nucleic Acids Res.* | 일치 (31251810). 같은 제목의 정정 공고(31612919)가 따로 있으니 원 논문을 인용할 것 |
| 4 | Holly 2013 *Mech. Ageing Dev.* | 일치 (23747814) |
| 5 | Donega 2026 *Mol. Cell. Biol.* | 서지 일치 (41772759). **저자가 4명이라 "et al."을 쓸 수 없음**: Donega, S., Gorospe, M., Harries, L. W. & Ferrucci, L. |
| 6 | Gorgoulis 2019 *Cell* | 일치 (31675495). 인용 위치의 주장을 직접 뒷받침하지 않음 |
| 7 | Ocampo 2016 *Cell* | 일치 (27984723). 위와 같음 |
| 8 | Latorre 2017 *BMC Cell Biol.* | 일치 (29041897) |
| 9 | Zhang 2026 bioRxiv | 일치 (42427714). DOI 접두어 10.64898은 openRxiv의 새 접두어로 확인, 2026-06-29 게시, 저자 5명 |
| 10 | Whitfield 2002 *Mol. Biol. Cell* | 일치 (12058064) |
| 11 | Koh 2015 *Nature* | 일치 (25970242) |
| 12 | Hsu 2015 *Nature* | 일치 (26331541) |
| 13 | Kwon 2021 *FASEB J.* | 일치 (33337569) |
| 14 | Sturm 2023 *Commun. Biol.* | 일치 (36635485). 자료 기술 논문 Sturm 2022 *Sci. Data*(36463290)를 함께 인용할 것 |
| 15 | Fleischer 2018 *Genome Biol.* | 일치 (30567591) |
| 16 | Wilks 2021 *Genome Biol.* | 일치 (34844637) |
| 17 | GTEx Consortium 2020 *Science* | 일치 (32913098) |
| 18 | Liu 2020 *Nature* | 일치 (32939092) |
| 19 | Shim et al. *J. Tissue Eng.* | **규정 위반.** 심사 중인 원고(JTE-Aug-26-0256, 2026-08-19 제출)는 번호 목록에 넣을 수 없음. 본문에 저자와 함께 "manuscript submitted"로 쓰고 사본을 편집자에게 제출 |
| 20 | Endicott 2022 *Nat. Commun.* | 일치 (36347867) |
| 21 | Sturm 2019 *Epigenetics* | 일치 (31156022) |
| 22 | Venet 2011 *PLoS Comput. Biol.* | 일치 (22028643) |
| 23 | Chatsirisupachai 2019 *Aging Cell* | 일치 (31560156) |
| 24 | Robinson & Oshlack 2010 *Genome Biol.* | 일치 (20196867) |
| 25 | Robinson, McCarthy & Smyth 2010 *Bioinformatics* | 일치 (19910308) |
| 26 | Ragueneau 2026 *Nucleic Acids Res.* | 일치 (41251150) |
| 27 | Benjamini & Hochberg 1995 *J. R. Stat. Soc. B* | 일치 (Crossref, 57(1), 289–300) |
| 28 | Ritchie 2015 *Nucleic Acids Res.* | 일치 (25605792) |
| 29 | Dobin 2013 *Bioinformatics* | 일치 (23104886) |
| 30 | Danecek 2021 *GigaScience* | 일치 (33590861) |

## 추가해야 할 인용

**GEO 자료의 원 논문.** 대비 64개를 만든 23개 series 중 어느 원 논문도 인용되지 않았다. 규정상 자료는 accession과 함께 참고문헌에 넣을 수 있다. NCBI GEO에 연결된 PubMed 기록은 다음과 같다.

| Series | 원 논문 |
|---|---|
| GSE93535 | Lämmermann 2018 *npj Aging Mech. Dis.* (29675264) |
| GSE109700 | De Cecco 2019 *Nature* (30728521) |
| GSE139563 | López-Antona 2022 *Aging Cell* (35266275) |
| GSE191055 | Hasegawa 2023 *Cell* (37001502) |
| GSE222414 | Yang 2024 *MedComm* (38919335) |
| GSE251807 | Lei 2025 *Sci. Rep.* (40461530) |
| GSE266052 | Brizio 2024 *Stem Cell Res. Ther.* (39334258) |
| GSE268248 | Imani 2024 *Cells* (39120336) |
| GSE279804 | Yao 2025 *J. Cell. Mol. Med.* (41058030) |
| GSE282054 | Abdelmohsen 2026 *Aging Cell* (41521531) |
| GSE293186 | Yuan 2026 *J. Invest. Dermatol.* (41161638) |
| GSE297233 | Lu 2025 *Cell* (40816266) |
| GSE306957 | Victorelli 2025 *Nat. Commun.* (41398033) |
| GSE307377 | Adams 2025 *Research Square* 프리프린트 (41333395), 연결된 두 번째 기록 41000743은 확인 필요 |
| GSE326951 | Zhang 2026 *Pharm. Biol.* (42165632) |
| GSE329475 | Zhong 2026 *Cells* (42645226) |
| GSE116968, GSE134533, GSE240486, GSE306748 | 연결된 논문 없음. GEO accession 자체를 자료 인용으로 |
| GSE113957, GSE149694, GSE179848 | 이미 인용됨(15, 18, 14) |

**본문 주장에 필요한 문헌.**
- López-Otín 2013 *Cell* 153, 1194–1217 (23746838), 2023 *Cell* 186, 243–278 (36599349): Fig. 7의 hallmark 선택
- Wang 2018 *Sci. Rep.* 8, 10929 (30026530): GTEx 조직의 연령 연관 splicing
- Ubaida-Mohien 2019 *eLife* 8, e49874 (31642809): 노화 인간 근육의 spliceosome 단백질 변화
- Eriksson 2003 *Nature* 423, 293–298 (12714972): HGPS 시료를 제외하는 근거로 언급할 경우
- Saul 2022 *Nat. Commun.* 13, 4827 (35974106): SenMayo로 확인할 경우

모두 넣으면 약 55개로 권장 상한 60개 안이다. 자료 원 논문은 보충 표 1의 참고문헌으로 옮겨 본문 목록을 35개 안팎으로 유지하는 방법도 있다.

# 투고 규정 대조

npj Aging의 Content types와 Submission guidelines 페이지(2026-09-15 확인) 기준이다. 규정에 따르면 첫 투고에서는 형식을 완벽히 맞추지 않아도 되고 수락 시점에 적용된다. 그래도 초록 길이와 필수 항목은 첫 인상에 직결된다.

| 항목 | 규정 | 현재 원고 | 조치 |
|---|---|---|---|
| 제목 | 15단어 이하, 구두점·관용구·말장난 없이 | 14단어, 하이픈 1개, "but not in skin"이 결과와 충돌 | 하이픈 없는 15단어 이하로 |
| 초록 | 150단어 이하, 소제목 없음 | **305단어** | 절반으로 |
| 서론·고찰 | 소제목 없음. 고찰에 한계·결론 절 금지 | 준수(한계는 문단) | 유지 |
| Methods | 소제목. 인체 유래 시료·자료는 윤리 승인과 동의 필수. Statistics and reproducibility 절 | 윤리 항목 없음 | 로컬 섬유아세포의 출처·승인 기재 |
| 그림 설명 | 제목 + 350단어 이하, 패널별 설명, 기호 대신 말 | 최대 224단어(Fig. 8) | 준수 |
| 그림 | 소문자 굵은 a, b. Arial/Helvetica. 그리스 문자는 symbol 서체. 적록 대비 피하기. 300 dpi 이상 | 대문자, "rho" 텍스트, Fig. 8F 적록, 600 dpi | 수정 |
| 참고문헌 | 60개 권장. 제출 중 원고는 목록 제외. 5명 초과만 et al. | 30개. 19번 제출 중. 5번 저자 4명에 et al. | 수정 |
| Data availability | 필수, Methods 뒤 별도 절 | "Data and code availability"로 묶임, 수락 후 공개 예정. **로컬 RNA-seq 미기탁** | GEO 기탁 후 accession 기재 |
| Code availability | 핵심 코드면 Methods 안에 필수, 편집자·리뷰어가 접근 가능해야 함 | 수락 후 공개 예정 | Zenodo·GitHub에 투고 전 공개 |
| Author contributions | 필수, 이니셜 | 없음 | 추가 |
| Competing interests | 필수 | 없음 | 추가 |
| Acknowledgements | 연구비는 여기에 | 없음 | 추가 |
| Supplementary Information | 단일 PDF. 큰 표는 "Supplementary Data N". Supplementary Methods 불가 | 보충 표 2·3이 프로젝트 파일 경로로만 지정됨. **보충 표 1의 대비 설명이 잘림**("DEX vs Contro", "RSeT-D13 vs Fi", "(both TGF-b;") | 파일로 만들고 명칭 변경, 잘린 칸 복원 |
| 관련 원고 | 겹치는 저자의 심사 중 원고 사본 제출, 커버레터에 명시 | JTE 원고 | 사본과 커버레터 |
| Reporting Summary · editorial policy checklist | Reporting Summary는 수정본에서 필수이고 투고 시 권장. editorial policy checklist는 심사로 보내기 전에 필수 | 없음 | 둘 다 작성 |
| 사전등록 문서 | (리뷰어 관점) 리뷰어가 읽을 수 있어야 함 | 한국어 | 영어판을 보충 자료로 |

# 고치는 순서 제안

앞 단계의 결정이 뒤 단계의 수치를 바꾸므로 이 순서를 권한다.

1. **코호트 정의를 확정한다.** GSE113957에서 HGPS 9명을 뺀다. 주 분석을 성인(20세 이상)으로 할지, 전체로 하고 성인을 민감도로 둘지 정한다. 세포주 출처(AG/GM)와 시퀀서를 공변량에 넣는다. 그 뒤 Fig. 3, 4, 5, 8E와 Fig. 1A 문구를 다시 계산한다.
2. **두 번째 기둥을 결정한다.** 권하는 방향은 (a)다.
   - (a) splicing outcome과 *COL1A2*를 핵심 주장에서 내린다. Fig. 4는 "기계 발현의 연령 효과 대부분이 증식으로 설명된다"만 남긴다. outcome 결과와 83세 이상 군집은 한 패널이나 보충 자료에 사실대로 보인다. 원고의 메시지는 "배양에서 splicing factor의 양은 증식 상태의 읽기값이고, 이 결합은 조직으로 옮겨지지 않는다"가 되고, Fig. 2, 6, 7, 8, 9가 이를 받친다.
   - (b) GTEx 배양 섬유아세포 492명의 junction 자료로 splicing outcome의 연령 효과를 재현해 본다. recount3 GTEx junction 파일을 새로 받아야 하므로, 용량과 목적을 먼저 확인하고 승인을 받은 뒤에만 진행한다.
3. **원고 문장을 고친다.** 39개 유전자 권고 삭제. Fig. 7의 "모두 0.62 이상" 수정과 초록·고찰 표현 조정. 2-deoxyglucose 명칭. 반대 방향 해리(3% O₂, β-hydroxybutyrate, mitoNUITs) 보고. "precision", "almost deterministically" 표현. 단측·양측. 작은 n의 정확 *P*. 177개 집합의 이름 통일.
4. **통계 보강을 넣는다.** Fig. 2는 혼합 모형과 세포주별 ρ. Fig. 6은 무작위 집합 비교와 senescence 제외 결과. 주요 ρ에 95% 신뢰구간.
5. **리뷰어 3의 조직 쪽 요구를 반영한다.** 증식 점수 분산·기울기, 구성 보정 후 결합, Fig. 8D의 양성 대조. 모두 이미 받은 GTEx 자료로 가능하다.
6. **그림을 한꺼번에 고친다.** 소문자 패널 표기, ρ·*R*²·minus, *P* 표기, 라벨 겹침, Fig. 7B 버그, Fig. 8F 색, 점수 이름 통일, Fig. 10을 보충 자료로.
7. **참고문헌을 고친다.** 5번 저자 목록, 19번을 본문 인용으로, GEO 원 논문·hallmark·조직 splicing 문헌 추가, 6·7번 교체 검토.
8. **투고 필수 항목을 채운다.** 150단어 초록과 제목, 윤리, 로컬 RNA-seq GEO 기탁과 코드 공개, 저자 기여·이해상충·감사, 보충 자료 단일 PDF와 Supplementary Data 명칭, 사전등록 영어판, Reporting Summary, 관련 원고를 밝힌 커버레터.
9. **한글판과 아티팩트를 맞춘다.**

# 부록: 확인 범위와 재현

## 이번에 확인하지 않은 것

- Fig. 5·7의 모든 수치를 원자료에서 하나씩 다시 계산하지는 않았다. 코호트 문제와 관련된 수치, 본문의 대표 수치만 대조했다.
- GSE307377에 연결된 두 번째 PubMed 기록(41000743)과 GSE297233의 "O4YRSK" 대비의 성격은 원 논문으로 확인해야 한다.
- recount3의 junction 수에 다중 정렬 read가 포함되는지는 recount3 문서로 확인해야 한다.
- 로컬 섬유아세포의 윤리·출처 정보는 JTE 원고에서 옮겨야 하며 여기서는 확인하지 않았다.

## 재현 스크립트 (`scripts/presubmission_checks/`)

| 파일 | 내용 |
|---|---|
| `hgps_sensitivity.R` | HGPS 제외·성인 한정 시 기계 발현, splicing outcome, 유전자 수준 결과 |
| `batch_sensitivity.R` | 세포주 출처·시퀀서 공변량, 83세 미만, AG 계열 한정 분석 |
| `gene_batch.R` | 조건별 "감소가 남는 유전자"와 공통 유전자 |
| `psi_batch.R` | PSI event와 *COL1A2* event의 조건별 결과 |
| `diag_pseudorep.R` | 328개 라이브러리의 세포주별 상관과 혼합 모형 |
| `diag_contrast_null.R` | 64개 대비의 발현량 일치 무작위 집합 비교, senescence 제외, 연구별 평균 |
| `refcheck2.py`, `lithelp.py` | 참고문헌 PubMed 대조 |
| `gse113957_status.tsv`, `gse113957_meta2.tsv` | SRA 시료 속성에서 뽑은 질병, 세포주 출처, 시퀀서 |
