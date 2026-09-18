---
title: "노화 프로그램 감사(Fig. 7) 재계산"
subtitle: "코호트 재정의와 cell-cycle 유전자 제거 · 2026년 9월 15일"
---

# 요약

Fig. 7은 두 가지를 다시 계산했다. 하나는 공여자 연령 쪽으로, 제출 전 모의 심사에서 확정한 코호트(GSE113957 정상 공여자 20세 이상 107명, HGPS 제외, 세포주 출처·시퀀서·성별 공변량)에서 recount3 gene sums로 다시 점수를 냈다. 다른 하나는 분열 속도 쪽으로, 각 프로그램에서 Reactome cell-cycle 유전자를 모두 뺀 뒤 측정된 분열 속도와의 상관을 다시 냈다(리뷰어 1 major 7a).

결론 세 가지.

- **splicing 쪽 결합은 프로그램 안의 cell-cycle 유전자 때문이 아니다.** pre-mRNA processing은 ρ 0.65 → 0.67, mRNA splicing은 0.63 → 0.64로 오히려 그대로다. 원고의 핵심 결합은 이 점검을 통과한다.
- **여덟 개 중 셋은 통과하지 못한다.** cellular senescence 0.67 → 0.38, base excision repair 0.56 → 0.35, interferon signalling 0.60 → 0.49로 내려간다. telomere maintenance는 82개 유전자 전부가 cell-cycle 집합에 들어 있어 판정 자체가 불가능하다.
- **collagen formation의 연령 신호는 세포주 출처를 공변량에 넣으면 사라진다**(*P* 0.006 → 0.21). 원고가 "분열과 무관하지만 연령과 강하게 연관된 프로그램"의 대표로 쓰던 항목이 새 코호트에서 유지되지 않는다. unfolded protein response는 약하게 남는다(*P* = 0.044, FDR = 0.07).

# 1. 분열 속도: cell-cycle 유전자를 뺀 뒤

측정된 분열 속도 328개 라이브러리(GSE179848), 각 프로그램에서 Reactome `Cell Cycle`, `Cell Cycle, Mitotic`, `DNA Replication`, `M Phase`의 합집합(668개 유전자)에 속한 유전자를 뺀 뒤 같은 점수를 다시 계산했다.

| 프로그램 | 유전자 | 제거 | 남음 | ρ (원래) | ρ (제거 후) | 판정 |
|---|---:|---:|---:|---:|---:|---|
| chromatin organisation | 312 | 38 | 274 | 0.69 | **0.63** | 유지 |
| pre-mRNA processing | 287 | 46 | 241 | 0.65 | **0.67** | 유지 |
| mRNA splicing | 198 | 12 | 186 | 0.63 | **0.64** | 유지 |
| DNA double-strand repair | 180 | 124 | 56 | 0.59 | **0.55** | 유지 |
| interferon signalling | 279 | 67 | 212 | 0.60 | 0.49 | 약화 |
| cellular senescence | 150 | 82 | 68 | 0.67 | 0.38 | 약화 |
| base excision repair | 63 | 44 | 19 | 0.56 | 0.35 | 약화 |
| telomere maintenance | 82 | 82 | 0 | 0.58 | — | 판정 불가 |
| cell cycle (양성 대조) | 499 | 499 | 0 | 0.62 | — | 판정 불가 |
| rRNA processing | 187 | 9 | 178 | 0.46 | 0.45 | |
| mitochondrial translation | 100 | 0 | 100 | 0.45 | 0.45 | |
| translation | 348 | 39 | 309 | 0.43 | 0.42 | |
| extracellular matrix | 235 | 1 | 234 | −0.33 | −0.34 | |
| collagen formation | 71 | 0 | 71 | −0.02 | −0.02 | |
| unfolded protein response | 89 | 3 | 86 | 0.05 | 0.05 | |

전체 24개 프로그램 표는 `public_data_tierA/derived/hallmark_revised/division_rate_cellcycle_removed.tsv`에 있다.

제거 후에도 |ρ| > 0.5를 지키는 프로그램은 chromatin organisation, pre-mRNA processing, mRNA splicing, DNA double-strand repair 넷이다.

## 본문 수정 문장

원고 본문의 "All eight equal or exceed the cell-cycle gene set itself (0.62), which is the natural ceiling for a division readout."는 사실과 다르다(interferon 0.60, DSB repair 0.59, telomere 0.58, BER 0.56). 다음으로 바꾼다.

> Four of the eight reach or exceed the cell-cycle gene set itself (0.62), and the remaining four sit just below it (0.56 to 0.60). The coupling is not carried by cell-cycle genes inside these sets: removing every Reactome cell-cycle gene from each programme leaves pre-mRNA processing at ρ = 0.67 and mRNA splicing at 0.64, and chromatin organisation (0.63) and double-strand break repair (0.55) above 0.5, while cellular senescence falls to 0.38, interferon signalling to 0.49 and base excision repair to 0.35; telomere maintenance cannot be assessed this way because all 82 of its measured genes belong to the cell-cycle union.

# 2. 공여자 연령: 새 코호트

정상 성인 107명(20–96세), recount3 gene sums, 공변량은 log depth, 세포주 출처(AG 대 그 외), 시퀀서, 성별이다. 발표판은 FPKM 97명(22–89세)에 성별만 공변량으로 넣었다.

| 프로그램 | ρ 분열 | ρ 연령(발표판) | FDR(발표판) | ρ 연령(새 코호트) | 10년당 계수 [95% CI] | FDR |
|---|---:|---:|---:|---:|---|---:|
| chromatin organisation | 0.69 | −0.24 | 0.030 | −0.51 | −0.038 [−0.055, −0.021] | 8 × 10⁻⁵ |
| cellular senescence | 0.67 | −0.25 | 0.051 | −0.37 | −0.026 [−0.042, −0.009] | 0.004 |
| pre-mRNA processing | 0.65 | −0.31 | 0.005 | −0.49 | −0.096 [−0.131, −0.061] | 5 × 10⁻⁶ |
| mRNA splicing | 0.63 | −0.33 | 0.004 | −0.51 | −0.098 [−0.132, −0.063] | 5 × 10⁻⁶ |
| cell cycle (양성 대조) | 0.62 | −0.33 | 0.004 | −0.47 | −0.088 [−0.128, −0.048] | 1 × 10⁻⁴ |
| interferon signalling | 0.60 | −0.25 | 0.048 | −0.30 | −0.027 [−0.044, −0.010] | 0.004 |
| DNA double-strand repair | 0.59 | −0.32 | 0.005 | −0.47 | −0.099 [−0.141, −0.057] | 7 × 10⁻⁵ |
| telomere maintenance | 0.58 | −0.34 | 0.005 | −0.49 | −0.096 [−0.138, −0.054] | 8 × 10⁻⁵ |
| base excision repair | 0.56 | −0.33 | 0.004 | −0.50 | −0.097 [−0.140, −0.054] | 9 × 10⁻⁵ |
| unfolded protein response | 0.05 | −0.39 | 0.004 | −0.43 | −0.020 [−0.039, −0.001] | **0.070** |
| collagen formation | −0.02 | −0.44 | 0.0002 | −0.44 | −0.015 [−0.039, +0.009] | **0.247** |
| extracellular matrix | −0.33 | −0.22 | 0.051 | −0.18 | +0.013 [−0.005, +0.030] | 0.209 |
| NF-kB / TNF | −0.14 | +0.03 | 0.86 | +0.13 | +0.028 [+0.010, +0.046] | 0.004 |
| insulin / IGF signalling | −0.27 | +0.16 | 0.11 | +0.08 | +0.020 [+0.002, +0.037] | 0.044 |

전체 표는 `public_data_tierA/derived/hallmark_revised/donor_age_primary_cohort.tsv`에 있고, 증식 보정 후 계수(`beta_adj`, `p_adj`, `FDR_adj`)도 같은 파일에 있다.

## collagen formation이 달라진 이유

같은 107명에서 공변량만 바꿔 보면 원인이 분명하다.

| 프로그램 | log depth만 | + 세포주 출처·시퀀서·성별 |
|---|---|---|
| collagen formation | −0.033 (*P* = 0.006) | −0.015 (*P* = 0.21) |
| unfolded protein response | −0.022 (*P* = 0.013) | −0.020 (*P* = 0.044) |
| mRNA splicing | −0.073 (*P* = 4 × 10⁻⁵) | −0.098 (*P* = 2 × 10⁻⁷) |
| chromatin organisation | −0.033 (*P* = 7 × 10⁻⁵) | −0.038 (*P* = 2 × 10⁻⁵) |

collagen formation의 연령 신호는 세포주 출처(AG 대 그 외)와 겹쳐 있어 공변량을 넣으면 절반으로 줄고 유의하지 않다. 반대로 splicing과 chromatin 쪽은 공변량을 넣으면 오히려 강해진다. 즉 이 코호트에서 연령과 세포주 출처의 혼동은 collagen 쪽 주장을 부풀리고 splicing 쪽 주장을 가리고 있었다.

이는 이미 기록된 음성 결과(GTEx 배양에서 collagen 연령 신호 불재현, CLAUDE.md §8)와 같은 방향이다. 원고는 collagen formation을 "분열과 무관하지만 연령과 연관된" 대표 사례로 더 이상 쓸 수 없다.

## 증식 보정 후

새 코호트에서 프로그램 점수에 증식 점수를 공변량으로 넣으면, 24개 중 연령 연관이 FDR < 0.05로 남는 것은 pre-mRNA processing(FDR = 0.040)과 mRNA splicing(0.040) 둘뿐이다. chromatin organisation은 *P* = 0.037이지만 FDR 0.25다. 이는 시료 수준 기계 발현에서 보정 후 −0.040/10년(*P* = 0.005)이 남는 것과 같은 방향이다.

# 3. senescence 집합

분열 속도 쪽 해부는 그대로다. Reactome Cellular Senescence 150개 중 82개가 cell-cycle 합집합에 속하고, 이들의 유전자별 ρ 중앙값은 0.26, 나머지 68개는 0.07이다(Wilcoxon *P* = 0.003). cell-cycle 유전자를 뺀 68개만으로 낸 집합 점수는 ρ = +0.38로, 여전히 양수다. 즉 "빠르게 분열하는 배양에서 senescence 점수가 높다"는 결과는 cell-cycle 유전자를 빼도 방향이 바뀌지 않지만 크기는 절반 가까이 줄어든다.

**SenMayo는 프로젝트에 없다.** `public_data_tierA/network/`를 포함해 프로젝트 전체를 검색했고 Reactome GMT와 fgsea 패키지의 mouse reactome GMT뿐이다. 새로 받지 않았다. 두 번째 senescence 패널로 확인하려면 사용자의 다운로드 승인이 필요하다(Saul 2022 *Nat. Commun.* 보충 자료 또는 MSigDB `SAUL_SEN_MAYO`, 유전자 125개 수준의 작은 목록).

# 4. 산출물

| 파일 | 내용 |
|---|---|
| `scripts/revision/hallmark_extensions.R` | 이 문서의 모든 계산 |
| `public_data_tierA/derived/hallmark_revised/division_rate_cellcycle_removed.tsv` | 24개 프로그램, cell-cycle 유전자 제거 전후 ρ |
| `public_data_tierA/derived/hallmark_revised/donor_age_primary_cohort.tsv` | 새 코호트 연령 연관, 발표판 수치 병기, 증식 보정 후 포함 |
| `public_data_tierA/derived/hallmark_revised/senescence_dissection_cellcycle_union.tsv` | senescence 150개 유전자별 ρ와 cell-cycle 소속 |

발표판 산출물(`public_data_tierA/derived/hallmark_audit/`)과 그림 스크립트는 건드리지 않았다.
