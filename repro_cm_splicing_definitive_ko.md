# splicing 질문: 어떤 데이터로 무엇에 답할 수 있는가 (확정판)

작성일: 2026-09-11
선행 문서: `repro_cm_phase2_splicing_results_ko.md` (§1–5의 결론을 **이 문서가 대체**한다)
스크립트: `run_phase2_H1_GSE109700.R`, `run_phase2_H2_quiescence.R`, `run_phase2_H1b_splicing_age.R`

> **이전 결론("splicing 결과는 이 논문의 기둥이 될 수 없다")은 데이터셋 선택과 통계 선택의 실패였다.** 적절한 자료로 다시 하면 senescence는 splicing을 분명히 바꾸고, 그 변화의 일부는 증식정지와 분리된다.

---

## 1. 질문을 네 개로 쪼갠다

splicing은 하나의 질문이 아니다. 네 개이고, 요구하는 데이터가 각각 다르다.

| | 질문 | 필요한 설계 |
|---|---|---|
| **Q1** | senescence가 splicing 결과를 바꾸는가 | senescent vs proliferating, 충분한 junction 깊이, n≥3 |
| **Q2** | 그 변화가 **증식정지와 분리**되는가 | 동일 PD의 **quiescent 대조** |
| **Q3** | 공여자 연령이 splicing 결과를 바꾸는가 | 연령 코호트, n 충분, 깊이 충분 |
| **Q4** | **Repro-CM**이 splicing 결과를 바꾸는가 | 반복이 있는 Repro-CM 처리 |

---

## 2. Q1 — 답: **예. 확실하다.**

**자료: GSE109700 (SRP131506, recount3, 이미 디스크에 있음).** LF1 폐 섬유아세포, **paired-end ribodepleted total RNA**, proliferating/early/deep senescent 각 n=3.

깊이가 결정적으로 다르다.

| | layout | 총 read | junction read 비율 | junction read |
|---|---|---:|---:|---:|
| **GSE109700** | paired, **ribodepleted** | 70–89 M | **35.3%** | **20–37 M** |
| GSE93535 | single, polyA | 30–50 M | 10.7% | 2.8–5.8 M |

전 시료를 19.46 M으로 깊이 맞춘 뒤 limma moderated t:

| 대비 | testable junction | FDR<0.05 | \|ΔPSI\|≥0.10 |
|---|---:|---:|---:|
| **deep senescence vs proliferating** | 97,563 | **651** | 168 |
| **early senescence vs proliferating** | 91,868 | **772** | 231 |
| deep vs early (둘 다 senescent) | 86,137 | **0** | 0 |

세 겹의 대조를 모두 통과한다.

1. **깊이 교락 제거.** deep 시료 3개가 원래 깊이 하위 3개와 정확히 일치했으나, 전 시료 downsampling 후에도 651개가 남는다. early arm은 애초에 proliferating과 깊이가 맞고(29.4 vs 30.2 M) 오히려 더 많은 event를 낸다.
2. **내부 음성 대조.** 두 senescent 상태끼리는 **0개**다.
3. **label permutation.** early-vs-proliferating 깊이보정판에서 참 분할 772개, 나머지 **9개 분할 전부 0개**.

> **senescence는 대체 splice-site 사용을 분명하게 바꾼다.**

---

## 3. Q2 — 답: **예. 다만 작다.**

**자료: GSE93535 (SRP096629, recount3, 이미 디스크에 있음).** HDF161, **전 시료 population doubling 15**. quiescent(Q) 3, Q+1201 3, SIPS 3, SIPS+1201 3.

n=3 vs 3으로는 검정력이 없다. 그러나 화합물 arm까지 넣어 **2×2 설계(state × compound)로 state 주효과**를 보면 n=6 vs 6이 된다.

| 필터 | testable junction | FDR<0.05 | \|ΔPSI\|≥0.10 | minP |
|---|---:|---:|---:|---:|
| mincov 20, 전 시료 | 3,441 | 23 | 12 | 6.0e−07 |
| mincov 40, 전 시료 | 1,729 | 14 | 8 | 5.8e−07 |
| 얕은 run 2개 제외, mincov 20 | 14,797 | **21** | 7 | — |
| 얕은 run 2개 제외, mincov 40 | 8,130 | **25** | 8 | — |

**label permutation (compound 층 내에서 state만 치환, 298회):**

| 필터 | 관측 | 순열 median | 순열 q95 | 순열 max | p |
|---|---:|---:|---:|---:|---:|
| mincov 20 | **21** | 1 | 5 | 11 | **0.0033** |
| mincov 40 | **25** | 0 | 4 | 6 | **0.0033** |

관측값이 **298개 순열 전부를 초과**한다(p는 순열 바닥값).

### 크기를 정직하게 비교한다

| 기준 대비 | 유의 event / testable | 비율 |
|---|---|---:|
| senescent vs **proliferating** (GSE109700) | 651 / 97,563 | **0.67%** |
| senescent vs **quiescent, 동일 PD** (GSE93535) | 21 / 14,797 | **0.14%** |

> **senescence의 splicing 변화 대부분은 증식정지와 공유된다. 그러나 약 5배 작은 senescence 특이 성분이 남고, 그것은 순열 대조를 통과한다.**

---

## 4. Q3 — 답: **event 수준은 아니오, 전역 noise는 약한 예. 단 깊이 한계가 있다.**

**자료: GSE113957 (SRP144355, recount3), 건강 공여자 104명, 22–96세.** junction read 중앙값 **4.03 M** — GSE109700의 1/5이다.

| 분석 | 보정 전 FDR<0.05 | 증식 보정 후 |
|---|---:|---:|
| 5'-anchored junction PSI (30,878) | 42 | **0** |
| 3'-anchored junction PSI (31,131) | 32 | **0** |

전역 지표(unannotated junction read fraction, 탐색적)는 살아남는다: 전체 보정(증식+platform+read length+mapping%) 후 β=+6.06e−05/10년, **p=0.016**; 단일 platform(n=91)에서도 p=0.027. 효과크기는 10년당 약 0.75% 상대 증가로 매우 작다.

**해석의 경계:** 증식 점수와 연령의 상관이 r=−0.47이고, 배양 섬유아세포에서 증식은 교란이 아니라 **매개**일 수 있다. Q1·Q2에서 senescence가 splicing을 바꾼다는 것이 확인된 지금, 연령→증식능 저하→splicing 변화라는 매개 경로가 오히려 자연스럽다. **증식 보정 후 0개라는 사실을 "연령 효과가 없다"로 읽으면 안 된다.**

---

## 5. Q4 — 답: **이 자료로는 불가능하다.**

Repro-CM은 조건당 **n=1**이다. 세 각도 모두 막힌다.

- junction ΔPSI: 두 비교자 일치 후보의 경험적 FDR 0.71–0.74 (`direction_probe` P6)
- intron retention: Repro-CM이 세 24h 조건 중 최고이고 발현과 ρ=−0.25로 교락 (P7)
- 반복이 없으므로 어떤 통계도 처치효과로 해석할 수 없다

**필요한 것:** Repro-CM·HDF-CM·iPSC-CM 각 n≥3의 새 RNA-seq. 가능하면 paired-end ribodepleted(= junction 수율 3배). 이것 없이는 "Repro-CM이 splicing을 바꾼다"를 쓸 수 없다.

---

## 6. 더 받아도 소용없는 것 / 받으면 도움 되는 것

### 소용없음

- **GSE93535의 FASTQ (17.36 GiB).** recount3가 이미 전체 read를 처리했다. SIPS_2는 애초에 4.5 M read짜리 얕은 run이므로 재정렬해도 junction이 늘지 않는다.

### 도움 될 수 있음 (필요해지면)

| 목적 | 자료 | 비용 | 한계 |
|---|---|---:|---|
| Q2 강화 | GSE179848 contact inhibition(21) vs late-passage control, FASTQ subset | 약 190 GiB (40 run) | polyA single-read라 junction 수율 ~10%. n이 크지만 시료당 깊이는 GSE93535 수준 |
| Q1·Q3 확장 | 다른 ribodepleted paired-end senescence/aging 연구 | 미탐색 | 별도 검색·승인 필요 |

---

## 7. 이 결과가 논문에서 갖는 위치

**바뀐 것:** splicing은 "검정했고 지지되지 않음"이 아니라 **"참조 축에 대해 확립됨"**이다.

- Q1·Q2의 결과는 **senescence 축의 성질**에 대한 독립적 검증이다. D 축이 포착하는 senescence 상태가 전사량뿐 아니라 **splicing 결과에서도** proliferating·quiescent와 구분된다는 뜻이다.
- 특히 Q2는 compendium의 contact-inhibition 결과(정지상태의 D=+0.168 vs senescence −0.5~−0.7)와 **같은 결론에 독립적으로 도달한다.** 하나는 유전자 발현, 하나는 splicing 결과다.

**바뀌지 않은 것:** Repro-CM에 대해서는 여전히 아무 말도 할 수 없다(Q4). 따라서 논문에서 splicing은

- ✅ 참조 축(senescence)의 성질로 보고한다 — Results
- ❌ Repro-CM의 작용기전으로 쓰지 않는다 — 불가
- ✅ 후속 실험 설계의 근거로 쓴다 — Discussion (n≥3, ribodepleted paired-end)
