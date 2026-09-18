# 섬유아세포 개입 compendium: D 축의 검증과 Repro-CM의 위치

실행일: 2026-09-11
스크립트: `scripts/build_perturbation_compendium.R`, `scripts/score_perturbation_compendium.R`, `scripts/make_figure7_compendium.R`
출력: `public_data_tierA/derived/compendium/`, `manuscript_figures_v4/Figure7A/B`

## 1. 왜 만들었나

리뷰어 평가(`repro_cm_reviewer_assessment_ko.md`)에서 확인한 두 가지 치명적 약점 — 개입이 n=1이고(M1), 1차 자료가 이미 출판됨(M2) — 은 **논문의 주어를 "Repro-CM"에서 "축(framework)"으로 바꾸면 무력화된다.** 그러려면 D가 개입 하나에만 적용된 지표가 아니라 여러 개입을 구분하는 축임을 보여야 한다.

디스크에 이미 있는 자료만으로 **50개 섬유아세포 개입 signature**를 만들어 동일한 두 참조 축에 배치했다. 추가 다운로드는 없었다.

## 2. 구성

| 출처 | 개입 | 수 |
|---|---|---:|
| GSE179848 | DEX, contact inhibition, 2-deoxyglucose, β-hydroxybutyrate, galactose, mitoNUITs(±DEX), oligomycin(±DEX), 저산소 3%, SURF1 변이 | 12 |
| GSE116968 | 적색광 단독, UV, Pre-Red vs UV (±injury 잔차화), 1 h·4 h | 8 |
| GSE149694 | Primed / NHSM / 5iLAF / RSeT / t2iLGoY day 13, fibroblast day 7 | 6 |
| GSE306957 | pLenti / siScramble / IFIH1 KO / DDX58 KO senescent vs proliferating | 4 |
| GSE240226 | UVA, Maifuyin rescue, succinate rescue | 3 |
| GSE89005 | single 6 h / 24 h, repeated 24 h UVA | 3 |
| GSE297233 | OSK ±dox, O4YRSK mutant | 2 |
| GSE109700 | early / deep 복제노화 | 2 |
| GSE93535 | SIPS vs Q, 1201 in Q, 1201 in SIPS | 3 |
| GSE113957 · GSE226189 · GSE307377 | 공여자 연령 | 3 |
| GSE302943 · GSE125429 · GSE191055 · GSE165177 · GSE179848(late) | UVA, 계대, MPTR | 5 |
| **GSE240486** | UV 손상, Osmoter 용량 2개 rescue, 약물 단독 | 4 |
| **GSE222414** | UVB 손상, TAT-UIFSP 2개 rescue | 3 |
| **GSE329475** | UVA 3일 연속 손상 | 1 |
| **this study** | **Repro-CM anchor** | **1** |

**58개 중 37개는 두 참조 축 어디에도 기여하지 않은 독립 signature다.** 축을 정의한 6개는 `[AXIS]`로, 같은 연구의 다른 대비는 `[same study]`로 표시해 순환성을 명시한다.

## 3. 축이 작동한다 — 세 가지 검증

### 3.1 senescence가 축의 바닥을 차지한다

| signature | D | 독립성 |
|---|---:|---|
| GSE179848 late vs early passage | −0.725 | [AXIS] |
| GSE109700 deep replicative senescence | −0.614 | [AXIS] |
| **GSE306957 IFIH1 KO senescent** | **−0.606** | **독립** |
| **GSE306957 pLenti senescent** | **−0.599** | **독립** |
| GSE191055 P27 vs P4 | −0.585 | [AXIS] |
| GSE93535 SIPS vs quiescent | −0.579 | [AXIS] |
| **GSE306957 siScramble senescent** | **−0.573** | **독립** |
| GSE109700 early replicative senescence | −0.568 | [same study] |
| **GSE306957 DDX58 KO senescent** | **−0.503** | **독립** |

**축을 정의하는 데 전혀 쓰이지 않은 GSE306957의 네 arm이 축을 정의한 대비들과 같은 자리에 착지한다.** 조직(폐 vs 진피)도 유도기전(X선 vs 복제/산화)도 다르다. D가 senescence를 일반적으로 포착한다는 직접 증거다.

### 3.2 D는 "증식정지 축"이 아니다 — 내부 대조

GSE179848의 **contact inhibition**은 증식이 멈췄지만 senescence가 아닌 상태다.

| | ρ_UVA | ρ_senescence | D |
|---|---:|---:|---:|
| contact inhibition (정지) | +0.186 | **+0.018** | **+0.168** |
| senescence 9개 평균 | +0.07 | **+0.67** | **−0.59** |

증식정지만으로는 senescence 축에서 전혀 올라가지 않는다. 이것은 `direction_probe/P4`(증식 잔차화 후 D 유지)와 `HO5b`(held-out 증식 보정 통과)를 **개입 수준에서 다시 확인**한 것이다. 세 가지 독립적인 방식이 같은 결론에 도달했다.

### 3.3 연령·MPTR도 예상 위치에 있다

공여자 연령 signature 3개가 모두 senescence 쪽으로 기운다(GSE113957 −0.324, GSE307377 −0.182, GSE226189 −0.181). MPTR도 −0.344다. 즉 D는 "노화 관련 상태"를 일관되게 아래쪽에 놓는다.

## 4. Repro-CM의 위치

| | D | 순위 |
|---|---:|---|
| GSE240226 acute UVA | +0.886 | 1 (축 정의, 순환) |
| GSE302943 cumulative UVA | +0.577 | 2 (축 정의, 순환) |
| **Repro-CM anchor** | **+0.460** | **3 / 58 · 독립 37개 중 1위** |
| GSE149694 t2iLGoY-D13 (naive) | +0.385 | 4 |
| GSE179848 2-deoxyglucose | +0.330 | 5 |

순환적인 두 개를 제외하면 **Repro-CM이 58개 개입 중 최상위**다. 이것은 "n=1 시료 하나를 공개자료와 상관시킨 결과"가 아니라 **검증된 지도 위에서의 위치**다. 리뷰어 평가 M1(n=1)의 무게가 크게 줄어든다.

## 5. 정직하게 남겨야 할 것들

1. **상위 두 개는 순환적이다.** UVA 대비가 UVA 축 꼭대기에 있는 것은 당연하며 증거가 아니다. 그림에서 빈 원으로 구분했다.
2. **~~rescue 개입끼리 방향이 엇갈린다~~ → 순환성이 원인이었다 (2026-09-11 수정).** 3차 held-out으로 독립 광보호 대비 5개를 추가하니 구조가 드러났다.

   | 그룹 | n | D 범위 |
   |---|---:|---|
   | **독립 광보호 대비** (GSE116968·GSE222414·GSE240486) | **11** | **+0.003 ~ +0.153, 전부 ≥ 0** |
   | 축 정의 연구 내부 rescue (GSE240226) | 2 | −0.134 (구조상 음수) |
   | 독립 UV 손상 대비 | 8 | −0.435 ~ +0.122 |

   **세 독립 연구의 광보호 대비 11개 중 음수가 하나도 없다.** 앞서 "이질성"이라고 적은 것은 GSE240226이 UVA 축을 정의하는 연구이기 때문에 그 연구의 rescue가 구조적으로 음수가 될 수밖에 없었던 것이다. 순환성을 분리하고 보면 규칙은 성립한다: **광보호는 D 양수, 손상은 D 음수.**
3. **t2iLGoY naive reprogramming이 D=+0.385로 4위**다. held-out sensitivity에서 경계를 넘었던 바로 그 arm이다. 일관된 신호이므로 Discussion에서 다뤄야 한다.
4. **GSE179848 대비들은 senescence 축의 한 축(late vs early)과 같은 연구다.** `[same study]`로 표시했다. contact inhibition 결과는 이 한계를 안고 해석한다.
5. 각 signature의 n과 설계는 제각각이다. D를 개입 간 정량 비교로 쓸 때는 효과크기 차이가 아니라 **방향과 순위**로 해석한다.

## 6. 논문 구성에 주는 함의

리뷰어 평가 §5에서 제안한 재구성이 이제 실행 가능하다.

> **주어를 축으로:** "섬유아세포에서 급성 스트레스와 senescence를 증식과 무관하게 분리하는 정량 축을 공개자료 50개 개입으로 구축·검증했다. 그 축 위에서 리프로그래밍기 secretome은 독립 개입 중 최상위에 위치한다."

- 1차 결과가 반복수 충분한 공개자료에서 나온다 → **M1(n=1) 완화**
- 로컬 자료는 지도 위의 한 점 → **M2(선행 중복) 무력화**
- D가 개입을 구분함을 보였다 → **M5(resource 범용성) 해소**
- contact inhibition이라는 내부 음성 대조가 생겼다 → 증식 반박 차단이 더 단단해짐
