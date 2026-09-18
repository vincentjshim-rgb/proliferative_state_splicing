# 사전등록 3차: 독립 광보호 연구와 UVA 축 재검정

문서 고정일: 2026-09-11
상태: **GSE240486 · GSE222414 · GSE329475를 아직 다운로드하지 않았고 어떤 결과도 보지 않은 시점에 작성됨**
선행: `repro_cm_preregistration_phase2_ko.md` (`41e4cf4b`), `repro_cm_preregistration_heldout2_senescence_ko.md` (`6e5a7703`)
근거 문서: `repro_cm_axis_characterisation_ko.md`

## 0. 왜 이 검정을 추가하는가

`repro_cm_axis_characterisation_ko.md`에서 두 가지가 드러났다.

1. **광보호 정렬이 사실상 한 연구에 의존한다.** 광보호 대비 6개 중 4개가 GSE116968 하나에서 나왔다. 독립 연구 수는 실질적으로 1개다.
2. **UVA 축이 독립 UV 손상으로 전이되지 않는다.** 독립 손상 6개의 median이 −0.017이다. 이로부터 "축이 포착하는 것은 손상이 아니라 초기 광적응 반응"이라는 **사후 해석**을 세웠다. 사후 해석은 새 데이터로 검정해야 한다.

이번 사전등록은 그 둘을 각각 검정한다. **특히 HO7은 내가 세운 사후 해석을 반증할 수 있는 검정이다.**

## 1. 다중검정 가족

1차(6개)와 2차(2개) 가족은 닫혔다. 재보정하지 않는다. **HO6·HO7은 세 번째 독립 가족**이며 이 가족 내부 5개에서만 BH 보정한다. 논문에서 세 가족을 분리해 보고하고, 분리 사실을 명시한다.

## 2. 데이터와 대비

| dataset | 세포 | 설계 | n |
|---|---|---|---|
| **GSE240486** | normal human dermal fibroblast | Control / 0.02% Osmoter / 0.06% Osmoter / UV / UV+0.02% / UV+0.06% | 군당 3 |
| **GSE222414** | WS1 human skin fibroblast | 무처리 / UVB / UVB+TAT-UIFSP5 / UVB+TAT-UIFSP6 | 군당 4 |
| **GSE329475** | human dermal fibroblast | NC / UVA 3일 연속 | 군당 3 |

세 연구 모두 섬유아세포다. 각질형성세포·in vivo 피부 연구는 세포 조성 교란 때문에 의도적으로 제외했다(GSE320415, GSE310833, GSE22083).

## 3. Primary tests — 정확히 5개

검정 통계량은 기존과 동일하다: `ρ_Spearman(local anchor score, contrast logFC)`, 발현량 decile 계층화 순열 10,000회, 양측. 로컬 anchor 정의는 변경하지 않는다.

### HO6 — 광보호 정렬의 독립 재현 (4개)

rescue signature에서 **injury 성분을 회귀 제거한 residual**을 사용한다. GSE240226·GSE116968과 동일한 절차다.

| ID | contrast | injury 축 |
|---|---|---|
| HO6a | GSE240486 · UV+0.02% Osmoter vs UV | UV vs Control |
| HO6b | GSE240486 · UV+0.06% Osmoter vs UV | UV vs Control |
| HO6c | GSE222414 · UVB+TAT-UIFSP5 vs UVB | UVB vs 무처리 |
| HO6d | GSE222414 · UVB+TAT-UIFSP6 vs UVB | UVB vs 무처리 |

- **예측:** 네 건 모두 ρ > 0. 기존 값(Maifuyin +0.127, succinate +0.094, 적색광 +0.122/+0.098)과 같은 범위를 기대한다.
- **성공:** 4개 중 **최소 2개**에서 ρ > 0 이고 BH-FDR < 0.05, **그리고 두 연구 각각에서 최소 1개씩**.
- **반증:** 3개 이상에서 ρ < 0 이고 BH-FDR < 0.05 → repair-aligned component 주장을 철회한다.
- **미지지:** 그 외.

### HO7 — UVA 축 일반화 실패의 재검정 (1개)

| ID | contrast | 성격 |
|---|---|---|
| HO7 | GSE329475 · UVA 3일 연속 vs NC | 독립 **UVA 손상** 대비 |

`axis_characterisation`의 사후 해석은 "축은 초기 광적응 반응을 포착하며, 누적 UV 손상은 로드되지 않는다"였다. 3일 연속 조사는 누적 노출이다.

- **예측:** `|ρ_UVA| < 0.15`. 축을 정의한 두 연구(+0.798, +0.816)와 전혀 다른 값이어야 한다.
- **성공(해석 지지):** `|ρ_UVA| < 0.15`
- **반증(해석 철회):** `ρ_UVA > 0.40` → 축이 누적 UVA 손상도 포착한다는 뜻이므로 "초기 광적응 축" 해석을 철회하고 원래 "UVA 축" 서술로 되돌린다.
- **중간:** 0.15 ≤ ρ_UVA ≤ 0.40 → 어느 쪽도 아님. 그대로 보고하고 해석을 약화한다.

**이 검정은 내 사후 해석에 불리하게 나올 수 있다. 그 가능성을 열어 두는 것이 이 사전등록의 목적이다.**

## 4. Sensitivity (BH 가족에 넣지 않음)

1. **약물 단독 특이성 대조.** GSE240486의 `0.02% Osmoter vs Control`, `0.06% Osmoter vs Control` (UV 없음). 광보호 정렬이 UV 맥락에 특이적이라면 약물 단독 대비는 rescue 대비보다 약해야 한다. 이 예측도 미리 적는다.
2. **injury 대비 자체.** GSE240486 `UV vs Control`, GSE222414 `UVB vs 무처리` — HO7과 같은 논리로 UVA 축 로딩을 본다.
3. **잔차화하지 않은 rescue 원값**도 함께 보고한다.
4. 각 신규 대비의 **compendium D 점수**를 계산해 기존 50개 지도에 올린다.

## 5. 제외 규칙 (사전 고정)

1. 발현 유전자(filterByExpr 통과) 기준 gene-symbol 매핑률 < 70%이면 primary에서 제외한다. 전체 행 기준도 함께 보고한다.
2. metadata로 사전 정의한 contrast를 구성할 수 없으면 NOT EVALUABLE. 실패로 계산하지 않는다.
3. 로컬 anchor의 comparator-consistent 집합과의 교집합 < 1,000개면 NOT EVALUABLE.
4. processed matrix가 이미 log 변환돼 있으면 재변환하지 않는다. 판정 근거를 기록한다.
5. 동일 연구 내 복수 arm은 독립 연구로 중복 계수하지 않는다. HO6a/b는 한 연구, HO6c/d는 다른 한 연구다.

## 6. 해석 경계 (결과와 무관하게 미리 고정)

- WS1은 **태아 유래 피부 섬유아세포주**이고 NHDF는 **초대 배양 진피 섬유아세포**다. 둘 다 섬유아세포지만 동일하지 않다. HO6c/d가 HO6a/b와 다르게 나오면 세포 출처 차이를 대안 설명으로 명시한다.
- Osmoter(사해수 미네랄)와 TAT-UIFSP(펩타이드 융합 단백질)는 기전이 전혀 다르다. **둘 다 정렬되면 기전 비특이적 공통 성분**이라는 뜻이고, 이는 기전 주장이 아니라 관찰이다.
- GSE329475의 UVA는 3일 연속 조사다. 시간·누적 노출이 축 로딩에 미치는 영향을 단일 연구로 분리할 수 없다.
- 이 검정들은 **transcriptional alignment**이며 phenotype 보호의 증거가 아니다.

## 7. 변경 이력

- 2026-09-11: 최초 고정. 이후 변경은 아래에만 append한다.
