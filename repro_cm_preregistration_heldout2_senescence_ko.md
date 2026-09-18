# 사전등록 2차: held-out senescence (GSE306957)

문서 고정일: 2026-09-11
상태: **GSE306957을 아직 다운로드하지 않았고 어떤 결과도 보지 않은 시점에 작성됨**
선행 문서: `repro_cm_preregistration_phase2_ko.md` (SHA256 `41e4cf4b548f74904db7ee8d6fefcf9e18cd801972bcdc1cec95de59a8788bd7`)
결과 문서: `repro_cm_heldout_results_ko.md`

## 0. 왜 이 검정을 추가하는가

1차 held-out(6개 primary test)은 `repro_cm_heldout_results_ko.md` §5에 기록한 대로 **설계 결손**이 있었다. 핵심 주장의 절반인 anti-senescence를 검정하는 데이터셋이 6개 중 하나도 없었다. 기존 senescence 참조 4개(GSE109700·GSE93535·GSE179848·GSE191055)는 discovery meta-signature에 이미 사용됐으므로 held-out이 아니다.

GSE306957은 이 결손을 메운다. `CLAUDE.md` §10에서 primary contrast가 이미 고정돼 있다.

## 1. 다중검정 가족(family) 처리 — 먼저 못박는다

1차 사전등록의 6개 primary test는 이미 실행되어 BH 보정이 끝난 **닫힌 가족**이다. 여기에 새 검정을 끼워 넣어 기존 FDR을 다시 계산하지 않는다. 그렇게 하면 이미 확정된 판정을 사후에 바꾸는 것이 된다.

따라서 **HO5는 독립된 두 번째 가족**으로 취급하고, 이 가족 내부(HO5a, HO5b 2개)에서만 BH 보정한다. 논문에서는 두 가족을 분리해 보고하고, 가족을 나눈 사실 자체를 명시한다.

## 2. 데이터와 대비

| 항목 | 내용 |
|---|---|
| dataset | GSE306957 |
| 세포 | MRC5 human **lung** fibroblast (기존 참조는 대부분 dermal → **조직 간 일반화 검정**이기도 하다) |
| primary contrast | **senescent control vs proliferating control** (`CLAUDE.md` §10에서 사전 고정) |
| sensitivity | RIG-I / MDA5 / MAVS knockout arm을 별도 검정. BH 가족에 넣지 않는다 |
| 파일 | `GSE306957_RAW.tar` 약 141 MB (processed) |

## 3. 가설과 성공 기준

검정 통계량은 1차와 동일하다: `ρ_Spearman(local anchor score, held-out contrast logFC)`, 발현량 decile 계층화 순열 10,000회, 양측. 로컬 anchor 정의는 변경하지 않는다.

### HO5a (primary)

- **예측:** ρ < 0. Repro-CM anchor는 senescence 방향에 반대한다.
- **성공:** ρ < 0 이고 BH-FDR < 0.05
- **반증:** ρ > 0 이고 BH-FDR < 0.05 → **핵심 주장의 anti-senescence 절반을 철회한다**
- **미지지:** 그 외

### HO5b (co-primary) — 증식 보정판

senescent vs proliferating 대비는 **정의상 증식 대비**이다. 따라서 보정 없는 ρ < 0만으로는 "senescence에 반대"와 "증식에 찬성"을 구분할 수 없다. 이것이 이 데이터셋에서 진짜로 위험을 감수하는 검정이다.

- **방법:** local anchor score와 held-out logFC를 각각 `direction_probe/P3_proliferation_loading_GSE179848.tsv`의 per-gene proliferation loading에 대해 선형 잔차화한 뒤 상관을 구한다. null은 잔차화한 anchor를 발현량 decile 내에서 섞어 만든다.
- **예측:** ρ_adjusted < 0. 즉 증식을 제거해도 anti-senescence 방향이 남는다.
- **성공:** ρ_adjusted < 0 이고 BH-FDR < 0.05
- **미지지:** ρ_adjusted ≥ 0 또는 비유의 → anti-senescence는 증식 상태와 분리되지 않는다. 이 경우 `direction_probe/P4`(증식 무관 유전자에서 D = 0.546)와 충돌하므로 **P4의 해석을 재검토해야 한다는 사실을 본문에 명시한다.**

### 종합 판정

- **CONFIRMED:** HO5a 성공 **그리고** HO5b 성공
- **PARTIAL:** HO5a만 성공
- **NOT SUPPORTED:** 둘 다 미지지
- **REFUTED:** HO5a 반증
- **NOT EVALUABLE:** §4 제외 규칙에 해당

## 4. 제외 규칙 (사전 고정)

1차 사전등록 §2.4와 동일하되, 1차 변경 이력에서 명확화한 내용을 처음부터 반영한다.

1. **발현 유전자(filterByExpr 통과) 기준** gene-symbol 매핑률 < 70%이면 primary에서 제외한다. 전체 행 기준 매핑률도 함께 보고한다.
2. metadata로 사전 정의한 contrast를 구성할 수 없으면 NOT EVALUABLE로 처리하고 실패로 계산하지 않는다.
3. 유전자 집합은 로컬 anchor의 comparator-consistent 집합(2,237개)과의 교집합으로 하고, 교집합 < 1,000개면 NOT EVALUABLE.
4. processed matrix가 이미 log 변환돼 있으면 재변환하지 않는다. 판정 근거를 기록한다.
5. 동일 연구 내 복수 arm은 독립 연구로 중복 계수하지 않는다.

## 5. 해석 경계 (결과와 무관하게 미리 고정)

- MRC5는 **폐** 섬유아세포다. 기존 로컬·참조 자료는 대부분 **진피** 섬유아세포다. 따라서 HO5가 성공하면 조직을 넘는 일반화 근거가 되지만, 실패하면 **주장 반증과 조직 특이성을 구분할 수 없다.** 실패 시 이 모호성을 그대로 보고하고, 조직 차이를 사후 변명으로 사용하지 않는다.
- 이 데이터셋의 senescence는 특정 모델(RIG-I/MDA5/MAVS·mtRNA-SASP 맥락)에서 유도된 것이다. 단일 모델 결과를 senescence 일반으로 확대하지 않는다.
- knockout arm 결과는 기전 주장이 아니라 sensitivity로만 사용한다.

## 6. 변경 이력

- 2026-09-11: 최초 고정. 이후 변경은 아래에만 append한다.

- 2026-09-11 (metadata 확인 후, **결과를 보기 전** 기록):
  1. **control arm 특정.** GSE306957의 MRC5 실험은 두 하위 실험으로 구성된다 — CRISPR 세트(pLenti 대조, IFIH1 KO, DDX58 KO)와 siRNA 세트(siScramble 대조, siMAVS). 사전등록 문서는 "senescent control vs proliferating control"만 지정했다. **두 비표적 대조(pLenti, siScramble)를 합쳐 primary로 하고, 하위 실험(CRISPR vs siRNA)을 batch covariate로 넣는다.** 결과 n = proliferating 6, senescent 6. 하위 실험 각각은 sensitivity로 별도 보고한다(pLenti 3 vs 2, siScramble 3 vs 4).
  2. **마우스 간 시료 제외.** GSE306957에는 MRC5 27개 외에 Bak/Bax 마우스 간 7개(ENSMUSG)가 포함돼 있다. 종·조직이 달라 primary·sensitivity 어디에도 사용하지 않는다.
  3. senescence 유도는 **20 Gy X선 조사**다(복제노화가 아님). §5의 "단일 모델을 senescence 일반으로 확대하지 않는다"가 그대로 적용된다.
