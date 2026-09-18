# Repro-CM 후속 연구 사전등록 (Phase 2 + held-out generalization)

문서 고정일: 2026-09-11
상태: **아래 데이터의 결과를 아직 열지 않은 시점에 작성됨**
근거 문서: `repro_cm_research_direction_v4_ko.md`, `CLAUDE.md` §10, §15-5
재현 스크립트: `scripts/run_direction_probe.R`

이 문서는 결과를 보기 전에 가설·대비·성공기준·제외규칙·다중검정 계획을 고정한다. 이후 어떤 결과가 나오든 이 문서의 기준을 바꾸지 않는다. 기준을 바꿔야 할 사유가 생기면 변경 이력을 이 파일 하단에 append하고, 변경 전후를 모두 보고한다.

---

## 0. 고정된 입력 (변경 금지)

로컬 anchor 정의는 `CLAUDE.md` §5를 그대로 사용한다. 새 데이터에 맞춰 재조정하지 않는다.

- STAR reverse-stranded gene counts → TMM → 2개 이상 시료에서 CPM ≥ 1 → pseudocount 0.5
- `dH = log2CPM(Repro) − log2CPM(HDF)`, `dI = log2CPM(Repro) − log2CPM(iPSC)`
- `S = sign(dH) · min(|dH|, |dI|)` (부호 불일치 시 0)
- 고정 산출물: `public_data_tierA/derived/decoupling_validation/primary_gene_axis_matrix.tsv`

증식 loading은 `public_data_tierA/derived/direction_probe/P3_proliferation_loading_GSE179848.tsv`를 그대로 사용한다(GSE179848 345 시료 기반, 이미 고정됨).

모든 공개 데이터와의 연관 검정은 **발현량 계층화(decile) 순열 null**을 사용한다. 순열 횟수 10,000회. 유전자는 생물학적 반복이 아니므로, 산출되는 p/FDR는 **signature association에 대한 유의성**이지 처치효과의 replicate p-value가 아니다. 이 문구를 모든 figure legend에 넣는다.

---

## 1. Phase 2 — splicing을 발현이 아니라 결과로 검정

### 배경 (이미 확정된 음성 결과)

`direction_probe/P5`에서 splicing factor **발현** 모듈은 증식 모듈과 r = 0.923이고, 연령 기울기는 증식 보정 후 p = 0.71로 소실됐다. 따라서 splicing factor 발현을 노화 축으로 사용하는 경로는 **닫혔다**. `direction_probe/P6`에서 로컬 n=1 junction ΔPSI 후보는 경험적 FDR 0.71–0.74로 증거가 될 수 없다.

Phase 2는 남은 단 하나의 열린 질문을 검정한다.

> PSI와 intron retention은 유전자 내부 비율이므로 발현량처럼 증식에 자동 비례하지 않는다. **splicing 결과는 증식과 분리 가능한 senescence 축인가?**

### 가설

**H1 (primary).** 동일 population doubling에서 senescent 섬유아세포는 quiescent 섬유아세포와 splicing 결과가 다르다.
- 자료: GSE93535 (HDF161, 전부 PD15) — `Q` (n=3) vs `SIPS` (n=3)
- 지표 (사전 정의, 우선순위 순):
  1. **IR burden**: 분석 가능한 intron 전체에서 `intronic coverage / flanking exonic coverage`의 시료별 중앙값
  2. **전역 |ΔPSI| 분포**: 검정 가능한 alternative junction의 |ΔPSI| 중앙값
  3. **차등 splicing event 수**: FDR < 0.05 & |ΔPSI| ≥ 0.10
- **성공:** 지표 1 또는 2에서 Q vs SIPS가 양측 p < 0.05 (시료 수준 검정, n=3 vs 3)이고, 같은 방향이 GSE191055 P4 vs P27에서도 재현
- **실패:** Q vs SIPS에 차이가 없다 → splicing 결과 역시 증식정지와 분리되지 않는다. **이 음성 결과도 그대로 보고한다.**
- **평가 불가:** 어느 한 arm의 usable run이 2개 미만

**H2.** senescence splicing 축은 UVA 급성 손상 splicing 축과 구별된다.
- 자료: GSE240226 UVA vs control (paired-end)
- 검정: H1에서 정의한 event 집합에서 `ρ(ΔPSI_senescence, ΔPSI_UVA)`
- **성공:** ρ가 0과 유의하게 다르지 않거나 음수 → 두 축이 분리됨 (gene-level decoupling의 splicing 수준 대응)
- **실패:** ρ > 0.3 이고 유의 → splicing 수준에서는 두 축이 같다

**H3 (application, 검정 아님).** 로컬 Repro-CM 시료를 H1 축 위에 배치한다.
- **로컬은 n=1이므로 이것은 통계 검정이 아니라 참조 분포 위의 한 점 배치(placement)다.** 이 문구를 figure와 본문에 그대로 쓴다.
- 참조 분포: H1/H2 공개 시료의 IR burden 및 splicing 축 점수
- 보고: 로컬 4개 시료(HDF-CM, Repro-CM, iPSC-CM, 15J_0hr)의 위치와 참조 분포의 사분위수
- **어떤 경우에도 p-value를 붙이지 않는다.**

### 데이터 확보 순서

1. recount3 / Snaptron에 GSE113957·GSE93535·GSE191055·GSE240226의 junction count가 있는지 먼저 확인한다. 있으면 FASTQ 대신 사용한다.
   - 특히 GSE113957(공여자 143명, 1–96세) junction이 확보되면 **n=143으로 "splicing 결과의 연령 변화가 증식과 분리되는가"를 직접 검정**할 수 있다. 이 경우 별도 가설 H1b를 추가한다: *IR burden은 공여자 연령과 연관되며 이 연관은 증식 모듈 보정 후에도 유지된다.*
2. recount3에 없는 것만 Tier B FASTQ로 받는다 (GSE93535 17.36 GiB, GSE191055 2.06 GiB, GSE240226 26.04 GiB).
3. 정렬은 로컬 BAM과 동일한 GENCODE v41 GRCh38 index, STAR 2.7.3a, `--outFilterType BySJout`를 사용한다.

---

## 2. Held-out generalization — 결과를 보기 전에 고정

### 2.1 데이터와 대비 (CLAUDE.md §10에서 이미 고정, 변경 없음)

| ID | dataset | primary contrast | sensitivity |
|---|---|---|---|
| HO1 | GSE297233 | OSK **+dox vs −dox, day 4** (96세 fibroblast) | mutant arm |
| HO2 | GSE307377 | **old vs young** primary dermal fibroblast, sex를 covariate | 없음 |
| HO3a | GSE116968 | **Pre-Red vs UV, 1 h** | — |
| HO3b | GSE116968 | **Pre-Red vs UV, 4 h** | — |
| HO4a | GSE149694 | **day 7 vs day 3** | — |
| HO4b | GSE149694 | **day 13 vs day 3** | — |

**primary test는 정확히 6개다.** BH 보정은 이 6개에 대해 수행한다. sensitivity·secondary 분석은 보정 집합에 넣지 않고 별도로 보고한다.

### 2.2 각 primary test의 검정 통계량과 사전 예측

모두 `ρ_Spearman(local anchor score, held-out contrast logFC)`, 발현량 계층화 순열 10,000회, 양측.

| ID | 사전 예측 | 성공 기준 | 실패(= 모델 반증) |
|---|---|---|---|
| HO1 | **양의 정렬이 없다** (외인성 secretome ≠ 세포 내재적 OSK) | \|ρ\| < 0.10 이거나 BH-FDR ≥ 0.05 | ρ > 0.15 이고 BH-FDR < 0.05 → "extrinsic ≠ intrinsic" 경계 주장 철회 |
| HO2 | **ρ < 0** (anchor가 노화 방향에 반대) | ρ < 0 이고 BH-FDR < 0.05 | ρ > 0 이고 BH-FDR < 0.05 → 핵심 주장 반증 |
| HO3a/b | UV injury를 잔차화한 rescue signature와 **ρ > 0** | 두 시점 중 최소 1개에서 ρ > 0, BH-FDR < 0.05 | 두 시점 모두 ρ < 0 이고 유의 → repair-aligned component 철회 |
| HO4a/b | **양의 정렬이 없다** (HO1과 동일 논리) | \|ρ\| < 0.10 이거나 BH-FDR ≥ 0.05 | ρ > 0.15 이고 BH-FDR < 0.05 |

HO2는 **증식 보정 전후를 모두** 보고한다. `P3_proliferation_loading_GSE179848.tsv`로 양 축을 잔차화한 뒤의 ρ를 함께 제시한다. 사전 예측: 보정 후에도 ρ < 0 유지.

HO3는 GSE240226 rescue 분석과 **동일한 잔차화 절차**를 쓴다: rescue contrast에서 injury contrast(UV vs control)를 회귀 제거한 residual을 사용한다.

### 2.3 전체 판정 기준 (사전 고정)

- **CONFIRMED:** HO2 성공 **그리고** HO3 성공 **그리고** HO1·HO4 어느 것도 반증 조건에 해당하지 않음
- **PARTIAL:** HO2 또는 HO3 중 하나만 성공
- **REFUTED:** HO2가 반증 조건 (ρ > 0, BH-FDR < 0.05)
- **NOT EVALUABLE:** 해당 dataset에서 사전 정의한 contrast를 metadata로 구성할 수 없음

판정이 PARTIAL 또는 REFUTED여도 **모든 결과를 그대로 보고하고 논문에 싣는다.** 성공한 subset만 고르지 않는다.

### 2.4 제외 규칙 (사전 고정)

1. gene symbol 매핑률 < 70%인 데이터셋은 primary에서 제외하고 사유를 명시한다.
2. metadata로 사전 정의한 contrast를 구성할 수 없으면 NOT EVALUABLE로 처리하고 실패로 계산하지 않는다.
3. 유전자 집합은 로컬 anchor의 comparator-consistent 집합(2,237개)과 해당 데이터셋의 교집합으로 한다. 교집합 < 1,000개면 NOT EVALUABLE.
4. processed matrix가 이미 log 변환돼 있으면 재변환하지 않는다. 변환 여부는 값 분포로 판정하고 기록한다.
5. 동일 연구 내 복수 contrast는 독립 데이터셋으로 중복 계수하지 않는다 (HO3a/b, HO4a/b는 각각 한 연구 내부의 두 시점이며, BH 보정 집합에는 포함하되 "독립 연구 수"로는 1개로 센다).

### 2.5 보고 규칙

- 6개 primary test 전부를 표로 제시한다. 성공·실패·평가불가를 모두 표시한다.
- forest plot의 가로 막대는 **matched-null 95% 구간**이며 confidence interval이 아님을 명시한다.
- pathway 분석은 secondary이며 BH 보정 집합에 넣지 않는다.
- HO1·HO4는 "null을 예측한 검정"이므로, 유의하지 않다는 결과를 **증거로 과대 해석하지 않는다.** 검정력 한계를 함께 보고한다.

---

## 3. 변경 이력

(결과를 본 뒤 기준을 바꿔야 할 경우 여기에만 append한다. 상단 본문은 수정하지 않는다.)

- 2026-09-11: 최초 고정.

- 2026-09-11 (분석 실행 직전, **결과를 보기 전** 기록):
  1. **§2.4-1 제외규칙의 측정 대상 명확화.** 규칙은 "gene symbol 매핑률 < 70%"라고만 썼으나, 이는 *ID 매핑 실패*를 걸러내려는 의도였다. 원시 ENSG 행렬은 lncRNA·pseudogene 때문에 전체 행 기준 매핑률이 정상적으로 60% 안팎이 된다. 따라서 판정은 **발현 유전자(filterByExpr 통과) 기준 매핑률**로 하고, 전체 행 기준 값도 함께 보고한다. 실측: GSE297233 전체 59.4% / 발현 91.2%, GSE149694 전체 75.4% / 발현 88.7%, GSE307377 100%. 세 데이터셋 모두 발현 기준으로 70%를 넘으므로 제외 없음. 이 변경은 "null을 예측한" HO1·HO4를 살려 두므로 우리 모델을 **더** 반증에 노출시키는 방향이다.
  2. **HO4b의 day 13 조건 특정.** GSE149694 day 13에는 배지별로 Primed-D13(n=3), NHSM-D13(n=2), 5iLAF-D13(n=2), RSeT-D13(n=2), t2iLGoY-D13(n=3)이 존재한다. 사전등록 문서는 "day 13"만 지정했다. 통상적 primed reprogramming 경로인 **Primed-D13을 primary**로 하고 나머지 네 배지는 sensitivity로 보고한다.
  3. **HO2의 sex covariate.** GSE307377 processed matrix의 열 이름과 GEO sample metadata를 연결하는 대응표가 GEO에 없다. 따라서 sex를 XIST 및 Y 염색체 유전자(RPS4Y1/DDX3Y/UTY/KDM5D/EIF1AY) 발현으로 추정해 covariate로 사용한다. 추정 결과는 `GSE307377_inferred_sex.tsv`에 기록한다.
     - 사후 검증: 추정된 성비는 전체 5 female / 4 male, young 2F·2M, old 3F·2M로 **GEO metadata와 완전히 일치**했다.

- 2026-09-11 (Phase 2 실행 시, 결과 산출 **전** 기록):
  1. **H1 지표 1(IR burden)은 recount3 자료로 계산할 수 없다.** recount3가 제공하는 것은 exon–exon junction count이며 intronic 염기 커버리지가 없다. IR을 계산하려면 base-level bigwig을 받아야 하는데 시료당 수백 MB로 이 단계의 비용 범위를 벗어난다. 따라서 **H1 판정은 사전등록 지표 2(전역 ΔPSI 분포)와 지표 3(차등 event 수)으로 한다.** 지표 1은 로컬 BAM에서만 계산했고 별도로 보고한다(`direction_probe` P7).
  2. **탐색 지표 추가.** junction 자료에서 계산 가능한 splicing 충실도 지표로 **unannotated(비주석) junction read fraction**을 추가한다. 사전등록 지표가 아니므로 **판정에 사용하지 않고 exploratory로만 보고한다.**
  3. **H1b의 공변량.** unannotated fraction은 depth·read length·mapping률에 민감하므로 sex·log(depth) 외에 platform·avg read length·unique mapping %를 민감도 모형에 넣는다.
