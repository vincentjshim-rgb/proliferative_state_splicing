# Held-out generalization 결과

실행일: 2026-09-11
사전등록: `repro_cm_preregistration_phase2_ko.md` (SHA256 `41e4cf4b548f74904db7ee8d6fefcf9e18cd801972bcdc1cec95de59a8788bd7`)
스크립트: `scripts/run_heldout_validation.R` · `scripts/make_figure6_heldout.R`
출력: `public_data_tierA/derived/heldout_validation/` · `public_data_tierA/derived/manuscript_figures_v4/Figure6*.png`

**대비·예측·성공기준·제외규칙은 데이터를 열기 전에 고정했고, 실행 후 변경하지 않았다.** 두 건의 사전 명확화(매핑률 측정 대상, day-13 배지 특정)는 결과를 보기 전에 변경 이력에 기록했다.

---

## 1. 총평

**종합 판정: PARTIAL.** 6개 primary test 중 2개 성공, 0개 반증, 3개 boundary 유지, 1개 미지지.

| ID | dataset | contrast | 사전 예측 | ρ | 증식 보정 ρ | BH-FDR | 판정 |
|---|---|---|---|---:|---:|---:|---|
| HO1 | GSE297233 | OSK +dox vs −dox, day 4 (96세) | 양의 정렬 없음 | +0.007 | −0.011 | 0.910 | BOUNDARY HELD |
| HO2 | GSE307377 | old vs young 진피 섬유아세포 | ρ < 0 | −0.032 | +0.008 | 0.197 | **NOT SUPPORTED** |
| HO3a | GSE116968 | Pre-Red vs UV 1 h, injury 잔차화 | ρ > 0 | **+0.122** | +0.131 | **0.0003** | **SUCCESS** |
| HO3b | GSE116968 | Pre-Red vs UV 4 h, injury 잔차화 | ρ > 0 | **+0.098** | +0.116 | **0.0003** | **SUCCESS** |
| HO4a | GSE149694 | fibroblast day 7 vs day 3 | 양의 정렬 없음 | −0.015 | −0.016 | 0.883 | BOUNDARY HELD |
| HO4b | GSE149694 | primed day 13 vs day 3 | 양의 정렬 없음 | −0.050 | −0.053 | 0.063 | BOUNDARY HELD |

공유 유전자는 모든 검정에서 2,196–2,236개였다. gene-symbol 매핑률(발현 유전자 기준)은 GSE297233 91.2%, GSE307377 100%, GSE149694 88.7%로 모두 제외 기준을 통과했다.

---

## 2. 가장 강한 결과: repair-aligned component의 독립 재현

v3에서 `B+ / 독립 rescue 연구 필요`로 유보됐던 항목이 **재현됐다.**

| 연구 | 개입 | 손상 | injury 제거 후 ρ | 단계 |
|---|---|---|---:|---|
| GSE240226 | Maifuyin 추출물 | UVA | +0.127 | discovery |
| GSE240226 | succinate | UVA | +0.094 | discovery |
| **GSE116968** | **적색광 전처치, 1 h** | **UVB** | **+0.122** | **held-out** |
| **GSE116968** | **적색광 전처치, 4 h** | **UVB** | **+0.098** | **held-out** |

서로 무관한 네 가지 광보호 개입(생약 추출물·대사산물·광생물조절 2시점), 두 독립 연구, 서로 다른 UV 파장에서 효과크기가 **+0.094 ~ +0.127로 사실상 동일**하다. 증식 보정 후에는 오히려 커진다(+0.116, +0.131).

### 결정적 보강: UVB injury 축과는 음의 상관

GSE116968에서 로컬 anchor는 **UVB injury 축 자체와 음의 상관**이다.

- UV vs Control, 1 h: ρ = −0.094 (p = 5 × 10⁻⁴)
- UV vs Control, 4 h: ρ = −0.169 (p = 1 × 10⁻⁴)

즉 이 연구에서는 "잔존 손상이라서 rescue와 같은 방향으로 보인다"는 대안설명이 **구조적으로 불가능하다.** injury와는 반대, rescue와는 같은 방향이다. 잔차화 전 원값도 이미 양수였다(+0.130, +0.105).

동시에 이것은 claim boundary를 좁힌다. **양의 stress 정렬은 UVA에 특이적이며(GSE240226 +0.255, GSE302943 +0.313), 급성 UVB 손상에는 일반화되지 않는다.** 본문에서 "acute photostress"를 "UVA"로 한정해야 한다.

---

## 3. 실패한 검정: 연대기적 노화 (HO2)

GSE307377에서 ρ = −0.032로 방향은 예측과 같았으나 유의하지 않았고(BH-FDR = 0.197), **증식 보정 후 +0.008로 소멸**했다.

| cohort | n | ρ |
|---|---:|---:|
| GSE113957 | 104 | −0.075 |
| GSE226189 | 82 | −0.037 |
| **GSE307377 (held-out)** | **9 (young 4, old 5)** | **−0.032, 보정 후 +0.008** |

세 cohort 모두 약한 음의 방향으로 일치하지만, 사전등록한 성공기준을 held-out에서 충족하지 못했다. 공여자 9명이라 검정력이 낮은 것도 사실이나, **검정력 부족을 이유로 성공으로 재해석하지 않는다.**

증식 보정 후 소멸한다는 점은 `direction_probe/P5`(splicing factor 발현의 연령 감소가 증식 보정 후 소멸, p = 0.71)와 같은 방향의 결과다. 따라서 다음을 명시적으로 결론에 넣는다.

> **Repro-CM anchor는 증식 상태를 보정하면 연대기적 노화 서명을 되돌리지 않는다.**

이는 v3의 `chronological-age reversal: C−, DO NOT CLAIM` 판정과 일치하며, held-out에서 이를 확정했다. senescence 궤적에 대한 반대 방향(D = 0.460, 증식 보정 후 0.468)과는 **별개의 축**이므로 핵심 주장을 반증하지는 않는다. 다만 두 축을 논문에서 혼동해 쓰면 안 된다.

성별은 GEO metadata에 열 대응표가 없어 XIST·Y 염색체 유전자로 추정했고, 추정 성비(전체 5F/4M, young 2F·2M, old 3F·2M)가 **GEO metadata와 완전히 일치**해 검증됐다.

---

## 4. Reprogramming 경계 (HO1, HO4)

primary 두 건 모두 사전 정의한 |ρ| = 0.15 경계 안에 있다. OSK 유도 자체(+0.007)도, primed reprogramming day 13(−0.050)도 로컬 anchor와 양의 정렬이 없다. mutant O4YRSK sensitivity도 −0.001이다.

따라서 **외인성 secretome 반응은 세포 내재적 부분 리프로그래밍의 전사체 모사가 아니다.** GSE165177 MPTR(ρ ≈ 0)에 이어 두 독립 연구에서 같은 결론이다.

### 보고해야 할 경계 위반 (sensitivity)

| day-13 배지 | ρ | p |
|---|---:|---:|
| Primed-D13 (**primary**) | −0.050 | 0.032 |
| NHSM-D13 | +0.007 | 0.864 |
| 5iLAF-D13 | +0.001 | 0.997 |
| RSeT-D13 | −0.036 | 0.692 |
| **t2iLGoY-D13** | **+0.180** | **1 × 10⁻⁴** |

naive 상태 유도 배지 t2iLGoY의 day 13은 사전 정의한 경계(+0.15)를 넘었다. 사전등록상 sensitivity arm은 판정을 바꾸지 않지만 **삭제하지 않고 그대로 보고한다.** 해석은 열어 둔다. naive 전환 특이적 프로그램과 부분적으로 겹칠 가능성은 후속 검정 대상이지, 현재 데이터로 결론 낼 사안이 아니다.

---

## 5. 이번 held-out 설계 자체의 결손 → **해소됨**

정직하게 기록한다. 이 1차 사전등록의 6개 primary test에는 **senescence 데이터셋이 하나도 없었다.** 핵심 주장의 절반(anti-senescence)이 held-out으로 검정되지 않은 상태였다. GSE109700·GSE93535·GSE179848·GSE191055는 이미 discovery meta-signature에 쓰였으므로 held-out이 아니다.

이 결손은 **2차 사전등록(GSE306957)으로 메웠고 결과는 CONFIRMED이다.** `repro_cm_heldout2_senescence_results_ko.md` 참조.

- HO5a senescent vs proliferating: ρ = −0.127, BH-FDR = 1 × 10⁻⁴
- HO5b 동일 대비, 증식 잔차화: ρ = −0.101, BH-FDR = 1 × 10⁻⁴
- 새 조직(폐 섬유아세포), 새 유도기전(20 Gy X선), 효과크기 축소 없음

2차는 **독립된 BH 가족**으로 보정했다. 1차 6개 가족은 닫힌 상태로 두고 재보정하지 않았다.

---

## 6. 결과를 반영한 핵심 주장 수정

수정 전(v4 §1)과 수정 후를 함께 남긴다.

**수정 후:**

> 부분 리프로그래밍기 secretome은 UVA 손상 인간 진피 섬유아세포에서 (i) UVA 특이적 스트레스 반응 및 **독립적으로 재현되는 광보호·repair 프로그램**과 정렬되고, (ii) 증식 상태와 무관하게 안정적 senescence 궤적에 반대하며, (iii) 세포 내재적 부분 리프로그래밍의 전사체 모사가 아닌 상태를 만든다. **연대기적 노화 서명은 되돌리지 않는다.**

영문:

> A reprogramming-phase secretome moves UVA-damaged human dermal fibroblasts into a state that is aligned with UVA-specific stress and with independently replicated photoprotective repair programs, opposes stable senescence trajectories independently of proliferative state, and is not a transcriptional mimic of cell-intrinsic partial reprogramming. It does not reverse a chronological-ageing signature.

변경점 세 가지: `acute photostress` → `UVA-specific`; repair 정렬에 `independently replicated` 추가; 연대기적 노화 비역전을 **주장문 안에** 명시.

---

## 7. 증거 등급 갱신

| 항목 | v3/v4 등급 | 갱신 | 근거 |
|---|---|---|---|
| repair-aligned component | B+ | **A−** | 2개 독립 연구·4개 개입·파장 무관, injury와는 음의 상관 |
| stress–senescence decoupling | A− | **A−** (유지) | 증식 잔차화 통과(D 0.461→0.468, 증식 무관 유전자 0.546) |
| pathway-level decoupling | C (UVA global 비유의) | **A−** | 합성축 검정 ρ = 0.476, p = 1.3 × 10⁻³² |
| extrinsic ≠ intrinsic reprogramming | D (주장 불가) | **B+ (음성 결과로서)** | 3개 독립 연구에서 일관, 단 naive arm 1건 예외 |
| chronological-age reversal | C− | **D / 명시적 배제** | held-out 미지지, 증식 보정 후 소멸 |
| anti-senescence의 held-out 재현 | — | **A−** | GSE306957에서 ρ=−0.127, 증식 보정 후 −0.101, 효과크기 축소 없음 |
| 증식 독립성 | A− (내부 민감도만) | **A** | held-out에서 co-primary로 사전 지정해 통과 |
