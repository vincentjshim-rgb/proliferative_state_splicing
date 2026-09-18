# Repro-CM 후속 연구: 확정 연구논리 v4

작성일: 2026-09-11
선행 문서: `repro_cm_core_logic_and_analysis_audit_v3_ko.md`, `CLAUDE.md`
새 산출물: `scripts/run_direction_probe.R`, `public_data_tierA/derived/direction_probe/`

이 문서는 "reprogramming / longevity / splicing / secretome" 네 키워드를 어떻게 **하나의 연구논리**로 묶을지에 대한 결정을 기록한다. 결정은 이번 세션에서 새로 수행한 6개 탐색 분석(P1–P6)의 결과에 근거한다. 추측이 아니라 수치가 방향을 정했다.

---

## 1. 확정 연구논리 (한 문장)

> **부분 리프로그래밍기 secretome은 UVA 손상 인간 진피 섬유아세포를 (i) 초기 광적응 반응 및 독립적으로 재현된 광보호 프로그램과 정렬되면서(누적·만성 UV 손상과는 정렬되지 않는다), (ii) 증식 상태와 무관하게 안정적 senescence 궤적에 반대하고, (iii) 세포 내재적 부분 리프로그래밍의 전사체 모사는 아닌 상태로 이동시킨다. 이 decoupling의 두 절반은 서로 다른 모듈이 담당한다 — 핵 내 pre-mRNA processing이 anti-senescence 쪽을, HSF1/UPR-PERK proteostatic stress response가 UVA 정렬 쪽을 담당한다. 연대기적 노화 서명은 되돌리지 않는다.**

영문 working claim:

> **A reprogramming-phase secretome moves UVA-damaged human dermal fibroblasts into a state that is aligned with an early photo-adaptive response and with independently replicated photoprotective programs, but not with accumulated or chronic UV injury, opposes stable senescence trajectories independently of proliferative state, and is not a transcriptional mimic of cell-intrinsic partial reprogramming. The two halves of this decoupling are carried by different module families — nuclear pre-mRNA processing on the anti-senescence side, and the HSF1/UPR-PERK proteostatic stress response on the UVA-aligned side. It does not reverse a chronological-ageing signature.**

v3 대비 달라진 점은 네 가지이며, 모두 새로 수행한 검정의 결과다.

1. **"증식 상태와 무관"이 검정을 통과했다**(§2 P4). senescence 전사체 연구 최대의 반박 — "그냥 세포가 더 잘 자라는 것 아니냐" — 을 정면으로 막는다.
2. **담당 모듈이 특정됐고, 하나가 아니라 둘이다**(§2 P1, P2). 553 pathway 전체에서 정량 순위를 매겼다.
3. **repair 정렬이 독립 재현됐다**(`repro_cm_heldout_results_ko.md`). 사전등록 held-out에서 두 시점 모두 성공.
4. **연대기적 노화 비역전을 주장문 안에 명시**했다. held-out에서 미지지였고 증식 보정 후 소멸했다.

> 참조 축이 어디까지 일반화되는지는 **`repro_cm_axis_characterisation_ko.md`**에 있다. senescence 축은 독립 연구로 거의 손실 없이 전이되지만(median +0.727 → +0.656), UVA 축은 독립 UV 손상으로 전이되지 않고(median −0.017) 광보호 대비에만 로드된다(+0.281).

> 사전등록 held-out 검정의 전체 결과·판정·증거등급 갱신은 **`repro_cm_heldout_results_ko.md`**에 있다. 종합 판정은 **PARTIAL**(6개 중 2 성공, 0 반증)이며, 그 문서의 §6 주장문이 최신 기준이다.

## 2. 이번 세션에서 새로 확인한 사실

전부 `scripts/run_direction_probe.R` 하나로 재현된다. 출력은 `public_data_tierA/derived/direction_probe/`.

### P1. pathway 수준 decoupling — 새로운 primary test

v3는 `local vs senescence` (ρ=-0.528)와 `local vs UVA` (ρ=+0.079, 비유의)를 **따로** 보고했고, 후자가 비유의하다는 점이 약점이었다. 그러나 decoupling 가설이 예측하는 것은 두 상관이 각각 유의한 것이 아니라, **local NES가 (UVA − senescence) 합성축과 정렬되는 것**이다. 이것을 직접 검정하면:

| 검정 | ρ | p | n |
|---|---:|---:|---:|
| local NES vs (UVA − senescence) decoupling 축 | **+0.476** | **1.28 × 10⁻³²** | 553 pathway |
| local NES vs UVA meta NES | +0.079 | 0.063 | 553 |
| local NES vs senescence meta NES | −0.528 | 4.3 × 10⁻⁴¹ | 553 |

첫 줄이 pathway 수준의 올바른 primary endpoint다. 이 한 줄로 "pathway 수준에서는 UVA 정렬이 유의하지 않다"는 v3의 약점이 해소된다. 단, **UVA 단독 상관이 비유의하다는 사실은 계속 명시해서 보고한다.**

### P2. 어떤 모듈이 decoupling을 담당하는가 (module × axis z-map)

발현량 계층화 순열로 각 모듈의 좌표 이동을 z로 환산했다. 양수 = 그 축에서 모듈이 함께 상승.

| module | local | UVA | senescence | MPTR reprog | fibroblast age | rodent max lifespan¹ | ITP max lifespan |
|---|---:|---:|---:|---:|---:|---:|---:|
| **splicing** | **+2.93** | +0.60 | **−7.40** | −0.18 | −4.06 | −1.90 | **+4.48** |
| **pre-mRNA processing** | **+3.51** | +0.75 | **−7.43** | −0.64 | −4.00 | −2.09 | **+5.13** |
| **mRNA export** | **+2.54** | +0.50 | **−4.59** | −1.01 | −2.50 | −1.91 | +2.92 |
| **3'-end processing** | **+2.07** | +0.63 | **−6.50** | −0.63 | −2.69 | −1.64 | +4.07 |
| splicing-factor regulon (curated) | +2.59 | +1.35 | −6.66 | −0.71 | −1.61 | −0.82 | +2.80 |
| **HSF1 heat-shock response** | **+2.71** | **+1.99** | −1.01 | −0.61 | −1.36 | −1.52 | +0.68 |
| **UPR / PERK** | **+2.65** | **+3.17** | −0.13 | +2.24 | −0.89 | −1.95 | +1.41 |
| ROS detoxification / NRF2 | +1.88 | +3.75 | **+2.20** | +1.57 | −1.11 | −1.98 | +1.39 |
| autophagy | +1.87 | +1.65 | +1.54 | +0.76 | +0.43 | −1.35 | +1.33 |
| NMD | +4.22 | **−4.05** | −6.84 | −0.01 | −5.09 | −1.18 | +4.23 |
| translation | +4.01 | **−4.12** | −6.94 | −0.37 | −5.30 | −2.16 | +4.78 |
| rRNA processing | +5.40 | **−2.95** | −7.40 | −0.58 | −6.03 | −1.41 | +5.46 |
| cell cycle / DNA replication | +3.96 | **−7.85** | −14.22 | −4.92 | −18.12 | −5.60 | +2.18 |
| ECM organization | −1.97 | −3.95 | +4.34 | +0.60 | +1.00 | −1.18 | −2.15 |
| interferon signaling | +5.73 | +0.10 | −4.70 | −3.03 | −1.15 | −6.73 | +0.30 |

¹ Tyshkovskiy 2026의 rodent 시트에는 발현량 열(`logCPM`)이 없어 발현량-매칭 null을 정의할 수 없다. 이 두 열만 **매칭 없는 순열**이고 나머지는 전부 발현량-매칭이다. `P2_module_axis_zmap.tsv`의 `null_type` 열에 기록돼 있다.

핵심 판독 세 가지:

1. **decoupling을 만들 수 있는 모듈은 두 종류뿐이다.**
   - **핵 내 pre-mRNA processing**(splicing, 3'-end, export): senescence에 강하게 반대(−4.6 ~ −7.4)하지만 UVA 정렬은 약하다(+0.5 ~ +1.4, 비유의).
   - **HSF1 / UPR-PERK proteostatic stress response**: UVA에 강하게 정렬(+2.0 ~ +3.2)하지만 anti-senescence는 약하다(−0.1 ~ −1.0).

   즉 decoupling의 두 절반을 서로 다른 모듈이 담당한다. "RNA processing이 전부"라고 쓰면 부정확하다.

2. **local에서 상승하는 나머지 모듈은 전부 UVA 축에서 음수다.** translation(−4.12), NMD(−4.05), rRNA(−2.95), cell cycle(−7.85). local anchor에서 신호 크기는 이쪽이 더 크지만 decoupling에는 기여할 수 없다. 이것이 P1에서 UVA 단독 pathway 상관이 0에 가까웠던 이유다.

3. **ROS detox는 decoupling 모듈이 아니다.** UVA 정렬은 가장 강하지만(+3.75) senescence에서도 상승한다(+2.20). 광보호 서사에 쓰고 싶어지는 모듈이지만 데이터가 지지하지 않는다.

### P3–P4. decoupling은 증식 신호가 아니다 — 핵심 신규 결과

GSE179848의 345개 섬유아세포 시료로 각 유전자의 **proliferation loading**(cell-cycle meta-gene과의 상관)을 추정한 뒤, decoupling index D를 다시 계산했다.

| 분석 | n genes | ρ_UVA | ρ_senescence | **D** |
|---|---:|---:|---:|---:|
| 보정 없음 | 2,216 | 0.289 | −0.172 | 0.461 |
| 세 축 모두 proliferation loading으로 잔차화 | 2,216 | 0.320 | −0.147 | **0.468** |
| 증식 무관 유전자만 (\|loading\| < 0.2) | 715 | 0.384 | −0.162 | **0.546** |
| 증식 연관 유전자만 (\|loading\| ≥ 0.4) | 861 | 0.176 | −0.223 | 0.400 |

또한 `ρ(local anchor score, proliferation loading) = +0.108`로, local anchor 자체가 증식 서명이 아니다.

**D는 증식을 제거해도 유지되고, 증식 무관 유전자에서 오히려 더 강하다(0.546).** v3는 Reactome cell-cycle 유전자 *목록*을 제거하는 방식만 시험했는데(D=0.475), 이번에는 데이터 기반 증식 loading으로 **연속적으로** 제거했다. 이것이 훨씬 강한 검정이며, 이 결과가 §1 claim에 "증식 상태와 독립적"을 넣을 수 있게 해준다.

### P5. splicing factor 발현과 증식은 이 자료에서 분리되지 않는다 (**2026-09-11 수정**)

GSE113957(건강 공여자 104명, 22–96세) 섬유아세포에서 모듈 점수의 연령 기울기를 증식 보정 전후로 비교했다.

| module | r(proliferation) | β/10yr (보정 전) | p | β/10yr (증식 보정) | p |
|---|---:|---:|---:|---:|---:|
| splicing | **+0.891** | −0.101 | 5.1e−04 | −0.016 | 0.36 |
| pre-mRNA processing | +0.902 | −0.096 | 8.6e−04 | −0.009 | 0.62 |
| rRNA processing | +0.864 | −0.101 | 4.1e−04 | −0.021 | 0.24 |
| translation | +0.569 | −0.131 | 1.8e−04 | **−0.079** | **0.012** |
| ECM organization | −0.081 | −0.018 | 0.017 | **−0.023** | **0.009** |

> **최초 판독은 틀렸다.** 나는 이것을 "splicing factor 발현은 노화 축으로 쓸 수 없다"는 결정적 음성으로 적었다. 그렇지 않다.

Reactome 1,169개 pathway 전체로 같은 검정을 돌려 보면, 증식 보정 후 연령 연관의 **생존 여부는 생물학이 아니라 증식과의 공선성이 결정한다.**

| 증식과의 \|r\| | 생존 | 소멸 | 생존률 |
|---|---:|---:|---:|
| < 0.5 | 57 | 17 | **77%** |
| 0.5–0.7 | 23 | 27 | 46% |
| 0.7–0.8 | 7 | 46 | 13% |
| 0.8–0.9 | 8 | 141 | **5%** |
| > 0.9 | 30 | 294 | **9%** |

연령 연관 pathway 650개 중 125개(19.2%)만 보정을 통과한다. splicing 모듈은 r≈0.89–0.90 구간에 있고, **그 구간의 어떤 모듈이든 5–9%만 살아남는다.** 즉 splicing의 소멸은 그 공선성 수준에서 기대되는 그대로이며 splicing 특이적 정보를 담고 있지 않다. r=0.92일 때 VIF는 6.5로, 보정이 사실상 퇴화(degenerate)한다.

게다가 배양 섬유아세포에서 증식은 **교란(confounder)이 아니라 매개(mediator)일 수 있다.** 공여자 연령 → 증식능 저하 → splicing 포함 생합성 모듈의 협응적 감소라면, 증식을 보정하는 것은 측정하려는 인과 경로 자체를 제거하는 over-adjustment다.

**올바른 결론:**

> 이 자료는 splicing 모듈 발현과 증식 상태를 **분리할 수 없다.** 따라서 splicing이 독립적인 노화 축이라는 주장도, 아니라는 주장도 이 자료로는 할 수 없다. Harries/Latorre 계열 문헌을 반박하는 근거가 아니다.

이 한계는 §2 P4(decoupling의 증식 독립성)에는 영향을 주지 않는다. 거기서는 유전자 수준 proliferation loading을 쓰고 **anchor와 loading의 상관이 0.108에 불과**해 공선성 문제가 없다. 두 분석은 regime이 다르다.

### P6. 로컬 n=1로는 junction 수준 splicing을 주장할 수 없다

로컬 STAR `SJ.out.tab`(HDF/REP/IPS)에서 alternative junction의 PSI를 계산하고, "세 시료 간 진짜 차이가 없다"는 binomial null과 비교했다.

| group | min coverage | n testable | ρ 관측 | ρ binomial null | \|ΔPSI\|≥0.1 후보 | noise 기대 | 경험적 FDR |
|---|---:|---:|---:|---:|---:|---:|---:|
| 5'-anchored | 20 | 21,554 | 0.498 | 0.467 | 606 | 446 | **0.74** |
| 5'-anchored | 100 | 5,389 | 0.473 | 0.454 | 15 | 5 | 0.35 |
| 3'-anchored | 20 | 21,538 | 0.509 | 0.465 | 653 | 466 | **0.71** |
| 3'-anchored | 200 | 2,349 | 0.501 | 0.445 | 3 | 0 | 0.02 |

두 비교자 사이의 ΔPSI 상관 ρ≈0.50은 **대부분 psi_REP를 공유하는 데서 오는 구조적 인공물**이다(null ρ≈0.46). 실제 초과분은 0.03–0.05에 불과하다. 낮은 coverage에서 후보의 70% 이상이 잡음이다.

**결론: `repro_cm_splicing_ubuntu` 파이프라인을 지금 4개 시료에 돌려 얻는 rMATS 후보 목록은 논문의 증거가 될 수 없다.** 높은 coverage 구간에서만 소수(3–17개)가 신뢰 가능하다. 로컬 splicing은 *증거*가 아니라 *적용 예시*로만 쓸 수 있다.

### P7. 로컬 intron retention도 splicing 개선을 지지하지 않는다

GENCODE v41로 동일 파라미터 재정렬한 BAM 4개에서, 다른 유전자의 exon·유전자 span과 겹치지 않는 intronic/exonic 블록의 read 밀도를 계산했다(`-split` 적용, spliced read가 intron을 가로지르는 오계수 방지). 네 라이브러리 모두에서 충분히 덮인 6,749개 유전자 기준.

| 시료 | intronic read fraction | 길이보정 IR | 유전자별 IR ratio 중앙값 |
|---|---:|---:|---:|
| Repro-CM | **5.56%** | 0.00627 | 0.0152 |
| HDF-CM | 5.19% | 0.00583 | 0.0136 |
| iPSC-CM | 4.93% | 0.00552 | 0.0122 |
| UVA 15J, 0 h | 2.86% | 0.00313 | 0.0065 |

- **24 h 세 조건 중 Repro-CM의 intron retention이 가장 높다.** splicing 효율 개선과 반대 방향이다.
- 0 h 시료는 약 2배 낮지만 경과시간·library batch와 완전히 교락돼 있어 해석할 수 없다.
- 유전자별 IR 변화는 발현 anchor score와 ρ = −0.25로 상관한다. 즉 이 지표는 발현량과 독립적이지 않다.
- **조건당 n=1이므로 기술적 기술(descriptive)이며 p-value를 붙이지 않는다.**

P6(junction), P7(intron retention) 두 각도에서 **로컬 자료는 splicing 개선 주장을 지지하지 않는다.** P5(발현)는 지지도 반박도 하지 못한다(공선성). splicing을 논문에 쓰려면 Phase 2에서 반복수가 있는 공개 자료로 검정해야 한다.

---

## 3. 버린 방향과 그 이유

| 버린 논리 | 근거 | 판정 |
|---|---|---|
| "Repro-CM이 부분 리프로그래밍을 모사한다" | MPTR 축에서 전체 ρ≈0, 모든 RNA processing 모듈 z ≈ −0.2~−1.0 | 사용 불가 |
| "Repro-CM이 longevity program을 켠다" | Tyshkovskiy 2026 rodent max-lifespan ρ=−0.055(반대); 2019 개입 signature에서 splicing NES 0.68–0.84(비유의); **같은 출처의 두 max-lifespan signature가 서로 부호가 반대**(rodent aggregate에서 splicing z=−1.90, ITP에서 z=+4.48) | 사용 불가 — 참조 signature 간 불일치 |
| "splicing factor 발현 회복이 노화를 되돌린다" | P5: 증식과 r≈0.89로 공선성이 커서 **분리 불가**. 문헌 반박 근거가 아니며, 지지 근거도 아니다 | **판정 보류** |
| "로컬 rMATS 후보가 기전이다" | P6: FDR 0.7 | 증거로 사용 불가 |
| "Repro-CM이 splicing 효율을 개선한다" | P7: intron retention이 세 조건 중 가장 높음, 발현과 ρ=−0.25로 교락 | **로컬 자료로는 사용 불가** |
| "miR-302/367이 매개한다" | 표적 GSEA NES=−1.05, FDR=0.399 | Supplementary identity로만 |

이 목록은 논문 Limitations가 아니라 **Results에 명시적으로 보고한다.** CLAUDE.md §11-11 원칙이며, 동시에 claim boundary를 만드는 자산이다.

---

## 4. 네 키워드의 최종 위치

| 키워드 | 논문 내 역할 | 증거 등급 |
|---|---|---|
| **secretome** | 개입 그 자체. 유일무이한 로컬 자산 | 로컬 기술적 A / 생물학적 n=1 |
| **splicing** | decoupling을 담당하는 모듈로 특정됨. **단, 발현이 아니라 splicing 결과(PSI/IR)로 검정해야 함** | 현재 B− → Phase 2에서 결정 |
| **reprogramming** | (1) secretome의 출처, (2) **경계 검정**: 외인성 secretome ≠ 세포 내재적 OSK. 음성이 오히려 기여 | 음성 결과로 A− |
| **longevity** | 선행 논문 *C. elegans* 결과를 배경으로 인용 + **포유류 lifespan signature끼리 서로 불일치한다는 사실**을 명시. 기전 가설은 Discussion 한정이며, splicing보다는 HSF1/UPR 스트레스 저항성 축이 무척추동물 수명과 더 잘 연결되는 보존 노드다 | Discussion 수준 |

longevity를 Results의 주장으로 올리면 자기 데이터와 충돌한다. Discussion에서 "왜 worm lifespan은 연장됐는데 포유류 lifespan signature와는 정렬되지 않는가"를 **열린 질문으로 제시**하는 것이 훨씬 강하다.

---

## 5. Phase 2: splicing을 발현이 아니라 결과로 검정한다

P5가 splicing **발현** 경로를 닫았지만, 동시에 진짜 질문을 만들었다.

> **PSI와 intron retention은 유전자 내부 비율이므로 발현량처럼 증식에 자동으로 비례하지 않는다. 그렇다면 splicing 결과는 증식과 분리 가능한 senescence 축인가?**

이 질문은 아직 아무도 이 방식으로 검정하지 않았고, 답이 어느 쪽이든 보고 가치가 있다. 그리고 이것이 **npj Aging 수준과 그 이하를 가르는 지점**이다.

### 사전등록 가설 (결과를 보기 전에 고정)

- **H1 (primary):** senescent 섬유아세포는 동일 population doubling의 quiescent 섬유아세포와 비교했을 때에도 splicing 결과(전역 ΔPSI 분포, intron retention burden)가 다르다. → 증식정지만으로 설명되지 않는다.
- **H2:** H1에서 정의한 senescence splicing 축은 UVA 급성 손상 축과 구별된다(decoupling이 splicing 수준에서도 성립).
- **H3 (application):** 로컬 Repro-CM 시료는 H1 축에서 anti-senescence 방향에 위치한다. *로컬 n=1이므로 이것은 검정이 아니라 배치(placement)이며, 사전에 정의한 참조 분포 위에 한 점을 찍는 것이다.*

성공/실패 기준, exclusion, BH 보정 계획은 데이터를 열기 전에 별도 파일로 고정한다.

### 필요 데이터와 비용 — recount3로 대부분 해결됨

FASTQ를 받기 전에 recount3/Snaptron 가용성을 먼저 확인했고, **senescence 축에 필요한 세 연구가 모두 junction count로 존재했다.** 이미 받았다(`public_data_tierA/recount3/`, SHA256 기록, 총 124 MB).

| 연구 | 설계상 역할 | recount3 | 확보 |
|---|---|---|---|
| GSE113957 (SRP144355) | 공여자 143명, 1–96세 → **H1b를 n=143으로 검정** | ✅ | junction + gene_sums + QC |
| GSE93535 (SRP096629) | **Q(PD15) vs SIPS(PD15)** → H1 핵심 | ✅ | junction + gene_sums + QC |
| GSE109700 (SRP131506) | deep RS 복제 | ✅ | junction + gene_sums + QC |
| GSE191055 / GSE240226 / GSE302943 / GSE179848 / GSE226189 / GSE165177 | 2021년 이후 | ❌ 404 | FASTQ 필요 |

recount3 freeze가 2019년경이라 최신 연구는 없다. 따라서 **H2(UVA splicing 축)만 FASTQ가 필요하다**: GSE240226 26.04 GiB (paired-end, UVA + 2개 rescue arm). H1/H1b는 추가 다운로드 없이 진행 가능하다.

recount3 junction 좌표계(0-based/1-based)와 로컬 STAR `SJ.out.tab` 좌표계를 실제 데이터로 대조한 뒤 결합한다. 이 확인 전에는 두 자료를 합치지 않는다.

### 이미 확보된 자산 (추가 비용 0)

- **GENCODE v41 동일 파라미터로 재정렬한 BAM 4개**: `public_data_tierA/derived/local_bam_v41/` (HDF, REP, IPS, **UVA0 = 15J_0hr**). 2022년 BAM과 달리 네 시료가 동일 조건이라 직접 비교 가능하다.
- **IR용 annotation**: `public_data_tierA/derived/ir_annotation_v41/` — 다른 유전자의 exon·유전자 span과 겹치지 않는 exonic 216,772 블록(19,183 genes) / intronic 146,529 블록(12,288 genes).
- 인간 STAR index(v41), STAR 2.7.3a, samtools, bedtools, trim_galore, sratoolkit, conda.
- 시스템: 62 GB RAM, 24 threads, 디스크 여유 약 950 GB.

## 6. Figure 1–6 재설계

각 figure는 검정하는 주장 한 문장과 증거 등급을 가진다.

### Figure 1 — 문제 정의: UVA 축과 senescence 축은 같은 축이 아니다
공개 UVA 2연구·senescence 4연구의 내부 재현성과 두 축의 분리. **Repro-CM 효과의 재현이 아니라 해석용 reference의 재현임을 legend에 명시.** 등급 B.

### Figure 2 — 로컬 anchor의 기술적 신뢰성
QC, 두 comparator 일치(8,921/12,319), strand/threshold/thinning 안정성. **n=1을 패널 안에 표기.** 등급: 기술 A / 생물 F.

### Figure 3 — 핵심: stress–senescence decoupling이 증식과 무관하다
A: 6개 연구 forest plot (가로 막대는 CI가 아니라 matched-null 95% 구간임을 명시)
B: D = 0.460과 10,000회 발현량-매칭 순열
C: **증식 잔차화 후 D = 0.468, 증식 무관 유전자에서 D = 0.546** ← 신규, 이 논문의 방법론적 핵심
D: UVA 잔차화 후 senescence 관계 ρ = −0.195
등급 **A−**.

### Figure 4 — 어떤 모듈이 decoupling을 담당하는가
A: 553 pathway에서 local NES vs (UVA − senescence), **ρ = 0.476, p = 1.3 × 10⁻³²** ← 새 primary pathway test
B: module × axis z-map (§2 P2 표)
C: RNA processing만 UVA 축에서 음수가 아님을 보이는 패널
D: nuisance 유전자군 제거 민감도
**UVA 단독 pathway 상관이 비유의(p=0.063)라는 사실을 같은 figure에 표시한다.** 등급 A− (anti-senescence) / B (UVA).

### Figure 5 — splicing: 발현이 아니라 결과로
A: **P5 표 — splicing factor 발현은 증식과 r=0.92, 보정 후 연령효과 소실** (음성 결과를 먼저 보여준다)
B–D: Phase 2 결과 (Q vs SIPS splicing 결과, IR burden, decoupling의 splicing 수준 재현)
Phase 2 미수행 시 이 figure는 A만 남기고 나머지는 Limitation으로 강등한다. 등급 미정.

### Figure 6 — 사전등록 held-out 일반화 **[완료]**
GSE297233·GSE307377·GSE116968·GSE149694, 6개 primary test, 사전등록 SHA256 `41e4cf4b`. 종합 판정 **PARTIAL**(2 성공, 0 반증, 3 boundary 유지, 1 미지지). 패널 A/B/C 제작 완료. 전체 결과는 `repro_cm_heldout_results_ko.md`.
- **A**: 6개 primary test forest plot. 막대는 matched-null 95% 구간(CI 아님).
- **B**: repair 정렬의 독립 재현 — 4개 광보호 개입, 2개 연구, 2개 파장에서 ρ = +0.094 ~ +0.127. **이 논문에서 가장 강한 단일 결과.**
- **C**: reprogramming 경계. primary 2건 모두 경계 내. naive 배지 sensitivity arm 1건(t2iLGoY-D13, +0.180)이 경계를 넘었고 삭제하지 않고 보고.
등급 **A−**(B 패널) / B(A·C 패널). 결손: **held-out 6개 중 senescence 데이터셋이 없다** — GSE306957 추가 필요.

---

## 7. 저널 판단 (held-out 결과 반영)

| 구성 | 상태 |
|---|---|
| Figure 3C 증식 독립성 (D 0.461→0.468, 증식 무관 유전자 0.546, 10개 decile 전부 양수) | **완료** |
| Figure 4A 새 pathway primary test (ρ = 0.476, p = 1.3 × 10⁻³²) | **완료** |
| Figure 4B 두 모듈 패밀리 분해 | **완료** |
| Figure 6 사전등록 held-out (PARTIAL, repair 재현 성공) | **완료** |
| Figure 6D held-out senescence (GSE306957, CONFIRMED) | **완료** |
| Figure 5B 로컬 intron retention (음성) | **완료** |
| Phase 2 splicing 결과 수준 검정 (H1/H1b) | 미수행, recount3로 추가 비용 없이 가능 |

핵심 주장의 **두 축이 모두 사전등록 held-out에서 지지**를 받았다 — repair 정렬(HO3a/b)과 anti-senescence(HO5a/b). 증식 독립성도 held-out에서 co-primary로 사전 지정해 통과했다. 사전등록·독립 재현·증식 독립성·명시적 음성 결과가 모두 갖춰졌으므로 **npj Aging 투고 가능** 수준이다.

남은 선택지는 하나다.

- **Phase 2 H1/H1b** — splicing 결과(PSI/IR)가 증식과 분리되는지를 반복수가 있는 공개 자료로 검정. recount3 자료가 이미 있어 추가 다운로드가 필요 없다. 성공하면 정량 framework 하나가 더 생기고, 실패해도 분야에 유용한 음성 결과다. 로컬 자료는 P5·P6·P7 세 각도 모두에서 splicing 개선을 지지하지 않으므로, **splicing을 논문에 넣으려면 이 검정이 필수**다.

Aging Cell은 기능 검증 없이는 여전히 어렵다.

논문의 독창성 네 가지. 선행 논문과 중복되지 않는다.

1. **증식 독립성을 검정한 decoupling metric** — 재사용 가능한 정량 framework
2. **모듈 분해로 담당 축을 특정** — 두 모듈 패밀리, 553 pathway 전체 보고
3. **사전등록 held-out 일반화** — 무관한 네 광보호 개입에서 repair 정렬 재현
4. **명시적 음성 결과 집합** — longevity, reprogramming mimicry, miRNA, 연대기적 노화, splicing factor 발현

## 8. 다음 작업 순서

완료된 항목은 취소선 대신 **[완료]**로 표시한다.

1. **[완료]** recount3/Snaptron 가용성 확인 → senescence 축 3개 연구 확보(124 MB)
2. **[완료]** Phase 2 + held-out 사전등록 문서 고정 (SHA256 `41e4cf4b`)
3. **[완료]** held-out 4개 데이터셋 다운로드 및 6개 primary test 실행 → PARTIAL
4. **[완료]** 로컬 4개 시료(15J_0hr 포함) 동일 파라미터 재정렬
5. **[완료]** Figure 3C · 4A · 4B · 6A/B/C 패널 제작
6. **[완료]** IR burden 산출 → splicing 개선 미지지 (P7)
7. **[완료]** GSE306957 2차 사전등록 및 HO5 검정 → **CONFIRMED**
8. Phase 2 H1/H1b: recount3 junction으로 splicing 결과가 증식과 분리되는지 검정 ← 다음 우선순위
9. H2(UVA splicing 축)가 필요하면 GSE240226 FASTQ(26 GiB)만 추가
10. Figure 1·2 최종 조판, 각 패널의 증거등급 표기
11. 최종 한 문장 주장이 모든 패널과 일치하는지 재점검
