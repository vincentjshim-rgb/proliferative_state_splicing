# Held-out senescence 결과 (HO5, GSE306957)

실행일: 2026-09-11
사전등록: `repro_cm_preregistration_heldout2_senescence_ko.md` (SHA256 `6e5a7703cc7a5ac6c695b4bd10ef53e872c709371e18ef9135cfde33f7be5ff4`)
스크립트: `scripts/run_heldout2_senescence.R`
출력: `public_data_tierA/derived/heldout2_senescence/`

## 1. 총평: **CONFIRMED**

1차 held-out(`repro_cm_heldout_results_ko.md` §5)에서 스스로 기록한 설계 결손 — *핵심 주장의 절반인 anti-senescence가 held-out 검정을 받지 않았다* — 가 해소됐다.

| ID | 검정 | 사전 예측 | ρ | BH-FDR | 판정 |
|---|---|---|---:|---:|---|
| HO5a | senescent vs proliferating control (pLenti+siScramble 통합, 하위실험 batch 보정) | ρ < 0 | **−0.127** | 1 × 10⁻⁴ | **SUCCESS** |
| HO5b | 동일 대비, 두 축 모두 proliferation loading으로 잔차화 | ρ < 0 | **−0.101** | 1 × 10⁻⁴ | **SUCCESS** |

BH 보정은 이 두 개만의 가족에서 수행했다. **1차 held-out 6개 가족은 닫힌 상태로 두고 재보정하지 않았다.**

gene-symbol 매핑률: 전체 행 61.3%, 발현 유전자 89.7% → 제외 기준 통과. 공유 유전자 2,219 / 2,211개.

## 2. 왜 이 결과가 강한가

### 2.1 증식 보정을 통과했다 — 이 데이터셋에서 가장 위험한 검정

`senescent vs proliferating`은 **정의상 증식 대비**다. 따라서 보정 없는 음의 상관만으로는 "senescence에 반대"와 "증식에 찬성"을 구분할 수 없다. 사전등록에서 증식 잔차화판을 co-primary로 지정한 이유다.

ρ가 −0.127에서 −0.101로, 크기의 **80%가 남았다.** 증식 성분은 이 연관의 소수 부분에 불과하다.

이것은 독립적으로 `direction_probe/P4`를 뒷받침한다. 거기서 decoupling index D는 증식 잔차화 후 0.461 → 0.468로 유지됐고, 증식 무관 유전자에서 0.546으로 오히려 커졌다. 서로 다른 데이터·다른 통계량이 같은 결론에 도달했다.

### 2.2 조직을 넘는다

MRC5는 **폐** 섬유아세포다. 로컬 시료와 discovery senescence 참조 4개는 대부분 **진피** 섬유아세포다. 사전등록 §5에서 "실패하면 주장 반증과 조직 특이성을 구분할 수 없다"고 미리 적어 두었는데, 성공했으므로 그 모호성은 발생하지 않았고 **조직을 넘는 일반화 근거**가 됐다.

### 2.3 senescence 유도 방식이 다르다

GSE306957은 **20 Gy X선 조사**로 유도한 senescence다. discovery 참조는 복제노화 3개(GSE179848·GSE109700·GSE191055)와 산화스트레스 SIPS 1개(GSE93535)였다. 즉 **새로운 유도 기전에서 재현**됐다.

### 2.4 효과크기가 discovery 범위와 일치한다

| 참조 | ρ | 단계 |
|---|---:|---|
| GSE93535 SIPS | −0.163 | discovery |
| GSE109700 deep RS | −0.134 | discovery |
| GSE179848 longitudinal RS | −0.100 | discovery |
| GSE191055 P27 vs P4 | −0.044 | discovery |
| **GSE306957 X선 유도 (held-out)** | **−0.127** | **held-out** |

held-out 값이 discovery 4개의 범위(−0.044 ~ −0.163) 한가운데에 있다. 축소(shrinkage)가 관찰되지 않았다.

### 2.5 RIG-I/MDA5 경로와 무관하다

| arm | ρ | 증식 보정 ρ |
|---|---:|---:|
| pLenti (CRISPR 대조) | −0.136 | −0.106 |
| siScramble (siRNA 대조) | −0.117 | −0.096 |
| IFIH1 KO (MDA5) | −0.119 | −0.093 |
| DDX58 KO (RIG-I) | −0.125 | −0.093 |

네 arm 전부, 보정 전후 모두 예측 방향이며 전부 p = 1 × 10⁻⁴ 수준이다. 이 연관은 해당 연구의 mtRNA–SASP 기전에 의존하지 않는다. **다만 이는 sensitivity 분석이며 기전 주장으로 쓰지 않는다.**

## 3. 해석 경계 (사전 고정된 대로 유지)

- 단일 senescence 모델(X선 유도) 결과이며, senescence 일반으로 확대하지 않는다.
- knockout arm 결과는 기전 증거가 아니라 sensitivity다.
- 이 검정은 **transcriptional alignment**이지 phenotype reversal이 아니다. 유전자는 생물학적 반복이 아니므로 여기의 p는 signature association 유의성이지 처치효과의 replicate p-value가 아니다.
- 로컬 anchor는 여전히 조건당 n=1에서 만들어졌다. 이 사실은 변하지 않는다.

## 4. 증거 등급 갱신

| 항목 | 이전 | 갱신 | 근거 |
|---|---|---|---|
| anti-senescence의 held-out 재현 | **미검정 (설계 결손)** | **A−** | 새 조직·새 유도기전에서 ρ=−0.127, 효과크기 축소 없음 |
| 증식 독립성 | A− (내부 민감도만) | **A** | held-out에서 co-primary로 사전 지정해 통과(−0.101) |
| stress–senescence decoupling | A− | **A−** (유지) | 양 축 모두 held-out 지지를 받았으나 D 자체를 held-out에서 재계산한 것은 아님 |

## 5. 두 held-out 가족 종합

| 가족 | 검정 수 | 성공 | 반증 | 종합 |
|---|---:|---:|---:|---|
| 1차 (GSE297233·GSE307377·GSE116968·GSE149694) | 6 | 2 | 0 | PARTIAL |
| 2차 (GSE306957) | 2 | 2 | 0 | **CONFIRMED** |

두 가족을 합치면 **8개 사전등록 검정 중 4개 성공, 0개 반증**이다. 성공한 4개는 핵심 주장의 두 축에 정확히 대응한다 — repair 정렬 2개(HO3a/b), anti-senescence 2개(HO5a/b). 미지지 1개(HO2, 연대기적 노화)는 주장문에서 **명시적으로 배제한** 항목이고, 나머지 3개(HO1, HO4a/b)는 null을 예측한 경계 검정이다.

**가족을 분리해 보정했다는 사실을 논문에 명시한다.** 합쳐서 8개로 재보정하면 1차 가족의 확정된 FDR이 사후에 바뀐다.
