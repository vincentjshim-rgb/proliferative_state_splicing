# 사전등록: secretome 개입은 재현 가능한 전사 class를 이루는가

작성 시점: 데이터 다운로드 **이전**. 이 문서를 SHA-256으로 잠근 뒤 내려받는다.

## 1. 배경과 동기

60개 교란의 연구 간 일치도를 계산한 결과(`public_data_tierA/derived/compendium_structure/`),
class마다 재현성이 극단적으로 달랐다.

| class | n | within-class 중앙값 r | 95% CI | 연구 간(different GEO series) |
|---|---:|---:|---|---:|
| senescence | 9 | +0.506 | 0.335–0.660 | +0.387 |
| donor age | 4 | +0.220 | 0.122–0.383 | — |
| photoprotection | 16 | +0.081 | 0.031–0.148 | — |
| **UV injury** | 9 | **+0.008** | −0.024–0.023 | **−0.007** |

UV 손상은 RNA-seq끼리만 봐도 −0.007이므로 플랫폼 artifact가 아니다.
노화 분야에서 secretome 기반 지표(proteomic aging clock 등)가 급증하고 있으나,
**secretome 개입이 수용 세포에서 재현 가능한 전사 상태를 만드는지는 검정된 바 없다.**

## 2. 데이터와 대비 — 결과를 보기 전에 고정

수용 세포가 정상 인간 섬유아세포이고, 처리가 무세포 secretome·conditioned medium·
extracellular vesicle 제제이며, 대응하는 대조 arm이 존재하는 연구만 포함한다.

| ID | dataset | 수용 세포 | 대비(처리 vs 대조) |
|---|---|---|---|
| SC1 | GSE139563 | IMR90 | 노화세포 CM day 4 vs day 0 |
| SC2 | GSE139563 | IMR90 | 노화세포 CM day 10 vs day 0 (민감도) |
| SC3 | GSE251807 | 정상 HDF | 골수 MSC-CM vs DMEM |
| SC4 | GSE266052 | HCA2 hTERT 섬유아세포 | cytokine-primed hMSC-CM vs 무처리 (둘 다 TGF-β 하) |
| SC5 | GSE266052 | HCA2 hTERT 섬유아세포 | naive hMSC-CM vs 무처리 (둘 다 TGF-β 하) |
| SC6 | GSE279804 | 복부 피부 섬유아세포 | ADSC nanovesicle vs 대조 |
| SC7 | GSE282054 | WI-38 | hTSC-CM vs ES-CM |
| SC8 | GSE293186 | 섬유아세포 | 내피세포 EV 72 h vs 대조 72 h |
| SC9 | GSE306748 | 초대 HDF | 불멸화 MSC secretome vs 대조 |
| SC10 | GSE306748 | 초대 HDF | 야생형 MSC secretome vs 대조 |
| SC11 | GSE268248 | **치은** 섬유아세포 | PRF serum vs 무처리 — 경계 대비, 진피가 아님을 명시 |

보너스(secretome class 아님, 기존 class 보강):
- GSE134533의 HDF UVB vs 무조사 → **UV injury class의 10번째 대비**
- GSE212873은 count matrix 없이 bigwig 7.35 GB뿐이므로 제외한다.

## 3. Primary endpoint — 정확히 3개

기존 분석과 동일한 통계량을 쓴다: 공통 유전자에서의 Spearman ρ, 쌍당 공유 유전자 1,000개 이상.
같은 연구에서 나온 대비 쌍은 연구 간 추정에서 제외한다.

### HS1 (주 검정) — secretome class의 연구 간 일치도
- 통계량: 서로 다른 GEO series에서 나온 secretome 대비 쌍의 **중앙값 Spearman r**
- **예측: |median r| < 0.15** (UV injury처럼 class를 이루지 못한다)
- **반증 기준: median r > 0.35** (senescence 수준으로 일관된다)
- 그 사이 값(0.15–0.35)은 "부분적 일관성"으로 보고하고 어느 쪽으로도 해석하지 않는다

### HS2 (양성 대조) — 파이프라인이 실제 secretome 효과를 잡는가
- GSE139563의 노화세포 CM(SC1)이 senescence 참조 축에 **r > +0.20**으로 적재될 것
- 실패 시 HS1의 음성 결과를 "효과를 못 잡은 것"과 구분할 수 없으므로 HS1 해석을 보류한다

### HS3 — Repro-CM의 위치
- Repro-CM과 **재생 지향 secretome**(SC3–SC10, MSC·hTSC·ADSC·내피 EV)의 중앙값 r이
  Repro-CM과 **노화 전파 CM**(SC1, SC2)의 중앙값 r보다 클 것
- 방향만 예측하며, 대비 수가 적으므로 유의성은 요구하지 않는다

## 4. 제외 규칙 (사전 고정)

1. 발현 유전자 중 gene symbol 매핑률 < 70%인 데이터셋은 제외
2. 대비를 구성했을 때 공유 유전자 < 1,000개면 해당 쌍을 결측 처리
3. 처리·대조 arm이 각각 n < 2이면 제외
4. 수용 세포가 섬유아세포가 아니거나(각질형성세포, MSC, 종양세포) 암 연관 섬유아세포(CAF)
   유도 설계이면 제외
5. 결과를 본 뒤 대비를 추가하거나 교체하지 않는다

## 5. 다중검정

HS1–HS3은 **다섯 번째 가족**이다. 기존 네 가족은 닫혔으므로 재보정하지 않는다.
HS1이 주 검정이고 HS2는 대조, HS3은 방향 예측이므로 BH 보정 대상은 HS1 하나다(p = FDR).

## 6. 미리 적어 두는 경계

- secretome의 **출처가 모두 다르다**(노화 섬유아세포, 골수 MSC, 지방 ADSC, 영양막 hTSC,
  내피세포, 혈소판). 일관성이 없더라도 그것이 "secretome이라는 범주가 무의미하다"는 뜻인지
  "출처가 다르면 다른 것을 만든다"는 뜻인지 이 설계로는 구분할 수 없다. 이 한계를 결과와 함께 쓴다.
- SC11(치은)은 진피가 아니므로 주 추정에서 빼고 민감도로만 쓴다.
- Repro-CM은 조건당 n=1이므로 HS3에 p-value를 붙이지 않는다.
