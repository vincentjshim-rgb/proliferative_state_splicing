# 사전등록: 배양 대 조직 경계 검정 (GTEx v8 via recount3)

- 작성일: 2026-09-13
- 상태: **GTEx 발현 자료를 내려받거나 열기 전에 작성·고정**
- 고정 방법: 이 파일의 SHA-256을 `public_data_tierA/logs/preregistration_sha256.txt`에 기록한 뒤 분석 스크립트를 실행한다.

## 1. 목적

원고의 현재 주장은 다음과 같다. 배양 인간 섬유아세포에서 splicing factor 발현과 8개 노화 프로그램(chromatin organisation, cellular senescence, pre-mRNA processing, mRNA splicing, interferon signalling, DNA double-strand break repair, telomere maintenance, base excision repair)은 분열 속도의 읽기값이며, collagen formation과 unfolded protein response는 그렇지 않다. 한계 절은 "분열 후 조직에 대해서는 말하지 않는다"고 적었다. 이 검정은 그 경계를 가정이 아니라 자료로 확인한다. 같은 공여자 집단에서 **배양 섬유아세포**, **피부 조직**, **골격근(분열 후 조직)**을 같은 정의로 비교할 수 있는 공개 자료는 GTEx가 유일하다.

## 2. 자료 (processed만)

| 파일 | 내용 | 크기 |
|---|---|---:|
| `gtex.gene_sums.SKIN.G026.gz` + metadata 3종 | Cells - Cultured fibroblasts, Skin - Sun Exposed (Lower leg), Skin - Not Sun Exposed (Suprapubic) | 약 196 MB |
| `gtex.gene_sums.MUSCLE.G026.gz` + metadata 3종 | Muscle - Skeletal | 약 81 MB |

합계 약 277 MB. FASTQ와 junction 파일은 받지 않는다. 따라서 GTEx에서는 splicing **outcome**(주석 없는 junction)은 계산하지 않으며, 프로그램 수준 발현 분석만 수행한다.

## 3. 결과를 보기 전에 고정한 정의

- **유전자 ID → symbol**: 프로젝트에서 이미 쓰는 매핑(`GSE282054_raw_counts.txt.gz`의 Geneid → Gene_name)을 그대로 사용.
- **정규화**: 조직 유형별로 따로. edgeR `filterByExpr` → TMM → log₂ CPM → 유전자별 *z*.
- **프로그램 점수**: `scripts/run_hallmark_audit.R`의 24개 Reactome 집합 목록 그대로. 해당 조직에서 발현된 유전자의 평균 *z*, 15개 미만이면 제외.
- **splicing 집합**: `conserved_core/age_down_splicing_core.txt`의 96개 유전자.
- **증식 점수**: `scripts/run_outcome_vs_expression.R`의 20개 marker(MKI67, CCNB1, CCNA2, CDK1, TOP2A, BIRC5, BUB1, PLK1, AURKA, TYMS, RRM2, PCNA, MCM2–MCM7, TK1, UBE2C) 평균 *z*.
- **순환성 방지**: 어떤 프로그램을 증식 점수와 비교하거나 증식으로 보정할 때, 그 프로그램에 포함된 증식 marker 20개는 프로그램 점수에서 먼저 제거한다. mitotic cell cycle 집합(양성 대조)은 보정 결과를 해석하지 않는다.
- **연령**: GTEx 공개 연령 구간(20–29, …, 70–79)의 중간값(25, 35, …, 75).
- **기술 공변량**: RIN, 허혈 시간, log₁₀ library size. recount3 metadata에 없는 필드는 쓰지 않고 그 사실을 기록한다.
- **생물 공변량**: 성별, 사망 양상 Hardy scale(결측은 별도 범주).
- **표본 제외**: 위 네 조직 유형이 아닌 표본; 연령 결측; 조직 유형별 library size 하위 1%. 한 공여자가 한 조직 유형에 표본이 여럿이면 library size가 가장 큰 1개.
- **상관의 정의**: 1차는 기술 공변량을 두 변수에서 모두 회귀로 제거한 뒤의 **partial Spearman ρ**(GTEx는 사후·수술 채취 자료로 기술 변동이 커서). 원시 Spearman ρ는 민감도 분석.
- **연령 모형**: `score ~ age + sex + Hardy + RIN + 허혈 시간 + log10 library size`, 증식 보정 모형은 여기에 증식 점수를 더한다. 감쇠율 = 1 − β_adj/β.

## 4. 1차 가설과 판정 기준

| | 가설 | 성공 | 실패 |
|---|---|---|---|
| **H1** 배양 재현 | GTEx 배양 섬유아세포에서 splicing 96 점수와 증식 점수의 partial ρ | ρ ≥ 0.5 | ρ < 0.3 (0.3–0.5는 부분 재현) |
| **H2** 조직에서 결합 약화 | lower-leg 피부 조직의 splicing 96–증식 ρ가 배양 섬유아세포보다 작다 | Fisher *z* 차이, 단측 *P* < 0.05 | 차이 없음 또는 반대 |
| **H3** 분열 후 조직 경계 | 골격근의 splicing 96–증식 partial ρ | ρ < 0.3 | ρ ≥ 0.5 |
| **H4** 배양의 연령 효과 감쇠 | 배양 섬유아세포에서 splicing 96의 연령 효과가 증식 보정 후 감쇠 | 감쇠 ≥ 50% | 감쇠 < 25%. 보정 전 연령 *P* ≥ 0.05면 **평가 불가**로 보고 |
| **H5** 배양의 collagen 해리 | 배양 섬유아세포에서 collagen formation–증식 partial ρ | \|ρ\| < 0.3 | \|ρ\| ≥ 0.5 |
| **H6** 조직 노화 양성 대조 | lower-leg 피부 조직에서 collagen formation이 연령에 따라 감소 | β < 0, *P* < 0.05 | β ≥ 0 또는 *P* ≥ 0.05 (이 경우 조직 연령 분석 전체를 검출력 부족으로 해석) |

H2 주석: 피부 표피에는 증식하는 basal keratinocyte가 있으므로, 조직에서도 결합이 유지되어 H2가 실패할 수 있다. 그 경우 결론은 "피부 조직은 분열 후 조직이 아니다"로 제한하고, 경계 판단은 H3에 맡긴다.

## 5. 탐색적 분석 (성공 기준 없음)

- **S1** 24개 프로그램의 증식 partial ρ를 배양·두 피부·골격근에서 나란히 비교 (Fig. 7A의 조직판).
- **S2** 일광 노출(lower leg) 대 비노출(suprapubic) 피부에서 collagen formation 연령 기울기 비교.
- **S3** 같은 공여자의 배양 섬유아세포와 피부 조직 사이, 프로그램 점수의 공여자 수준 일치도.
- **S4** 피부 조직 구성 민감도: keratinocyte marker(KRT5, KRT14, KRT1, KRT10)와 dermal fibroblast marker(COL1A1, COL1A2, DCN, LUM, PDGFRA) 점수를 공변량으로 추가.
- **S5** 24개 프로그램 × 조직 유형의 연령 효과 전체표, 증식 보정 전후.

## 6. 다중검정

- H1–H6: 성공 판정은 위의 효과크기 기준으로 한다. 각 가설의 *P* 값은 6개 1차 가설 안에서 Benjamini–Hochberg 보정 후 함께 보고한다.
- S1, S5: 조직 유형별로 증식 계열과 연령 계열을 각각 따로 BH 보정한다(기존 감사와 동일).

## 7. 보고 원칙

- 성공과 실패를 모두 보고한다. 실패한 가설을 원고에서 빼지 않는다.
- 결과를 본 뒤 정의, 임계값, 표본 제외 규칙을 바꾸지 않는다. 불가피하게 바꾸면 아래 수정 이력에 날짜와 이유를 적고, **원래 정의의 결과도 함께** 보고한다.
- 분석 대상은 재처리된 공개 자료이며 인과를 주장하지 않는다. GTEx 연령은 10년 구간이므로 연령 효과는 구간 해상도에서만 해석한다.

## 수정 이력

(없음)
