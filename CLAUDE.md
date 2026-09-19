# Claude 연구 인계 지침

이 파일은 이 프로젝트를 Claude에게 인계할 때 가장 먼저 읽는 지속 지침이다. 과거 작업 요약이 아니라, 이후 분석·해석·그림·원고 작성에서 지킬 연구 설계와 판단 기준을 기록한다.

- 최종 갱신: 2026-09-20 (원고 v13. **Fig. 7을 graphical abstract로 재작도(패널 b 삭제), 근육을 사전등록 결산 문단으로 재강등, 조직 패턴을 "배양에서 결합된 19개가 19개 모두 피부에서 약해진다 + 복제 4개만 71–75% 유지"로 확정**(§24). 그 앞은 2026-09-19 v12: 리뷰어 7명·검증자 7명 심사에서 확정된 52건 반영(§23). 그 앞은 v11. **개념 요약 Fig. 7 신설, 근육 강등, Results 절 사이 예측 문장으로 논리 사슬 보강**(§22). 그 앞은 2026-09-18 v10. **재심사 패널 4명 + 편집자 결정(Major revision)을 반영하기 시작했다**: 두 연구 역분산 결합 추정치 철회, GTEx 잔여를 suppression으로 명시, 피부 내부 양성 대조를 Fig. 6 패널 e로 승격, 대리 지표 품질을 Fig. 2 패널 f로 신설, 사전등록 스코어카드를 보충 표 9로 이동, 제목을 *weakens in skin*으로 확정(§19). 그 앞은 v9. **로컬 자료를 전부 내리고 공개 자료 재분석으로 확정**했고, 선행 연구는 고찰의 인용으로 잇는다(§4·§18). 그 앞은 v8의 **본문 그림 8 → 6장 일직선 재편**: 관찰 → 측정된 분열 속도 → 증식이 설명하는 몫 → 63개 대비 → 23개 프로그램 → 사전등록 조직 경계. 곁가지 분석은 전문을 보충 Results 1–8로 옮기고 본문에는 결과를 밝힌 짧은 포인터만 남겼다. 제목 변경. 재편 과정에서 나온 오류 다수 수정(§8 끝의 '2026-09-17 수정 기록'))
- 직전 상태(v7, 그림 8장): `archive/2026-09-17_before_six_figure_restructure/`
- 직전 판(v5, 2026-09-14): `archive/2026-09-15_before_cohort_redefinition/CLAUDE_2026-09-15_v5.md`. 그림과 원고 사본도 같은 디렉터리에 있다.
- 그 이전 판(Repro-CM stress–senescence decoupling 방향, 2026-09-11): `archive/CLAUDE_2026-09-11_repro_cm_decoupling.md`.

## 1. 프로젝트의 목적과 범위

선행 논문 `JTE-Aug-26-0256_Proof_hi.pdf`의 후속 in silico 연구다. 선행 연구는 부분 리프로그래밍 단계에서 얻은 conditioned medium(secretome)이 UVA 손상 인간 진피 섬유아세포의 표현형을 개선하고 *C. elegans* 수명을 연장했음을 보였다.

**연구 방향의 변화.** 2026-09-11까지는 Repro-CM decoupling 가설이었고, 2026-09-12–13에 "배양 섬유아세포의 노화 signature는 증식 상태의 읽기값"으로 바꿨다. 2026-09-15 제출 전 모의 심사 뒤 **방향 (a)** 를 택했다(사용자 결정): splicing outcome(주석 없는 junction 비율)과 *COL1A2* event를 핵심 주장에서 내리고, 원고는 "배양에서 splicing factor의 양은 증식 상태의 읽기값이며 이 결합은 조직으로 옮겨지지 않는다" 하나만 주장한다. outcome은 음성 결과로 사실대로 보고한다(Fig. 4, 보충 Fig. S1).

다음 원칙은 고정한다.

- AlphaGenome은 사용하지 않는다.
- 공간 오믹스는 사용하지 않는다.
- 신규 wet-lab 검증을 전제로 스토리를 만들지 않는다. 투고 시점과 목표 모두 wet-lab 자료 없이 진행한다.
- 공개 데이터가 뒷받침하는 범위 안에서만 aging, reprogramming, longevity를 논한다.
- 분석 결과가 약하거나 반대 방향이면 숨기지 않는다.
- 사용자가 다운로드를 명시적으로 승인하기 전에는 새 공개 데이터를 다운로드하지 않는다.

## 2. 현재 핵심 주장 (GeroScience 원고 v8, 2026-09-17)

- 목표 학술지: **GeroScience (Original Article)**. npj Aging은 사전문의 단계에서 보류
- 영문 제목: *Splicing-factor expression in cultured human fibroblasts reports proliferative state and the coupling weakens in skin*
- 한글 제목: 배양 인간 섬유아세포의 splicing factor 발현은 증식 상태를 읽으며, 이 결합은 피부에서 약해진다
- 이전 제목(*Ageing signatures and transcriptomic age predictors…*)은 2026-09-17 독립 심사 패널 5명 전원이 부적합 판정: 나이 예측기는 사후 분석 패널 하나이고, 정작 주제인 splicing factor가 제목에 없었다. 제목에 "not in skin"이라고 쓰지 않는다(피부는 0.23으로 유의). "rather than age"도 쓰지 않는다(잔여 감소가 남는다)

한 문장 주장:

> In cultured human fibroblasts the abundance of splicing-factor transcripts is a readout of proliferative state; in skin from the same donor pool the coupling weakens while the same score still recovers the replication programmes. (Skeletal muscle cannot test it and is reported as a preregistered outcome only, §22.)

> 배양 인간 섬유아세포에서 splicing factor 전사체의 양은 증식 상태의 읽기값이며, 같은 공여자 집단의 피부에서는 같은 점수가 복제 프로그램은 여전히 회수하는데도 이 결합은 약해진다. (골격근으로는 검정할 수 없어 사전등록 결과로만 보고한다, §22.)

다음과 같이 더 강한 문장으로 바꾸지 않는다.

- "splicing factor는 노화와 무관하다" — 증식 보정 후에도 잔여 감소가 남는다(성인 107명 −0.040/10년 *P* = 0.005, GTEx 배양 −0.053/10년 *P* = 0.002, 신뢰구간 겹침).
- "노화에서 splicing은 중요하지 않다" — 조직 outcome 근거(Zhang et al. bioRxiv 2026)는 별개의 측정이며 이 연구는 그것을 지지도 반박도 하지 않는다.
- **"splicing outcome이 나이에 따라 나빠진다" 또는 "*COL1A2*가 나이에 따라 점진적으로 재분배된다"** — 방향 (a)에서 내린 주장이다. 성인 코호트에서 outcome의 연령 효과는 없다(*P* = 0.51).
- "배양 결과가 조직에서도 성립한다" — GTEx에서 피부는 약화(0.23, 구성 보정 후 0.10), 근육은 소실(−0.12).
- "collagen 연령 신호" 또는 "proteostasis 프로그램은 분열과 무관하다" — GTEx에서 재현되지 않았고, 재정의 코호트에서도 collagen은 공변량 보정 후 비유의(FDR 0.25).
- GTEx 조직의 연령 효과를 근거로 한 모든 주장 — 사전등록 양성 대조(H6)가 실패했다.
- Repro-CM에 대해: "광노화를 역전시킨다", "부분 리프로그래밍을 직접 모사한다", "전신 rejuvenation 또는 longevity program을 활성화한다", "특정 miRNA나 TIMP2가 효과를 유발한다".
- "문헌이 거명하는 splicing factor"를 출처 없이 나열하지 않는다(§3의 목록 규칙).
- **"결합은 배양에 속한다" / "세포의 성질이 아니라 배양의 성질"** — bulk 조직은 대부분 분열하지 않는 세포라 세포 안에서 결합이 같아도 공여자 간 비교에서는 사라질 수 있다. 조직 결과는 "이 측정이 무엇을 보고하는가"에 관한 진술로만 쓴다.
- **"세어서 얻은 분열 속도의 읽기값"** — 처리 섭동에서 splicing 점수는 세어서 얻은 속도보다 전사체 증식 프로그램을 더 가깝게 따라간다. 주장은 "증식 상태(proliferative state)"로 쓴다.
- **methylation 시계에 대한 어떤 주장도** — 사전등록은 시계가 같은 변수를 따라갈 것이라고 *예측*했고 그 예측이 실패했다. 이를 "주장의 경계를 세우는 대조"로 포장하지 않는다. 본문은 "예측했고 충족되지 않았다, 시계에 대해서는 진술하지 않는다"까지다.

## 3. 연구 재료의 정확한 생물학적 의미

**로컬 mRNA-seq은 현재 원고(v9)에서 쓰지 않는다.** 2026-09-17에 전부 내렸다(§18). 원고는 공개 자료만 분석하며, 선행 연구는 고찰에서 `(Shim, V. *et al.*, manuscript submitted)`로 인용한다. 아래 §4·§5의 로컬 자료 규칙은 다시 쓸 일이 생길 때를 위해 보존한다.

**공개 자료의 역할**

| 자료 | 현재 역할 |
|---|---|
| GSE179848 (Cellular Lifespan Study) | **측정된 분열 속도**(계대마다 세포 계수), 328개 라이브러리. 배양 시간 종단 자료(HC1–HC4) |
| GSE179847 (같은 연구의 DNA methylation) | EPIC array 479개. GSE179848과 **같은 세포주·같은 계대**라 시료 수준 짝짓기가 된다. 역할은 **경계 대조 하나**: 전사체 읽기값이 분열 속도를 보고한다는 주장이 노화 바이오마커 전반으로 확대되지 않음을 보인다. 사전등록 고정 뒤 내려받았다(§10). 절대 DNAm age는 근사값이므로 순위·잔차만 해석한다 |
| GSE113957 (Fleischer 2018) | 공여자 연령. **주 코호트는 정상 성인 107명(20–96세)**: 기탁 143개에서 HGPS 10개와 20세 미만 26개를 뺀 것(`scripts/revision/cohort_gse113957.R`). recount3 SRP144355 gene·junction sums를 쓰고, 세포주 출처(AG/기타)·시퀀서·성별·depth를 공변량에 넣는다. 민감도: 정상 133명, 성인 20–82세 76명, 기탁 142개 |
| GSE226189 | **연구에서 제외**(2026-09-14). 처리 batch가 연령 층과 겹친다(batch별 연령 Kruskal–Wallis *P* < 0.001). 남은 흔적은 GTEx 사전등록이 고정한 96개 판 splicing 집합의 정의뿐이다(Fig. 8, §9). 제외 근거를 점검하는 `run_core_depth2.R`만 이 자료를 읽는다 |
| 23개 GEO series, 63개 대비 | senescence, 공여자 연령, UV 손상, 광보호, reprogramming, 대사·배양, secretome. 대비 수준 분석만. 보충 표 1은 28행이다(종류가 다른 대비를 가진 GSE116968·GSE179848·GSE222414·GSE240486은 여러 행). 이전 원고의 "29개 연구"는 표의 행 수를 센 오류였고 당시 실제 series는 24개였다 |
| GTEx v8 (recount3) | 배양 섬유아세포, 일광 노출 하지 피부, 비노출 치골상부 피부, 골격근. 사전등록한 배양 대 조직 경계 검정 |

**용어 구분.** "측정된 분열 속도"는 GSE179848에만 쓴다. 다른 모든 자료의 증식 측정은 20개 marker(MKI67, CCNB1, CCNA2, CDK1, TOP2A, BIRC5, BUB1, PLK1, AURKA, TYMS, RRM2, PCNA, MCM2–MCM7, TK1, UBE2C)의 "증식 점수"다. 두 측정은 모든 프로그램에서 일치하지 않는다(§7 Fig. 8).

## 4. 로컬 원본 데이터와 사용 규칙 (현재 원고에서는 쓰지 않음, 2026-09-17)

외부 원본 프로젝트는 읽기 전용으로 취급한다: `/home/shim/Downloads/project/project_shim/`

주요 RNA-seq 원본은 `/home/shim/Downloads/project/project_shim/RNA_seq/00.data/`에 있다. iPSC-CM R1은 반드시 `/home/shim/Downloads/project/project_shim/RNA_seq/00.data/iPSC_CM_24h.R1.fastq.gz`를 쓴다.

- 검증된 SHA256: `7040ea05c002460f22ef884071a6c9faa7359fb1cd95fd833aef2eb8c10e26b6`
- workspace에 있던 SHA256 `74351c62...` 사본은 손상되어 사용 금지

검증된 trimmed reads는 `/home/shim/Downloads/project/project_shim/RNA_seq/01.trim/`에 있다(HDF, Repro, iPSC의 paired `*_val_1.fq.gz`, `*_val_2.fq.gz`; raw/trimmed 12개 파일 gzip 무결성 통과).

STAR gene counts: `public_data_tierA/derived/local_repro_anchor/star_genecounts/{HDF,REP,IPS}/` — GENCODE v41 GRCh38 STAR index, reverse-stranded column 4가 주 분석값. 과거 RSEM 결과(실제 reads의 2–3%만 계수)는 쓰지 않는다.

| 조건 | input reads | uniquely mapped | reverse-strand assignment | mismatch rate |
|---|---:|---:|---:|---:|
| HDF-CM | 16,676,510 | 94.48% | 85.69% | 0.26% |
| Repro-CM | 15,342,652 | 93.67% | 84.62% | 0.26% |
| iPSC-CM | 16,843,629 | 94.66% | 86.21% | 0.26% |

**read 수준 coverage용 네 라이브러리**(UVA 15 J cm⁻² 0 h, 그리고 같은 세포의 HDF-CM·Repro-CM·iPSC-CM 24 h): `scripts/align_local_four_samples.sh`, 출력 `public_data_tierA/derived/local_bam_v41/`. STAR 2.7.3a. samtools 1.3.1은 이 BAM을 색인하지 못해 qiime2 환경의 samtools를 쓴다(원고 Methods에 1.21로 기재).

가장 중요한 제한은 조건당 생물학적 반복이 `n=1`이라는 점이다. 다음을 금지한다.

- 로컬 조건 간 group-level DEG p-value 또는 FDR 계산
- 로컬 샘플을 이용한 통계적 유의성 주장
- 유전자를 biological replicate처럼 취급하는 통계

로컬 자료는 effect-size 기반 기술(descriptive) 용도로만 쓴다. Fig. 10은 측정을 기록할 뿐 처리 효과를 검정하지 않는다.

## 5. 고정된 로컬 anchor 정의 (현재 원고에서는 쓰지 않음, 2026-09-17)

1. STAR reverse-stranded gene counts에 TMM normalization
2. 적어도 2개 샘플에서 CPM ≥ 1인 유전자만
3. pseudocount 0.5

```text
dH = log2(CPM_Repro + 0.5) - log2(CPM_HDF + 0.5)
dI = log2(CPM_Repro + 0.5) - log2(CPM_iPSC + 0.5)
if sign(dH) == sign(dI) and both are non-zero:
    S = sign(dH) * min(abs(dH), abs(dI))
else:
    S = 0
```

필터 통과 12,319 genes; 방향 일치 8,921 (72.4%); `|S| ≥ log2(1.25)` 2,204; `|S| ≥ log2(1.5)` 580; `|S| ≥ 1` 60. reverse-stranded 대 unstranded 점수 상관 ρ = 0.9869.

`scripts/run_local_repro_anchor.R`의 `Repro_specific_score`가 S다. DEG 통계량이 아니라 기술적 effect-size anchor다. 대비 표에 `LOCAL_ReproCM`을 다시 넣으려면 `run_stats_supplements.R`의 `INCLUDE_LOCAL`을 TRUE로 바꾸고 Fig. 4·S6·보충 표 4를 다시 만들어야 한다. 그때는 §18의 3,118개 유전자 표 문제부터 고쳐야 한다.

## 6. 공개 데이터 현황과 출처 기록

- `public_data_tierA/` 원자료는 derived를 빼고 약 1.8 GB다. 그중 `logs/` 473 MB의 대부분은 격리된 불완전 파일이다.
- **격리 파일** `public_data_tierA/logs/rejected/GPL16956_family.soft.gz.incomplete`(약 472.6 MiB)는 분석에 쓰거나 임의 삭제하지 않는다.
- 주요 디렉터리: `aging/`, `senescence/`, `uva/`, `longevity/`, `network/`(Reactome GMT), `secretome/`, `heldout/`, `heldout2/`, `heldout3/`, `recount3/`, `gtex/`.
- **출처 기록.** 2026-09-13에 `heldout/`(GSE297233, GSE307377, GSE116968, GSE149694), `heldout2/`(GSE306957), `heldout3/`(GSE240486, GSE222414, GSE329475, GSE326951)의 파일 24개를 중앙 목록 `download_manifest.tsv`(71행)에 추가했다. GEO 파일 크기와 일치한 것만 넣었고, 추가 전 사본은 `logs/download_manifest_before_heldout_20260913.tsv`다. `secretome/` 12개 파일은 아직 `secretome/logs/sha256.txt`에만 기록되어 있다. 모든 파일은 각 디렉터리 sha256 대조를 통과했다. 새로 받는 파일은 중앙 목록과 디렉터리 sha256 양쪽에 기록한다.
- **`heldout*/` 자료의 현재 역할.** 이전 방향에서는 결과를 보기 전에 예측을 고정한 held-out 검증 자료였다. 현재 원고에서는 그 역할이 없고, 63개 대비 중 31개(23개 GEO series 중 9개)를 공급하는 일반 대비다: reprogramming 7/7, 광보호 13/13, senescence 4/9, UV 손상 5/6, 공여자 연령 1/2. 현재 원고의 검증은 GTEx 사전등록 검정(Fig. 8)이 맡는다.
- 이전 판의 "다운로드 전 대기 목록" 다섯 데이터셋은 모두 받았다.
- GTEx: `gtex/`에 SKIN·MUSCLE gene sums와 metadata 8개 파일(265 MB), `gtex_sha256.txt`, 중앙 목록에도 기록.

자료별 주의사항:

- GSE93535: fold change의 39%가 ±4.056에 포화(Cufflinks/FPKM 재구성)되어 GSEA에서 제외했고, 시료별 GTF가 `CUFF.*` ID라 시료 수준 점수화가 불가능하다.
- GSE165177: GEO의 PubMed ID `5390271`은 오기이며 실제는 35390271(Gill et al., eLife 2022)이다.
- GSE113957 연령은 열 이름에 두 가지 형식이 섞여 있다(`101_19yr_...`, `108_31_female_...`). 정규식 `^[0-9]+_([0-9]+)(yr|YR)?_.*$`로 성인 97명을 복원한다.
- recount3 gene sums는 base coverage 합이다. GSE113957과 GTEx 모두 이를 count로 쓰는 같은 규칙을 적용했다.
- 826,224 × 143 junction 행렬은 dense로 만들면 멈춘다. `Matrix::sparseMatrix`로 읽는다.

## 7. 현재 핵심 결과 (원고 v8 수치, 2026-09-17)

본문은 여섯 절, 그림 여섯 장이다. 곁가지 분석은 보충 Results(SR) 1–8에 전문이 있고 본문에는 결과를 밝힌 포인터만 있다.

| Fig. | 결과 |
|---|---|
| 1 | 세 자료에서 모두 감소한 유전자의 1위 과정은 mitotic cell cycle(283/470, BH-FDR 3 × 10⁻⁶⁶), 다음이 capped intron-containing pre-mRNA processing(177/280, 8 × 10⁻⁴⁵). 이 177개가 splicing 집합(v2). 성인 107명 연령 *r* = −0.50(*P* = 4 × 10⁻⁸). 배양 시간 HC1 ρ = −0.89, HC2 −0.83, HC3 −0.96(정확 *P* = 0.003), HC4 +0.77(정확 *P* = 0.10) |
| 2 | 측정된 분열 속도와 ρ = 0.71 [0.65, 0.76](328개). 혼합모형 β = 0.69(*P* = 3 × 10⁻⁴⁶), 세포주별 ρ 0.32–0.88(7개 중 6개 유의), 41개 계열 ρ = 0.54. 거명 인자 20개(Holly 2013) 중앙값 0.56. **처리 섭동 포인터(SR1, 보충 Fig. S1, 표 2; 사후)**: 177개 점수로 Δsplicing 대 Δ분열속도 ρ = 0.49, 대 Δ증식점수 **ρ = 0.79**; contact inhibition 속도 −0.35/day·점수 −0.73; 측정 속도 고정 시 9개 중 **5개**가 잔여(DEX +0.53, oligomycin +0.56, mitoNUITs+DEX +0.45, oligo+DEX +0.41, contact inhibition −0.17), 증식 점수 고정 시 9개 중 7개가 잔여 없음. **어느 측정으로도 설명 안 되는 잔여는 둘**: oligomycin(속도 고정 +0.56, 증식 고정 +0.28)과 3% 산소(속도 고정 −0.21, 증식 고정 −0.33). **2-DG는 셋째가 아니다** — 세어서 얻은 속도를 고정하면 +0.10 [−0.04, 0.24]로 설명되고, −0.32는 2-DG가 splicing은 그대로 둔 채 증식 점수만 올려서 생기는 값이다. '잔여를 남긴다'의 기준은 BH-FDR < 0.05로 통일한다(95% CI 기준으로는 galactose가 여섯째가 된다) |
| 3 (신규, 옛 3·4 병합) | 증식과 *r* = 0.89. 연령 효과 −0.133 → −0.040/10년(70% 제거), 네 코호트 정의에서 63–73%. depth만 공변량이면 보정 후 −0.022(*P* = 0.11). QC 공변량 추가해도 −0.110 → −0.039. **무작위 집합 1,000개**: 205개가 자체 연령 연관(대부분 증가), 보정이 중앙값 56% 제거, 70% 이상 감쇠는 27%; 크기(|β| ≥ 0.133)에 도달한 집합 0개. 유전자 137/177 연령 연관 → 19개 유지, 거명 인자 11개 중 0개. 코호트별 40/56/19/0, 네 정의 공통 0(넷째가 0개라 강제된 결과; 나머지 셋 공통 10개). **포인터**: outcome 지표 β = +1.5 × 10⁻⁵ [−3.0, +6.0 × 10⁻⁵](SR3, 보충 Fig. S3·S4), 나이 예측기 *r* = 0.72·증식과 −0.74·증식점수 단독 0.46(SR4, 보충 Fig. S5; 사후) |
| 4 (옛 5) | **63개 대비** ρ = 0.86 [0.78, 0.92], *R*² = 0.92(공유 23개 제거 0.85/0.91). 무작위 집합 최대 **0.56**. senescence 제외 0.79/0.66(**n = 54**), 연구당 평균 0.79/0.94(**n = 29; secretome 제제를 각각 따로 센 단위**). 종류별 잔차 모두 비유의. **개입 포인터(SR5, 보충 Fig. S6; 사후)**: 예측 대상 종류를 뺀 적합에서 reprogramming 평균 잔차 +0.02 [−0.14, +0.18](대비 7개·series 2개, 검정력 거의 없음), secretome **−0.02** [−0.07, +0.02]; 개별 조건은 양쪽으로 벗어남(naive +0.34, OSK −0.18). **특이성 포인터(SR6, 보충 Fig. S7, 표 5; 사후)**: 세어서 얻은 분열 속도에 대해서는 선택되지 않은 Reactome 정의에서도 비대칭 유지(66–91%), **GTEx 배양의 전사체 증식 점수에 대해서는 재현되지 않음**(translation이 더 가깝다) |
| 5 (옛 6) | hallmark 23개 중 8개가 |ρ| > 0.5. cell-cycle 유전자 제거 후 pre-mRNA 0.67, mRNA splicing 0.64, chromatin 0.63, DSB repair 0.55만 유지. 보정 후 연령 연관이 남는 것은 splicing 두 프로그램(FDR_adj 0.04)이지만 **이를 특이성의 근거로 쓰지 않는다**: rRNA processing(−0.038), translation(−0.032), chromatin organisation(−0.015)도 보정 후 명목 *P* = 0.04–0.05이고 다중검정 보정에서만 떨어진다(FDR_adj 0.25). UPR의 0.07은 *보정 전* FDR이므로 '근소하게 놓쳤다'는 식으로 쓰지 않는다(보정 후 FDR 0.35). senescence 집합 +0.67(150개 중 82개 세포주기, 제거 후 +0.38), SenMayo −0.10. **패널 포인터(SR7, 보충 Fig. S8, 표 7; 사후)**: 7개 패널에 걸쳐 ρ = 0.71이지만 **정확 *P* = 0.088, n = 7**; CellAge 억제 패널은 예외(21%인데 +0.69); **연령 연관이 있던 5개 패널 중** 증식 보정 후 유지되는 것 0개(66–100% 손실). **methylation 포인터(보충 Note 2, Fig. S9, 표 8)**: M1 충족(ρ = 0.36, n = 86), M2 ρ = −0.01 [−0.10, 0.08], M3 ρ = 0.16 [−0.03, 0.34], M4 불충족. 시계에 대해 진술하지 않는다 |
| 6 (옛 7, 패널 d 제외) | 사전등록 GTEx(96개 판): 배양 492명 partial ρ = 0.74, 하지 피부 0.23, 치골상부 0.06, 근육 −0.12. 근육에서 양의 결합 최대는 **base excision repair 0.30, proteasome은 −0.37**(v7 본문의 '0.37, proteasome'은 부호 오류였다). 가설 6개 중 4개 지지, H6 불지지, H4 평가 불가. 사후: 기울기 0.55/0.17/0.07/−0.13, 분산 맞춘 부분표본 0.77–0.87, 구성 보정 후 피부 0.10·근육 −0.04·배양 0.72. 공여자 일치 0/25 + 양성 대조 0.75(SR8, 보충 Fig. S10c). **잔여 연령 효과(사전등록 아님)**: GSE113957 −0.128 → −0.033(ρ(age,prolif) = −0.46), GTEx −0.025 → −0.053(**ρ = +0.10, 모형에 쓴 492명 기준**), 결합 −0.042 [−0.065, −0.019](599명, 사후; *I*²는 연구 2개라 적지 않는다), Reactome 프로그램 −0.041 |

보충 그림: S1 섭동, S2 유전자 수준(옛 Fig. 3a–d; 패널 d 단위는 **10년당 log₂ CPM** — 옛 그림의 'z per decade'는 단위·크기 오류였다), S3 outcome 지표(옛 Fig. 4a·c·e), S4 event 수준(옛 S1), S5 나이 예측기(옛 Fig. 4f + 모든 팔), S6 개입(옛 Fig. 8, 표본 외 예측), S7 특이성(선택·비선택 정의), S8 senescence 패널, S9 methylation, S10 GTEx 사후 + 공여자 일치(옛 Fig. 7d).

보충 표(새 본문의 첫 인용 순): 1 코호트 민감도, 2 섭동(+산소 대비), 3 공여자 민감도 분석, 4 자료와 대비(+개입 표본 외 예측), 5 특이성과 결합 추정, 6 프로그램 감사, 7 senescence 패널, 8 methylation 검정, 9 GTEx 프로그램. 옛 번호 → 새 번호: 1→1, 7→2, 5→3, 3→4, 8→5, 4→6, 9→7, 6→8, 2→9.

## 8. 음성·불재현 결과 (삭제 금지)

현재 방향:

- **H6 실패**: 하지 피부 collagen formation은 연령에 따라 증가(+0.027/10년, *P* = 0.018).
- **collagen 연령 신호 불재현**: GSE113957 배양의 ρ = −0.44가 GTEx 배양에서는 +0.017/10년(FDR 0.33).
- **측정 분열 속도와 증식 점수의 불일치**: GTEx 배양에서 UPR 0.70(Lifespan 0.05), TCA 0.74(0.29), translation 0.75(0.43), chromatin organisation 0.43(0.69).
- **H4 평가 불가**: GTEx 배양에 보정 전 연령 효과 없음(*P* = 0.33).
- **배양 이탈 표본 32개**(GTEx, 사후): 특정 추출·시퀀싱 batch에 몰려 있고 RIN이 낮다. 제외하면 원시 ρ 0.72 → 0.89. 사전등록 추정치는 이들을 포함한다.
- GSE226189 제외. HC4는 감소하지 않음(ρ = +0.77, *P* = 0.07, 55일까지만 추적).
- 증식 보정 후 연관이 양의 방향으로 나타나는 6개 유전자(*NDC1*, *NUP43*, *NUP50*, *NUP155*, *SF1*, *SNRNP48*). 보정 후 연관 유전자의 결합도가 나머지보다 낮다는 차이는 v2에서 비유의(*P* = 0.12).
- **GSE113957 코호트 문제 — 2026-09-15에 해결하고 원고에 반영했다.** 기탁된 143개 라이브러리에는 HGPS 10명과 20세 미만 26명이 들어 있고(분석 표에 들어간 142개 집합에는 HGPS 9명 — 한 라이브러리에 junction 지표가 없다), 83세 이상 31명은 전원 Coriell AG 계열이다. 주 코호트를 정상 성인 107명으로 바꾸고 세포주 출처·시퀀서·성별을 공변량에 넣었다. 결과:
  - 기계 발현의 증식 보정 감소는 모든 코호트 정의에서 유지된다(63–73%). 성인 107명: −0.133 → −0.040/10년(*P* = 6 × 10⁻⁷ → 0.005).
  - **splicing outcome의 연령 효과는 사라진다**(성인 *P* = 0.51; 20–82세는 반대 방향 *P* = 0.012). 발표판의 *P* = 0.002는 코호트 구성이 만든 것이다.
  - **보정 후 감소 유지 유전자는 코호트마다 40/56/19/0개이고 공통 0개**다.
  - *COL1A2* 원위 junction(chr7:94,406,304–94,418,498)은 GENCODE v41 미주석이다(로컬 정렬에서 고유 read 0개·다중정렬 2–3개였다는 근거는 v9에서 원고에서 내렸다). PSI는 83세 이상에서만 계단형으로 올라간다(0.035–0.049 → 0.072).
  - 상세: `review_presubmission_npjaging_ko.md`, 스크립트 `scripts/presubmission_checks/`, `scripts/revision/`.
- **증식 보정의 감쇠는 splicing에 특이적이지 않다 (2026-09-16 음성 대조).** 발현량을 맞춘 무작위 유전자 집합 1,000개 중 205개가 연령 연관을 가졌고, 같은 보정이 그 중앙값 56%를 제거했다. 70% 감쇠는 특별하지 않다. 특이한 것은 연관의 크기다(보정 전 *P* ≤ 6 × 10⁻⁷, 보정 후 *P* ≤ 0.005, |β| ≥ 0.133에 도달한 무작위 집합 0개). 원고 결과·고찰에 그대로 적었다.
- **SenMayo로 senescence 결과를 재확인 (2026-09-16).** 측정 88개 중 세포주기 유전자 1개뿐이고 분열 속도와 ρ = −0.10(*P* = 0.07). Reactome Cellular Senescence의 +0.67은 구성(150개 중 82개가 세포주기)이 만든 것이다. 목록은 `public_data_tierA/senescence/senmayo/`(중앙 목록 73행, sha256 기록).
- **senescence 패널 일곱 개 중 증식 보정 후 연령 연관을 유지하는 것은 0개다 (2026-09-16, 보충 Fig. S7).** 배양 섬유아세포에서 어떤 senescence 패널도 증식과 독립된 연령 신호를 주지 않는다. 이것은 음성 결과로 보고하며, 특정 패널을 노화 지표로 권하지 않는다.
- **outcome 지표는 라이브러리 품질과 함께 움직인다 (2026-09-16).** 고유 정렬 비율 ρ = −0.43, 다중정렬 +0.42, intron 신호 +0.31. QC 공변량을 넣으면 남아 있던 연령 연관도 사라진다(*P* = 0.96). 기계 발현은 QC 보정 후에도 남는다(−0.110/10년, 보정 후 −0.039).
- **거명 splicing factor 목록의 출처 부재(2026-09-15).** 이전 그림의 28개·25개 목록 중 상당수가 인용 문헌 어디에서도 사람 연령 연관으로 거명되지 않았다(*U2AF1*, *PRPF8*, *SNRPA/B/D2/F*, *SF3B5*, *RBM39*, *PRPF19/31*, *DDX39B*, *CELF1*, *MBNL1*; *SF1*은 *C. elegans*만; *ESRP1*·*PTBP1*은 오히려 증가 보고). Holly 2013 표 3의 20개로 교체했다.
- **collagen formation의 연령 신호는 공변량 보정 후 사라진다**(성인 107명 ρ = −0.44이지만 FDR 0.25). GTEx 배양에서 재현되지 않은 것과 같은 방향이다.
- **cell-cycle 유전자를 빼면 결합이 유지되는 것은 8개 중 4개뿐**이다(senescence 0.67 → 0.38, interferon 0.60 → 0.49, BER 0.56 → 0.35; telomere는 판정 불가). splicing 두 프로그램은 유지되거나 상승한다.
- **H4 평가 불가, H6 불지지**는 그대로다.

- **2026-09-17 수정 기록 (재편 중 독립 검증이 잡은 오류).** 다음은 v7까지 원고·산출물에 있던 오류이고 v8에서 고쳤다. 같은 실수를 반복하지 않도록 남긴다.
  - **처리 섭동 분석이 177개가 아니라 96개 점수로 돌았다.** `hallmark_audit/audit_inputs.rds`의 `splice` 열은 177개 집합이 생기기 전(2026-09-13)에 쓴 96개 점수다. `run_treatment_perturbation.R`는 이제 `division_rate/sample_division_rate.tsv`의 Fig. 2 점수를 읽고 ρ가 0.71이 아니면 멈춘다. **`audit_inputs.rds`의 `splice` 열을 쓰지 않는다.**
  - **3% 산소 문장('속도 −0.02/day, splicing −0.40')은 대조만의 짝 비교가 아니었다.** 3% 라이브러리 35개(그중 20개 contact inhibition) 대 21% 229개(그중 166개 처리)를 배양 시간 맞춤 없이 합친 값이었다. 대조 배양만 짝지으면 속도 −0.10/day, splicing −0.51이고 증식 점수를 고정해도 −0.33이 남는다(증식으로 설명되지 않는 잔여).
  - **특이성 분석의 GTEx 수치(0.70 대 0.54)는 증식 점수가 아니라 Reactome Cell Cycle 프로그램 점수에 대해 계산된 것이었다.** 20개 marker 증식 점수로는 방향이 뒤집힌다(0.29 대 0.58). 또한 177개 집합은 이 자료에서 *선택된* 집합이므로 선택되지 않은 Reactome 정의로 함께 돌려야 한다.
  - **senescence 패널 문장 '일곱 모두 49–100% 손실'은 틀렸다.** SenMayo와 CellAge 유도 패널은 처음부터 연령 연관이 없었고 CellAge 유도의 손실률은 −107%다. 맞는 문장은 '연령 연관이 있던 다섯 중 유지되는 것 없음'. 일곱 패널 ρ = 0.71은 정확 *P* = 0.088이다.
  - **개입의 '전형적 효과의 3분의 2'는 표본 내 적합값이었다.** 종류를 뺀 적합으로는 0.59(secretome 0.50, reprogramming 0.13). 'reprogramming 배지가 가장 밀접하다'는 series 2개·대비 7개이고 naive 배지에 의존한다.
  - **근육의 '최대 결합 0.37, proteasome'은 부호 오류**(실제 −0.37; 양의 최대는 base excision repair 0.30).
  - **GTEx 증식–연령 ρ는 +0.11(511명)이 아니라 +0.10(모형에 쓴 492명).**
  - **옛 Fig. 3d는 유전자 수준 계수에 10을 곱하고 'z per decade'라고 적었다**(실제는 10년당 log₂ CPM).
  - **Fig. 4 범례에 패널 f가 없었고 본문 패널 호출이 어긋나 있었다.** '발표판 분석', 'published list' 같은 표현은 독자가 볼 수 없는 이전 초고를 가리키므로 'deposited 142-sample set'으로 바꿨다.
  - **보충 Note 1 항목 6과 `results_gtex_posthoc_ko.md`가 첫 사후 실행값(0.111/−0.052/−0.021)에 머물러 있었다.** 같은 날 섬유아세포 marker에 PDGFRA를 되돌린 뒤의 현재값은 0.10/−0.058/−0.040이다. 두 기록을 현재값으로 고치고 변경 사실을 적었다.
  - **한글 원고에 이전 삽입 스크립트가 남긴 깨진 문장 2개**('표준편차 0. [다른 문장] .04', '0. [다른 문장] .49')와 엉뚱한 절에 들어간 '이탈 표본 32개' 문장이 있었다. **정규식으로 문장을 끼워 넣을 때 '0.' 뒤를 문장 끝으로 오인하지 않는다.**
  - **progeria 공여자는 9명이 아니라 10명이다.** 기탁 143개 = HGPS 10 + 정상 133(20세 미만 26 + 성인 107). '기탁 142개 시료 집합'에 HGPS가 9명인 것은 HGPS 라이브러리 하나(SRR7093948)에 junction 지표가 없어 분석 표에 들어가지 않았기 때문이다. 본문은 '10명 제외'로 쓰고 Methods에 142개인 이유를 적는다.
  - **GTEx 배양 이탈 표본을 뺀 원시 상관은 0.73이 아니라 0.72 → 0.89**(Spearman 0.7246; 0.74는 Pearson).
  - **커버레터의 '시료 3,197개'는 어느 표에도 없는 수였다.** 조직별 공여자 수의 합은 2,715다.
  - **'공변량 중 세포주 출처와 시퀀서 항이 잔여를 0에서 떨어뜨린다'는 틀렸다.** 세포주 출처 항 하나가 그렇게 하고(−0.038, *P* = 0.004), 시퀀서 항은 오히려 0 쪽으로 움직인다. 근거 표 `reviewer_sensitivities/covariate_attribution.tsv`.
  - *HNRNPK*의 분열 속도 상관은 0.67이 아니라 0.66(0.6649)이다.
  - Fig. 2e에는 혼합모형 추정치가 없다(β = 0.69는 패널 a). '사전등록 2차 분석 S3'이라는 표현은 사전등록에 없다 — S1–S5는 '성공 기준 없는 탐색 분석'이다.
  - **해소됨(2026-09-17)**: 로컬 대비 `LOCAL_ReproCM`이 §5의 12,319개 anchor가 아니라 3,118개 유전자 표에서 계산되던 문제는, 로컬 자료를 전부 내리면서 사라졌다(§18).

- **세포 구성 보정은 증식 축을 눕히지 않는다. splicing만, 피부에서만 지운다 (2026-09-18, 사후).** 심사 지적(편집자 BC 2)에 따라 `gtex_posthoc.R`의 구성 보정을 splicing 두 집합에서 **전체 프로그램으로 확장**했다(`composition_adjusted_all_programmes.tsv`). 결과:

| 조직 | cell cycle 양성 대조 | splicing 96 |
|---|---|---|
| 배양 | 0.941 → **0.952** | 0.739 → **0.722** |
| 하지 피부(일광) | 0.670 → **0.666**(−0.004) | 0.225 → **0.101**(−0.124) |
| 치골상부(비노출) | 0.597 → **0.586**(−0.011) | 0.056 → **−0.058**(−0.114) |
| 근육 | 0.031 → 0.157 | −0.115 → −0.040 |

  피부에서 구성 보정은 **양성 대조를 건드리지 않는다**(0.004·0.011 하락). 같은 보정이 splicing의 결합은 절반 이상 지운다. 선택되지 않은 정의도 같은 방향이다(mRNA splicing 0.134 → 0.063, pre-mRNA processing 0.145 → 0.082). collagen formation은 0.030 → 0.027로 움직이지 않는다. 배양에서는 아무것도 움직이지 않는다. **따라서 "피부에 남은 결합"은 세포 구성이 만든 것이고, 범위 제한·검정력 반론은 양성 대조가 보정 뒤에도 0.67·0.59로 살아남는 것으로 기각된다.** 사전등록 추정치는 보정 전 0.225이므로 본문은 그것을 앞세우고 이 값은 사후로 밝힌다.
- **일광 노출 대 비노출 피부 차이는 유의하다 (2026-09-18, 사후).** splicing 96은 0.225 대 0.056(Fisher *z* = 3.09, *P* = 0.0020), 양성 대조는 0.670 대 0.597(*z* = 2.19, *P* = 0.029). **양성 대조도 차이가 나므로 "대조는 같은데 splicing만 다르다"고 쓰지 않는다.** 쓸 수 있는 것은 기울기 차이다(대조 11% 하락, splicing 4배 하락). 치골상부는 사전등록 H5이고 지지되었다.

이전 방향(상세는 보관본):

- miRNA target GSEA: miR-302/367 NES = −1.05, FDR = 0.399. miRNA 인과를 지지하지 않는다.
- longevity: Tyshkovskiy maximum-lifespan signature 상관 약 −0.055(기대와 반대). GSE165177 MPTR ρ = 0.0004.
- secretome class 사전등록: HS2 실패, HS1 해석 보류(`repro_cm_secretome_class_results_ko.md`).

## 9. 현재 Figure 구성과 생성 스크립트

출력 `public_data_tierA/derived/figures_sciadv/`, 테마 `scripts/sciadv_theme.R`(2026-09-15 재작성: 소문자 패널, `rp()`·`pfmt()`·`sci()`·`num_axis()`, ρ·*R*²·수학 minus, 기본 8 pt, 최소 7 pt, 판정색은 PASS 파랑·FAIL 주황). 스크립트 이름은 역사적 이유로 그림 번호와 다르다.

| Fig. | 질문 | 스크립트 |
|---|---|---|
| 1 | splicing 기계가 공여자 연령·배양 시간에 따라 감소하는가 | `make_sa_fig1.R` |
| 2 | 측정된 분열 속도를 따라가는가(반복 측정 보정, **패널 e: 20개 marker 대리 지표의 품질**) | `make_bio_fig2.R` |
| 3 | 증식이 연령 연관의 얼마를 설명하고 무엇이 남는가(무작위 집합 대조 포함) | `make_sa_fig3merged.R` |
| 4 | 63개 대비에서 결합이 일반화되는가 | `make_sa_fig3.R` |
| 5 | 어떤 노화 프로그램이 분열 읽기값인가 | `make_sa_fighall.R` |
| 6 | 결합은 조직에서도 성립하는가(사전등록 GTEx, **패널 a: 배양·하지 피부·치골상부 피부**(근육은 b의 heatmap에만), **패널 d: 구성 보정의 프로그램별 비교**) | `make_sa_figgtex.R` (Fig6.png·FigS8.png·보충 표 9 스코어카드를 함께 쓴다) |
| 7 | graphical abstract: 배양 접시·공여자·피부 그림 + 각 대상에서 한 측정의 실제 산점도(2026-09-20 재작도, §24) | `make_sa_fig7summary.R` |
| S1 | 처리로 증식을 섭동해도 결합이 유지되는가 | `make_sa_figpert.R` |
| S2 | 어떤 유전자가 연령 연관을 유지하는가 | `make_sa_figS_gene.R` |
| S3 | splicing outcome 지표와 event 수준 splicing(2026-09-18 병합, 5패널) | `make_sa_figS_outcome.R` |
| S4 | 개입(표본 외 예측) | `make_bio_fig6.R` |
| S5 | 결합이 splicing에 특이적인가 | `make_sa_figspec.R` |
| S6 | 어떤 senescence 패널이 분열을 보고하는가 | `make_sa_figsen.R` |
| S7 | methylation 시계(사전등록) | `make_sa_figmeth.R` |
| S8 | GTEx 사후 분석 + 공여자 일치 | `make_sa_figgtex.R` |

**보충 그림 10 → 8 (2026-09-18, 사용자 결정).** 옛 S3(outcome 지표)과 S4(event 수준)는 같은 음성 결과 하나를 나눠 지고 있어 **한 그림 5패널로 병합**했고(`make_splicing_fig.R` 은퇴), 옛 S5(전사체 나이 예측기)는 사후 분석이고 수치가 이미 보충 표 3에 있어 **그림을 내리고 표만 남겼다**(`make_sa_figS_agepred.R` 은퇴). 나머지는 두 칸씩 당겼다(옛 S6–S10 → 현 S4–S8). 본문 그림 6개는 그대로다.

**덮어쓰기 주의 (2026-09-18 실행 차단 완료).** 은퇴 스크립트 **26개**가 살아있는 그림 파일명에 쓴다. 특히 `make_sa_fig4.R`은 **Fig6.png**(GTEx 경계)에 쓰면서 `LOCAL_ReproCM` 행이 남은 `conserved_core/splicing_vs_proliferation.tsv`를 읽으므로, 한 번만 잘못 실행해도 내린 비공개 대비가 게재 그림으로 돌아간다. 다른 예: `make_journal_fig6.R`·`make_pub_fig6.R`(→ Fig6.png), `make_sa_fig5gene.R`·`make_bm_fig3.R`·`make_journal_fig3.R`·`make_pub_fig3.R`(→ Fig3.png), `make_sa_fig4new.R`(→ Fig4.png), `make_sashimi_local.R`(→ FigS3.png).

**2026-09-18에 26개 전부에 하드 가드를 넣었다.** 각 파일 맨 위에 `FIGURES_WRITTEN` 선언과 `stop()`이 있어 그냥 실행하면 멈춘다. 일부러 돌려야 하면 `ALLOW_RETIRED_FIGURE_SCRIPT=1`을 준다. 가드를 지우지 않는다. 그림을 다시 만들 때는 위 표의 스크립트 15개만 쓴다.

**그림 규칙(리뷰 반영).** 패널 문자는 소문자 굵은 a/b/c. `rho`·`R2`·`9e-04`·`P = < 0.001`·하이픈 음수는 쓰지 않는다. 모든 상관은 n과 단위(공여자/시료/유전자/대비)를 밝힌다. 177개 점수의 축 이름은 "pre-mRNA processing score" 하나로 통일한다. Fig. 7만 사전등록 96개 판을 쓰고 범례에 밝힌다. 렌더 뒤 PNG를 직접 열어 겹침·잘림·보이지 않는 점을 고친다.

**거명 splicing factor 목록.** `scripts/revision/named_splicing_factors.txt`(20개) 하나만 쓴다. 출처는 Holly 2013 표 3(배양 사람 섬유아세포에서 측정된 패널)이며 근거는 `scripts/revision/named_splicing_factors_source.md`와 `scripts/revision/references_v6.md` §7에 있다. 빨강 강조는 "인용 문헌 2편 이상에서 측정으로 연령 연관이 보고된 4개(*SRSF1*, *SRSF2*, *SRSF6*, *HNRNPK*)"라는 기준으로만 쓴다. *U2AF1*·*PRPF8*·*SNRPA/B/D2/F*·*SF3B5*·*RBM39* 등은 사람 연령 연관 출처가 없으므로 "문헌이 거명한다"고 쓰지 않는다.

**개정 분석 스크립트(2026-09-15).** `scripts/revision/`: `cohort_gse113957.R`(코호트 정의), `run_cohort_core.R`(Fig. 1c·3·4 수치), `run_stats_supplements.R`(Fig. 2 혼합모형·Fig. 1d 정확검정·Fig. 5 대비와 무작위 집합·집합 재정의 민감도), `run_psi_cohorts.R`(보충 S1), `run_interventions.R`(Fig. 8), `gtex_posthoc.R`(Fig. S2), `hallmark_extensions.R`(Fig. 6), `make_supp_tables.R`(보충 표 1·4), `references_v6.md`, `FIGURE_SPEC.md`.
## 10. 사전등록 기록

규칙: 가설, 정의, 제외 기준, 성공 기준을 자료를 받기 전에 Markdown으로 고정하고 SHA-256을 기록한다. 분석 스크립트는 실행 전에 해시를 확인한다. 고정 뒤에는 문서 끝의 변경 이력에만 추가하고, 바꾼 경우 원래 정의의 결과도 함께 보고한다.

| 문서 | 방향 | 고정 해시 | 현재 파일 |
|---|---|---|---|
| `preregistration_gtex_culture_tissue_boundary_ko.md` | 현재 | `61972b99…` | 동일. 결과 `results_gtex_culture_tissue_boundary_ko.md`, 영문 번역 `supplementary_note1_preregistration_en.md` |
| `preregistration_methylation_clock_ko.md` | 현재 | `be51a755…` | 동일(변경 없음). 결과 `public_data_tierA/derived/methylation_clock/preregistered_tests.tsv`, 영문 번역·결과·이탈 기록 `supplementary_note2_preregistration_methylation_en.md`. `scripts/revision/run_methylation_clock.R`가 실행 전 해시 확인 |
| `repro_cm_preregistration_secretome_class_ko.md` | 이전 | `55bf016d…` | 동일 |
| `repro_cm_preregistration_heldout3_photoprotection_ko.md` | 이전 | `315e68ea…` | 동일 |
| `repro_cm_preregistration_heldout4_sauchinone_ko.md` | 이전 | `c542d32f…` | 동일 |
| `repro_cm_preregistration_phase2_ko.md` | 이전 | `41e4cf4b…` | `c02865ed…`. 앞 130줄이 고정 해시를 재현한다(원문 보존). 변경 이력 2회 추가. 기록 시점 문제는 아래 |
| `repro_cm_preregistration_heldout2_senescence_ko.md` | 이전 | `2b2f298e…` (11:43:03) | `6e5a7703…`. 앞 76줄이 고정 해시를 재현한다(원문 보존). 자료 다운로드 뒤 metadata를 보고 11:44:27에 변경 이력을 추가했고, 첫 분석 출력은 11:45:27이다. `.sha256` 파일만 갱신되지 않았을 뿐 **문제 없음** |

**heldout 사전등록 점검 (2026-09-13, 세션 기록과 파일 시각 대조, KST).** phase2의 첫 변경 이력(11:11:33)은 "결과를 보기 전 기록"이라고 적혀 있지만, 11:09:45 실행이 이미 primary test 결과와 판정(OVERALL: PARTIAL)을 출력한 뒤다. 그중 HO4b의 Primed-D13 primary와 HO2의 성별 추정은 결과 출력 전(11:07:46)에 작성된 스크립트에 이미 들어 있었다. 매핑률 제외 규칙을 발현 유전자 기준으로 해석한 것은 결과를 출력한 같은 실행에서 도입되었고, 원래 문구(전체 행 기준 < 70%)대로라면 GSE297233(59.4%)은 제외되었을 것이다. 둘째 변경 이력(19:16:19, "결과 산출 전")도 Phase 2 결과가 출력된 19:09·19:12 뒤에 기록되었다(IR burden 계산 불가, unannotated junction fraction을 탐색 지표로 추가). 따라서 이전 문서(`manuscript_benchmark_v1.md`, `repro_cm_reviewer_assessment_ko.md`)의 "결과 확인 전 고정" 서술은 phase2에 대해 정확하지 않다. unannotated junction fraction(현재 원고의 splicing outcome)은 이때 사전 지정되지 않은 탐색 지표로 처음 등장했으므로, 원고에서 사전 지정 지표로 표현하지 않는다. heldout3, heldout4, secretome class는 고정 뒤 변경이 없다.

GTEx 사전등록 판정: H1 지지(0.74), H2 지지, H3 지지(−0.12), H4 평가 불가, H5 지지(0.05), H6 불지지. 사전등록 규칙에 따라 GTEx 조직의 연령 효과는 어떤 주장의 근거로도 쓰지 않는다.

methylation 사전등록 판정: M1 지지, M2·M3·M4 불지지(§7 S4). 사전등록이 정한 보고 한계(문단 하나, 패널 한 세트, 절대 나이 해석 금지)를 지킨다. M3의 신뢰구간이 +0.34까지 걸치므로 **"연관 없음의 증거"로 쓰지 않는다**. 이탈 두 건(Hannum 시계 추가, clock acceleration의 공통 기울기 적합)은 보충 Note 2에 적었다.

## 11. 논리적 오류와 금지 표현

1. 로컬 `n=1` 자료에 DEG/FDR 언어를 쓰지 않는다.
2. 유전자를 biological replicate처럼 써서 작은 p-value를 만들지 않는다.
3. pathway enrichment p-value를 시료 수준 처리 유의성으로 해석하지 않는다.
4. 선택한 pathway만 보여주고 전체 결과를 숨기지 않는다.
5. "측정된 분열 속도"를 GSE179848 밖의 자료에 쓰지 않는다.
6. 배양 결과를 조직으로 일반화하지 않는다.
7. 보정 전 연령 효과를 코호트끼리 비교하지 않는다. 비교는 보정 후 추정치로 한다.
8. GTEx 조직의 연령 효과, collagen 연령 신호, proteostasis의 분열 비의존성을 주장 근거로 쓰지 않는다.
9. 사후 분석을 validation이나 held-out 결과로 포장하지 않는다. GTEx 사후 3종(범위 제한·구성 보정·양성 대조)은 그림과 본문 모두에서 post hoc이라고 밝힌다.
10. 상관과 인과, 전사 정렬과 표현형 역전을 구분한다.
11. Repro-CM에 대해 miRNA transfer·인과, TIMP2 인과, broad age reversal, reprogramming mimicry, longevity program을 주장하지 않는다.
12. 선행 논문과의 중복을 무시하지 않는다. JTE 원고는 심사 중이므로 번호 참고문헌에 넣지 않고 본문에 `(Shim, V. *et al.*, manuscript submitted)`로 쓴다.
13. 음성 결과를 삭제하지 않는다.
14. `heldout*/` 자료를 held-out 검증이라고 쓰지 않고, 이전 방향의 사전등록(phase2, heldout2–4, secretome class)을 현재 원고에서 인용하지 않는다.
15. **splicing outcome(주석 없는 junction 비율)과 *COL1A2*에 대해**: 연령 증가, 점진적 재분배, 증식 보정에 견디는 성질을 주장하지 않는다. 성인 코호트에서 연령 효과는 없다. 이 지표는 사전 지정된 적이 없는 탐색 지표라고 Methods에 적는다.
16. **"보정 후 감소가 남는 유전자"를 대안 지표로 권하지 않는다.** 코호트 정의에 따라 40/56/19/0개로 바뀌고 공통 유전자가 없다.
17. **GSE113957을 쓸 때는 반드시 코호트를 밝힌다.** HGPS 10명(142개 집합 안에서는 9명)과 20세 미만 26명을 빼고, 세포주 출처·시퀀서·성별을 공변량에 넣는다(`scripts/revision/cohort_gse113957.R`). 83세 이상 31명은 전원 AG 계열이므로 공변량으로 분리되지 않는다는 한계를 함께 적는다.
18. **거명 splicing factor는 출처 있는 목록만 쓴다**(§9).
## 12. 원고·통계·그림 작성 규칙 (사용자 지정)

**한글판.** 영문 원고나 보고서를 만들거나 고칠 때마다 한글판을 함께 만든다(md, docx, 읽기용 아티팩트). 수치·그림은 동일하게, 유전자명·데이터셋 번호·통계 용어는 영문 유지.

**통계.** 전통적인 biology 저널(npj Aging, Aging Cell) 방식을 쓴다: Spearman/Pearson과 *P*, 공변량을 넣은 선형 모형, Fisher *z* 신뢰구간, BH-FDR, Wilcoxon, 평균과 95% 신뢰구간의 forest plot. ML·통계 저널식 지표, 임의 지수("module score +0.27", "framework", "benchmark" 식 표현), 통계 지표 나열을 피한다. 모든 정량 주장은 이름 있는 생물학적 양으로 쓰고, 생물학(유전자, locus, pathway)을 통계보다 앞에 둔다.

**그림.** Kerepesi et al., *Sci. Adv.* 2021(doi 10.1126/sciadv.abg6082)의 형식을 따른다: 대문자 패널, 패널 위 r/ρ와 *P*, small multiples, 단순한 도식. 각 그림은 생물학적 이유와 내용(이름 붙인 유전자, sashimi, 공여자 산점도, 대비 이름)을 담는다. 라벨 없는 heatmap이나 AI티 나는 dot plot은 쓰지 않는다. 실패·음성 결과는 한 패널로 보여주면 충분하다. 렌더 후 PNG를 직접 보고 겹침·잘림을 고친다.

**글쓰기.** 제목은 흔한 과학 단어로 된 서술문으로 쓴다("전사 축" 같은 조어 금지). 논문 안에 자기비판 문단을 넣지 않고, 실패는 사실만 짧게 적는다. "our own data"라고 쓰지 않는다. contact inhibition 같은 낡은 개념으로 틀을 짜지 않는다.

## 13. 분석 품질 점검

**자료와 전처리:** 파일 무결성과 시료 정체, library size·mapping·RNA 무결성, 유전자 ID 매핑과 중복 처리, metadata 수작업 검증(연령 형식, batch와 연령의 혼동).

**통계:** 효과 크기와 방향, 코호트 간 이질성, 증식 보정 전후 비교와 증식–연령 방향 보고, 기술 공변량(RIN, 허혈 시간, library size), 다중검정 보정, 같은 연구의 여러 대비를 독립 자료처럼 중복 계수하지 않기.

**일반화:** 결과 확인 전에 고정한 가설과 기준, 독립 코호트에서의 재현, 음성 대조와 경계 조건, 전체 결과 공개.

**그림마다 답할 질문:** 이 패널이 검정하는 주장은 무엇인가? biological n은? 효과 크기와 방향이 보이는가? 불확실성이나 적절한 null이 표시되어 있는가? *P*/FDR의 단위(유전자, 시료, 대비)가 분명한가? 상관인지 인과인지 구분되는가?

## 14. 투고 가능성에 대한 현재 판단

목표는 GeroScience(Original Article)다. 2026-09-17 독립 심사 패널(편집자·생물학 리뷰어·통계 리뷰어·스토리 에디터·최소주의자)의 만장일치 권고에 따라 본문 그림을 6장으로 줄였다. 새 기여는 셋으로 정의한다: (1) 추정이 아니라 *센* 분열 속도를 앵커로 쓴 것, (2) 증식이 설명하는 몫과 두 코호트에서 재현되는 작은 잔여의 정량, (3) 사전등록한 배양–조직 경계.

강점:

- 세포 계수로 측정한 분열 속도(328개 라이브러리, 7개 세포주), 반복 측정을 고려한 혼합모형에서도 유지
- 코호트를 정확히 정의한 공여자 분석(정상 성인 107명, 세포주 출처·시퀀서 보정)과 네 가지 정의의 민감도 공개
- 63개 대비(23개 GEO series) 일반화 + 발현량을 맞춘 무작위 집합 1,000개 대조
- 사전등록한 GTEx 배양 대 조직 검정(배양 511/492, 피부 739·626, 근육 839명)과 사후 반론 검토(범위 제한·세포 구성·양성 대조)
- 실패·불재현·철회한 주장을 모두 그림과 본문에 남김

위험:

- in silico만으로 구성. 우리가 분열 속도를 조작하지는 않았지만, 기탁 자료의 처리 9종이 조작한 섭동에서 결합이 유지된다(보충 Fig. S5). 다만 contact inhibition 외에는 속도 변화 폭이 작고, 처리 2개는 증식으로 설명되지 않는 잔여를 남긴다
- 결합 자체는 알려져 있어 '새 발견'으로 읽히면 기각 위험. 기여는 크기·경계·잔존으로 정의한다(§1의 위치 잡기 규칙)
- GSE179848 밖에서는 증식 측정이 전사체 기반
- 조직은 프로그램 발현 수준만 검토했고 조직 splicing outcome은 측정하지 않음. Zhang/Gladyshev bioRxiv 2026과의 관계를 고찰에서 분명히 한다
- 두 번째 기둥(splicing outcome)을 철회했으므로 논문의 주장이 하나로 줄었다. 그만큼 "알려진 결합의 정량"이라는 성격이 강해졌다
- methylation 대조의 핵심 검정(M3)은 짝지은 배양이 109개뿐이라 검정력이 부족하다. 신뢰구간이 +0.34까지 걸치므로 경계 진술로만 쓰고 "연관 없음"으로 쓰지 않는다
- 사전등록이 공개 등록이 아니라 내부 기록이다(보충 Note 1에 명시)
- **투고 필수 항목 중 사용자만 채울 수 있는 것이 남아 있다**(§17)

## 15. 먼저 읽을 파일

- 원고 v8: `manuscript_splicing_v2.md`, `manuscript_splicing_v2_ko.md` (+ 같은 이름의 `.docx`, `.html`). **보충 Results 전문: `supplementary_results_en.md`, `supplementary_results_ko.md`** — 본문에서 옮긴 분석 8건이 여기에 있다. 본문 포인터를 고치면 이 파일도 같이 고친다
- 덜어내기 판단 근거(2026-09-17): 내용 단위 목록 아티팩트 https://claude.ai/artifact/2D9nKSKmpsBDmeNxaxwHdk , 심사 패널 결과 `…/scratchpad/panel_result.json`은 세션 임시 파일이므로 요지는 `revision_response_ko.md` 9차 처리에 옮겨 적었다
- 제출 전 모의 심사(2026-09-15): `review_presubmission_npjaging_ko.md` (+ `.html`, 아티팩트 https://claude.ai/code/artifact/aa94de80-764a-46e9-af8b-edf1eedd2a43 ). 필수 6개 항목과 리뷰어 3명 의견. **모든 항목의 처리 결과는 §7·§8과 `scripts/revision/`에 있다.**
- 요청사항 재점검: `review_audit_ko.md` (항목별 판정, 남은 취약점 4가지, 아티팩트 https://claude.ai/artifact/WnpBWqXL26bYxQyqpUdPQe )
- 투고 서류(2026-09-16 작성): `supplementary_information.pdf`(+ `.docx`, 빌드 `bash scripts/artifact_build/build_supplementary_pdf.sh`) — 보충 Results 1–8, 보충 그림 S1–S8, 표 1–10, Note 1–2, Supplementary Data 목록을 단일 PDF로(86쪽). `cover_letter_geroscience.md`(+ `.docx`, 한글판 `_ko.md`)
- 투고 규정 대조: `submission_checklist_ko.md` (2026-09-16 저널 페이지 직접 대조, 남은 항목 7개)
- 대응 기록: `revision_response_ko.md` — 모의 심사 지적별로 무엇을 어떻게 고쳤는지와 남은 항목
- 합본 읽기 페이지: `revision_bundle_ko.html` (1부 대응 기록, 2부 원고 v6 한글 전문과 그림 11개, 3부 모의 심사 전문). 빌드는 `bash scripts/artifact_build/build_review_bundle.sh`, 아티팩트 https://claude.ai/artifact/4pfbcSCkBZFyc5twWiAVot (2026-09-19 새 계정으로 재게시)
- 개정 설계서: `scripts/revision/FIGURE_SPEC.md`(그림 번호 변경, 새 수치, 타이포그래피 규칙), `scripts/revision/references_v6.md`(참고문헌 36개와 자료 분류 확인, §7에 거명 유전자 출처 판정)
- 사전등록: `preregistration_gtex_culture_tissue_boundary_ko.md`, 결과 `results_gtex_culture_tissue_boundary_ko.md`, 영문판 `supplementary_note1_preregistration_en.md`, 사후 분석 `results_gtex_posthoc_ko.md`, 프로그램 재계산 `results_hallmark_revised_ko.md`. **methylation**: `preregistration_methylation_clock_ko.md`, 영문판·결과 `supplementary_note2_preregistration_methylation_en.md`, 산출물 `public_data_tierA/derived/methylation_clock/`
- 읽기용 아티팩트(**2026-09-19 새 계정 k0611c@gmail.com으로 재게시** — 이전 계정의 `PtwCBB7Y…`/`5hYrp2ow…` 링크는 이 계정에서 갱신 불가): 영문 https://claude.ai/artifact/TwvN7zBRxMxXpEg4ZrChjY , 한글 https://claude.ai/artifact/KH59jVRBduohNoJtqqAw93 . 모의 심사 https://claude.ai/artifact/N4ikBartGYSpgcN2Fs6R38 . 다시 게시할 때는 이 URL을 `url`로 넘긴다(`/code/artifact/` 형식 URL은 fetch가 실패한다)
- **쉬운 설명 길잡이(2026-09-20 신설)**: `review_pages/story_guide_ko.html`, 아티팩트 https://claude.ai/artifact/EwGgVtCfoRxNyYF9oJYYFP . 한 문장 주장 → 왜 이 질문인가 → 논리 사슬 6단계 → **그림 15장이 각각 왜 있는지(없으면 무엇이 무너지는지)** → 하지 않는 주장 → v10–v13 변경 → 남은 일. 원고를 고치면 이 페이지의 수치도 같이 본다.
- **리뷰 페이지 로컬 사본**(git 제외, `review_pages/`): 편집자 결정문 `editorial_decision_20260918.html`, 리퍼리 4명 요약 `referee_reports_20260918.md`, 스토리 지도 `story_map_20260919.html`(아티팩트 https://claude.ai/artifact/6juTkLDKRePKVjuHU5aMP4 ), 합본 `revision_bundle_ko.html`. 계정이 바뀌면 아티팩트 링크는 죽으므로 이 사본이 기준이다.
- 결과 디렉터리(`public_data_tierA/derived/`): `cohort_revised/`, `revision_stats/`, `hallmark_revised/`, `gtex_posthoc/`, `psi_cohorts/`, `interventions_revised/`, `methylation_clock/`, `treatment_perturbation/`, `machinery_specificity/`, `residual_meta/`, `senescence_panels/`, `interventions_revised/`(표본 외 예측 포함), `division_rate/`, `conserved_core/`, `gtex_boundary/`, `figures_sciadv/`
## 16. 작업 순서

1. 이 파일과 원고 v8, `supplementary_results_en.md`, `revision_response_ko.md`, GTEx·methylation 결과 기록을 읽는다.
2. 새 공개 데이터가 필요하면 후보, 크기, 디스크 여유, 목적, 대비를 먼저 제시하고 명시적 승인을 받는다. processed 파일만 받고 중앙 목록과 디렉터리 sha256 양쪽에 기록한다.
3. 받기 전에 사전등록을 고정한다(§10).
4. 성공과 실패를 모두 보고하고, 결과를 본 뒤 정의를 바꾸지 않는다.
5. GSE113957을 쓰는 분석은 `scripts/revision/cohort_gse113957.R`의 코호트 정의를 통해서만 한다(§11.17).
6. 그림은 `sciadv_theme.R`로 만들고 렌더 결과 PNG를 직접 열어 확인한다. 번호를 바꾸면 `save_fig(..., "FigN.png")`와 원고 범례를 함께 바꾼다.
7. 원고를 고치면 한글판을 같이 고치고 `bash scripts/artifact_build/build_reading_artifacts.sh`로 docx와 HTML을 다시 만든 뒤, 두 아티팩트를 같은 URL에 다시 게시한다. 제목이나 요약 수치가 바뀌면 `scripts/artifact_build/tpl_v5*.html`의 머리말과 요약 띠를 먼저 고친다.
8. 한 문장 주장(§2)이 모든 그림과 맞는지 다시 점검하고, 본문 수치를 산출 파일과 대조한다.
## 17. 완료 기준

- Fig. 1–6과 보충 S1–S8이 각각 하나의 질문과 생물학적 읽을거리를 가진다.
- 모든 주장이 §2의 경계 안에 있다.
- 사전등록 결과가 H4·H6을 포함해 모두 보고되어 있고, 영문판(보충 Note 1)이 있다.
- 제출 전 모의 심사의 필수 6개 항목이 모두 해결되어 있다(2026-09-15 완료).
- 한글판이 영문판과 수치·그림·구조에서 일치한다.
- 저널 형식을 충족한다(GeroScience: 초록 250단어 이하 — 현재 250; 제목 15단어; 참고문헌 **41개**).
- **남은 것(사용자만 채울 수 있음): 저자 목록·소속·교신저자, 코드·자료 저장소 DOI, 저자 기여·이해상충·감사, Reporting Summary와 editorial policy checklist.** 윤리 승인 문장·로컬 GEO 기탁·JTE 사본 동봉은 2026-09-17에 해소됐다(§18).

이 조건을 충족하지 못하면 "투고 준비 완료"라고 보고하지 않는다. 현재 증거 수준, 가장 큰 결손, 추가 자료로 해소할 수 있는 부분을 구분해 사용자에게 제시한다.
## 18. 로컬 자료 제거 (2026-09-17, 사용자 결정)

**결정.** "로컬자료는 빼고 진행하고 대신 논문을 citation해서 선행연구가 진행됐다는 식으로 마무리한다." 이 연구는 **공개 자료 재분석**이고, 선행 연구(JTE-Aug-26-0256)는 자료가 아니라 고찰의 인용으로 잇는다.

**내린 것**: `LOCAL_ReproCM` 대비 1개(Fig. 4·보충 Fig. S6·보충 표 4), 로컬 정렬 4개 라이브러리의 read 수준 패널(보충 Fig. S4c), Methods "Read-level junction support" 절, 그 절에서만 인용하던 STAR·samtools 참고문헌 2개, 윤리·데이터 가용성의 placeholder 두 개, 통계 절의 "반복이 없는 단일 대비" 단서.

**수치**: 64 → **63개 대비**(series는 23개 그대로). ρ 0.865 → **0.864 [0.78, 0.92]**, *R*² 0.923 그대로, 공유 유전자 제거 0.851/0.911 그대로. 무작위 집합 최대 0.58 → **0.56**. senescence 제외 n = 55 → **54**, 연구당 한 값 n = 30 → **29**. 개입 22 → **21**개, 종류 제외 적합 중앙값 |잔차|/|변화| 0.048/0.116 → **0.045/0.111**, secretome 종류 평균 −0.03 → **−0.02**. 참고문헌 43 → **41개**. **판정은 하나도 바뀌지 않았다.**

**더한 것**: 고찰 마지막 문단 하나. 선행 연구를 `(Shim, V. *et al.*, manuscript submitted)`로 인용하고(§11.12에 따라 번호 목록에는 넣지 않는다), 자료를 공유하지 않으며 그 결과를 검정하지 않는다는 것, 이 분석이 내놓는 것은 측정 자체의 읽기값이라는 것을 적었다. **Repro-CM에 대한 주장은 하나도 하지 않는다**(§2·§11.11 그대로).

**투고 부담**: 윤리 승인·시료 출처 문장, 로컬 RNA-seq GEO 기탁, JTE 원고 사본 동봉 세 항목이 사라졌다. 커버레터의 관련 원고 선언은 유지한다.

**검증에서 함께 잡은 오래된 오류 4개** (로컬 제거와 무관): Metabolism of RNA ρ 0.77 → **0.76**, translation 0.58 → **0.57**, ECM organisation −0.59 → **−0.60**(`fig6_other_programmes.tsv`), 종류 잔차 원시 *P* 0.25–0.94/BH 0.94 → **0.22–0.85/BH 0.85**(`fig6_class_residuals.tsv`). 마지막 것은 Fig. 4c 패널 주석에도 0.94로 박혀 있어 파일을 읽도록 고쳤다(`make_sa_fig3.R`). 그림 주석에 통계값을 문자열로 박지 않는다.

**되돌리려면**: `scripts/revision/run_stats_supplements.R`의 `INCLUDE_LOCAL <- FALSE`를 TRUE로 바꾸고 `run_interventions.R`·`make_supp_tables.R`·`make_sa_fig3.R`·`make_bio_fig6.R`·`make_splicing_fig.R`을 되돌린 뒤 그림과 표를 다시 만든다. 그 전에 §8의 3,118개 유전자 표 문제를 먼저 고쳐야 한다.

## 22. 개념 요약 그림, 근육 강등, 논리 사슬 (2026-09-19, 사용자 요청; 원고 v11)

**Fig. 7 신설** (`scripts/make_sa_fig7summary.R` → `Fig7.png`; 고찰 첫 문단 끝과 신설된 둘째 문단에서 인용). 사용자가 Sturm et al. 2023 *Commun. Biol.* (s42003-022-04303-x) Fig. 9 같은 연구 전체 개념도를 요청했다. (a) 배양 / 배양의 연령 효과 / bulk 피부 세 열의 도식에 원고 수치를 실었고, (b) 자료 7행의 결합을 Results 순서로 늘어놓은 막대 + Fisher *z* 95% CI: 센 속도 0.71, 비선택 Reactome mRNA splicing 0.63, 공여자 107명 **Spearman 0.91**(본문 Fig. 3a는 Pearson 0.89 — 범례에 밝혔다), 63개 대비 0.86, GTEx 배양 0.74, 피부 0.23·0.06(+ 구성 보정 후 0.10·−0.06을 빈 기호로). **(b)의 값은 전부 산출 표에서 계산하고, (a)에 적은 수치는 `chk()`로 같은 표와 대조한 뒤에야 그림을 쓴다.** 근육은 그리지 않는다(범례에 이유). §21의 "자료 없는 flowchart 금지"의 유일한 예외이며, 자료 패널 (b)와 짝지어 둔다. 값 파일 `figures_sciadv/fig7_coupling_strip.tsv`.

**근육 강등(사용자 결정: "근육은 검정할 수 없는데 굳이 언급 안 해도").** 잠긴 H3 판정과 보고 의무는 그대로다(보충 표 9·10, Fig. 6b heatmap). 내린 곳: 초록, 서론 마지막 문장, 6절 제목, Fig. 6 제목과 **패널 a(셋째 facet을 골격근 → 비노출 피부로 교체; `SC`/`LB`/`STRIP`)**, 기울기 목록과 구성 보정 문장의 근육 값, 고찰 첫 문단과 조직 문단. 남긴 곳: 6절 둘째 문단 끝의 한 문장(H3 값 −0.12 + 증식 점수가 어떤 프로그램도 회수 못 함 + 근거로 쓰지 않음), 고찰 조직 문단의 한 문장(refs 29·30을 유지하려고), 한계 문단("피부와 골격근만"), Methods. 커버레터의 근육 한 문장은 편집자 BC2 대응이므로 유지.

**논리 사슬 보강(사용자 질문: "다 독립된 자료라 각 연결점이 약하다").** 해법은 각 절 첫머리에 **앞 절이 내놓는 예측**과 **넘어가는 도구**를 적는 것이다. 2절: Fig. 1a에서 세포주기가 앞선다 → 센 속도를 따라가야 한다. 3절: Fig. 2의 결합이 Fig. 1b의 연령 효과를 만든다면 보정이 대부분을 없애야 하며, 20개 marker 점수는 Fig. 2e에서 계수에 대해 보정한 도구다. 4절: 연령·계대만이 아니라 어떤 섭동이든 움직여야 한다. 5절: Fig. 4가 이미 기계 전체를 보였다 → 다른 프로그램도 읽기값일 수 있다. 6절: 분열하는 세포의 성질이면 조직에서 약해지고 복제 프로그램은 회수되며 collagen은 무관해야 한다(복제 프로그램 회수는 **사전등록 아님, 내부 대조**라고 명시). 고찰 둘째 문단(신설)이 "시료를 공유하지 않는 자료 + 두 도구(177개 집합, 보정된 증식 점수) + 예측 사슬"을 정리하고 Fig. 7b를 가리킨다. 사슬 값: 모든 배양 자료에서 ρ 0.63–0.91, 조직에서만 약해진다.

**함께 고친 것.** GTEx 배양의 CI "0.69 to 0.78"은 어떤 산출물에도 없는 값이었다 → Methods가 말하는 Fisher *z*(완전 사례 n)로 **0.70–0.78**, 비노출 피부 **−0.02–0.13**(하지 피부 0.15–0.29는 그대로). 한글판의 묵은 인용 3건(Fig. 6d → 6c 두 곳, Methods의 보충 Fig. S8 → S6)과 깨진 문장 1건("공여자의 그러나"), 영문 4절의 빠진 마침표, 한글 4절의 중첩 괄호. 초록 245단어. `xref.py`의 H3/H4/H6·"보충 표 6·10" 불일치 보고는 정규식 한계(한글 조사, "Tables 6 and 10")이지 실제 불일치가 아니다.

## 24. Fig. 7을 graphical abstract로, 근육 재강등, 조직 패턴 문장 (2026-09-20; 원고 v13)

**사용자 판단 세 가지에 대한 처리.**

1. **"25개 전부 같은 패턴이어야 의미 있는 것 아닌가."** 맞다. 그래서 확인했더니 **의미 있는 형태로 전부 성립한다**: 배양에서 증식을 양으로 따라가는 프로그램 **19개(ρ ≥ 0.3)는 19개 모두** 두 피부 부위에서 덜 따라간다. 순서를 벗어나는 4개(chaperone 0.17, collagen 0.05, ECM −0.21, NF-κB −0.57)는 애초에 배양에서 증식과 함께 오르지 않는 것들이다. 게다가 **유지율이 두 무리를 가른다**: 복제에 묶인 4개는 배양 결합의 71–75%를 유지하고(telomere 75.0, BER 72.3, DSB 71.8, cell cycle 71.2), 19개 중 나머지는 **54%를 넘지 못한다**(5위 mitochondrial translation 53.8). 이것이 내부 대조의 정량판이라 본문에 넣었다. 앞선 "Across all programmes the pattern is the same"(틀림)이나 "표를 보라"(무의미)보다 강하다. 계산은 `programme_by_tissue.tsv`의 `rho_prolif`만 쓰면 재현된다.
2. **"근육은 삭제? 강등?"** → **삭제 불가, 재강등**. 삭제할 수 없는 이유: 잠긴 사전등록 H3은 보고 의무가 있고(§10·§11.13, 편집자 구속 조건), 고찰의 근육 문장 하나가 ref 29(Ubaida-Mohien 근육 proteomics)를 유일하게 인용하므로 지우면 참고문헌 41개가 재번호된다. 그래서 **결과 2문단(결합 서술)에서는 완전히 들어내고 3문단(사전등록 결산)으로 옮겼다**. 이제 결합 이야기는 배양 → 피부만으로 흐르고, 근육은 "H3 기준은 충족했으나 그곳에서는 증식 점수가 세포주기조차 회수하지 못하므로 판정만 보고하고 해석하지 않으며 어디에서도 근거로 쓰지 않는다(사후 판단, 잠긴 판정 불변)" 한 문장으로 끝난다. 남은 곳: Fig. 6b heatmap, 보충 표 9·10, 고찰 한 문장, 한계 문단, Methods.
3. **"Fig. 7은 글 상자뿐이고 7b는 왜 있는지 모르겠다."** → **패널 b 삭제, 전체를 graphical abstract로 재작도**(`make_sa_fig7summary.R`). 세 열 각각이 *연구 대상의 그림 + 그 대상에서 한 측정의 실제 산점도*다: ① 배양 접시(분열 중인 세포 포함) + 센 속도 대 splicing 점수(328개, ρ = 0.71), ② 젊은/고령 공여자 접시 + 연령 대 splicing(107명, r = −0.50) + 70%/잔여 막대, ③ 피부 단면(표피 keratinocyte·진피 섬유아세포·collagen) + **같은 크기·같은 축의 두 상자**(세포주기 0.67 위, splicing 0.23 아래). 셀·접시·조직은 ggplot polygon으로 직접 그렸다(ggforce 없음).
   - **피부 두 상자는 기술 공변량의 표준화 잔차로 그린다.** 원자료로 그리면 두 프로그램의 *기울기*가 0.167 대 0.256으로 비슷해 보여(splicing 점수의 산포가 더 커서) 그림이 본문 수치와 어긋난다. 표준화 잔차로 그리면 **그린 선의 기울기 = 인용한 편상관**이 되고 렌더 전에 `chk()`로 대조한다. 범례에 명시.
   - §21의 "자료 없는 flowchart 금지"에 대한 예외는 이제 필요 없다. Fig. 7의 모든 패널에 실제 점과 선이 있다.
   - 7b(모든 자료의 결합 띠)를 지우면서 거기서만 살던 수치 두 개(비선택 Reactome 0.63은 Fig. 5a에 있음, 공여자 코호트 Spearman 0.91은 본문이 Pearson 0.89를 쓰므로 불필요)가 사라졌다. 고찰 문장은 논문이 실제로 보고하는 네 값(0.71 / 0.89 / 0.86 / 0.74)을 열거하는 형태로 바꿨다. `fig7_coupling_strip.tsv`는 더 이상 쓰지 않는다.
   - **`coord_fixed(ratio = 1)` + 캔버스 0–180 × 0–112 + `save_fig(..., 182, 114)`** 로 원이 원으로 보인다. 7 pt에서 한 줄은 **44자**를 넘기지 않는다(52 mm 열 폭). 넘기면 옆 열을 침범한다.

## 23. v11 심사 워크플로와 v12 수정 (2026-09-19, 사용자 요청 "리뷰어로서 한번 더 비판적으로")

리뷰어 7명(논리·그림·수치·한영·주장 경계·교차참조·심사위원) + 반박 검증자 7명. 61건 중 52건 확정, 9건 판단 필요, 반박 0건. 전부 v12에 반영했다. 되풀이하지 말아야 할 것만 적는다.

- **v11에서 새로 넣은 문단이 가장 큰 오류원이었다.** 고찰 "독립 자료" 문단의 '시료를 공유하지 않는다', 'batch·코호트·파이프라인을 공유하지 않는다'는 본문 Methods와 모순(GSE179848이 Fig. 1c·2·4의 12개 대비·5a에 들어가고, 공여자 코호트는 Fig. 1b·3·4·5b·6c에, GSE113957과 GTEx는 둘 다 recount3). 지금은 '겹침을 밝힌 자료 + 두 도구 + 기대의 사슬 + 전부에 공통인 batch는 없음'으로 쓴다. "cannot be an artefact common to them all" 같은 문장은 쓰지 않는다.
- **'calibrated'를 쓰지 않는다.** 20개 marker 점수는 계수에 대해 *비교*(ρ = 0.63)했을 뿐 재척도·재가중하지 않았다. 'compared with / checked against counting'. 한글은 '계수와 비교한'. 그리고 오차를 안은 공변량이 잔여를 남긴다는 문장(측정 오차 → 보정 불완전 → 잔여는 그 자체로 근거가 아님)을 3절에 넣었다.
- **"not short of power"·"specific to the programme"(피부 내부 대조)은 한계 문단과 모순이었다.** 복제 프로그램은 점수가 지표로 삼는 분열 세포에서 발현되므로, 회수된다는 것은 분열 세포 비율 차이가 분해된다는 뜻이지 모든 세포에서 발현되는 프로그램에 대한 검정력 증명이 아니다. 그렇게 쓴다.
- **구성 보정은 '전체 프로그램'이 아니라 6개 집합에만 돌렸다**(`composition_adjusted_all_programmes.tsv` 24행 = 6 × 4). 'specifically' 주장 금지. 섬유아세포 marker에 COL1A1·COL1A2가 들어 있어 collagen formation은 독립 점검이 아니고, 배양에서는 같은 보정이 collagen을 0.05 → 0.42로 움직인다(원고에 '배양의 splicing 결합은 불변'으로 한정).
- **"recovers no programme"(근육)은 표와 모순** — BER 0.30, collagen 0.24, ECM 0.22가 FDR < 0.05. 맞는 문장은 '유사분열 세포주기 프로그램을 회수하지 못한다(0.03, P = 0.38)'. H3 '지지' 판정을 본문에 명시하고 근육 강등은 사후 판단이라고 적었다.
- 'Across all programmes the pattern is the same'은 25개 중 21개만 맞다 → 표를 가리키는 문장으로. 'weakens only in tissue' 금지(근육은 검정 불가). '63 to 78%'는 Fig. 3b·표 1이 63–73%이고 78%는 depth-only(표 3). 하지 피부 CI는 plain Fisher z로 **0.16–0.29**(0.15는 1.06/(n−3−k) 분산).
- **묵은 인용이 또 남아 있었다**: 'panels c and d'(Fig. 1은 a–c), 두 번째 'Fig. 2d'(→ 2e), 보충 Results EN/KO의 S4·S5a,c·S3a,c, Note 1의 Fig. 6d·S7b, Note 2의 S9 ×3, 커버레터의 S3c·S4a, 빌더 머리말 'S1 to S10, Tables 1 to 9', 빌더 표 10의 Fig. 6c/6e/S10c. **그림·패널 번호를 바꾸면 원고 2 + 보충 Results 2 + Note 1·2 + 빌더 + 커버레터 2 = 9개 파일을 돈다.**
- **철회한 pooled 추정치가 Methods 제목·Note 1 항목 9·빌더 표 5(d)·한글 커버레터에 남아 있었다.** 이제 Note 1에 철회 문장을 append(append-only), 표 5(d)는 두 코호트만(combined 행 필터), 제목은 '두 코호트의 잔여 연령 효과'.
- **한글판에 v9 잔재 '재현되는 것은 보정 후 추정치다'가 살아 있었다**(바로 다음 문장과 모순). 영문에서 지운 주장은 한글에서 grep으로 확인한다.
- **Supplementary Data 목록의 설명 8/10이 실제 파일과 달랐다**(Data 5–10이 서로 다른 표를 가리킴). 설명을 파일 헤더대로 다시 썼고, 보충 표 8의 'Data 6' 인용은 저장소 경로로 바꿨다. 파일은 그대로다.
- Fig. 7: 7a의 '(H5)'는 배양 검정이므로 피부 열에서 빼고 배양 열에 적었다, 'growth-promoting media'·'MYC targets' 삭제, reprogramming은 '어느 쪽으로도', 7b에 연구당 한 값(◆ 0.79, n = 29)과 축 라벨(행마다 비교 대상이 다름) 추가. **Fig. 7b의 행들은 하나의 추정량이 아니다**라고 범례에 적었다.
- Fig. 1c(71개)와 2c(70개)의 라이브러리 수 차이는 HC2 164일째 라이브러리(Sample_328)에 속도가 없어서다 — 범례에 적었다.
- Methods에 빠져 있던 것 추가: 대비 점수 정의, 63개 대비의 Spearman·정의 series 제외 재적합·무작위 집합 1,000개·일표본 t, 41개 계열 집계, Thorndike case II, 'Results 2–6절 첫머리의 기대는 저자의 추론'이라는 공개 문장.
- `pdfFonts()` 등록은 명명 인자여야 한다(`do.call(pdfFonts, setNames(list(Type1Font(...)), FONT))`). 앞 판의 등록은 조용히 실패해 스크립트가 시작하자마자 멈췄고, 렌더 확인을 md5 불변으로 오판했다. **스크립트를 고친 뒤에는 PNG의 mtime이 바뀌었는지 본다.**

## 21. 그림 형식 정리 (2026-09-19, 사용자 요청)

**왜.** 사용자 지적: Fig. 6d·S2c 같은 형태는 "AI가 만드는 대표적 바 형태"라 그 자체로 reject 사유가 되고, 다른 논문에서 안 쓰는 형태가 너무 많다. splicing·노화 논문(Harries 2011 *Aging Cell*, Holly 2013 *MAD*, Lee 2016 *Aging Cell*, Latorre 2017, naked mole-rat *GeroScience* 2019)의 그림 어휘로 맞췄다: **산점도+회귀선, 그룹별 box, 오차막대 막대그래프, 시료 주석 heatmap, 계대별 line, volcano.**

**바꾼 것(58 → 52패널).**
- **삭제 5**: flowchart 도식 Fig. 1a·6a·S3a(자료 없음), 동어반복 meta-scatter Fig. 2c(ρ 대 ρ, R1 지적)·S2c(유지 비율 대 결합도, R2 "독립 검정 아님").
- **텍스트 컬럼 forest → 오차막대 막대** (13): Fig. 2e(옛 e)·3b·4c·6c(옛 d)·S1c·S3b·S5a·S5b·S5c·S6b·S7d·S8a. 추정치·*P* 컬럼은 전부 보충 표로 보냈다.
- **롤리팝/덤벨 → 막대·heatmap** (6): Fig. 2b(거명 인자, † = 후기 계대에서 증가 보고)·5a(23×2 **heatmap**, 값 인쇄, 행 이름 색 = 분류)·5d·6d(옛 e)·S2c(옛 d, 짝 막대)·S6a(막대, 세포주기 비율은 이름에)·S8b.
- 유의성 표기 통일: 채운/빈 기호 대신 **막대 색(파랑/빨강) 또는 40% 불투명도**.

**패널 번호 변경.** Fig. 1: b→a, c→b, d→c. Fig. 2: d→c, e→d, f→e. Fig. 6: b→a, c→b, d→c, e→d. S2: d→c. S3: b→a, c→b, d→c, e→d. 본문·범례·보충 Results 영/한 전부 단일 패스로 재번호했고 검증 통과(삭제 패널 인용 0, 범례 문자 연속, 인용 패널이 범례에 존재).

**함께 잡은 것.** `supplementary_results_en/ko.md`가 **어제(9/18)의 보충 그림 10 → 8 재번호를 받지 않은 채** 옛 번호(S4b, S5a·b, S6, S7, S8, S10c)를 인용하고 있었다. PDF 빌더가 이 파일을 그대로 넣으므로 보충 PDF의 Results 절 그림 번호가 전부 틀려 있었다. 오늘 고쳤다. **그림 번호를 바꿀 때는 원고 2개 + 보충 Results 2개, 네 파일을 함께 돌린다.**

**규칙.** 새 패널을 만들 때 롤리팝, 덤벨, 추정치/*P* 텍스트 컬럼이 붙은 forest, 자료 없는 flowchart는 쓰지 않는다. 그런 정보는 막대+95% CI로 그리고 수치는 보충 표에 둔다.

## 20. 두 표를 위치로 join하지 않는다 (2026-09-18, 실수 기록)

정의 series 3개를 뺀 검증값을 처음에 **ρ = 0.881, *R*² = 0.940**으로 계산해 원고에 넣었다. 틀렸다. 맞는 값은 **ρ = 0.841 [0.733, 0.907], *R*² = 0.929**(n = 49, 20개 series)다.

원인: `revision_stats/fig6_contrasts_revised.tsv`(id 순 정렬)와 `figures_sciadv/supp_table1.tsv`(class·series·id 순 정렬)를 **행 위치로 zip** 해서 대비마다 엉뚱한 accession이 붙었다. 두 파일은 정렬이 다르다.

**규칙.** 대비를 series나 accession에 붙일 때는 `make_supp_tables.R`이 쓰는 경로를 따른다: `id`에서 `^(GSE[0-9]+)`를 뽑고, 그게 없는 secretome 대비는 `secretome_class/secretome_signatures.tsv`의 `dataset` 열로 `id`를 match한다. `supp_table1.tsv`에는 id 열이 없으므로 그 파일로는 join하지 않는다.

리뷰어 2가 같은 검정에서 0.840을 보고했는데 내가 0.881을 내놓고 "방향은 같으니 재확인 필요"로 넘겼다. **독립 계산이 어긋나면 어느 쪽이 맞는지 정하기 전에는 원고에 넣지 않는다.**

## 19. 재심사 패널과 편집자 결정 (2026-09-18)

리뷰어 4명 전원 **Major revision**, reject 없음. 사전등록 해시·그림 주석·헤드라인 수치 재계산에서 오류 0건. 결정문 아티팩트 https://claude.ai/artifact/MNLXonfYkHcHbq1rpt953K (2026-09-19 재게시; 옛 `7H63hPbz…`는 이전 계정 소유)

**핵심 판정.** 기여 (ii)("두 코호트에서 재현되는 잔여 연령 효과")는 **철회**한다. GTEx 배양의 −0.053은 사전등록이 "평가 불가"로 선언한 모형에서 나오고(감쇠 −1.098), 25개 프로그램 중 15개가 같은 보정 후 음의 연령 효과를 얻으며, **세포주기 양성 대조 자신이 +0.0064(P = 0.748) → −0.0208(P = 0.0018)**이 된다. 사전등록 Note 1이 이미 해석 불가로 정해 둔 비교다. 세 기여 중 (i) 유지, (ii) 삭제, (iii) 좁혀서 유지.

**반영 완료(BC 1·2 + 리뷰어 1 권고).** pooled 추정치 전면 삭제(초록·Results·고찰·Methods·Fig. 6d·커버레터), suppression 명시와 분모 보고, specification 범위(−0.040 → −0.016~−0.022), 피부 내부 양성 대조 승격과 Fig. 6 패널 e 신설, 근육 판별 불가 재분류(**잠긴 H3 판정은 건드리지 않는다**), Fig. 2 패널 f 신설, 스코어카드를 보충 표 9로 이동(Fig. 6 243 → 196 mm), 보충 표 9·10 번호 교환.

**남은 구속 조건.** 3 선택·재정의 이력 공개(정의 series 3개 뺀 ρ = 0.84 [0.73, 0.91], R² = 0.93 포함 — 저자에게 유리), 4 대비를 series 단위로, 5 산출물과 모순되는 문장 3개, 6 **코드 저장소와 Supplementary Data 실제 제작**.

**편집자가 명시적으로 면제한 것(하지 않는다).** GTEx 혈액·고회전 상피 추가(새 다운로드 필요), epiTOC2/RepliTali, 20개 marker 점수의 split-half 신뢰도, GTEx에서의 무작위 집합 1,000개, 근육 분산 맞춤, GSE179848 밖에서 더 나은 증식 측정치 도출, Cochran *Q*·*I*².

**제목 확정:** *Splicing-factor expression in cultured human fibroblasts reports proliferative state and the coupling weakens in **skin**.* "in tissue"는 판별 가능한 조직이 하나뿐이라 지지되지 않는다.
