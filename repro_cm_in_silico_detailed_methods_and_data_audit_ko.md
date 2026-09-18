# Repro-CM 후속 in silico 연구: 상세 Methods와 데이터 감사표

작성일: 2026-09-10  
상태: 이 문서는 2026-09-10의 승인 전 설계본이다. Tier A 실행 후 확정된 Methods는 `repro_cm_in_silico_methods_v2_ko.md`, 실제 결과와 스토리 수정은 `repro_cm_tierA_execution_and_story_revision_ko.md`, Tier B 크기 재감사는 `tierB_raw_download_decision_ko.md`를 우선 참조한다.  
대상 Figure: 새 후속 논문의 Figure 1–6

## 1. 분석 전제와 primary estimand

### 1.1 연구 단위

로컬 mRNA-seq의 네 조건은 각각 biological replicate n=1이다.

- UVA 직후 0 h
- UVA 후 HDF-CM 24 h
- UVA 후 Repro-CM 24 h
- UVA 후 iPSC-CM 24 h

따라서 로컬 자료만으로 biological dispersion, p-value, confidence interval 또는 FDR을 추정할 수 없다. 로컬 자료의 역할은 Repro-CM 연관 effect-rank와 splicing candidate를 만드는 것이다. 통계적 추론과 일반화는 biological replicate가 있는 독립 공개자료에서 수행한다.

### 1.2 primary comparison

주 비교는 같은 24 h 시점의 다음 두 비교다.

- Repro-CM 대 HDF-CM
- Repro-CM 대 iPSC-CM

15J_0hr는 time-zero context로만 사용한다. 0 h와 24 h의 차이는 시간 경과, UVA 후 자연 회복, 배지 노출 및 CM 효과가 모두 섞여 있으므로 treatment effect가 아니다.

### 1.3 primary outcome

primary outcome은 “Repro-CM 연관 gene-rank가 독립 UVA injury signature와 반대 방향인지”이다. secondary outcomes는 다음이다.

- UVA protective/rescue signature와 같은 방향인지
- chronological age 및 senescence signature와 반대인지
- fibroblast identity를 보존하는지
- splicing/RBP module이 독립 자료에서 재현되는지
- extracellular network model이 held-out response를 예측하는지
- mammalian mortality/longevity signature와 분자적으로 수렴하는지

## 2. 로컬 원자료 감사

### 2.1 FASTQ inventory

모든 파일은 paired-end 2×76 bp로 확인되었다.

| sample | R1 크기, byte | R2 크기, byte | mate당 record 수 | full gzip integrity | 판정 |
|---|---:|---:|---:|---|---|
| 15J_0hr | 654,998,717 | 667,798,085 | 14,978,682 | R1 pass, R2 pass | 사용 가능 |
| HDF_CM_24h | 735,267,786 | 752,453,493 | 16,715,203 | R1 pass, R2 pass | 사용 가능 |
| Rep_CM_24h | 678,923,176 | 692,243,768 | 15,442,649 | R1 pass, R2 pass | 사용 가능 |
| iPSC_CM_24h | 743,929,659 | 756,626,916 | payload상 16,911,181 | R1 fail, R2 pass | 현재 사용 금지 |

중요: iPSC_CM_24h.R1.fastq.gz는 끝부분의 gzip CRC 검증에서 실패했다. decompress되는 payload의 line 수가 R2와 같더라도 CRC failure가 있으면 파일이 원본과 동일하다고 보장할 수 없다. 이 파일은 원래 sequencing 전달본에서 정확히 다시 복사하거나 재전송받고, 제공자 checksum과 비교해야 한다.

### 2.2 기존 preflight의 사각지대

현재 repro_cm_splicing_ubuntu 설정은 각 FASTQ의 앞 10,000 record만 읽어 구조와 pair를 점검하고 전체 파일의 SHA-256을 계산한다. SHA-256은 현재 파일을 식별할 뿐, 그 파일이 정상 gzip stream인지 판정하지 않는다. 뒤쪽 CRC 오류는 앞 10,000 record scan으로 잡히지 않는다.

따라서 workflow 시작 전에 반드시 다음을 별도 필수 gate로 추가한다.

1. 모든 FASTQ 전체 gzip 또는 pigz integrity test
2. R1/R2 전체 record 수 일치
3. 가능하면 전체 read ID의 mate pairing 일치
4. 제공자 checksum과 로컬 SHA-256 비교
5. 파일 교체 후 immutable manifest 재생성

### 2.3 library metadata의 미확정 항목

JTE proof에는 TruSeq RNA sample preparation과 mRNA purification은 기재되어 있으나 strandedness는 명확하지 않다.

- sequencing core의 원래 sample sheet 또는 delivery report가 최우선
- 없으면 RSeQC infer_experiment로 추정
- 불명확한 경우 gene-level은 unstranded로 처리
- splicing junction count는 stranded와 unstranded 설정의 민감도 분석

추정 결과를 확정 metadata처럼 서술하지 않는다.

## 3. 내부에서 먼저 찾아야 하는 자료

공개자료보다 아래 내부자료 확인이 우선이다.

| 필요한 내부자료 | 이유 | 없을 때 영향 |
|---|---|---|
| 정상 iPSC_CM_24h R1 원본과 제공자 checksum | time-matched comparator 완성 | Repro-specific signature와 splicing 분석 중단 |
| mRNA-seq sample sheet, lane, strandedness, insert size | 정확한 정렬과 batch 추적 | splicing 신뢰도 하락 |
| mRNA-seq GEO/SRA/ENA accession | 후속 논문의 투명성 및 저널 요건 | Aging Cell 투고에 치명적 |
| extracellular miRNA raw FASTQ 또는 raw count matrix | Figure 4 signed miRNA network | miRNA network를 exploratory로 축소 |
| miRNA sample별 biological replicate 및 batch | 차등 abundance 추론 | p-value/FDR 사용 불가 |
| antibody-array 원 spot intensity, blank, duplicate, exposure 정보 | 단백질 seed 가중치와 QC | TIMP2 포함 unweighted list만 가능 |
| CM batch, collection day, pooling metadata | biological unit와 source-state 재현성 | CM specificity 해석 약화 |
| phenotype 원수치와 replicate ID | 선행논문 효과 크기와 omics 연계 | 새 주분석에는 필수 아님, 재현성·보조분석 제한 |
| 최종 JTE Data Availability 문구 | 데이터 재사용 disclosure | 중복 출판 및 repository 문제 |

PDF 그림에서 spot intensity를 역추출하거나 heatmap 색을 정량값으로 취급하지 않는다.

## 4. 공개 데이터 다운로드 후보: 승인 전 최종 목록

### 4.1 Tier A — 먼저 받을 소용량 processed 자료

Tier A는 gene-level 방향과 연구 설계의 go/no-go를 확인하기 위한 최소 세트다. 현재는 목록만 확정했으며 다운로드하지 않았다.

| 우선 | accession/resource | 설계와 표본 | Figure | 받을 파일 | 예상 압축 용량 | 핵심 한계 |
|---:|---|---|---|---|---:|---|
| 1 | [GSE125429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125429) | 4 HDF strain, chronic UVA1 대 paired control, 8 arrays | 1 | log2/linear processed matrix | 1.9 MB | array, chronic high-dose 조건 |
| 2 | [GSE240226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240226) | HDF control/UVA/UVA+Maifuyin/UVA+succinic acid, 총 12 RNA-seq | 1 | sample raw-count/RPKM tables와 metadata | 약 6–20 MB | Repro-CM이 아닌 별도 보호 처치 |
| 3 | [GSE302943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302943) | UVA 대 control 각 n=5 RNA-seq | 1, held-out | 10개 gene-count 파일 | 약 1 MB | BJ-5ta hTERT 불멸화, 누적 UVA |
| 4 | [GSE89005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89005) | sham n=3, UVA/UVB 각 n=2, single/repeated, 6/24 h | 1 | normalized sample tables; 필요 시 RAW array archive | matrix는 소용량, RAW 49.5 MB | array, irradiated arm n=2 |
| 5 | [GSE109700](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109700) | proliferating/early/deep senescence 각 n=3 | 2 | count matrix | 0.43 MB | lung LF1, dermal HDF 아님 |
| 6 | [GSE191055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191055) | HDF P4 n=4, P27 n=4 | 2 | P27 대 P4 gene table | 1.4 MB | 별도 young/old skin 각 n=1은 제외 |
| 7 | [GSE226189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189) | 22–89세 primary skin fibroblast 82명 | 2, 4 | 82개 geneCOUNT와 metadata | 약 18 MB | hg19, 연속 연령·배양 covariate 점검 필요 |
| 8 | [GSE113957](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113957) | 건강인 133명, 1–94세; HGPS 10명 | 2 | processed FPKM matrix와 phenotype | 7.2 MB | single-end stranded; HGPS는 별도 분석 |
| 9 | [GSE179848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179848) | 345 longitudinal fibroblast samples, control donor 4명과 SURF1 donor 3명 | 2, 4 | raw-count matrix 52.0 MB와 sample metadata | 약 52 MB | sample 345개가 독립 n이 아님; donor random effect 필수 |
| 10 | [GSE165177](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165177) | 3 donor와 여러 시점/대조의 MPTR | 2 | 두 processed log2-RPM table | 15.3 MB | 복잡한 sample labeling, subset 사전 정의 필요 |
| 11 | [GSE93535](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93535) 및 논문 supplement | quiescent/SIPS HDF, 각 조건 ± plant extract | 2, held-out | 우선 논문의 differential-expression supplementary table과 metadata | 수 MB 수준, 승인 직전 확인 | GEO는 raw 중심이며 rescue가 특정 추출물에 한정 |
| 12 | [Cell Metabolism 2019 supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC6907080/) | 17 lifespan intervention, 77 comparison | 5 | intervention effect/signature tables | 소용량 | 일부 신규 실험 포함, 피부 특이 아님 |
| 13 | [Nature 2026 supplements](https://www.nature.com/articles/s41586-026-10542-3) | 11,000+ transcriptomes, 25+ tissues, 4 mammals | 5 | aging/mortality/lifespan 계수와 module table | 소용량 표 중심 | 원 연구 tissue composition 반영 |
| 14 | curated network tables | OmniPath, Reactome, miRTarBase, Ensembl ortholog, ENCODE eCLIP summary | 3–5 | versioned TSV/bed 또는 API snapshot | 승인 직전 계산 | source별 license와 version 기록 필요 |

Tier A의 예상 총량은 대략 150–250 MB 범위다. GSE93535 및 network snapshot의 실제 파일 구조를 metadata 단계에서 확정한 뒤 정확한 download manifest와 byte 수를 다시 제시한다.

### 4.2 Tier B — Tier A 결과가 통과한 뒤에만 받을 raw RNA-seq

raw FASTQ는 exact splicing과 동일 pipeline 재처리에 필요하지만 대용량이다. 먼저 SRA RunInfo만 확인해 paired/single-end, read length, spot 수, 총 byte를 확정한다.

| 순서 | accession | raw를 받는 이유 | 다운로드 범위 | 현재 추정/주의 |
|---:|---|---|---|---|
| 1 | GSE191055 / PRJNA789540 | P4 대 P27 splicing, 각 n=4 | HDF 8 sample만; skin n=1 두 개 제외 | layout와 총 byte 승인 직전 확인 |
| 2 | GSE109700 / PRJNA431791 | early/deep senescence splicing, 각 n=3 | 9 sample | layout와 read length 확인 |
| 3 | GSE302943 / PRJNA1292871 | n=5+5 UVA splicing held-out | 10 sample | 불멸화 세포이므로 보조 검증 |
| 4 | GSE240226 / PRJNA1003002 | UVA injury와 두 rescue의 splicing | 12 sample | 약 24 GB로 추정, processed signal 후 진행 |
| 5 | GSE165177 / PRJNA693500 | MPTR time-course splicing | 사전 선택한 donor/time-point subset 우선 | 전체 raw 약 45.7 GB 수준, 95 sample 전체는 후순위 |
| 6 | GSE179848 / PRJNA745334 | longitudinal replicative-aging splicing | healthy control donor의 선택 시점 우선 | 345개 2×150 bp 전체는 매우 큼 |
| 7 | GSE226189 | 연령 연속축의 splicing | 82 sample 전체 또는 극단 연령 subset | raw 대용량; gene-level 결과가 강할 때만 |
| 8 | GSE307377 | 4 young 대 5 aged donor held-out | 9 sample | GEO series에 processed matrix가 불명확하여 raw 가능성 높음 |
| 9 | GSE93535 | SIPS 및 rescue splicing | design 확인 후 관련 12 sample | library layout와 파일 구조 확인 필요 |

Tier B는 한꺼번에 받지 않는다. 각 단계에서 exact event 또는 module validation 가능성을 확인한 뒤 다음 dataset으로 넘어간다.

### 4.3 Tier C — 선택적 또는 제외

| 자료 | 판정 | 이유 |
|---|---|---|
| GSE274955 | 핵심 분석 제외 | 연령별 donor 한 명, 두 부위가 같은 donor; cell 수를 biological replicate로 세면 pseudoreplication |
| GSE203034 | exploratory만 | 한 남성 donor를 36.2세와 72.8세에 측정 |
| GSE130973 | 보조 sensitivity만 | young donor 2, old donor 3의 skin scRNA-seq; spatial 분석은 하지 않고 donor-level pseudobulk만 가능 |
| GSE125429 | splicing 제외 | microarray |
| GSE89005 | splicing 제외 | microarray |
| GSE63577 | 낮은 우선순위 | 50 bp single-end 다수 sample; exact isoform 해상도 낮음 |
| GTEx skin | 보조만 | tissue cell composition이 cultured HDF와 다름 |
| 전체 LINCS L1000 | 초기 다운로드 제외 | 대용량이며 후보가 정해지기 전에는 탐색 자유도만 증가 |
| 임의 aging GWAS/AlphaGenome | 제외 | CM response 질문과 직접 맞지 않음 |
| spatial transcriptomics | 제외 | 단층 HDF/공개 bulk 중심 질문과 불일치 |

## 5. 사전등록과 discovery–validation 분리

분석 시작 전에 protocol 문서를 tag하고 변경 이력을 남긴다.

### 5.1 discovery set

- UVA: GSE125429, GSE89005, GSE240226의 control 대 UVA injury arm
- chronological age: GSE113957, GSE226189
- senescence: GSE109700, GSE191055
- partial reprogramming: GSE165177의 사전 정의된 donor와 time course
- network training: GSE226189 cross-sectional coexpression

### 5.2 held-out set

- UVA: GSE302943
- senescence rescue: GSE93535
- aging: GSE307377
- network robustness: GSE179848에서 donor/culture-day effect를 제거한 network
- partial reprogramming: metadata 검토 후 GSE237269 또는 GSE297233 중 하나

held-out dataset의 outcome은 후보 ranking 및 parameter tuning에 사용하지 않는다.

### 5.3 protocol deviation

dataset 교체, cutoff 변경, sample 제외는 이유와 시간을 기록한다. 결과를 본 뒤 바꾼 분석은 confirmatory가 아니라 sensitivity/exploratory로 표시한다.

## 6. M1 — FASTQ QC 및 reference

### 6.1 immutable inputs

각 파일에 대해 다음 manifest를 만든다.

- absolute source path
- file size
- SHA-256
- gzip integrity
- record 수
- read length distribution
- first/last read ID hash
- sample, condition, mate, lane
- acquisition source와 날짜

원본 FASTQ는 read-only로 두고 workflow output과 분리한다.

### 6.2 reference

- genome: GRCh38 primary assembly
- annotation: [GENCODE human release 50](https://www.gencodegenes.org/human/), GRCh38.p14
- FASTA와 GTF의 exact URL, release, MD5/SHA-256 기록
- ERCC 사용 여부는 sample sheet에서 확인
- gene identifier는 Ensembl stable ID를 기준으로 하고 version suffix를 별도 보존

선행논문은 GRCh37/GENCODE v19를 사용했지만, 새 분석은 모든 raw validation cohort를 같은 GRCh38/GENCODE v50으로 재처리한다. 선행 결과와의 차이는 reference update sensitivity로 보고한다.

### 6.3 read QC와 trimming

- FastQC 및 MultiQC
- fastp로 adapter와 low-quality tail 처리
- 76 bp read이므로 aggressive trimming 금지
- trimming 전후 read 수, Q30, adapter fraction, length distribution 보고
- read 길이가 너무 짧아진 pair는 사전 설정 길이 아래에서 제거

### 6.4 alignment와 quantification

- STAR two-pass alignment
- coordinate-sorted BAM, splice junction table, gene count 생성
- featureCounts 또는 STAR GeneCounts로 gene count
- Salmon selective alignment를 transcript abundance sensitivity로 병행
- Picard/RSeQC로 insert size, duplication, gene-body coverage, junction saturation, rRNA fraction, strandedness 평가

QC threshold는 자동 삭제 기준과 경고 기준을 분리한다.

- gzip integrity 또는 mate mismatch: 자동 중단
- uniquely mapped read가 매우 낮음, rRNA가 높음, severe 3-prime bias: 원인 조사 후 exclusion 여부를 blind하게 결정
- sample n=1이므로 불리한 결과라는 이유로 sample을 제거하지 않음

## 7. M2 — 로컬 Repro-CM 연관 gene-rank

### 7.1 표현량

TMM-normalized CPM을 사용하고 다음처럼 변환한다.

\[
x_{g,s} = \log_2(CPM_{g,s} + 0.5)
\]

gene g가 24 h 세 sample 중 최소 두 sample에서 CPM 1 이상일 때 주분석에 포함한다. CPM 0.5와 2를 sensitivity cutoff로 사용한다.

### 7.2 두 time-matched contrast

\[
d_H(g)=x_{g,Repro}-x_{g,HDF}
\]

\[
d_I(g)=x_{g,Repro}-x_{g,iPSC}
\]

두 차이의 부호가 같을 때만 Repro-specific score를 부여한다.

\[
S(g)=
\begin{cases}
sign(d_H(g)) \times \min(|d_H(g)|, |d_I(g)|), & sign(d_H)=sign(d_I) \\
0, & otherwise
\end{cases}
\]

이 점수는 한 comparator와의 큰 차이만으로 후보가 상위에 오는 것을 막는다. p-value는 계산하지 않는다.

### 7.3 robustness

- pseudocount 0.1, 0.5, 1.0
- expression filter CPM 0.5, 1, 2
- TMM-CPM 대 Salmon TPM
- effect cutoff fold change 1.25, 1.5, 2.0
- top 100, 250, 500, 1000 gene rank
- read down-sampling 및 Salmon bootstrap

bootstrap은 sequencing/quantification uncertainty만 반영한다. biological confidence interval로 부르지 않는다.

### 7.4 0 h sample 사용

15J_0hr는 PCA/trajectory 그림에서 위치를 보여주는 데 사용할 수 있지만, 24 h untreated UVA가 없으므로 natural recovery와 CM effect를 분리하지 못한다. 0 h와의 차이를 치료 효과나 회춘 효과로 사용하지 않는다.

## 8. M3 — 공개 코호트 전처리와 study-level effect

### 8.1 절대 원칙

서로 다른 연구의 normalized expression matrix를 한 행렬로 합쳐 case/control DE를 수행하지 않는다. 플랫폼과 batch가 biological group과 완전히 섞일 수 있기 때문이다.

각 연구 안에서 effect를 먼저 추정하고, gene별 effect 또는 signed statistic을 이후 meta-analysis한다.

### 8.2 RNA-seq

- raw count가 있으면 edgeR quasi-likelihood, DESeq2 또는 limma-voom 중 사전 선택한 하나를 primary로 사용
- paired donor는 design에 donor term 포함
- batch, sex, passage, treatment, time이 있으면 사전 정의된 covariate 사용
- low-count filtering은 group label을 보지 않고 수행
- shrinkage log fold change와 signed z 또는 moderated t를 저장

예시 UVA model:

\[
expression \sim donor + UVA
\]

예시 rescue model:

\[
expression \sim treatment
\]

필요 contrasts:

- UVA 대 control
- UVA+rescue 대 UVA
- rescue-only가 있으면 rescue-only 대 control을 별도 평가

### 8.3 microarray

- 원 probe intensity가 있으면 platform-appropriate background correction과 quantile/RMA 계열 정규화
- limma로 paired 또는 unpaired model
- 여러 probe가 한 gene에 매핑되면 outcome을 보지 않고 median effect를 사용
- uniquely mapped probe만 쓰는 sensitivity 분석

### 8.4 chronological age

가능한 기본식:

\[
expression_g \sim age + sex + batch + passage
\]

metadata가 없는 covariate를 임의로 대입하지 않는다. age coefficient를 10년당 effect 또는 standardized age effect로 변환한다.

GSE179848의 primary analysis는 healthy control donor의 longitudinal trajectory다.

\[
expression_g \sim culture\_day + (1 + culture\_day \mid donor)
\]

기술 반복과 treatment arm을 biological replicate로 세지 않는다. SURF1 arm은 별도 interaction 분석으로 둔다.

### 8.5 senescence

- GSE191055: P27 대 P4, n=4+4
- GSE109700: early 또는 deep senescence 대 proliferating, 각 n=3
- deep와 early를 합치지 않고 progression을 별도 contrast로 둠
- GSE93535: SIPS 대 quiescent, extract rescue 대 SIPS

## 9. M4 — Figure 1 외부 UVA concordance

### 9.1 primary test

각 독립 UVA 연구의 모든 gene을 signed statistic으로 정렬하고 로컬 S(g)의 up/down set에 대해 preranked GSEA를 수행한다.

성공 방향:

- Repro-up gene은 UVA-down 쪽 또는 rescue-up 쪽에 enrichment
- Repro-down gene은 UVA-up 쪽 또는 rescue-down 쪽에 enrichment

두 방향을 하나의 directional concordance score로 요약하되 각각의 NES와 FDR도 공개한다.

### 9.2 보조 test

- Spearman correlation of gene-level signed effects
- RRHO 또는 rank–rank overlap
- top-k overlap의 hypergeometric test
- pathway-level NES correlation

top-k overlap만 유의하고 full-rank 결과가 약하면 강한 replication으로 해석하지 않는다.

### 9.3 null model

- expression decile, gene length, GC, number of annotated transcripts를 맞춘 random gene set
- Repro rank label permutation
- unrelated perturbation signature
- HDF-CM 및 iPSC-CM comparator-derived rank

최소 10,000 permutation으로 empirical p-value를 구하고 같은 signature family 안에서 BH FDR을 적용한다.

### 9.4 heterogeneity

각 study의 NES 또는 standardized concordance를 random-effects REML meta-analysis한다.

- pooled effect와 95% CI
- prediction interval
- Cochran Q
- I-squared
- leave-one-study-out

UVA dose, time, primary 대 immortalized cell을 moderator로 탐색하되 작은 연구 수에서 moderator p-value를 과신하지 않는다.

## 10. M5 — Figure 2 다축 상태 분해

### 10.1 gene-set freeze

분석 전에 다음 gene set의 source와 version을 고정한다.

- chronological fibroblast aging
- replicative senescence
- stress-induced senescence
- SenMayo/SASP
- fibroblast identity 및 ECM
- pluripotency
- E2F/G2M proliferation
- apoptosis
- p53/DNA damage
- interferon/inflammation
- oxidative stress
- autophagy/proteostasis

한 결과에 맞춰 gene set을 교체하지 않는다. Hallmark/Reactome 등 외부 gene set의 version과 이용 조건을 기록한다.

### 10.2 sample-level score

cross-platform 비교에는 rank-based singscore를 primary로 사용한다. GSVA/ssGSEA는 sensitivity로 사용한다. score는 각 연구 안에서 control 중심화한다.

\[
Z_{score} = \frac{score - mean(score_{control})}{SD(score_{control})}
\]

control n이 너무 작으면 z-score 대신 raw rank score와 exact permutation을 사용한다.

### 10.3 MPTR model

GSE165177의 sample name에서 donor, experiment, treatment success, intermediate/final state, day를 구조화한다. 분석 전에 sample inclusion table을 사람이 검수한다.

\[
score \sim treatment \times day + (1 \mid donor) + (1 \mid experiment)
\]

핵심은 Repro-CM이 MPTR과 완전히 같다는 것을 보이는 것이 아니라 다음을 분리하는 것이다.

- aging/senescence reversal component
- fibroblast identity retention
- pluripotency/dedifferentiation component
- proliferation component

### 10.4 three-axis paracrine rejuvenation index

일반화 가능한 방법론을 만들려면 세 독립 축을 유지한다.

- A: aging/senescence-opposing score
- I: fibroblast identity-preservation score
- D: dedifferentiation/proliferation penalty

단일 합성 점수는 secondary summary로만 쓰고, 세 원래 축을 항상 함께 표시한다. 가중치는 discovery 결과를 보고 최적화하지 않고 동일 가중 또는 외부 benchmark에서 사전 결정한다.

“rejuvenation-like” gate:

- A가 사전 설정 방향으로 유의
- I가 control 대비 의미 있게 감소하지 않음
- D가 상승해 A 전체를 설명하지 않음

절대 transcriptomic age의 연 단위 변화는 보고하지 않는다.

## 11. M6 — Figure 3 alternative splicing과 RBP

### 11.1 local discovery

- STAR two-pass junction
- rMATS-turbo의 statoff 모드로 event count/PSI 산출
- event type: SE, A5SS, A3SS, MXE, RI
- MAJIQ 또는 LeafCutter를 annotation-independent sensitivity로 사용
- Salmon transcript usage를 gene-level DTU sensitivity로 사용

local candidate는 다음을 모두 만족한다.

- Repro 대 HDF-CM과 Repro 대 iPSC-CM에서 동일 방향
- 두 비교 중 작은 absolute delta PSI가 0.15 이상
- informative junction read 합 20 이상
- 서로 다른 quantifier 또는 bootstrap에서 방향 안정
- multi-mapping 또는 paralog 문제가 없는 locus

### 11.2 public statistical validation

biological replicate가 있는 raw dataset은 같은 GRCh38/GENCODE v50 pipeline으로 정렬한다.

- rMATS 또는 LeafCutter group test
- BH FDR 0.05
- absolute delta PSI 0.10 이상을 최소 biological effect로 설정
- exact coordinate 또는 동일 transcript consequence의 방향 일치
- 최소 두 독립 study에서 방향 일치해야 main text event

플랫폼 때문에 exact event가 관찰 불가능하면 그 study를 negative replication으로 세지 않는다.

### 11.3 RBP triangulation

각 candidate exon에 대해 세 증거를 분리한다.

1. motif: upstream intron, exon, downstream intron의 RBP motif enrichment
2. binding: ENCODE eCLIP high-confidence peaks/targets
3. activity: RBP expression, regulon activity 또는 known position-dependent effect

motif background는 exon length, flanking intron length, GC, splice-site strength, expression을 맞춘다. 최소 두 종류의 증거가 같은 RBP를 지지해야 main candidate로 올린다.

### 11.4 functional annotation

- coding frame 유지 여부
- protein domain 포함/제외
- premature termination codon과 NMD 가능성
- conserved exon 여부

이 annotation은 기능 예측이지 기능 검증이 아니다.

## 12. M7 — Figure 4 fibroblast-context extracellular network

### 12.1 context network

primary network:

- GSE226189의 82 donor gene expression
- age, sex, batch 등 이용 가능한 covariate 효과를 제거한 residual
- robust correlation 또는 WGCNA topological overlap
- OmniPath/Reactome의 curated signaling edge로 graph를 제한 또는 가중

robustness network:

- GSE179848에서 donor와 culture-day 구조를 고려한 residual network
- 두 network에서 edge/module preservation 평가

단일-cell의 cell 수를 network n으로 사용하지 않는다.

### 12.2 secreted protein seed

원 antibody-array 수치가 있을 때:

- local background subtraction
- duplicate spot aggregation
- exposure saturation 확인
- Repro-CM 대 HDF-CM/iPSC-CM/mTeSR의 descriptive log-ratio
- biological replicate가 없으면 p-value 없음

ligand–receptor mapping은 curated source와 version을 기록한다. 수용 HDF에서 receptor가 발현되지 않으면 primary seed에서 제외하고 sensitivity에만 둔다.

### 12.3 miRNA edge

- miRTarBase의 strong experimental evidence target을 primary
- TargetScan prediction은 sensitivity
- miRNA 증가 → target 억제의 negative edge
- seed abundance가 없는 miRNA는 네트워크에 임의 추가하지 않음

### 12.4 propagation

receptor seed에서 personalized PageRank 또는 random walk with restart를 수행한다.

\[
p_{t+1} = (1-r)W p_t + r p_0
\]

- W: degree-normalized fibroblast-context network
- p0: ligand abundance와 receptor expression으로 가중한 seed vector
- r: restart probability

r과 edge threshold는 discovery dataset의 nested resampling에서 정하고 held-out dataset에서 고정한다. signed edge는 positive/negative channel을 분리한 diffusion 또는 signed network sensitivity로 구현한다.

### 12.5 평가

- gene-level AUROC와 AUPRC
- top-ranked gene의 관찰 response enrichment
- pathway NES overlap
- degree-matched random receptor seed
- edge-rewired network
- ligand label permutation
- leave-one-ligand-out
- bootstrap network

AUPRC는 positive gene 비율에 민감하므로 random prevalence baseline과 함께 제시한다.

### 12.6 Figure 4 go/no-go

- 원 secretome/miRNA 수치 없음: weighted network 중단
- model이 degree-matched null을 이기지 못함: mechanism figure에서 제거
- discovery에서는 맞고 held-out에서 실패: exploratory supplement로 이동
- TIMP2 제거 시에만 모델이 무너짐: 흥미로운 민감도이지 TIMP2 necessity 증거는 아님

## 13. M8 — Figure 5 cross-species aging/longevity

### 13.1 ortholog

- Ensembl release를 고정
- human–mouse one-to-one ortholog primary
- human–worm one-to-one ortholog secondary
- one-to-many mapping은 제외하거나 orthogroup sensitivity로 별도 처리
- 매핑 전후 gene 수와 selection bias 보고

### 13.2 signature

- Cell Metabolism 2019의 intervention-specific effect와 lifespan-extension meta-signature
- Nature 2026의 age, expected mortality, maximum lifespan, module coefficients

각 source가 제공한 coefficient 방향을 그대로 보존하고 tissue/sex/study label을 유지한다.

### 13.3 test

- signed preranked GSEA
- gene-effect Spearman correlation
- pathway-level concordance
- tissue-stratified random-effects meta-analysis
- leave-one-intervention-out
- ortholog, expression, length를 맞춘 random null

Repro-CM local rank를 C. elegans 선행 phenotype과 결합해 p-value를 만들지 않는다. 그 phenotype은 같은 선행 연구의 결과다.

## 14. M9 — Figure 6 후보 통합과 held-out 검증

### 14.1 nonredundant evidence

각 후보에 대해 다음 rank를 만든다.

- local Repro specificity
- UVA anti-concordance
- aging/senescence anti-concordance
- identity/dedifferentiation safety
- splicing/RBP
- extracellular network
- longevity/mortality

서로 강하게 종속된 두 지표는 한 evidence family로 묶어 중복 가중을 막는다.

### 14.2 rank aggregation

primary는 equal-weight robust rank aggregation 또는 rank product를 사용한다. arbitrary hand-tuned weight는 사용하지 않는다.

- dataset bootstrap
- 한 dataset씩 제외
- cutoff multiverse
- network source 교체
- ortholog rule 교체

각 후보의 top-k 포함 빈도를 stability-selection probability로 보고한다.

### 14.3 held-out evaluation

후보와 parameter를 고정한 뒤 held-out dataset에서 다음만 평가한다.

- 예상 방향의 gene-set enrichment
- top-k precision
- AUROC/AUPRC
- pathway-level concordance
- predeclared primary endpoint

held-out 결과를 본 뒤 후보를 다시 정렬하면 그 결과는 새 discovery이며 같은 자료를 validation으로 부르지 않는다.

### 14.4 선택적 LINCS/CMap

최종 후보가 정해진 뒤에만 다음을 조회한다.

- candidate gene overexpression
- candidate gene knockdown
- receptor/RBP 관련 compound

동일 perturbation의 replicate는 MODZ 등 robust consensus로 합치고, cell-line별 direction과 dose/time을 보존한다. HDF가 아닌 암세포주에서의 일치는 보조 증거다. 약물 후보는 target evidence, 독성, 승인 상태를 별도 표로 분리하며 계산 일치를 치료 가능성으로 해석하지 않는다.

## 15. M10 — 다중검정, 불확실성, negative controls

### 15.1 다중검정

- gene-level public DE: BH FDR
- splicing event: event family 내 BH FDR
- gene-set family: tested set 전체에서 BH FDR
- network permutation: empirical p와 BH FDR
- primary endpoint는 한 개로 제한

### 15.2 effect size 우선

p-value와 함께 반드시 다음을 보고한다.

- log fold change 또는 standardized beta
- delta PSI
- NES
- AUROC/AUPRC
- 95% CI 또는 bootstrap interval
- heterogeneity와 prediction interval

### 15.3 negative controls

- UVA와 무관한 perturbation signature
- random gene set
- gene expression/length matched set
- random receptor seeds
- rewired network
- HDF-CM/iPSC-CM comparator signature
- shuffled ortholog mapping
- 실패 또는 neutral longevity intervention이 제공되면 함께 평가

### 15.4 donor가 분석 단위

single-cell 또는 longitudinal sample이 많아도 donor 수를 넘어 biological n을 부풀리지 않는다.

- scRNA: donor-level pseudobulk
- longitudinal: donor random effect
- technical replicate: 평균 또는 nested technical term
- 여러 tissue가 같은 donor에서 온 경우 donor cluster-robust inference

## 16. M11 — 재현 가능한 구현

### 16.1 workflow

- Snakemake로 raw-to-result DAG
- conda-lock 또는 container digest
- R renv lockfile
- software version과 command-line parameter 자동 기록
- reference 및 input checksum
- 고정 random seed
- figure data를 TSV로 함께 출력

### 16.2 권장 산출물 구조

- metadata: sample inclusion, exclusion, provenance
- qc: raw/trim/alignment/junction reports
- expression: per-study counts, normalized matrices, effect tables
- splicing: event counts, PSI, validation status
- signatures: frozen gene sets와 version
- network: edge provenance, seed weights, parameter
- meta: study effects, heterogeneity, leave-one-out
- figures: plot source tables와 scripts
- reports: HTML/Quarto analysis report

### 16.3 공개

- 로컬 sequencing은 GEO/SRA/ENA
- code는 GitHub와 Zenodo release
- environment lockfile과 container
- public dataset accession 및 sample inclusion table
- full Data Availability statement
- 선행논문의 같은 data를 재사용했음을 명시

[Aging Cell 저자 지침](https://onlinelibrary.wiley.com/page/journal/14749726/homepage/forauthors.html)은 supporting data의 public repository 보관과 repository link를 요구하고, 가능한 경우 script도 공개하도록 권고한다.

## 17. M12 — 단계별 go/no-go

| gate | 통과 조건 | 실패 시 조치 |
|---|---|---|
| Gate 0: input | iPSC R1 교체, 전 FASTQ integrity pass, accession/metadata 확보 | 분석 시작하지 않음 |
| Gate 1: local rank | 두 time-matched contrast에서 rank가 QC/normalization에 안정 | Repro-specific claim 중단 |
| Gate 2: UVA | 최소 두 독립 UVA injury와 반대, rescue와 같은 방향 | photoaging generalization 축소 |
| Gate 3: aging | 최소 두 aging/senescence 코호트에서 반대, identity 보존 | rejuvenation 용어 제거 |
| Gate 4: splicing | main event가 최소 두 독립 raw dataset에서 방향 재현 | splicing을 제목/main mechanism에서 제거 |
| Gate 5: network | degree-matched null보다 우수하고 held-out에서 성능 유지 | network를 exploratory supplement로 이동 |
| Gate 6: longevity | tissue/study-stratified 수렴성과 null 대비 우수 | lifespan 관련 문구 제거 |
| Gate 7: final | held-out endpoint 통과, code/data 공개 가능 | Aging Cell 대신 범위가 맞는 저널 재선정 |

## 18. Figure별 최종 Methods 대응표

| Figure | primary method | 필요한 데이터 | 가장 큰 bias control |
|---|---|---|---|
| 1 | within-study effect, directional GSEA, random-effects meta-analysis | local mRNA rank, 4 UVA cohorts | external concordance와 direct replication 구분 |
| 2 | rank-based multi-axis scoring, mixed models | age, senescence, MPTR cohorts | proliferation/identity/dedifferentiation 분리 |
| 3 | rMATS/LeafCutter/MAJIQ, eCLIP/motif | local 및 public raw RNA-seq | N=1 p-value 금지, independent event validation |
| 4 | context network, RWR/PPR, permutation | quantitative secretome/miRNA, fibroblast network | degree bias, circular tuning, seed uncertainty |
| 5 | one-to-one ortholog signed enrichment | longevity/mortality coefficient tables | species/tissue 차이와 ortholog null |
| 6 | robust rank aggregation, stability selection, held-out benchmark | 모든 evidence family | discovery–validation leakage 차단 |

## 19. 현재 다운로드 상태와 사용자의 다음 승인 단위

현재 공개 데이터 다운로드: 0 byte.

사용자가 다운로드를 지시할 때 권장되는 최소 승인 단위는 Tier A processed package다.

- UVA: GSE125429, GSE89005, GSE240226, GSE302943
- aging/senescence: GSE109700, GSE191055, GSE226189, GSE113957, GSE179848, GSE93535
- partial reprogramming: GSE165177
- longevity: Cell Metabolism 2019 및 Nature 2026 supplementary tables
- network: 소용량 curated snapshot만

Tier A 분석이 Gate 2–3을 통과한 뒤에만 SRA RunInfo와 정확한 byte 수를 다시 제시하고, 사용자의 별도 승인을 받아 Tier B raw FASTQ를 순차적으로 받는다.

## 20. 분석 시작 전 체크리스트

- [ ] iPSC_CM_24h.R1.fastq.gz 정상본 교체
- [ ] 모든 로컬 FASTQ full gzip integrity pass
- [ ] 제공자 checksum 비교
- [ ] strandedness/sample sheet 확인
- [ ] local mRNA-seq accession 확보
- [ ] miRNA raw/count와 replicate 확인
- [ ] antibody array 원 intensity 확인
- [ ] 선행논문과 후속논문의 중복·재사용 문구 확정
- [ ] discovery/held-out dataset table 동결
- [ ] primary endpoint 및 gene-set version 동결
- [ ] Tier A download manifest의 URL, byte, checksum 재확인
- [ ] 사용자 다운로드 승인

이 체크리스트가 완료되기 전에는 본 분석의 p-value, 후보 순위 또는 저널 적합성 결론을 확정하지 않는다.
