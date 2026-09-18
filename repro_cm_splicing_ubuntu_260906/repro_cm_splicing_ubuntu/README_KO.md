# Repro-CM alternative-splicing discovery workflow

이 폴더는 Ubuntu 분석컴에서 네 개의 paired-end mRNA FASTQ를 안전하고 재현 가능하게 처리하기 위한 Snakemake 워크플로입니다.

현재 조건마다 biological sample이 하나뿐이므로 rMATS 통계 검정은 의도적으로 비활성화되어 있습니다. 결과는 splice-event 존재, junction coverage, PSI 및 ΔPSI를 이용한 **후보 발굴용**이며 유의한 differential splicing으로 보고하면 안 됩니다.

## 입력 시료

입력 경로는 `config/samples.tsv`에 등록되어 있습니다.

- `15J_0hr`
- `HDF_CM_24h`
- `Rep_CM_24h`
- `iPSC_CM_24h`
- 각 시료는 `.R1.fastq.gz`와 `.R2.fastq.gz`로 구성

디렉터리 이름 `00.RNA-seq_samle`은 제공받은 경로를 그대로 사용했습니다. 분석컴의 실제 이름이 `00.RNA-seq_sample`이라면 원본 파일을 옮기지 말고 `config/samples.tsv`의 경로만 수정합니다.

## 워크플로 구성

```text
FASTQ preflight 및 SHA-256 inventory
  ├─ raw FastQC
  └─ fastp paired-end trimming
       └─ trimmed FastQC
            └─ STAR two-pass alignment (GRCh38/GENCODE v50)
                 ├─ BAM index 및 flagstat
                 ├─ featureCounts gene-count QC
                 └─ rMATS event counting
                      ├─ SE, A5SS, A3SS, MXE, RI
                      ├─ JC/JCEC PSI 및 ΔPSI
                      └─ Repro-CM 후보 catalog
```

사용 버전은 Conda 환경 파일에 고정되어 있고, 실제 실행 버전은 `results/metadata/software_versions.tsv`에 기록됩니다.

## 권장 분석컴 사양

- Ubuntu x86-64
- CPU 12–24 threads
- RAM 최소 40 GB, 권장 64 GB 이상
- 여유 디스크 최소 150 GB, FASTQ 크기가 크면 250 GB 이상
- reference와 Conda package를 받을 인터넷 연결
- `micromamba`, `mamba` 또는 `conda` 중 하나

GPU는 필요하지 않습니다. AlphaFold 단계는 이 워크플로에 포함하지 않았으며, RNA 수준에서 검증된 coding isoform이 확보된 뒤 별도 단계로 진행합니다.

## 최초 실행

워크플로 폴더로 이동한 후 실행합니다.

```bash
cd /home/shim/Downloads/develop_shim_project_260906/repro_cm_splicing_ubuntu
chmod +x setup_runner_env.sh setup_reference.sh system_check.sh run.sh run_with_env.sh
./system_check.sh
./setup_runner_env.sh
```

먼저 FASTQ만 점검합니다. 이 단계는 파일 크기에 따라 SHA-256 계산 시간이 걸립니다.

```bash
CORES=4 ./run_with_env.sh results/metadata/preflight.json
cat results/metadata/preflight.json
column -ts $'\t' results/metadata/fastq_profiles.tsv | less -S
```

`analysis_status`가 `ready_for_exploratory_analysis`인지 확인한 다음 GENCODE v50 reference를 준비합니다.

```bash
./setup_reference.sh
```

전체 DAG를 실행하지 않고 먼저 확인합니다.

```bash
CORES=16 ./run_with_env.sh --dry-run
```

오류가 없으면 전체 분석을 시작합니다.

```bash
CORES=16 ./run_with_env.sh
```

터미널 연결이 끊길 수 있는 환경에서는 `tmux` 안에서 실행합니다.

```bash
tmux new -s repro_splicing
CORES=16 ./run_with_env.sh 2>&1 | tee workflow_run.log
```

`tmux`에서 빠져나오기는 `Ctrl-b` 다음 `d`, 다시 접속하기는 다음 명령입니다.

```bash
tmux attach -t repro_splicing
```

## 주요 결과

| 결과 | 위치 | 의미 |
|---|---|---|
| 입력 감사 | `results/metadata/preflight.json` | 시료 수, read length, 오류와 경고 |
| 입력 checksum | `results/metadata/input_sha256.tsv` | 원본 FASTQ 식별 및 변경 감시 |
| 통합 QC | `results/qc/multiqc/multiqc_report.html` | FastQC, fastp, STAR, flagstat, counts QC |
| 정렬 BAM | `results/alignment/<sample>/Aligned.sortedByCoord.out.bam` | 좌표 정렬된 RNA-seq alignment |
| splice junction | `results/alignment/<sample>/SJ.out.tab` | STAR junction counts |
| gene counts | `results/counts/gene_counts.tsv` | 기존 gene-level 결과와의 방향성 확인용 |
| contrast별 event | `results/rmats/<contrast>/summary/all_events.JC.tsv` | 모든 junction-count event와 ΔPSI |
| contrast별 후보 | `results/rmats/<contrast>/summary/candidates.JC.tsv` | coverage/effect 기준 통과 후보 |
| 통합 후보표 | `results/rmats/candidate_catalog.ranked.tsv` | primary contrast 재현성 기반 순위 |
| 실행 버전 | `results/metadata/software_versions.tsv` | software provenance |

JC 결과를 주 분석으로 사용하고 JCEC는 보조 민감도 분석으로 사용합니다. `delta_psi_group1_minus_group2`는 반드시 `config/contrasts.tsv`의 group 순서와 함께 해석합니다. 예를 들어 `Rep_CM_24h_vs_HDF_CM_24h`에서 양수는 Repro-CM의 exon/inclusion-form 사용이 더 높다는 뜻입니다.

기본 후보 필터는 다음과 같습니다.

- `abs(ΔPSI) >= 0.20`
- 양쪽 group의 모든 시료에서 inclusion+skipping junction read 합계 `>=20`
- P value와 FDR는 계산하지 않음

이는 통계적 유의성 기준이 아니라 N-of-1 자료에서 검증 후보 수를 줄이기 위한 보수적 탐색 기준입니다.

## 결과 해석 제한

현재 자료에서 허용되는 표현:

> Repro-CM-treated HDFs showed candidate differences in splice-junction usage compared with HDF-CM- or iPSC-CM-treated HDFs.

현재 자료에서 사용하면 안 되는 표현:

> Repro-CM significantly restored UVA-induced alternative splicing.

그 이유는 다음과 같습니다.

- 조건당 biological replicate가 하나
- `15J_0hr`와 CM-treated 24 h 사이에 시간 차이가 있음
- matched non-UVA 24 h control이 없음
- matched UVA+vehicle 24 h control이 없음
- treatment, donor 및 CM batch 효과가 분리되지 않음

상위 후보는 IGV/Sashimi plot 육안 확인 후, independent biological experiments `n>=4`에서 junction RT-PCR 또는 ddPCR로 검증해야 합니다.

## Library strandedness

기본값은 `fr-unstranded`입니다. 원래 library preparation 기록에서 stranded protocol이 확인된 경우에만 `config/config.yaml`을 수정합니다.

```yaml
analysis:
  library_type: fr-firststrand
```

rMATS 기준으로 Salmon의 `ISR`은 `fr-firststrand`, `ISF`는 `fr-secondstrand`에 해당합니다. 추정만으로 바꾸지 않습니다.

## 중단 후 재실행

같은 명령을 다시 실행하면 완료된 output은 재사용하고 미완료 단계만 이어서 수행합니다.

```bash
CORES=16 ./run_with_env.sh
```

원본 FASTQ를 삭제·이동·수정하지 않습니다. `results/`만 워크플로가 생성합니다.

## replicate가 추가될 때

1. 각 biological replicate에 고유 `sample_id`를 부여합니다.
2. `config/samples.tsv`에 행을 추가하고 같은 condition 이름을 사용합니다.
3. 조건별 최소 3개 이상, 권장 4–6개의 독립 biological replicate를 확보합니다.
4. `config/config.yaml`의 `statistical_test`를 `true`로 바꿉니다.
5. 기존 결과와 섞이지 않도록 새 results directory를 지정합니다.

```yaml
results_dir: results_replicated
analysis:
  rmats:
    statistical_test: true
```

technical lane 또는 FASTQ 분할본은 biological replicate로 등록하면 안 됩니다.

## 자주 발생하는 문제

- `MissingInputException`: `config/samples.tsv`의 절대경로와 `samle` 철자를 확인합니다.
- STAR index 중 메모리 부족: 동시에 실행 중인 작업을 줄이고 RAM 40 GB 이상을 확보합니다.
- Conda solve 실패: `.snakemake/conda`를 임의 삭제하기 전에 로그를 보존하고 Codex CLI에 진단을 요청합니다.
- MultiQC가 끝나지 않음: 먼저 `results/logs/multiqc.log`를 확인합니다.
- rMATS output이 비어 있음: read length, library type, STAR mapping rate, junction coverage를 순서대로 확인합니다.

## 근거와 버전 선택

- rMATS v4.4.0은 `--statoff`로 통계 분석을 건너뛰고 event count를 생성할 수 있습니다: <https://github.com/Xinglab/rmats-turbo/blob/v4.4.0/README.md>
- GENCODE v50은 이 워크플로 작성 시점의 current GRCh38 release입니다: <https://www.gencodegenes.org/human/>
- Bioconda package 환경으로 실행 도구를 고정합니다: <https://bioconda.github.io/>

