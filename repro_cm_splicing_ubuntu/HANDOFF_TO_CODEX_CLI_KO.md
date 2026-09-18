# 분석컴 Codex CLI 인계문

아래 내용을 Ubuntu 분석컴의 Codex CLI에 그대로 전달할 수 있습니다.

---

나는 Repro-CM mRNA-seq의 alternative-splicing 탐색 분석을 하려고 한다. 프로젝트 폴더는 다음 위치에 있다.

`/home/shim/Downloads/develop_shim_project_260906/repro_cm_splicing_ubuntu`

먼저 이 폴더의 `README_KO.md`, `config/config.yaml`, `config/samples.tsv`, `config/contrasts.tsv`, `workflow/Snakefile`을 전부 읽고 현재 시스템을 점검해라. 원본 FASTQ는 다음 상위 디렉터리에 있으며 절대로 이동, 이름 변경, 수정 또는 삭제하지 마라.

`/home/shim/Downloads/develop_shim_project_260906/00.RNA-seq_samle`

주의: `samle` 철자는 오타처럼 보이지만 현재 제공된 실제 경로이므로, 존재 여부를 확인하기 전에는 `sample`로 고치지 마라.

연구 설계에는 `15J_0hr`, `HDF_CM_24h`, `Rep_CM_24h`, `iPSC_CM_24h`가 각각 한 biological sample만 있다. 따라서 다음 원칙을 지켜라.

1. rMATS 통계 검정, differential-splicing P value 및 FDR을 현재 자료에서 사용하지 마라.
2. `--statoff` 상태에서 PSI, ΔPSI, junction coverage 기반 exploratory candidate만 생성하라.
3. `15J_0hr`와 24 h 조건 비교는 treatment와 elapsed time이 confounded되었다고 표시하라.
4. matched healthy 24 h와 UVA+vehicle 24 h control이 없으므로 “UVA splicing rescue”라고 단정하지 마라.
5. primary contrast는 `Rep_CM_24h vs HDF_CM_24h`와 `Rep_CM_24h vs iPSC_CM_24h`이다.
6. 원본 데이터와 기존 결과를 덮어쓰지 마라. 워크플로가 지정한 `results/` 아래에만 출력하라.

다음 순서로 작업하고 각 단계의 실제 결과를 확인하라.

```bash
cd /home/shim/Downloads/develop_shim_project_260906/repro_cm_splicing_ubuntu
chmod +x setup_runner_env.sh setup_reference.sh system_check.sh run.sh run_with_env.sh
./system_check.sh
./setup_runner_env.sh
CORES=4 ./run_with_env.sh results/metadata/preflight.json
cat results/metadata/preflight.json
```

preflight가 실패하면 오류 원인만 진단하고, 원본 파일을 변경하지 않은 채 `config/samples.tsv`의 잘못된 경로만 수정해라. 성공하면 다음을 수행한다.

```bash
./setup_reference.sh
CORES=16 ./run_with_env.sh --dry-run
```

dry-run에서 생성 예정 job, reference 경로, sample 수 4개, contrast 수 5개를 확인해라. 분석컴의 RAM이 40 GB 미만이거나 여유 디스크가 150 GB 미만이면 전체 실행 전에 사용자에게 보고해라. 충분하면 다음 명령으로 실행하고 로그를 모니터링해라.

```bash
CORES=16 ./run_with_env.sh 2>&1 | tee workflow_run.log
```

실패하면 다음을 우선 확인해라.

- `results/logs/preflight.log`
- `results/logs/fastp/`
- `results/logs/star_index.log`
- `results/logs/star_align/`
- `results/logs/rmats/`
- `.snakemake/log/`

임의로 낮은 품질 기준을 적용하거나 통계 검정을 켜서 오류를 숨기지 마라. 수정이 필요하면 변경 전 원인, 변경 파일, 생물정보학적 영향을 설명하고 최소 변경만 수행해라.

완료 후 다음 내용을 한국어로 보고해라.

- 네 FASTQ pair의 크기, SHA-256, detected read length
- raw/trimmed read 수와 Q30
- sample별 uniquely mapped rate, multimapping rate, splice-junction 수
- BAM 및 BAI 존재 여부
- contrast별 SE/A5SS/A3SS/MXE/RI event 수
- 기본 필터를 통과한 candidate 수
- 두 primary contrast에서 같은 방향으로 지지된 event 수와 상위 후보 20개
- 각 후보의 gene, event type, ΔPSI, junction coverage
- `N=1 exploratory` 한계와 다음 RT-PCR 검증 우선순위

결과를 설명할 때 `significant`, `proved`, `restored`라는 표현은 독립 biological validation 전에는 사용하지 마라. AlphaFold 분석은 junction RT-PCR로 확인되고 coding sequence가 달라지는 후보가 선택될 때까지 시작하지 마라.

---

