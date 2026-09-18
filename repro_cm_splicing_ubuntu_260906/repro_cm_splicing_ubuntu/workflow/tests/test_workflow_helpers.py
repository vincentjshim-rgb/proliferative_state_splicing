#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import unittest
from contextlib import contextmanager
from pathlib import Path

import yaml


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
TEST_TMP_ROOT = Path(__file__).resolve().parents[2] / ".test_tmp"
TEST_TMP_ROOT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(SCRIPT_DIR))

from rmats_summary import COORDINATE_COLUMNS, EVENT_TYPES, summarize_rmats


@contextmanager
def workspace_test_dir(name: str):
    """Use a predictable directory so Windows sandbox ACLs remain inherited."""
    path = TEST_TMP_ROOT / name
    path.mkdir(parents=True, exist_ok=True)
    yield str(path)


class WorkflowHelperTests(unittest.TestCase):
    def test_rmats_no_stats_summary(self):
        with workspace_test_dir("rmats_summary") as temp_name:
            root = Path(temp_name)
            rmats_dir = root / "rmats"
            rmats_dir.mkdir(exist_ok=True)
            for idx, event_type in enumerate(EVENT_TYPES, start=1):
                coordinates = COORDINATE_COLUMNS[event_type]
                definition_fields = ["ID", "GeneID", "geneSymbol", "chr", "strand", *coordinates]
                definition_values = [
                    str(idx),
                    f"ENSG_TEST_{idx}",
                    f"GENE{idx}",
                    "chr1",
                    "+",
                    *[str(100 * idx + offset) for offset in range(len(coordinates))],
                ]
                with (rmats_dir / f"fromGTF.{event_type}.txt").open(
                    "w", newline="", encoding="utf-8"
                ) as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(definition_fields)
                    writer.writerow(definition_values)
                raw_fields = [
                    "ID",
                    "IJC_SAMPLE_1",
                    "SJC_SAMPLE_1",
                    "IJC_SAMPLE_2",
                    "SJC_SAMPLE_2",
                    "IncFormLen",
                    "SkipFormLen",
                ]
                for mode in ("JC", "JCEC"):
                    with (rmats_dir / f"{mode}.raw.input.{event_type}.txt").open(
                        "w", newline="", encoding="utf-8"
                    ) as handle:
                        writer = csv.writer(handle, delimiter="\t")
                        writer.writerow(raw_fields)
                        writer.writerow([str(idx), "80", "20", "20", "80", "100", "100"])

            summary_dir = root / "summary"
            counts = summarize_rmats(
                rmats_dir=rmats_dir,
                contrast="Rep_vs_HDF",
                group1="Rep",
                group2="HDF",
                summary_dir=summary_dir,
                min_reads=20,
                min_abs_delta=0.2,
            )
            self.assertEqual(counts["all_JC"], 5)
            self.assertEqual(counts["candidates_JC"], 5)
            with (summary_dir / "candidates.JC.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["delta_psi_group1_minus_group2"], "0.600000")
            self.assertEqual(rows[0]["statistical_test"], "not_run_no_biological_replicates")

    def test_preflight_paired_fastq(self):
        with workspace_test_dir("preflight") as temp_name:
            root = Path(temp_name)
            r1 = root / "test.R1.fastq.gz"
            r2 = root / "test.R2.fastq.gz"
            with gzip.open(r1, "wt", encoding="ascii") as h1, gzip.open(
                r2, "wt", encoding="ascii"
            ) as h2:
                for idx in range(3):
                    h1.write(f"@read{idx} 1:N:0:1\nACGTACGT\n+\nFFFFFFFF\n")
                    h2.write(f"@read{idx} 2:N:0:1\nTGCATGCA\n+\nFFFFFFFF\n")
            manifest = root / "samples.tsv"
            manifest.write_text(
                "sample_id\tcondition\tbiological_replicate\tr1\tr2\n"
                f"test\ttest_condition\t1\t{r1}\t{r2}\n",
                encoding="utf-8",
            )
            contrasts = root / "contrasts.tsv"
            contrasts.write_text(
                "contrast_id\tgroup1\tgroup2\tpriority\tinterpretation\n"
                "self\ttest_condition\ttest_condition\tprimary\ttest\n",
                encoding="utf-8",
            )
            config = root / "config.yaml"
            config.write_text(
                yaml.safe_dump(
                    {
                        "contrasts": str(contrasts),
                        "analysis": {
                            "fastq_scan_records": 100,
                            "rmats": {"statistical_test": False},
                        },
                    }
                ),
                encoding="utf-8",
            )
            json_out = root / "preflight.json"
            profiles = root / "profiles.tsv"
            checksums = root / "checksums.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_DIR / "preflight.py"),
                    "--manifest",
                    str(manifest),
                    "--config",
                    str(config),
                    "--json",
                    str(json_out),
                    "--profiles",
                    str(profiles),
                    "--checksums",
                    str(checksums),
                    "--skip-sha256",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads(json_out.read_text(encoding="utf-8"))
            self.assertEqual(report["nominal_read_length"], 8)
            self.assertEqual(report["analysis_status"], "ready_for_exploratory_analysis")


if __name__ == "__main__":
    unittest.main()
