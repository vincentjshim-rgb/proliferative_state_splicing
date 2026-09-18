#!/usr/bin/env python3
"""Conservatively map GPL16956 60-mer probes to GENCODE v50 genes."""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
PLATFORM = ROOT / "public_data_tierA/metadata/GPL16956_platform_full.soft"
GENOME = ROOT / "repro_cm_splicing_ubuntu/resources/gencode_v50/GRCh38.primary_assembly.genome.fa"
GTF = ROOT / "repro_cm_splicing_ubuntu/resources/gencode_v50/gencode.v50.primary_assembly.annotation.gtf.gz"
OUT = ROOT / "public_data_tierA/derived/GPL16956_probe_mapping"
OUT.mkdir(parents=True, exist_ok=True)

FASTA = OUT / "GPL16956_probes.fa"
SAM = OUT / "GPL16956_to_GRCh38.sam"
PROBE_BED = OUT / "GPL16956_high_confidence_alignments.bed"
GENE_BED = OUT / "GENCODE_v50_gene_spans.bed"
INTERSECT = OUT / "probe_gene_intersections.tsv"
MAPPING = OUT / "GPL16956_probe_to_GENCODE_v50.tsv"
METRICS = OUT / "mapping_metrics.tsv"


def parse_platform():
    probes = {}
    in_table = False
    header = None
    with PLATFORM.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
                continue
            row = line.split("\t")
            if header is None:
                header = row
                continue
            record = dict(zip(header, row))
            probe = record.get("ID", "")
            sequence = record.get("SEQUENCE", "").upper()
            if probe and re.fullmatch(r"[ACGT]+", sequence or ""):
                probes[probe] = {
                    "sequence": sequence,
                    "transcript_type": record.get("TRANSCRIPT_TYPE", ""),
                    "build": record.get("BUILD", ""),
                }
    return probes


def cigar_lengths(cigar):
    query = 0
    reference = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in "MI=X":
            query += length
        if op in "MDN=X":
            reference += length
    return query, reference


probes = parse_platform()
with FASTA.open("w", encoding="ascii") as handle:
    for probe, record in probes.items():
        handle.write(f">{probe}\n{record['sequence']}\n")

with SAM.open("w", encoding="utf-8") as stdout, (OUT / "minimap2.stderr.log").open("w", encoding="utf-8") as stderr:
    subprocess.run(
        ["minimap2", "-a", "-x", "sr", "--secondary=no", "-t", "4", str(GENOME), str(FASTA)],
        check=True,
        stdout=stdout,
        stderr=stderr,
    )

alignments = {}
with SAM.open("r", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith("@"):
            continue
        row = line.rstrip("\n").split("\t")
        probe, flag_text, chrom, pos_text, mapq_text, cigar = row[:6]
        flag = int(flag_text)
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        query_len, ref_len = cigar_lengths(cigar)
        nm = None
        for tag in row[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":")[-1])
                break
        start = int(pos_text) - 1
        alignments[probe] = {
            "chrom": chrom,
            "start": start,
            "end": start + ref_len,
            "mapq": int(mapq_text),
            "strand": "-" if flag & 0x10 else "+",
            "query_aligned": query_len,
            "nm": nm,
        }

passed = {}
for probe, aln in alignments.items():
    if aln["mapq"] >= 30 and aln["query_aligned"] >= 55 and aln["nm"] is not None and aln["nm"] <= 2:
        passed[probe] = aln

with PROBE_BED.open("w", encoding="utf-8") as handle:
    for probe, aln in sorted(passed.items(), key=lambda item: (item[1]["chrom"], item[1]["start"], item[0])):
        handle.write(
            f"{aln['chrom']}\t{aln['start']}\t{aln['end']}\t{probe}\t{aln['mapq']}\t{aln['strand']}\n"
        )

attr_re = re.compile(r'(\S+) "([^"]*)";')
with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as source, GENE_BED.open("w", encoding="utf-8") as dest:
    for line in source:
        if line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        if len(row) != 9 or row[2] != "gene":
            continue
        attrs = dict(attr_re.findall(row[8]))
        gene_id = attrs.get("gene_id", "")
        gene_name = attrs.get("gene_name", gene_id)
        dest.write(f"{row[0]}\t{int(row[3]) - 1}\t{row[4]}\t{gene_id}\t{gene_name}\t{row[6]}\n")

with INTERSECT.open("w", encoding="utf-8") as stdout:
    subprocess.run(
        ["bedtools", "intersect", "-a", str(PROBE_BED), "-b", str(GENE_BED), "-wa", "-wb"],
        check=True,
        stdout=stdout,
    )

probe_genes = defaultdict(set)
with INTERSECT.open("r", encoding="utf-8") as handle:
    for line in handle:
        row = line.rstrip("\n").split("\t")
        probe_genes[row[3]].add((row[9], row[10]))

status_counts = defaultdict(int)
with MAPPING.open("w", newline="", encoding="utf-8") as handle:
    fields = [
        "probe", "transcript_type", "submitted_build", "chrom", "start_0based", "end",
        "mapq", "query_aligned_nt", "NM", "gene_id", "gene_symbol", "status",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for probe, record in probes.items():
        aln = alignments.get(probe)
        genes = probe_genes.get(probe, set())
        status = "PASS"
        gene_id = ""
        gene_symbol = ""
        if aln is None:
            status = "no_primary_alignment"
        elif aln["mapq"] < 30:
            status = "low_mapq"
        elif aln["query_aligned"] < 55:
            status = "short_alignment"
        elif aln["nm"] is None or aln["nm"] > 2:
            status = "high_or_missing_NM"
        elif not genes:
            status = "no_overlapping_gene"
        else:
            symbols = {symbol for _, symbol in genes}
            if len(symbols) != 1:
                status = "ambiguous_gene"
            else:
                gene_symbol = next(iter(symbols))
                gene_ids = sorted(gene for gene, symbol in genes if symbol == gene_symbol)
                gene_id = ";".join(gene_ids)
        status_counts[status] += 1
        writer.writerow({
            "probe": probe,
            "transcript_type": record["transcript_type"],
            "submitted_build": record["build"],
            "chrom": "" if aln is None else aln["chrom"],
            "start_0based": "" if aln is None else aln["start"],
            "end": "" if aln is None else aln["end"],
            "mapq": "" if aln is None else aln["mapq"],
            "query_aligned_nt": "" if aln is None else aln["query_aligned"],
            "NM": "" if aln is None else aln["nm"],
            "gene_id": gene_id,
            "gene_symbol": gene_symbol,
            "status": status,
        })

with METRICS.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["valid_ACGT_probes", len(probes)])
    writer.writerow(["primary_alignments", len(alignments)])
    writer.writerow(["alignment_threshold_pass", len(passed)])
    for status, count in sorted(status_counts.items()):
        writer.writerow([f"status_{status}", count])
    writer.writerow(["pass_fraction", status_counts["PASS"] / len(probes) if probes else 0])

print(f"Mapped {status_counts['PASS']} of {len(probes)} valid probes to one GENCODE v50 gene.")
