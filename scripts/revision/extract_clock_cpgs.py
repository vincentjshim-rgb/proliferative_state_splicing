"""Pull the clock CpGs out of the 2 GB methylation matrix without loading it.

The matrix has one row per CpG and two columns per sample ("<sentrix> beta" and
"<sentrix> Detection Pval"). Only the ~1,300 CpGs used by the four clocks are kept,
with their detection P values, so that the analysis works on a small table.

Usage: python3 scripts/revision/extract_clock_cpgs.py
Output: public_data_tierA/methylation/clock_cpg_betas.tsv.gz
"""
import csv, gzip, pathlib, sys

WS = pathlib.Path(os.environ.get("PROJECT_ROOT", pathlib.Path(__file__).resolve().parents[2]))
MET = WS / "public_data_tierA/methylation"
COEF = MET / "clock_coefficients"
SRC = MET / "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz"
OUT = MET / "clock_cpg_betas.tsv.gz"

## the union of CpGs used by the clocks
want = set()
for f, col in (("Horvath2013_AdditionalFile3.csv", "CpGmarker"),
               ("Horvath2.csv", "CpGmarker"),
               ("PhenoAge.csv", "CpGmarker"),
               ("Hannum.csv", "CpGmarker")):
    with open(COEF / f, newline="") as fh:
        for row in csv.DictReader(fh):
            cg = (row.get(col) or "").strip()
            if cg.startswith("cg"):
                want.add(cg)
print(f"clock CpGs wanted: {len(want)}", flush=True)

kept = 0
with gzip.open(SRC, "rt", newline="") as fh, gzip.open(OUT, "wt", newline="") as out:
    reader = csv.reader(fh)
    header = next(reader)
    header[0] = "cpg"
    out.write("\t".join(header) + "\n")
    for row in reader:
        if row[0] in want:
            out.write("\t".join(row) + "\n")
            kept += 1
            if kept % 200 == 0:
                print(f"  {kept} of {len(want)}", flush=True)
print(f"kept {kept} CpGs of {len(want)}; columns: {len(header)}")
