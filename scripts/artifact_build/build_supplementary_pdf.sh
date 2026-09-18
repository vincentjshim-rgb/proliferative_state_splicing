#!/usr/bin/env bash
# Build the single Supplementary Information PDF the journal asks for.
# markdown (assembled by build_supplementary_pdf.py) -> docx (pandoc) -> pdf (LibreOffice)
set -euo pipefail
WS="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
B=$WS/scripts/artifact_build
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

python3 "$B/build_supplementary_pdf.py" "$T"

# web-size copies of the three supplementary figures keep the PDF manageable
python3 - "$WS" "$T" <<'PY'
import sys, pathlib, re
from PIL import Image
ws, t = map(pathlib.Path, sys.argv[1:])
md = (t / "supplementary.md").read_text()
(t / "fig").mkdir(exist_ok=True)
for p in sorted((ws / "public_data_tierA/derived/figures_sciadv").glob("FigS*.png")):
    im = Image.open(p); w, h = im.size
    im.convert("RGB").resize((1800, round(h * 1800 / w)), Image.Resampling.LANCZOS).save(t / "fig" / p.name, optimize=True)
md = md.replace(str(ws / "public_data_tierA/derived/figures_sciadv"), str(t / "fig"))
(t / "supplementary.md").write_text(md)
PY

pandoc "$T/supplementary.md" -o "$T/supplementary.docx"
soffice --headless --convert-to pdf --outdir "$T" "$T/supplementary.docx" >/dev/null 2>&1
cp "$T/supplementary.docx" "$WS/supplementary_information.docx"
cp "$T/supplementary.pdf"  "$WS/supplementary_information.pdf"
echo "built supplementary_information.pdf ($(stat -c%s "$WS/supplementary_information.pdf") bytes) and .docx"
