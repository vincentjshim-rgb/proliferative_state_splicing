#!/usr/bin/env bash
# Rebuild manuscript_splicing_v2{,_ko}.{docx,html}: web-size figures embedded in the reading pages.
# Masthead title and summary strip live in tpl_v5*.html; edit them when the title or headline numbers change.
set -euo pipefail
WS="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
B=$WS/scripts/artifact_build
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/figweb"
python3 - "$WS" "$T" <<'PY'
import sys, pathlib
from PIL import Image
ws, t = map(pathlib.Path, sys.argv[1:])
for p in sorted((ws / "public_data_tierA/derived/figures_sciadv").glob("Fig*.png")):
    im = Image.open(p); w, h = im.size
    im.convert("RGB").resize((1800, round(h * 1800 / w)), Image.Resampling.LANCZOS).save(t / "figweb" / p.name, optimize=True)
PY
cp "$B/tpl_v5.html" "$B/tpl_v5_ko.html" "$T/"
for L in "" "_ko"; do
  sed "s#public_data_tierA/derived/figures_sciadv/#$T/figweb/#g" "$WS/manuscript_splicing_v2${L}.md" > "$T/web${L}.md"
  pandoc "$T/web${L}.md" -s --embed-resources --standalone -o "$T/body${L}.html"
  python3 "$B/bld_v5.py" "$T" "$L"
  cp "$T/art_v5${L}.html" "$WS/manuscript_splicing_v2${L}.html"
done
cd "$WS"
pandoc manuscript_splicing_v2.md -o manuscript_splicing_v2.docx
pandoc manuscript_splicing_v2_ko.md -o manuscript_splicing_v2_ko.docx
echo "built manuscript_splicing_v2{,_ko}.{html,docx}"
