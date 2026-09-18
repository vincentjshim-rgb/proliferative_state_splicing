#!/usr/bin/env bash
# One Korean reading page that carries, in this order:
#   1부  요청사항 재점검 (review_audit_ko.md)
#   2부  모의 심사 지적과 대응 (revision_response_ko.md)
#   3부  원고 v8 한글 전문, 그림 6개와 보충 그림 10개 포함 (manuscript_splicing_v2_ko.md)
#   4부  투고 규정 대조표 (submission_checklist_ko.md)
#   5부  커버레터 초안 (cover_letter_geroscience_ko.md)
#   6부  모의 심사 전문 (review_presubmission_npjaging_ko.md)
# Output: revision_bundle_ko.html  (published to the review artifact URL)
set -euo pipefail
WS="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
B=$WS/scripts/artifact_build
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/figweb"

# web-size figures, as in build_reading_artifacts.sh
python3 - "$WS" "$T" <<'PY'
import sys, pathlib
from PIL import Image
ws, t = map(pathlib.Path, sys.argv[1:])
for p in sorted((ws / "public_data_tierA/derived/figures_sciadv").glob("Fig*.png")):
    im = Image.open(p); w, h = im.size
    im.convert("RGB").resize((1800, round(h * 1800 / w)), Image.Resampling.LANCZOS).save(t / "figweb" / p.name, optimize=True)
PY

# concatenate the three documents, demoting their headings one level under a part heading
python3 - "$WS" "$T" <<'PY'
import sys, pathlib, re
ws, t = map(pathlib.Path, sys.argv[1:])
def load(name, drop_frontmatter=True):
    s = (ws / name).read_text()
    if drop_frontmatter and s.startswith("---"):
        s = s.split("---", 2)[2].lstrip("\n")
    return s
def demote(s):
    # ###### stays put; everything else moves one level down
    return re.sub(r'^(#{1,5}) ', lambda m: "#" + m.group(1) + " ", s, flags=re.M)

parts = [
    ("# 1부 · 요청사항 재점검\n\n편집자와 리뷰어의 요청을 원고 v8과 항목별로 대조한 결과, 남은 항목, 그리고 이번 보완이 만들어낸 새 약점이다.\n\n",
     demote(load("review_audit_ko.md"))),
    ("\n\n# 2부 · 모의 심사 지적과 대응\n\n원고 v5의 모의 심사에서 나온 지적과, 원고 v8에서 그것을 어떻게 고쳤는지의 기록이다.\n\n",
     demote(load("revision_response_ko.md"))),
    ("\n\n# 3부 · 원고 v8 한글 전문\n\n본문 그림 6개와 보충 그림 10개를 포함한 전문이다. 영문판은 `manuscript_splicing_v2.md`이며 수치와 그림은 동일하다.\n\n",
     demote(load("manuscript_splicing_v2_ko.md"))),
    ("\n\n# 4부 · 투고 규정 대조표\n\n2026년 9월 16일에 저널의 Content types와 Submission guidelines를 직접 받아 대조한 결과다.\n\n",
     demote(load("submission_checklist_ko.md"))),
    ("\n\n# 5부 · 커버레터 초안 (한글)\n\n투고 시 제출할 커버레터의 한글 확인본이다. 영문 원본은 `cover_letter_geroscience.md`.\n\n",
     demote(load("cover_letter_geroscience_ko.md"))),
    ("\n\n# 6부 · 모의 심사 전문\n\n2026년 9월 15일에 npj Aging 편집자와 리뷰어 3명의 관점으로 수행한 제출 전 심사의 전문이다. 여기의 지적이 1부의 대응으로 이어진다.\n\n",
     demote(load("review_presubmission_npjaging_ko.md"))),
]
md = "".join(h + b for h, b in parts)
md = md.replace("public_data_tierA/derived/figures_sciadv/", str(t / "figweb") + "/")
(t / "bundle_ko.md").write_text(md)
print(f"bundle: {len(md):,} characters")
PY

pandoc "$T/bundle_ko.md" -s --embed-resources --standalone -o "$T/body_ko.html"

# template: masthead and summary strip for the bundle, plus the review severity chips
python3 - "$WS" "$T" <<'PY'
import sys, pathlib, re
ws, t = map(pathlib.Path, sys.argv[1:])
tpl = (ws / "scripts/artifact_build/tpl_v5_ko.html").read_text()
tpl = tpl.replace("<title>분열 속도와 노화 signature</title>", "<title>모의 심사와 원고 v8</title>")
tpl = re.sub(r'<p class="eyebrow">.*?</p>', '<p class="eyebrow">모의 심사 · 대응 · 원고 v8 · 한글판</p>', tpl, count=1)
tpl = re.sub(r'<h1 class="doc">.*?</h1>',
             '<h1 class="doc">모의 심사, 대응, 그리고 원고 v8</h1>', tpl, count=1, flags=re.S)
tpl = re.sub(r'<p class="sub">.*?</p>',
             '<p class="sub">지적과 대응 · 그림이 들어간 원고 전문 · 투고 규정 대조 · 심사 전문 · 2026년 9월 16일</p>',
             tpl, count=1)
strip = re.search(r'<div class="strip">.*?</div>\n</div>', tpl, re.S)
tpl = tpl[:strip.start()] + '''<div class="strip">
  <div><b>6 / 6</b><span>필수 항목 처리</span></div>
  <div><b>107명</b><span>정상 성인 공여자</span></div>
  <div><b>70%</b><span>증식으로 설명되는 연령 효과</span></div>
  <div><b>0.74 &rarr; 0.23 &rarr; &minus;0.12</b><span>배양, 피부, 근육 (GTEx)</span></div>
  <div><b>8 + 3</b><span>본문 그림과 보충 그림</span></div>
  <div><b>7</b><span>투고 전 사용자 입력이 남은 항목</span></div>
</div>
</div>''' + tpl[strip.end():]
chips = '''
/* ---- review: severity chips ---- */
.sev{display:inline-block;font-family:var(--sans);font-size:10.5px;font-weight:600;letter-spacing:.06em;
  line-height:1;padding:3px 6px 2px;border:1px solid;border-radius:2px;vertical-align:2px;margin-right:6px;white-space:nowrap}
.sev.c{color:var(--flag);border-color:var(--flag)}
.sev.m{color:var(--accent);border-color:var(--accent)}
.sev.n{color:var(--ink3);border-color:var(--rule)}
/* part dividers */
h1.s1{border-top:2px solid var(--ink);padding-top:14px}
'''
tpl = tpl.replace("</style>", chips + "</style>", 1)
(t / "tpl_v5_ko.html").write_text(tpl)
PY

cp "$B/bld_v5.py" "$T/"
python3 "$T/bld_v5.py" "$T" "_ko"
cp "$T/art_v5_ko.html" "$WS/revision_bundle_ko.html"
echo "built revision_bundle_ko.html ($(stat -c%s "$WS/revision_bundle_ko.html") bytes)"
