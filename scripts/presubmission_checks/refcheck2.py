import re, json, time, sys, urllib.parse, html as H
sys.path.insert(0, sys.argv[1])
from lithelp import req, B
md = open("manuscript_splicing_v2.md").read()
block = md.split("# References\n", 1)[1].split("\n# ", 1)[0]
refs = re.findall(r"^(\d+)\. (.+)$", block, re.M)
def clean(s): return H.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s or ""))).strip()
def norm(s): return re.sub(r"[^a-z0-9]", "", H.unescape(s).lower())
def expand_pages(pg):
    a = pg.split("-")
    if len(a) == 2 and a[1].isdigit() and a[0].isdigit() and len(a[1]) < len(a[0]): return a[0] + "-" + a[0][:len(a[0])-len(a[1])] + a[1]
    return pg
rows = []
for num, ref in refs:
    jm = re.search(r"\*([^*]+)\* \*\*([^*]+)\*\*, ([^ (]+) \((\d{4})\)", ref)
    if "et al.*" in ref.split("*" + (jm.group(1) if jm else "@@") + "*")[0]:
        auth, rest = ref.split("*et al.*", 1); auth = auth.strip() + " et al."
    elif ref.startswith("The GTEx Consortium."):
        auth, rest = "The GTEx Consortium", ref[len("The GTEx Consortium."):]
    else:
        am = re.match(r"^(.*?& [^,]+, (?:[A-Z]\. ?)+)", ref)
        auth, rest = (am.group(1), ref[am.end():]) if am else ("?", ref)
    if jm:
        title = rest.split("*" + jm.group(1) + "*")[0].strip().rstrip(".")
        journal, vol, pages, year = jm.groups()
    else:
        title = rest.split(". *")[0].strip(); journal = vol = pages = ""; year = (re.findall(r"\((\d{4})\)", ref) or [""])[-1]
    title_q = re.sub(r"\*", "", title)
    ids = []
    for term in [f'"{title_q}"[TI]', " AND ".join(w + "[TI]" for w in re.findall(r"[A-Za-z0-9]{5,}", title_q)[:8])]:
        d = req(B + "esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": 2, "retmode": "json", "tool": "claude-code"})); time.sleep(0.4)
        try: ids = json.loads(d)["esearchresult"]["idlist"]
        except Exception: ids = []
        if ids: break
    if not ids:
        rows.append((num, "NOT-IN-PUBMED", f"{auth[:28]} | {title_q[:80]} | {journal} {vol}, {pages} ({year})")); continue
    x = req(B + "efetch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ids[0], "retmode": "xml", "tool": "claude-code"})); time.sleep(0.4)
    g = lambda p: clean(re.search(p, x, re.S).group(1)) if re.search(p, x, re.S) else ""
    pt = g(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>").rstrip(".")
    jab = g(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>"); pvol = g(r"<Volume>(.*?)</Volume>")
    pg = g(r"<MedlinePgn>(.*?)</MedlinePgn>"); yr = g(r"<PubDate>.*?<Year>(\d{4})</Year>")
    eyr = g(r'<ArticleDate DateType="Electronic">.*?<Year>(\d{4})</Year>')
    doi = g(r'<ArticleId IdType="doi">(.*?)</ArticleId>')
    au = [(clean(l), clean(i)) for l, i in re.findall(r"<Author[^>]*>\s*<LastName>(.*?)</LastName>.*?<Initials>(.*?)</Initials>", x, re.S)]
    coll = g(r"<CollectiveName>(.*?)</CollectiveName>")
    issues = []
    if norm(pt) != norm(title_q): issues.append(f"TITLE: PubMed='{pt}'")
    if vol and pvol and vol != pvol: issues.append(f"VOL {vol}≠{pvol}")
    mp = pages.replace("–", "-")
    if pg and norm(mp) not in (norm(pg), norm(expand_pages(pg))): issues.append(f"PAGES {pages}≠{pg}")
    if year not in (yr, eyr): issues.append(f"YEAR {year}≠{yr} (epub {eyr})")
    if au and auth != "The GTEx Consortium":
        if norm(au[0][0]) != norm(auth.split(",")[0]): issues.append(f"1st AUTHOR {auth.split(',')[0]}≠{au[0][0]}")
        ini = re.match(r"^[^,]+, ((?:[A-Z]\.\s?-?)+)", auth)
        if ini and norm(ini.group(1)) != norm(au[0][1]): issues.append(f"INITIALS {ini.group(1).strip()}≠{au[0][1]}")
    n = len(au)
    if n > 5 and "et al." not in auth: issues.append(f"STYLE {n} authors → et al.")
    if 0 < n <= 5 and "et al." in auth: issues.append(f"STYLE only {n} authors → list all")
    if jab and norm(jab) != norm(journal): issues.append(f"JOURNAL abbrev '{journal}' vs NLM '{jab}'")
    rows.append((num, "OK" if not issues else "CHECK", f"PMID {ids[0]} | {jab} {pvol}:{pg} ({yr}/e{eyr}) | {doi} | {n} au" + (f" | coll={coll}" if coll else "") + ("" if not issues else "\n        ⚑ " + " ; ".join(issues))))
for r in rows: print(f"[{r[0]:>2}] {r[1]:<13} {r[2]}")
