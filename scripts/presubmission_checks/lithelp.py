import json, urllib.parse, urllib.request, time, re, textwrap
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def req(u):
    for a in range(3):
        try: return urllib.request.urlopen(u, timeout=40).read().decode()
        except Exception as e:
            if a == 2: return ""
            time.sleep(0.8)
def esearch(term, n=1):
    r = req(B + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "term": term, "retmax": n, "retmode": "json", "tool": "claude-code"}))
    time.sleep(0.4)
    try: return json.loads(r)["esearchresult"]["idlist"]
    except Exception: return []
def fetch(pmids):
    if not pmids: return {}
    x = req(B + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "xml", "tool": "claude-code"}))
    time.sleep(0.4)
    out = {}
    for art in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", x, re.S):
        pid = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
        ti = re.search(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", art, re.S)
        jr = re.search(r"<Title>(.*?)</Title>", art, re.S)
        yr = re.search(r"<Year>(\d{4})</Year>", art)
        ab = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.S))
        clean = lambda s: re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s or "")).strip()
        if pid: out[pid.group(1)] = (yr.group(1) if yr else "", clean(jr.group(1) if jr else ""),
                                     clean(ti.group(1) if ti else ""), clean(ab))
    return out
