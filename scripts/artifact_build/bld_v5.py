import re, html, pathlib, sys
S = pathlib.Path(sys.argv[1]); lang = sys.argv[2]           # "" or "_ko"
src = (S / f"body{lang}.html").read_text()
body = src[src.find("<body>")+6 : src.rfind("</body>")]
body = body[body.find("</header>")+9:]
body = re.sub(r'(<img[^>]*?)\s+style="width:[^"]*"', r'\1', body)

## a caption is "Fig. 3." in the main text and "Supplementary Fig. S1." / "보충 Fig. S1."
## in the supplementary section; both are laid out as plates, the supplementary ones
## keeping their S number so that the text references match
CAPHEAD = r'(?:Supplementary\s+Fig\.|보충\s+Fig\.|Fig\.)'

def figrepl(m):
    img, cap = m.group(1), m.group(2)
    num = re.search(CAPHEAD + r'\s*(S?\d+)\.', cap)
    n = num.group(1) if num else ''
    supp = n.startswith('S')
    cap = re.sub(r'^<strong>' + CAPHEAD + r'\s*S?\d+\.\s*(.*?)</strong>',
                 r'<span class="fig-t">\1</span>', cap, flags=re.S)
    label = ('Supplementary Fig. ' if supp and lang == '' else '보충 Fig. ' if supp else 'Fig. ') + n
    return (f'<figure class="plate" id="fig{n}">'
            f'<div class="plate-hd"><span class="fig-n">{label}</span></div>'
            f'<div class="plate-art">{img}</div>'
            f'<figcaption>{cap}</figcaption></figure>')
body, nfig = re.subn(r'<p>(<img[^>]*/?>)</p>\s*<p>(<strong>' + CAPHEAD + r'.*?)</p>', figrepl, body, flags=re.S)

# section -> figure number, read straight from the rendered order
FIGOF, order = {}, []
for m in re.finditer(r'<h2 id="([^"]+)"|<figure class="plate" id="fig(\d+)"', body):
    if m.group(1): order.append([m.group(1), None])
    elif order and order[-1][1] is None: order[-1][1] = m.group(2)
FIGOF = {h: f for h, f in order if f}

nav = []
for m in re.finditer(r'<h([12]) id="([^"]+)">(.*?)</h[12]>', body, flags=re.S):
    lvl, hid, txt = m.group(1), m.group(2), re.sub(r'<[^>]+>', '', m.group(3))
    nav.append((lvl, hid, txt, FIGOF.get(hid)))

def h_repl(m):
    lvl, hid, txt = m.group(1), m.group(2), m.group(3)
    f = FIGOF.get(hid)
    tag = f'<span class="h-fig">Fig.&nbsp;{f}</span>' if f else ''
    return f'<h{lvl} id="{hid}" class="s{lvl}">{tag}<span>{txt}</span></h{lvl}>'
body = re.sub(r'<h([12]) id="([^"]+)">(.*?)</h[12]>', h_repl, body, flags=re.S)
body = re.sub(r'<table', '<div class="tw"><table', body)
body = re.sub(r'</table>', '</table></div>', body)

navhtml = "".join(f'<a class="nv l{lvl}" href="#{hid}">{html.escape(txt)}'
                  f'{f"<em>{f}</em>" if f else ""}</a>' for lvl, hid, txt, f in nav)
TPL = (S / f"tpl_v5{lang}.html").read_text()
out = TPL.replace("<!--NAV-->", navhtml).replace("<!--BODY-->", body)
(S / f"art_v5{lang}.html").write_text(out)
print(f"lang{lang or '_en'}: figures {nfig}, nav {len(nav)}, bytes {len(out):,}")
