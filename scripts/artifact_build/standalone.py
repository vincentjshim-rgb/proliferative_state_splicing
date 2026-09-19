"""Wrap a page fragment as a complete HTML document.

The reading pages were written to be published as artifacts, where the host supplies
<!doctype>, <head> and <body>. Opened from disk that leaves them without a character-set
declaration, so a browser guesses the encoding and Korean text comes out as mojibake, and
without a doctype it lays the page out in quirks mode. This adds the missing skeleton so the
same file can be double-clicked and also published unchanged.
"""
import re

HEAD = """<!doctype html>
<html lang="{lang}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="color-scheme" content="light dark">
<style>html{{-webkit-text-size-adjust:100%}}body{{margin:0}}img{{max-width:100%}}
[hidden]{{display:none!important}}</style>
{head}</head>
<body>
"""
TAIL = "\n</body>\n</html>\n"
## <title>, <link> and <style> written at the top of a fragment belong in the head; a browser
## hoists them anyway, but styles left in the body can paint unstyled content first
LEADING = re.compile(r"\s*(<title>.*?</title>|<link\b[^>]*>|<meta\b[^>]*>|<style>.*?</style>)",
                     re.I | re.S)


def is_document(text):
    return bool(re.match(r"\s*<!doctype", text, re.I))


def wrap(text, lang="en"):
    """Return text as a complete document; already-complete input is returned unchanged."""
    if is_document(text):
        return text
    head, rest = [], text.lstrip("\n")
    while True:
        m = LEADING.match(rest)
        if not m:
            break
        head.append(m.group(1))
        rest = rest[m.end():]
    head = ("\n".join(head) + "\n") if head else ""
    return HEAD.format(lang=lang, head=head) + rest.lstrip("\n") + TAIL


def strip(text):
    """Inverse of wrap: the page fragment, for hosts that supply their own skeleton.

    The head keeps whatever the page itself declared (its title, fonts and styles) and loses
    only the skeleton this module added, so wrap(strip(x)) == x.
    """
    if not is_document(text):
        return text
    skeleton = HEAD.format(lang="", head="")
    injected = skeleton[skeleton.index("<meta charset"):skeleton.index("</head>")]
    head = text[text.index("<head>") + 6:text.index("</head>")].replace(injected, "").strip("\n")
    body = text[text.index("<body>") + 6:text.rindex("</body>")].strip("\n")
    return (head + "\n" + body if head else body) + "\n"


if __name__ == "__main__":
    import pathlib, sys
    mode, path = sys.argv[1], pathlib.Path(sys.argv[2])
    out = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else path
    lang = sys.argv[4] if len(sys.argv) > 4 else ("ko" if "_ko" in path.name or "ko" in path.stem else "en")
    t = path.read_text(encoding="utf-8")
    out.write_text(wrap(t, lang) if mode == "wrap" else strip(t), encoding="utf-8")
    print(f"{mode}: {out} ({out.stat().st_size:,} bytes)")
