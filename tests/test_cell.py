"""Cell overview and element cycles: rules and layout.

1. Every metabolic function named in the cell template exists in the KEGG-Decoder
   function map, and every KO the template uses directly is listed in refdata.
2. Layout (needs Playwright with Chromium; skipped otherwise): the cell is drawn for the
   test genome and for a synthetic worst case in which every component is present and the
   outer membrane is on. Labels must not overlap each other or the drawn components
   (rotation-aware), must stay inside the drawing, and tile text must fit its tile.
Run: python3 -I tests/test_cell.py
"""

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import decoder_core as dc  # noqa: E402
import refdata as rd  # noqa: E402

TEMPLATE = (ROOT / "report_template.html").read_text(encoding="utf-8")
CELL_JS = TEMPLATE[TEMPLATE.index("/* Cell overview:"):TEMPLATE.index("/* Element cycles:")]


def check_rules():
    fmap = dc.load_function_map()
    names = set(re.findall(r'fn: "([^"]+)"', CELL_JS))
    for best in re.findall(r"best: \[([^\]]+)\]", CELL_JS):
        names.update(re.findall(r'"([^"]+)"', best))
    unknown = sorted(n for n in names if n not in fmap)
    assert not unknown, f"functions not in the KEGG-Decoder map: {unknown}"
    kos = set(re.findall(r'"(K\d{5})"', CELL_JS))
    assert kos <= set(rd.CELL_EXTRA_KOS), sorted(kos - set(rd.CELL_EXTRA_KOS))
    return len(names), len(kos)


def worst_case(path):
    """Annotations in which every function of the map is complete, LPS markers are present
    and more transporter families than there are places."""
    fmap = dc.load_function_map()
    kos = sorted({k for ks in fmap.values() for k in ks} | set(rd.CELL_EXTRA_KOS))
    header = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
              "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
              "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]
    tcs = ["3.A.1.1.1", "2.A.1.1.1", "1.B.1.1.1", "2.A.6.2.1", "3.A.3.1.1", "2.A.3.1.1", "2.A.7.1.1", "1.A.8.1.1", "3.D.1.1.1"]
    with open(path, "w", encoding="utf-8") as out:
        out.write("#" + "\t".join(header) + "\n")
        for i, k in enumerate(kos):
            row = dict.fromkeys(header, "-")
            row.update({"query": f"g{i:05d}", "COG_category": "S", "Description": "synthetic", "KEGG_ko": "ko:" + k,
                        "KEGG_TC": tcs[i % len(tcs)] if i < 200 else "-"})
            out.write("\t".join(row[h] for h in header) + "\n")


def layout(report):
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        return None
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1300, "height": 900})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(Path(report).as_uri())
        page.wait_for_function("window.__cellQA !== undefined", timeout=20000)
        page.evaluate("document.fonts ? document.fonts.ready : null")
        res = page.evaluate("window.__cellQA()")
        res["cycles"] = page.evaluate("window.__cycleQA()")
        browser.close()
    assert not errors, errors
    return res


def offline(report):
    """The default (inline) report must make no network requests and use the embedded fonts."""
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        external = []
        page.on("request", lambda r: external.append(r.url) if not r.url.startswith(("file:", "data:", "blob:")) else None)
        page.goto(Path(report).as_uri())
        page.wait_for_function("window.__cellQA !== undefined", timeout=20000)
        page.evaluate("document.fonts.ready.then(() => 1)")
        fonts = page.evaluate("['IBM Plex Sans', 'IBM Plex Sans Condensed', 'IBM Plex Mono'].map(f => document.fonts.check('600 12px \"' + f + '\"') || document.fonts.check('12px \"' + f + '\"'))")
        loaded = page.evaluate("[...document.fonts].filter(f => f.status === 'loaded').map(f => f.family)")
        browser.close()
    assert not external, external
    assert all(fonts) and loaded, (fonts, loaded)
    return len(loaded)


def main():
    n_fn, n_ko = check_rules()
    msg = f"OK: cell rules use {n_fn} known functions and {n_ko} listed KOs"
    with tempfile.TemporaryDirectory() as t:
        worst = Path(t) / "worst.emapper.annotations"
        worst_case(worst)
        results = {}
        for name, src, lang in (("test genome", ROOT / "test-data" / "synthetic.emapper.annotations", "en"),
                                ("test genome (bg)", ROOT / "test-data" / "synthetic.emapper.annotations", "bg"),
                                ("worst case", worst, "en"), ("worst case (bg)", worst, "bg")):
            out = Path(t) / (re.sub(r"\W+", "_", name) + ".html")
            subprocess.run([sys.executable, "-I", str(ROOT / "eggsplorer.py"), "-i", str(src), "-o", str(out),
                            "--plotly", "cdn", "--lang", lang], check=True)
            res = layout(out)
            if res is None:
                print(msg + "; layout check skipped (Playwright not installed)")
                return
            assert res["overlaps"] == 0 and res["outside"] == 0 and res["tile_overflow"] == 0, (name, res)
            for cid, c in res["cycles"].items():
                assert c["overlaps"] == 0 and c["outside"] == 0, (name, cid, c)
            results[name] = res["labels"]
        inline = Path(t) / "inline.html"
        subprocess.run([sys.executable, "-I", str(ROOT / "eggsplorer.py"), "-i",
                        str(ROOT / "test-data" / "synthetic.emapper.annotations"), "-o", str(inline)], check=True)
        n_fonts = offline(inline)
    print(msg + "; layout clean for " + ", ".join(f"{k} ({v} labels)" for k, v in results.items())
          + f"; inline report makes no network requests and uses {n_fonts} embedded font faces")


if __name__ == "__main__":
    main()
