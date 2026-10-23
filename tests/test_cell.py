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
CYCLE_JS = TEMPLATE[TEMPLATE.index("/* Element cycles:"):TEMPLATE.index("/* CAZy by substrate")]
HEADER = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
          "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
          "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]


def check_rules():
    fmap = dc.load_function_map()
    names = set(re.findall(r'fn: "([^"]+)"', CELL_JS))
    for best in re.findall(r"best: \[([^\]]+)\]", CELL_JS):
        names.update(re.findall(r'"([^"]+)"', best))
    unknown = sorted(n for n in names if n not in fmap)
    assert not unknown, f"functions not in the KEGG-Decoder map: {unknown}"
    kos = set(re.findall(r'"(K\d{5})"', CELL_JS))
    known = set(rd.CELL_EXTRA_KOS) | set(rd.CYCLE_EXTRA_KOS)  # the cell also reads a few cycle KOs (C1 tile)
    assert kos <= known, sorted(kos - known)
    cnames = set(re.findall(r'S\("([^"]+)"\)', CYCLE_JS))
    for grp in re.findall(r"(?:best|allOf)\(\[([^\]]+)\]\)", CYCLE_JS):
        cnames.update(re.findall(r'"([^"]+)"', grp))
    unknown = sorted(n for n in cnames if n not in fmap)
    assert not unknown, f"cycle functions not in the KEGG-Decoder map: {unknown}"
    ckos = set(re.findall(r'"(K\d{5})"', CYCLE_JS))
    assert ckos == set(rd.CYCLE_EXTRA_KOS), sorted(ckos ^ set(rd.CYCLE_EXTRA_KOS))
    return len(names), len(kos)


def worst_case(path):
    """Annotations in which every function of the map is complete, LPS markers are present
    and more transporter families than there are places."""
    fmap = dc.load_function_map()
    kos = sorted({k for ks in fmap.values() for k in ks} | set(rd.CELL_EXTRA_KOS) | set(rd.CYCLE_EXTRA_KOS))
    header = HEADER
    tcs = ["3.A.1.1.1", "2.A.1.1.1", "1.B.1.1.1", "2.A.6.2.1", "3.A.3.1.1", "2.A.3.1.1", "2.A.7.1.1", "1.A.8.1.1", "3.D.1.1.1"]
    with open(path, "w", encoding="utf-8") as out:
        out.write("#" + "\t".join(header) + "\n")
        for i, k in enumerate(kos):
            row = dict.fromkeys(header, "-")
            row.update({"query": f"g{i:05d}", "COG_category": "S", "Description": "synthetic", "KEGG_ko": "ko:" + k,
                        "KEGG_TC": tcs[i % len(tcs)] if i < 200 else "-"})
            out.write("\t".join(row[h] for h in header) + "\n")


def write_genes(path, genes):
    """Annotations with one row per (gene, [KOs])."""
    with open(path, "w", encoding="utf-8") as out:
        out.write("#" + "\t".join(HEADER) + "\n")
        for g, ks in genes:
            row = dict.fromkeys(HEADER, "-")
            row.update({"query": g, "COG_category": "S", "Description": "synthetic", "KEGG_ko": ",".join("ko:" + k for k in ks)})
            out.write("\t".join(row[h] for h in HEADER) + "\n")


SULFUR_CASES = [
    # proteobacterial assimilation (cysDN, cysC, cysH, cysJI), TsdA, and one protein carrying both phsA and soeA KOs
    ([("d", ["K00957"]), ("n", ["K00956"]), ("c", ["K00860"]), ("h", ["K00390"]), ("j", ["K00380"]), ("i", ["K00381"]),
      ("t", ["K19713"]), ("x", ["K08352", "K21307"])],
     {"sulfate activation": 1, "APS / PAPS reduction": 1, "assimilatory reduction": 1, "to tetrathionate": 1,
      "thiosulfate oxidation": 0, "sulfate ⇄ APS": 0}, True),
    # sat only, no cysH: activation (and the dissimilatory sat arrow) light, reduction of APS does not;
    # sat serving both sulfate arrows is shared by design and not reported
    ([("s", ["K00958"])], {"sulfate activation": 1, "APS / PAPS reduction": 0, "sulfate ⇄ APS": 1}, False),
    # cysD without cysN: half an ATP sulfurylase; Sox complete goes to sulfate, not to tetrathionate
    ([("d", ["K00957"])] + [(f"sox{i}", [k]) for i, k in enumerate(["K17222", "K17223", "K17224", "K17225", "K17226", "K17227"])],
     {"sulfate activation": 0.5, "thiosulfate oxidation": 1, "to tetrathionate": 0}, False),
]


def sulfur_rules(tmp):
    from playwright.sync_api import sync_playwright
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        for i, (genes, want, shared) in enumerate(SULFUR_CASES):
            src, out = Path(tmp) / f"s{i}.annotations", Path(tmp) / f"s{i}.html"
            write_genes(src, genes)
            subprocess.run([sys.executable, "-I", str(ROOT / "eggsplorer.py"), "-i", str(src), "-o", str(out), "--plotly", "cdn"], check=True)
            page.goto(out.as_uri())
            page.wait_for_function("window.__cycleEdges !== undefined", timeout=20000)
            edges = {e["name"]: (e["s"] or 0) for e in page.evaluate("window.__cycleEdges['cy-s']")}
            for name, s in want.items():
                assert abs(edges[name] - s) < 1e-9, (i, name, edges[name], s)
            got_shared = page.evaluate("window.__cycleShared")
            assert bool(got_shared) == shared, (i, got_shared)
        browser.close()
    return len(SULFUR_CASES)


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
        n_sulfur = sulfur_rules(t)
    print(msg + "; layout clean for " + ", ".join(f"{k} ({v} labels)" for k, v in results.items())
          + f"; inline report makes no network requests and uses {n_fonts} embedded font faces"
          + f"; {n_sulfur} sulfur-diagram cases give the expected arrows")


if __name__ == "__main__":
    main()
