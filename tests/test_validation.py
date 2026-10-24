"""Cell overview and element cycles on twenty reference genomes with well-known biology.

The KO assignments in test-data/validation come from eggNOG-mapper runs on the RefSeq proteomes
(see test-data/validation/README.md). For each genome the test lists components that must be drawn,
components that must not be drawn, cycle steps that must light up and steps that must stay absent
or pale. It guards the rules in the report template against regressions; it is not a benchmark of
eggNOG-mapper. Needs Playwright with Chromium (skipped otherwise).
Run: python3 -I tests/test_validation.py
"""

import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "test-data" / "validation"
HEADER = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
          "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
          "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]

# draw / nodraw: cell components; lit: cycle steps drawn at >= 0.50 with their catalytic subunit;
# off: steps absent or pale; tile: text expected in a tile; shared: steps reported as resting on the same proteins
EXPECT = {
    "GCF_000005845.2": dict(  # Escherichia coli K-12 MG1655
        envelope="Diderm", draw=["NDH-1", "NDH-2", "bo3 oxidase", "bd oxidase", "F-type ATPase", "Tat", "flagellum"],
        nodraw=["bc1 / bcc", "aa3 oxidase", "cbb3 oxidase", "Na-NQR", "T3SS", "competence (com)"],
        lit=["NO₃⁻ reduction", "DNRA", "sulfate activation", "APS / PAPS reduction", "assimilatory reduction", "H₂ oxidation"],
        off=["N₂ fixation", "NO₂⁻ reduction", "N₂O reduction", "thiosulfate reduction", "methanogenesis", "CBB cycle", "RuMP pathway"],
        shared=["NO₂⁻ oxidation / NO₃⁻ reduction"]),
    "GCF_000009045.1": dict(  # Bacillus subtilis 168
        envelope="Monoderm", draw=["NDH-2", "bc1 / bcc", "aa3 oxidase", "aa3-600 oxidase", "bd oxidase", "flagellum", "competence (com)"],
        nodraw=["NDH-1", "cbb3 oxidase", "T3SS"],
        lit=["NO₃⁻ reduction", "NO₂⁻ assimilation", "sulfate activation", "assimilatory reduction"],
        off=["DNRA", "N₂ fixation", "N₂O reduction", "methanogenesis"]),
    "GCF_000006765.1": dict(  # Pseudomonas aeruginosa PAO1
        envelope="Diderm", draw=["NDH-1", "aa3 oxidase", "cbb3 oxidase", "bo3 oxidase", "bd oxidase", "T2SS", "T3SS", "T6SS", "flagellum"],
        nodraw=["aa3-600 oxidase"],
        lit=["NO₃⁻ reduction", "NO₂⁻ reduction", "NO reduction", "N₂O reduction", "NO₂⁻ assimilation"],
        off=["N₂ fixation", "DNRA", "methanogenesis", "CBB cycle"]),
    "GCF_000021045.1": dict(  # Azotobacter vinelandii DJ
        envelope="Diderm", draw=["NDH-1", "bd oxidase", "flagellum"],
        lit=["N₂ fixation", "NO₃⁻ assimilation", "NO₂⁻ assimilation"],
        off=["DNRA", "NO₂⁻ reduction", "N₂O reduction", "methanogenesis"]),
    "GCF_000009145.1": dict(  # Nitrosomonas europaea ATCC 19718
        envelope="Diderm", draw=["NDH-1", "bc1 / bcc"],
        lit=["NH₃ oxidation", "NH₂OH oxidation", "NO₂⁻ reduction", "CBB cycle"],
        off=["N₂ fixation", "methanogenesis"], shared=["NH₃ oxidation / CH₄ oxidation"]),
    "GCF_000012725.1": dict(  # Nitrobacter winogradskyi Nb-255
        envelope="Diderm", draw=["NDH-1", "bc1 / bcc", "aa3 oxidase"],
        lit=["NO₂⁻ oxidation", "CBB cycle"], off=["N₂ fixation", "NH₃ oxidation", "methanogenesis"],
        shared=["NO₂⁻ oxidation / NO₃⁻ reduction"]),
    "GCF_000195755.1": dict(  # Desulfovibrio vulgaris Hildenborough
        envelope="Diderm", draw=["flagellum", "aa3 oxidase"], nodraw=["NDH-1", "bo3 oxidase"],
        lit=["sulfate ⇄ APS", "APS ⇄ sulfite", "sulfite ⇄ sulfide", "DNRA", "H₂ oxidation"],
        off=["methanogenesis", "thiosulfate oxidation", "CBB cycle"],
        tile={"Hydrogenases": "NiFe"}, notile=["Methanogenesis"]),
    "GCF_000025485.1": dict(  # Allochromatium vinosum DSM 180
        envelope="Diderm", draw=["NDH-1", "flagellum"],
        lit=["thiosulfate oxidation", "to tetrathionate", "sulfite ⇄ sulfide", "APS ⇄ sulfite", "sulfide oxidation", "CBB cycle", "N₂ fixation"],
        off=["methanogenesis", "thiosulfate reduction"], tile={"Phototrophy": "RC-II"}),
    "GCF_000009725.1": dict(  # Synechocystis sp. PCC 6803
        envelope="Diderm", draw=["NDH-1", "b6f", "bd oxidase"], nodraw=["flagellum", "bc1 / bcc"],
        lit=["NO₃⁻ assimilation", "NO₂⁻ assimilation", "CBB cycle", "assimilatory reduction"],
        off=["N₂ fixation", "DNRA", "methanogenesis"], tile={"Phototrophy": "PSII", "C fixation and C1": "CBB"}),
    "GCF_000008325.1": dict(  # Methylococcus capsulatus Bath
        envelope="Diderm", draw=["NDH-1", "bc1 / bcc"], nodraw=["flagellum"],
        lit=["CH₄ oxidation", "CH₃OH oxidation", "RuMP pathway", "CBB cycle", "N₂ fixation"],
        off=["methanogenesis"], shared=["NH₃ oxidation / CH₄ oxidation"], tile={"C fixation and C1": "sMMO"}),
    "GCF_000006985.1": dict(  # Chlorobaculum tepidum TLS (green sulfur bacterium; RuBisCO-like protein, no PRK)
        envelope="Diderm", draw=["bd oxidase"], nodraw=["flagellum"],
        lit=["rTCA cycle", "N₂ fixation", "thiosulfate oxidation", "sulfide oxidation", "sulfite ⇄ sulfide"],
        off=["CBB cycle", "methanogenesis"], tile={"Phototrophy": "RC-I", "C fixation and C1": "rTCA"}),
    "GCF_000143685.1": dict(  # Clostridium ljungdahlii DSM 13528 (acetogen)
        envelope="Monoderm", draw=["Rnf", "F-type ATPase", "flagellum"], nodraw=["NDH-1", "bd oxidase", "aa3 oxidase"],
        lit=["Wood–Ljungdahl"], off=["methanogenesis", "CBB cycle"], tile={"C fixation and C1": "WL"}, notile=["Methanogenesis"]),
    "GCF_011066545.1": dict(  # Candidatus Kuenenia stuttgartiensis CSTR1 (anammox; hzs/hdh have no KO in eggNOG)
        envelope="Diderm", lit=["NO₂⁻ reduction"], off=["methanogenesis", "N₂ fixation"]),
    "GCF_001458695.1": dict(  # Candidatus Nitrospira inopinata (comammox)
        envelope="Diderm", lit=["NH₃ oxidation", "NH₂OH oxidation", "NO₂⁻ oxidation", "rTCA cycle"],
        off=["CBB cycle", "methanogenesis"], shared=["NH₃ oxidation / CH₄ oxidation", "NO₂⁻ oxidation / NO₃⁻ reduction"]),
    "GCF_000006885.1": dict(  # Streptococcus pneumoniae TIGR4 (no respiratory chain, competent)
        envelope="Monoderm", draw=["F-type ATPase", "competence (com)", "Sec-SRP"],
        nodraw=["NDH-1", "NDH-2", "bd oxidase", "aa3 oxidase", "flagellum", "T3SS"],
        off=["N₂ fixation", "NO₃⁻ reduction", "sulfate activation", "methanogenesis"]),
    "GCF_000091545.1": dict(  # Thermus thermophilus HB8 (V/A-type ATP synthase only)
        draw=["V-type ATPase", "NDH-1"], nodraw=["F-type ATPase", "flagellum"], off=["N₂ fixation", "methanogenesis"],
        envelope="Monoderm"),
    "GCF_000195955.2": dict(  # Mycobacterium tuberculosis H37Rv
        envelope="Monoderm", draw=["NDH-1", "NDH-2", "bc1 / bcc", "bd oxidase", "T7SS (ESX)"], nodraw=["flagellum", "T3SS"],
        lit=["NO₃⁻ reduction", "NO₂⁻ assimilation"], off=["DNRA", "N₂ fixation", "methanogenesis"]),
    "GCF_000027325.1": dict(  # Mycoplasma genitalium G37 (wall-less, minimal)
        envelope="No cell wall", draw=["F-type ATPase"], nodraw=["NDH-1", "NDH-2", "flagellum", "Tat", "bd oxidase"],
        off=["N₂ fixation", "NO₃⁻ reduction", "sulfate activation", "CBB cycle", "methanogenesis"]),
    "GCF_000011065.1": dict(  # Bacteroides thetaiotaomicron VPI-5482
        envelope="Diderm", draw=["Na-NQR", "Rnf", "SusCD (PUL)", "bd oxidase"], nodraw=["flagellum"],
        lit=["DNRA"], off=["N₂ fixation", "methanogenesis"]),
    "GCF_000009085.1": dict(  # Campylobacter jejuni NCTC 11168
        envelope="Diderm", draw=["NDH-1", "bc1 / bcc", "cbb3 oxidase", "flagellum"], nodraw=["Na-NQR", "T3SS"],
        lit=["NO₃⁻ reduction", "DNRA", "H₂ oxidation"], off=["N₂ fixation", "methanogenesis"]),
}

JS = """() => {
  const svg = document.getElementById('ce-svg');
  const drawn = [...svg.querySelectorAll('g > title')].map(t => t.textContent.split(':')[0]);
  const tiles = [...svg.querySelectorAll('.ce-tile')].map(g => [...g.querySelectorAll('text')].map(t => t.textContent));
  const edges = {}; for (const k in window.__cycleEdges) window.__cycleEdges[k].forEach(e => { edges[e.key] = [e.s, e.nocore]; });
  return {drawn, tiles, edges, shared: window.__cycleShared};
}"""


def annotations(acc, path):
    with gzip.open(VAL / f"{acc}.kos.tsv.gz", "rt", encoding="utf-8") as src, open(path, "w", encoding="utf-8") as out:
        next(src)
        out.write("#" + "\t".join(HEADER) + "\n")
        for line in src:
            q, kos = line.rstrip("\n").split("\t")
            row = dict.fromkeys(HEADER, "-")
            row.update({"query": q, "COG_category": "S", "Description": "-", "KEGG_ko": ",".join("ko:" + k for k in kos.split(","))})
            out.write("\t".join(row[h] for h in HEADER) + "\n")


def check(acc, res, want):
    errs = []
    env = [t for t in res["tiles"] if t[0] in ("Diderm envelope", "Monoderm envelope", "No cell wall")]
    if not env or not env[0][0].startswith(want["envelope"]):
        errs.append(f"envelope {env}")
    errs += [f"not drawn: {n}" for n in want.get("draw", []) if n not in res["drawn"]]
    errs += [f"drawn: {n}" for n in want.get("nodraw", []) if n in res["drawn"]]
    for n in want.get("lit", []):
        s, nocore = res["edges"].get(n, [0, ""])
        if not (s or 0) >= 0.5 or nocore:
            errs.append(f"step not lit: {n} {s} {nocore}")
    for n in want.get("off", []):
        s, nocore = res["edges"].get(n, [0, ""])
        if (s or 0) >= 0.5 and not nocore:
            errs.append(f"step lit: {n} {s}")
    for prefix in want.get("shared", []):
        if not any(x.startswith(prefix) for x in res["shared"]):
            errs.append(f"not reported as shared: {prefix}")
    tiles = {t[0]: t[1] if len(t) > 1 else "" for t in res["tiles"]}
    for name, text in want.get("tile", {}).items():
        if text not in tiles.get(name, ""):
            errs.append(f"tile {name}: {tiles.get(name)!r} lacks {text!r}")
    errs += [f"tile drawn: {n}" for n in want.get("notile", []) if n in tiles]
    assert not errs, (acc, errs)


def main():
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("validation skipped (Playwright not installed)")
        return
    n_checks = 0
    with tempfile.TemporaryDirectory() as t, sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        for acc, want in EXPECT.items():
            src, out = Path(t) / f"{acc}.annotations", Path(t) / f"{acc}.html"
            annotations(acc, src)
            subprocess.run([sys.executable, "-I", str(ROOT / "eggsplorer.py"), "-i", str(src), "-o", str(out), "--plotly", "cdn"],
                           check=True, stdout=subprocess.DEVNULL)
            page.goto(out.as_uri())
            page.wait_for_function("window.__cycleEdges !== undefined && window.__cellQA !== undefined", timeout=30000)
            check(acc, page.evaluate(JS), want)
            n_checks += sum(len(v) if isinstance(v, (list, dict)) else 1 for v in want.values())
        browser.close()
    assert not errors, errors
    print(f"OK: {len(EXPECT)} reference genomes, {n_checks} expectations on the cell overview and element cycles")


if __name__ == "__main__":
    main()
