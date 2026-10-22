"""Cross-check the report numbers against an independent pandas computation.

Run: python3 -I tests/test_report.py [annotations file]
Defaults to test-data/synthetic.emapper.annotations. Needs pandas.
"""

import json
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import decoder_core as dc  # noqa: E402
import eggsplorer as er  # noqa: E402


def load(path):
    import io
    text = "".join(l for l in open(path, encoding="utf-8") if not l.startswith("##"))
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def check(path):
    df = load(path)
    rows, _ = er.read_annotations(path)
    d = er.summarise(rows, "en", 25)

    assert d["counts"]["proteins"] == len(df)
    for c in d["coverage"]:
        assert c["n"] == int((df[c["col"]] != "-").sum()), c

    cog = Counter()
    for v in df["COG_category"]:
        if v != "-":
            cog.update(set(v))
    assert {c["letter"]: c["n"] for c in d["cog"] if c["n"]} == dict(cog)

    kos = set(re.findall(r"K\d{5}", " ".join(df["KEGG_ko"])))
    assert d["counts"]["kos"] == len(kos)

    maps = Counter()
    for v in df["KEGG_Pathway"]:
        if v != "-":
            maps.update({"map" + p[-5:] for p in v.split(",")})
    assert {p["id"]: p["genes"] for p in d["pathways"]} == dict(maps)

    import csv as _csv, gzip as _gz
    moved = {r["ec"]: r["replaced_by"] for r in _csv.DictReader(_gz.open(ROOT / "refdata_ec" / "ec_names.tsv.gz", "rt"), delimiter="\t")
             if r["obsolete"] == "1" and r["replaced_by"] and "," not in r["replaced_by"]}

    def now(e):
        while e in moved:
            e = moved[e]
        return e
    ec1 = Counter()
    for v in df["EC"]:
        if v != "-":
            ec1.update({now(e).split(".")[0] for e in v.split(",")})
    assert {k: n for k, n in d["ec"]["nodes"].items() if "." not in k} == dict(ec1)

    pf = Counter()
    for v in df["PFAMs"]:
        if v != "-":
            pf.update(set(v.split(",")))
    assert d["pfam"]["top"][:10] == pf.most_common(10) or \
        [n for _, n in d["pfam"]["top"][:10]] == [n for _, n in pf.most_common(10)]

    # CAZy by substrate: recount genes per main group straight from the bundled table
    sub = {r["family"]: r for r in _csv.DictReader(open(ROOT / "refdata_cazy" / "cazy_substrates.tsv"), delimiter="\t")
           if r["excluded"] == "0"}
    main_genes, any_genes = defaultdict(set), defaultdict(set)
    for q, v in zip(df["query"], df["CAZy"]):
        fams = {m.group(1) + m.group(2) for m in (re.match(r"^(CBM|GH|GT|PL|CE|AA)(\d+)", x) for x in v.split(",")) if m}
        for fam in fams if v != "-" else ():
            if fam in sub:
                if sub[fam]["main_group"]:
                    main_genes[sub[fam]["main_group"]].add(q)
                for g in sub[fam]["groups"].split(","):
                    any_genes[g].add(q)
    got = {c["key"]: c for c in d["cazy"]["substrates"]}
    assert set(got) == {g for g in any_genes if any_genes[g]}, (set(got), set(any_genes))
    for g, c in got.items():
        assert c["main"] == len(main_genes[g]) and c["main"] + c["possible"] == len(any_genes[g]), (g, c)

    scores = dc.score_kos(kos)
    assert {f: e["score"] for f, e in d["decoder"]["evidence"].items()} == scores
    return d["counts"]


def check_cli(path):
    with tempfile.TemporaryDirectory() as t:
        out = Path(t) / "r.html"
        subprocess.run([sys.executable, "-I", str(ROOT / "eggsplorer.py"), "-i", str(path),
                        "-o", str(out), "--tables", str(Path(t) / "tables"), "--plotly", "cdn"], check=True)
        html = out.read_text()
        payload = re.search(r'<script type="application/json" id="report-data">(.*?)</script>', html, re.S).group(1)
        data = json.loads(payload)
        assert html.startswith("<!doctype html>") and "<title>" in html.split("<body>")[0]
        assert data["d"]["counts"]["proteins"] > 0
        assert len(list((Path(t) / "tables").glob("*.tsv"))) == 10


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "test-data" / "synthetic.emapper.annotations"
    counts = check(target)
    check_cli(target)
    print("OK:", target.name, counts)
