"""Check GO counts with an independent, top-down computation.

go_core propagates each protein's terms up to their ancestors. Here we go the
other way: build each term's descendants from a child index and count proteins
with any annotation inside that set. Both must agree for every reported term,
and the filters must hold (no roots, no 'do not annotate' terms, depth limit,
no redundant parent kept next to its specific child).
Run: python3 -I tests/test_go.py [annotations file]
"""

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import eggsplorer as er  # noqa: E402
import go_core as gc  # noqa: E402


def main(path):
    rows, _ = er.read_annotations(path)
    gene_go = {r["query"]: er.split(r["GOs"]) for r in rows if er.is_set(r["GOs"])}
    res = gc.top_terms(gene_go, n=20)
    terms, alt = gc.load()

    children = defaultdict(set)
    for t, info in terms.items():
        if not info["obsolete"]:
            for p in info["parents"]:
                children[p].add(t)

    def descendants(t):
        out, stack = {t}, [t]
        while stack:
            for c in children[stack.pop()]:
                if c not in out:
                    out.add(c)
                    stack.append(c)
        return out

    mapped = {g: {gc.current(x, terms, alt) for x in ids} - {None} for g, ids in gene_go.items()}
    mapped = {g: s for g, s in mapped.items() if s}
    assert res["annotated"] == len(mapped)

    ancestors, depth = gc.make_helpers(terms)
    n = 0
    for ns, items in res["top"].items():
        assert len(items) <= 20
        assert [i["genes"] for i in items] == sorted((i["genes"] for i in items), reverse=True)
        for item in items:
            t = item["id"]
            down = descendants(t)
            expected = sum(1 for s in mapped.values() if s & down)
            assert item["genes"] == expected, (t, item["genes"], expected)
            assert terms[t]["ns"] == ns and t not in gc.ROOTS and not terms[t]["general"]
            assert depth(t) >= gc.MIN_DEPTH
            n += 1
        ids = {i["id"] for i in items}
        for item in items:  # a kept term has no kept-or-hidden descendant covering >= 90% of it
            for other in items:
                if other["id"] != item["id"] and item["id"] in ancestors(other["id"]):
                    assert other["genes"] < gc.REDUNDANT * item["genes"], (item["id"], other["id"])
        assert ids
    print(f"OK: {Path(path).name}: {res['annotated']} GO-annotated proteins, {n} top terms verified top-down")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ROOT / "test-data" / "synthetic.emapper.annotations")
