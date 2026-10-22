"""The bundled CAZy and TC names must match the reviewed names in refdata_checks/.

Each check file records, for every name, the evidence from its CAZy or TCDB
page, the verdict and the final label shipped in refdata.py.
Run: python3 -I tests/test_refdata.py
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import refdata as rd  # noqa: E402


def check(name, key, table):
    with open(ROOT / "refdata_checks" / name, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r[key] for r in rows} == set(table), f"{name}: codes differ from refdata"
    for r in rows:
        assert r["verdict"] in ("OK", "FIX"), r
        assert table[r[key]] == r["final_label"], (r[key], table[r[key]], r["final_label"])
    return len(rows)


def check_substrates():
    """Main substrate groups must be among the groups dbCAN lists for the family."""
    keys = [k for k, _, _ in rd.CAZY_SUBSTRATE_GROUPS]
    with open(ROOT / "refdata_cazy" / "cazy_substrates.tsv", encoding="utf-8") as fh:
        rows = {r["family"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for fam, main_group in rd.CAZY_MAIN_SUBSTRATE.items():
        assert rows[fam]["main_group"] == main_group, fam
        assert main_group in rows[fam]["groups"].split(","), (fam, main_group, rows[fam]["groups"])
    for fam, r in rows.items():
        assert not fam.startswith("GT"), fam
        assert r["groups"] and set(r["groups"].split(",")) <= set(keys), (fam, r["groups"])
        assert (r["excluded"] == "1") == (fam in rd.CAZY_SUBSTRATE_EXCLUDE), fam
    return len(rows)


def main():
    n1 = check("cazy_check_2026-10.tsv", "family", rd.CAZY_FAMILIES)
    n2 = check("tc_check_2026-10.tsv", "code", {**rd.TC_SUBCLASSES, **rd.TC_FAMILIES})
    n3 = check_substrates()
    print(f"OK: {n1} CAZy and {n2} TC names match the reviewed check files; {n3} CAZy substrate rows consistent")


if __name__ == "__main__":
    main()
