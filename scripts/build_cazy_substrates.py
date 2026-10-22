#!/usr/bin/env python3
"""Build the CAZy family -> substrate group table bundled with the tool.

Input: the dbCAN-sub family-substrate mapping (fam-substrate-mapping-08012023.tsv,
dbCAN3, Zheng et al. 2023, NAR 51:W115; https://bcb.unl.edu/dbCAN2/).
Usage: python build_cazy_substrates.py fam-substrate-mapping.tsv ../refdata_cazy/cazy_substrates.tsv

For every non-GT family it keeps the dbCAN high-level substrates, the display groups
they fall into, and the main group from refdata.CAZY_MAIN_SUBSTRATE. It stops if a
main group is not among the groups dbCAN lists for that family, or if a dbCAN
substrate term has no group (so new terms are reviewed, not silently dropped).
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import refdata as rd  # noqa: E402

TERM_GROUP = {
    "cellulose": "cellulose", "beta-glucan": "cellulose", "cellooligosaccharide": "cellulose",
    "glucooligosaccharide": "cellulose", "beta-glycan": "cellulose",
    "xylan": "hemicellulose", "xyloglucan": "hemicellulose", "beta-mannan": "hemicellulose", "mannan": "hemicellulose",
    "glucomannan": "hemicellulose", "galactomannan": "hemicellulose", "galactoxyloglucan": "hemicellulose",
    "pectin": "pectin", "arabinan": "pectin", "beta-galactan": "pectin", "arabinogalactan": "pectin",
    "arabinogalactan protein": "pectin",
    "starch": "starch", "glycogen": "starch", "alpha-glucan": "starch",
    "sucrose": "sucrose", "fructan": "sucrose", "palatinose": "sucrose",
    "trehalose": "trehalose",
    "chitin": "chitin", "chitosan": "chitin", "chitooligosaccharide": "chitin",
    "peptidoglycan": "peptidoglycan",
    "host glycan": "host", "sialic acid": "host", "human milk polysaccharide": "host", "fucose": "host",
    "alginate": "algal", "agarose": "algal", "carrageenan": "algal", "porphyran": "algal", "ulvan": "algal",
    "fucoidan": "algal",
    "lignin": "lignin",
    "polyphenol": "other", "glycolipid": "other", "exo-polysaccharide": "other", "galactosaminogalactan": "other",
    "gellan": "other", "xanthan": "other", "beta-glucuronan": "other", "abscisic acid": "other", "alkaloid": "other",
    "cephalosporin C": "other", "cutin": "other", "uric acid": "other", "glucosylglycerate": "other",
    "glucosylglycerol": "other", "alpha-mannan": "other", "alpha-galactan": "other", "raffinose": "other",
    "rhamnose": "other", "alpha-rhamnoside": "other", "beta-fucosides": "other", "glycosidases": "other",
    "lipopolysaccharide": "other",
}


def terms(text):
    for part in re.split(r",| and ", text):
        part = re.sub(r"^and ", "", part.strip())
        if part:
            yield part


def main(src, dst):
    substrates = defaultdict(set)
    with open(src, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) >= 3 and row[2].strip():
                substrates[row[2].strip()].update(terms(row[0]))
    unknown = sorted({t for v in substrates.values() for t in v} - set(TERM_GROUP))
    if unknown:
        sys.exit(f"dbCAN terms without a group: {unknown}")
    keys = {k for k, _, _ in rd.CAZY_SUBSTRATE_GROUPS}
    order = [k for k, _, _ in rd.CAZY_SUBSTRATE_GROUPS]
    rows, problems = [], []
    for fam in sorted(substrates, key=lambda f: (re.match(r"[A-Z]+", f).group(), int(re.search(r"\d+", f).group()))):
        if fam.startswith("GT"):
            continue
        groups = {TERM_GROUP[t] for t in substrates[fam]}
        main_group = rd.CAZY_MAIN_SUBSTRATE.get(fam, "")
        if main_group and (main_group not in keys or main_group not in groups):
            problems.append(f"{fam}: main group {main_group!r} not in dbCAN groups {sorted(groups)}")
        rows.append([fam, main_group, ",".join(g for g in order if g in groups), "; ".join(sorted(substrates[fam])),
                     "1" if fam in rd.CAZY_SUBSTRATE_EXCLUDE else "0"])
    missing = sorted(set(rd.CAZY_MAIN_SUBSTRATE) - set(substrates))
    if missing:
        problems.append(f"main group given for families not in dbCAN: {missing}")
    if problems:
        sys.exit("\n".join(problems))
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "main_group", "groups", "dbcan_substrates", "excluded"])
        w.writerows(rows)
    print(f"{len(rows)} families written to {dst}; {sum(1 for r in rows if r[1])} with a main group")


if __name__ == "__main__":
    main(*sys.argv[1:3])
