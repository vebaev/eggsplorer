#!/usr/bin/env python3
"""Write the synthetic test genome used by the tests and the Galaxy wrapper.

The file has the layout of an eggNOG-mapper 2.x annotations table, but every row is
invented: KOs come from the KEGG-Decoder function map (some functions complete, others
partial), GO terms, EC numbers, TC and CAZy families from the bundled reference tables.
It is not a real genome. Deterministic (fixed seed).
Usage: python scripts/make_test_data.py   (writes test-data/synthetic.*)
"""

import csv
import gzip
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import decoder_core as dc  # noqa: E402
import refdata as rd  # noqa: E402

N = 300
COMPLETE = ["glycolysis", "TCA Cycle", "F-type ATPase", "NADH-quinone oxidoreductase", "Cytochrome bd complex",
            "Flagellum", "Chemotaxis", "Sec-SRP", "Twin Arginine Targeting", "transporter: phosphate",
            "dissim nitrate reduction", "sulfur assimilation", "histidine", "tryptophan",
            "riboflavin biosynthesis", "MEP-DOXP pathway", "Copper transporter CopA"]
HEADER = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
          "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
          "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]
PFAMS = ["ABC_tran", "BPD_transp_1", "Response_reg", "HATPase_c", "HisKA", "Acetyltransf_1", "MFS_1", "AAA",
         "adh_short", "Methyltransf_11", "TetR_N", "GntR", "LysR_substrate", "HTH_1", "Glycos_transf_1"]
MAPS = ["00010", "00020", "00190", "00230", "00240", "00260", "00340", "00400", "00500", "00520", "00910",
        "00920", "01100", "01110", "01120", "01200", "01230", "02010", "02020", "02030", "02040", "03010", "03060"]
MODULES = ["M00001", "M00002", "M00009", "M00026", "M00023", "M00144", "M00157", "M00176", "M00178", "M00222"]
LEVELS = ["2|Bacteria"] * 2 + ["1239|Firmicutes"] * 5 + ["91061|Bacilli"] * 12 + ["186826|Lactobacillales"] * 3


def main():
    rng = random.Random(2026)
    fmap = dc.load_function_map()
    kos = []
    for f in COMPLETE:
        kos.extend(k for k in fmap[f] if k not in kos)
    others = sorted({k for f, ks in fmap.items() if f not in COMPLETE for k in ks} - set(kos))
    kos += rng.sample(others, 40)
    kos += rng.sample(kos, 15)  # a few second copies
    with gzip.open(ROOT / "refdata_go" / "go_terms.tsv.gz", "rt", encoding="utf-8") as fh:
        go = [r["go_id"] for r in csv.DictReader(fh, delimiter="\t") if r["obsolete"] == "0" and r["general"] == "0"]
    with gzip.open(ROOT / "refdata_ec" / "ec_names.tsv.gz", "rt", encoding="utf-8") as fh:
        ecs = [r["ec"] for r in csv.DictReader(fh, delimiter="\t") if r["obsolete"] == "0" and r["ec"].count(".") == 3
               and "-" not in r["ec"] and not r["ec"].split(".")[-1].startswith("n")]
    go = rng.sample(go, 400)
    tcs = [k + ".1.1" for k in sorted(rd.TC_FAMILIES)]
    cazy = sorted(k for k in rd.CAZY_FAMILIES)
    rows = []
    for i in range(N):
        q = f"SYN_{i + 1:05d}"
        r = dict.fromkeys(HEADER, "-")
        cog = rng.choice("CEGJKLMPSTVOHIF")
        r.update({"query": q, "seed_ortholog": f"0.{q}", "evalue": f"{rng.uniform(1, 9):.2f}e-{rng.randint(20, 150)}",
                  "score": f"{rng.uniform(60, 900):.1f}", "eggNOG_OGs": f"COG{rng.randint(1, 5000):04d}@1|root,COG{rng.randint(1, 5000):04d}@2|Bacteria",
                  "max_annot_lvl": rng.choice(LEVELS), "COG_category": cog,
                  "Description": "synthetic test protein", "PFAMs": ",".join(rng.sample(PFAMS, rng.randint(1, 2)))})
        if i < len(kos):
            r["KEGG_ko"] = "ko:" + kos[i]
            r["Description"] = f"synthetic protein for {kos[i]}"
            r["BRITE"] = "ko00000,ko01000"
            if rng.random() < 0.4:
                m = rng.choice(MAPS)
                r["KEGG_Pathway"] = f"ko{m},map{m}"
                r["KEGG_Module"] = rng.choice(MODULES)
                r["KEGG_Reaction"] = f"R{rng.randint(1, 9999):05d}"
                r["KEGG_rclass"] = f"RC{rng.randint(1, 3000):05d}"
            if rng.random() < 0.35:
                r["EC"] = rng.choice(ecs)
        elif rng.random() < 0.3:
            r["COG_category"] = "S"
            r["Description"] = "-"
            r["PFAMs"] = "-"
        if rng.random() < 0.45:
            r["GOs"] = ",".join(rng.sample(go, rng.randint(2, 6)))
        if rng.random() < 0.13:
            r["KEGG_TC"] = rng.choice(tcs)
        if rng.random() < 0.06:
            r["CAZy"] = rng.choice(cazy)
        rows.append(r)
    out = ROOT / "test-data" / "synthetic.emapper.annotations"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("## synthetic test data, not a real genome (scripts/make_test_data.py)\n## emapper-2.1.13\n##\n")
        fh.write("#" + "\t".join(HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(r[h] for h in HEADER) + "\n")
        fh.write("## 300 queries scanned\n")
    with open(ROOT / "test-data" / "synthetic.kofam.tsv", "w", encoding="utf-8") as fh:
        fh.write("#\tgene name\tKO\tthrshld\tscore\tE-value\tKO definition\n#\t---------\t------\t-------\t------\t---------\t-------------\n")
        for r in rows[:5]:
            fh.write(f"*\t{r['query']}\t{r['KEGG_ko'][3:]}\t100.00\t350.0\t1.1e-100\tsynthetic test entry\n")
    print(f"wrote {out.name}: {N} proteins, {len(kos)} with a KO")


if __name__ == "__main__":
    main()
