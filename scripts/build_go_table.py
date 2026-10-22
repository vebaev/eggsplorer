#!/usr/bin/env python3
"""Build the compact GO table bundled with the tool from a GO OBO file.

Usage: python build_go_table.py go.obo ../refdata_go/go_terms.tsv.gz

Keeps, for every term: id, name, namespace (BP/MF/CC), parents through is_a,
part_of and the genus of intersection_of (so terms defined only logically in
the editors' file still connect to the hierarchy), obsolete flag, replaced_by
alternative ids and whether GO flags the term as too general to annotate
(subset gocheck_do_not_annotate).
"""

import gzip
import sys

NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}


def parse(path):
    term = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                if term and term.get("id", "").startswith("GO:"):
                    yield term
                term = {"parents": [], "alt": [], "replaced_by": ""} if line == "[Term]" else None
                continue
            if term is None or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            value = value.split(" ! ")[0].strip()
            if key == "id":
                term["id"] = value
            elif key == "name":
                term["name"] = value
            elif key == "namespace":
                term["ns"] = NS.get(value, value)
            elif key == "is_a":
                term["parents"].append(value.split()[0])
            elif key == "relationship" and value.startswith("part_of "):
                term["parents"].append(value.split()[1])
            elif key == "intersection_of" and value.startswith("GO:"):
                term["parents"].append(value.split()[0])  # genus of a logical definition
            elif key == "is_obsolete" and value == "true":
                term["obsolete"] = True
            elif key == "replaced_by":
                term["replaced_by"] = value.split()[0]
            elif key == "alt_id":
                term["alt"].append(value)
            elif key == "subset" and value == "gocheck_do_not_annotate":
                term["general"] = True
    if term and term.get("id", "").startswith("GO:"):
        yield term


def main(src, dst):
    n = 0
    with gzip.open(dst, "wt", encoding="utf-8", compresslevel=9) as out:
        out.write("go_id\tname\tnamespace\tparents\tobsolete\treplaced_by\talt_ids\tgeneral\n")
        for t in parse(src):
            parents = sorted({p for p in t["parents"] if p.startswith("GO:") and p != t["id"]})
            out.write("\t".join([t["id"], t.get("name", ""), t.get("ns", ""), ",".join(parents),
                                 "1" if t.get("obsolete") else "0", t["replaced_by"], ",".join(t["alt"]),
                                 "1" if t.get("general") else "0"]) + "\n")
            n += 1
    print(f"{n} terms written to {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
