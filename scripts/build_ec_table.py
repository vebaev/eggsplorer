#!/usr/bin/env python3
"""Build the compact EC name table bundled with the tool.

Input: the ENZYME nomenclature (ExPASy, CC BY 4.0) as OBO, e.g. the PyOBO
export at https://github.com/biopragmatics/obo-db-ingest (export/ec/ec.obo).
Usage: python build_ec_table.py ec.obo ../refdata_ec/ec_names.tsv.gz

Keeps id (without the 'EC:' prefix), name, whether the entry is obsolete
(deleted or transferred) and, for transferred entries, the current EC number(s),
for classes, subclasses, sub-subclasses and enzymes.
"""

import gzip
import re
import sys


def unescape(text):
    return re.sub(r"\\(.)", r"\1", text)


def main(src, dst):
    rows, term = [], None
    for line in open(src, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("["):
            if term and term.get("id"):
                rows.append(term)
            term = {} if line == "[Term]" else None
            continue
        if term is None:
            continue
        if line.startswith("id: EC:"):
            term["id"] = line[7:].strip()
        elif line.startswith("name: "):
            term["name"] = unescape(line[6:].strip())
        elif line.startswith("is_obsolete: true"):
            term["obsolete"] = True
        elif line.startswith("replaced_by: EC:"):
            term.setdefault("replaced_by", []).append(line[16:].strip())
    if term and term.get("id"):
        rows.append(term)
    with gzip.open(dst, "wt", encoding="utf-8", compresslevel=9) as out:
        out.write("ec\tname\tobsolete\treplaced_by\n")
        for t in rows:
            out.write(f"{t['id']}\t{t.get('name', '')}\t{1 if t.get('obsolete') else 0}\t"
                      f"{','.join(t.get('replaced_by', []))}\n")
    print(f"{len(rows)} EC entries written to {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
