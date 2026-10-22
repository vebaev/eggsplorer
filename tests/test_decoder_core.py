"""Check decoder_core against the original KEGG-Decoder run as a script.

Builds synthetic eggNOG-mapper files from random subsets of the KOs that
KEGG-Decoder inspects, runs both implementations and requires identical
matrices. Run: python3 -I tests/test_decoder_core.py
"""

import csv
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import decoder_core as dc  # noqa: E402

HEADER = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\t"
          "Description\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\t"
          "KEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs")


def write_emapper(path, kos, rng):
    kos = sorted(kos)
    rng.shuffle(kos)
    with open(path, "w") as fh:
        fh.write("## emapper-2.1.12\n## synthetic test file\n" + HEADER + "\n")
        i = 0
        while i < len(kos):
            chunk = kos[i:i + rng.choice([1, 1, 2, 3])]  # some genes carry several KOs
            i += len(chunk)
            ko_field = ",".join("ko:" + k for k in chunk)
            fh.write(f"gene_{i}\t1.x\t1e-50\t200\t-\t-\tC\tdesc\t-\t-\t-\t{ko_field}"
                     + "\t-" * 9 + "\n")
        fh.write("gene_none\t1.x\t1e-50\t200\t-\t-\tS\t-\t-\t-\t-\t-" + "\t-" * 9 + "\n")
        fh.write("## 10 queries scanned\n")


def run_original(samples, tmp):
    koala = tmp / "koala.txt"
    with open(koala, "w") as fh:
        for name, kos in samples.items():
            for j, k in enumerate(sorted(kos)):
                fh.write(f"{name}_{j}\t{k}\n")
    out = tmp / "original.tsv"
    subprocess.run([sys.executable, "-I", str(dc.DECODER_PATH), "-i", str(koala), "-o", str(out)],
                   check=True, cwd=tmp)
    return read_matrix(out)


def read_matrix(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = [h.strip() for h in rows[0]]
    return {r[0].strip(): {h: float(v) for h, v in zip(header[1:], r[1:])} for r in rows[1:]}


def main():
    rng = random.Random(42)
    universe = sorted(set().union(*dc.SCORER_KOS.values()))
    extra = [f"K{n:05d}" for n in rng.sample(range(1, 26000), 300)]
    fractions = {"Empty": 0.0, "Sparse": 0.1, "Half": 0.5, "Dense": 0.85, "Full": 1.0}
    samples = {n: set(rng.sample(universe, int(len(universe) * f))) | set(extra[:50])
               for n, f in fractions.items()}

    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        files = []
        for name, kos in samples.items():
            f = tmp / f"{name}.emapper.annotations"
            write_emapper(f, kos, rng)
            assert dc.read_emapper_kos(f) == kos, f"parser lost KOs for {name}"
            files.append(f)

        ours_path = tmp / "ours.tsv"
        dc.main(["-i", *map(str, files), "-n", *samples, "-o", str(ours_path)])
        ours = read_matrix(ours_path)
        # The original cannot score a genome with no KOs, so compare the rest.
        original = run_original({k: v for k, v in samples.items() if v}, tmp)

    assert list(next(iter(ours.values()))) == dc.FUNCTION_ORDER
    for name in original:
        diff = {f: (ours[name][f], original[name][f])
                for f in dc.FUNCTION_ORDER if ours[name][f] != original[name][f]}
        assert not diff, f"{name}: {list(diff.items())[:5]}"
    print(f"OK: {len(dc.FUNCTION_ORDER)} functions, {len(dc.SCORERS)} scoring blocks, "
          f"{len(universe)} KOs; {len(original)} genomes identical to original KEGG-Decoder")


if __name__ == "__main__":
    main()
