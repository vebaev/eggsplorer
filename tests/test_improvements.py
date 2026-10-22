"""Tests for the additions on top of KEGG-Decoder.

Run: python3 -I tests/test_improvements.py
"""

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import decoder_core as dc  # noqa: E402

HEADER = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\t"
          "Description\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\t"
          "KEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs")


def emapper_row(gene, kos, evalue="1e-80", score="300"):
    ko = ",".join("ko:" + k for k in kos) if kos else "-"
    return f"{gene}\t1.x\t{evalue}\t{score}\t-\t-\tC\tdesc\t-\t-\t-\t{ko}" + "\t-" * 9 + "\n"


def write(path, text):
    Path(path).write_text(text)
    return path


def test_function_map():
    fmap = dc.build_function_map(seed=0)
    assert fmap == dc.build_function_map(seed=3), "map depends on the random seed"
    unassigned = set().union(*map(set, dc.unassigned_kos(fmap).values()))
    # only the alcohol dehydrogenases hit by the '==' bug are never counted
    assert unassigned == set(dc._ADH_KOS) - {"K04072"}, unassigned
    assert [f for f, k in fmap.items() if not k] == list(dc.NOT_EVALUATED)
    assert fmap["nitrogen fixation"] == ["K02586", "K02588", "K02591"]
    assert set(fmap["glycine"]) == {"K00600", "K01620"}


def test_evidence_and_copies(tmp):
    glyco = ["K01835", "K01810", "K01623", "K00927", "K01689", "K00134", "K01834", "K00873"]
    text = "## test\n" + HEADER + "\n"
    text += "".join(emapper_row(f"g{i}", [k]) for i, k in enumerate(glyco))
    text += emapper_row("pfkA", ["K00850"]) + emapper_row("pfkB", ["K00850"])  # two copies
    text += emapper_row("weak", ["K01647"], evalue="1e-3", score="20")         # citrate synthase
    f = write(tmp / "a.emapper.annotations", text)

    all_hits = dc.read_emapper(f)
    filtered = dc.read_emapper(f, max_evalue=1e-10, min_score=50)
    assert "K01647" in all_hits and "K01647" not in filtered

    ev = dc.function_evidence(filtered)["glycolysis"]
    assert ev["score"] == 1.0 and ev["genes"] == 10
    pfk = next(x for x in ev["found"] if x["ko"] == "K00850")
    assert pfk["copies"] == 2 and pfk["genes"] == ["pfkA", "pfkB"]
    assert "K00895" in ev["missing"] and ev["gaps"] == []  # alternative, not a real gap

    tca = dc.function_evidence(filtered)["TCA Cycle"]
    assert "K01647" in tca["gaps"]  # finding citrate synthase would raise the score


def test_kofam_and_modes(tmp):
    egg = write(tmp / "b.emapper.annotations", "## x\n" + HEADER + "\n"
                + emapper_row("g1", ["K00016"]) + emapper_row("g2", ["K00656"]))
    detail = write(tmp / "b.kofam.tsv",
                   "#\tgene name\tKO\tthrshld\tscore\tE-value\tKO definition\n"
                   "#\t---------\t------\t-------\t------\t---------\t-------------\n"
                   "*\tg1\tK00016\t100.0\t400.1\t1.2e-120\tL-lactate dehydrogenase\n"
                   "\tg2\tK00656\t500.0\t120.0\t3e-30\tformate C-acetyltransferase\n"
                   "*\tg3\tK00156\t300.0\t510.0\t2e-150\tpyruvate dehydrogenase (quinone)\n")
    mapper = write(tmp / "b.kofam.mapper", "g1\tK00016\ng3\tK00156\ng4\n")

    kof = dc.read_kofam(detail)
    assert kof == {"K00016": {"g1"}, "K00156": {"g3"}}  # g2 is below threshold
    assert dc.read_kofam(mapper) == kof
    assert "K00656" in dc.read_kofam(detail, significant_only=False)

    e = dc.read_emapper(egg)
    union, labels = dc.combine(e, kof, "union")
    assert set(union) == {"K00016", "K00656", "K00156"}
    assert labels == {"K00016": "both", "K00656": "eggNOG only", "K00156": "KOfam only"}
    assert set(dc.combine(e, kof, "consensus")[0]) == {"K00016"}
    assert set(dc.combine(e, kof, "eggnog")[0]) == {"K00016", "K00656"}
    assert set(dc.combine(e, kof, "kofam")[0]) == {"K00016", "K00156"}

    lactate = dc.function_evidence(*dc.combine(e, kof, "union"))["Mixed acid: Lactate"]
    assert lactate["found"][0]["source"] == "both"


def test_corrected_mode():
    adh = {"K00001"}
    assert dc.score_kos(adh)["Mixed acid: Ethanol, Acetylaldehyde to Ethanol"] == 0
    fixed = dc.score_kos(adh, corrected=True)
    assert fixed["Mixed acid: Ethanol, Acetylaldehyde to Ethanol"] == 1.0
    assert fixed["pectinesterase"] is None

    everything = set().union(*dc.SCORER_KOS.values())
    raw, cor = dc.score_kos(everything), dc.score_kos(everything, corrected=True)
    assert raw["nitrogen fixation"] == 0.99 and cor["nitrogen fixation"] == 1.0
    assert all(v == 1.0 for f, v in cor.items() if f not in dc.NOT_EVALUATED), \
        [f for f, v in cor.items() if v not in (1.0, None)]
    assert dc.same_kos_groups()["valine"] == ["isoleucine"]

    # Sec-SRP without SecDF (e.g. streptococci): original 0, corrected counts the other genes
    sec = {"K03070", "K03073", "K03075", "K03076", "K03106", "K03110", "K03210", "K03217"}
    assert dc.score_kos(sec)["Sec-SRP"] == 0
    ceiling = dc.load_ceilings().get("Sec-SRP") or 1.0
    assert dc.score_kos(sec, corrected=True)["Sec-SRP"] == round(min(1.0, 8 * 0.083 / ceiling), 2)
    # with secD present the original rule applies unchanged
    assert dc.score_kos(sec | {"K03072"}, corrected=True)["Sec-SRP"] == round(min(1.0, dc.score_kos(sec | {"K03072"})["Sec-SRP"] / ceiling), 2)


def test_cli(tmp):
    egg = write(tmp / "c.emapper.annotations", "## x\n" + HEADER + "\n"
                + emapper_row("g1", ["K00016"]) + emapper_row("g2", ["K00001"]))
    out, ev = tmp / "c.tsv", tmp / "c.json"
    dc.main(["-i", str(egg), "-n", "S_1", "-o", str(out), "-e", str(ev), "--corrected"])
    header, row = [l.split("\t") for l in out.read_text().splitlines()]
    values = dict(zip(header[1:], row[1:]))
    assert row[0] == "S_1"  # underscores in sample names are kept
    assert values["pectinesterase"] == "NA"
    assert values["Mixed acid: Ethanol, Acetylaldehyde to Ethanol"] == "1.0"
    data = json.loads(ev.read_text())
    assert data["corrected"] and "S_1" in data["samples"]
    assert data["samples"]["S_1"]["Mixed acid: Lactate"]["found"][0]["genes"] == ["g1"]


def main():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_function_map()
        test_evidence_and_copies(tmp)
        test_kofam_and_modes(tmp)
        test_corrected_mode()
        test_cli(tmp)
    print("OK: function map, evidence/copies, filters, KofamScan modes, corrected mode, CLI")


if __name__ == "__main__":
    main()
