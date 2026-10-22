"""Render the Galaxy <command> with Cheetah for the XML <tests> and check outputs.

A stand-in for `planemo test` where a Galaxy server cannot be started. It
fills the template the way Galaxy would for each test case, runs the command
and evaluates the test's assert_contents / collection checks.

Run with a Python that has Cheetah3: python render_galaxy_tests.py
"""

import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from Cheetah.Template import Template

ROOT = Path(__file__).resolve().parents[1]


class Data(str):
    """A dataset parameter: renders as its path, has element_identifier."""
    def __new__(cls, path, ident):
        obj = str.__new__(cls, path)
        obj.element_identifier = ident
        return obj


class NS(dict):
    __getattr__ = dict.__getitem__


def case_params(test, tmp):
    td = ROOT / "test-data"
    flat = {}
    for p in test.iter("param"):
        flat[p.get("name")] = p.get("value")
    ann = flat["annotations"]
    params = {
        "annotations": Data(str(td / ann), ann),
        "sample_name": flat.get("sample_name", ""),
        "kofam_cond": NS(use_kofam=flat.get("use_kofam", "no"),
                         kofam=Data(str(td / flat["kofam"]), flat["kofam"]) if "kofam" in flat else None,
                         mode=flat.get("mode", "union")),
        "filters": NS(max_evalue=flat.get("max_evalue", ""), min_score=flat.get("min_score", "")),
        "top_n": flat.get("top_n", "25"),
        "lang": flat.get("lang", "en"),
        "output_tables": flat.get("output_tables", "false") == "true",
        "report": str(tmp / "report.html"),
        "__tool_directory__": str(ROOT),
    }
    return params


def check_contents(text, asserts, label):
    for a in asserts:
        tag = a.tag
        if tag == "has_text":
            assert a.get("text") in text, f"{label}: missing {a.get('text')!r}"
        elif tag == "not_has_text":
            assert a.get("text") not in text, f"{label}: unexpected {a.get('text')!r}"
        elif tag == "has_n_lines":
            n = len(text.splitlines())
            assert n == int(a.get("n")), f"{label}: {n} lines, expected {a.get('n')}"
        elif tag == "has_line_matching":
            assert any(re.fullmatch(a.get("expression"), l) for l in text.splitlines()), \
                f"{label}: no line matches {a.get('expression')!r}"
        else:
            raise NotImplementedError(tag)


def main():
    tool = ET.parse(ROOT / "eggsplorer.xml").getroot()
    command = tool.find("command").text
    for i, test in enumerate(tool.find("tests").findall("test"), 1):
        tmp = Path(tempfile.mkdtemp())
        try:
            params = case_params(test, tmp)
            rendered = str(Template(command, searchList=[params]))
            cmd = " ".join(rendered.split())
            subprocess.run(cmd, shell=True, check=True, cwd=tmp)
            n_out = 1 + (1 if params["output_tables"] else 0)
            assert n_out == int(test.get("expect_num_outputs")), f"test {i}: outputs"
            for out in test.findall("output"):
                check_contents((tmp / "report.html").read_text(encoding="utf-8"),
                               out.find("assert_contents"), f"test {i} {out.get('name')}")
            for coll in test.findall("output_collection"):
                files = {p.stem: p for p in (tmp / "tables").glob("*.tsv")}
                assert len(files) == int(coll.get("count")), f"test {i}: {len(files)} tables"
                names = [el.get("name") for el in coll.findall("element")]
                ordered = sorted(files)  # Galaxy checks that listed elements follow the collection order
                assert [n for n in ordered if n in names] == names, f"test {i}: element order {names}"
                for el in coll.findall("element"):
                    check_contents(files[el.get("name")].read_text(encoding="utf-8"),
                                   el.find("assert_contents"), f"test {i} {el.get('name')}")
            print(f"test {i}: OK  ({cmd[:110]}...)")
        finally:
            shutil.rmtree(tmp)


if __name__ == "__main__":
    sys.exit(main())
