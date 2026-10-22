# Eggsplorer

Galaxy tool that turns an eggNOG-mapper `*.emapper.annotations` table into one
self-contained, interactive HTML report. It runs offline: plotly.js is embedded
and no database is queried.

## What the report shows

Overview and annotation coverage · taxonomic level of annotation transfer · COG categories (radar) ·
most frequent GO terms · **cell overview** (respiratory chain, uptake, secretion, transporter families,
motility and metabolic function groups on one fixed template) · metabolic function completeness
(KEGG-Decoder definitions, corrected, with per-function evidence) · nitrogen, sulfur and carbon cycle
diagrams · TC transporters (Sankey) · EC classes with ENZYME names (Sankey) · CAZy families and
CAZy by substrate · Pfam domains · searchable gene table. English or Bulgarian.

## Use in Galaxy

Install from the Galaxy Tool Shed (owner `vebaev`, repository `eggsplorer`) or add the folder to a
local `tool_conf.xml`. Galaxy removes JavaScript from HTML outputs by default, so an administrator
should add `eggsplorer` to the HTML allowlist (`sanitize_allowlist_file`, or *Admin → Manage
Allowlist*) for the report to work inside Galaxy. A downloaded report always works in a browser.

## Layout

| File | Role |
|---|---|
| `eggsplorer.xml`, `macros.xml` | Galaxy wrapper |
| `eggsplorer.py` | parses the annotations, computes summaries, writes HTML and TSV tables |
| `report_template.html` | page layout, styles and charts (Plotly) |
| `refdata.py` | static reference tables: COG categories, EC classes, TC subclasses and families, CAZy classes and families, layout groups |
| `refdata_checks/` | review of every bundled CAZy and TC name against its CAZy or TCDB page (October 2026): evidence, verdict, final label |
| `refdata_cazy/` | CAZy family → substrate groups derived from dbCAN-sub by `scripts/build_cazy_substrates.py`, with provenance |
| `go_core.py` | most frequent informative GO terms (mapping, propagation, filters) |
| `refdata_go/` | Gene Ontology table (CC BY 4.0) derived by `scripts/build_go_table.py`, with provenance |
| `decoder_core.py` | metabolic function scoring with the KEGG-Decoder logic, evidence, KofamScan support |
| `third_party/kegg_decoder/` | unmodified KEGG-Decoder v1.3 (MIT), its licence, provenance notes and the derived KO map |
| `test-data/` | synthetic 300-protein test genome and KofamScan file, written by `scripts/make_test_data.py` (not a real genome) |
| `tests/` | equivalence with KEGG-Decoder, checks of the additions, cross-check of report numbers with pandas, layout checks of the SVG figures |
| `.github/workflows/ci.yml` | runs the tests, `planemo lint` and `planemo test` on every push |

## Command line

```bash
python eggsplorer.py -i sample.emapper.annotations -n Sample1 -o report.html \
    [--kofam sample.kofam.tsv --mode consensus] [--max-evalue 1e-10] [--min-score 50] \
    [--tables tables/] [--lang en|bg] [--plotly inline|cdn]
```

Only the Python standard library is needed, plus the `plotly` package as the
source of `plotly.min.js` in the default inline mode.

## Tests

Test dependencies: `pip install plotly==5.24.1 pandas Cheetah3 playwright seaborn matplotlib scipy` and
`playwright install chromium` (seaborn, matplotlib and scipy are only needed to run the original KEGG-Decoder in
the equivalence test).

```bash
python -I tests/test_decoder_core.py   # identical scores to the original KEGG-Decoder
python -I tests/test_improvements.py   # KO map, evidence, filters, KofamScan, corrected mode
python -I tests/test_report.py         # report numbers vs an independent pandas computation
python -I tests/test_refdata.py        # bundled CAZy/TC names and CAZy substrate groups are consistent
python -I tests/test_go.py             # GO counts verified top-down, independent of go_core
python -I tests/test_cell.py           # cell overview rules; layout check in Chromium (test genome and worst case)
planemo test eggsplorer.xml
```

## Licences and citations

Tool code: MIT (see `LICENSE`). KEGG-Decoder: MIT, Copyright (c) 2019 Benjamin Tully.
Gene Ontology and ENZYME data: CC BY 4.0 (see the `NOTICE.md` files).
Cite eggNOG-mapper (Cantalapiedra et al. 2021), KEGG-Decoder (Graham et al. 2018)
and KEGG (Kanehisa & Goto 2000) when using the report; the CAZy-by-substrate chart uses dbCAN-sub (Zheng et al. 2023). KEGG names are not
bundled; KO identifiers link to the KEGG website.
