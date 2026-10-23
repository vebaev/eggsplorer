# Reference genomes for the validation test

KO assignments (`query`, `KEGG_ko`) of ten bacterial genomes whose cell envelope, respiratory chain,
secretion systems and nitrogen, sulfur and C1 metabolism are well described. They are used by
`tests/test_validation.py` to check the cell overview and element-cycle rules.

| Assembly | Organism |
|---|---|
| GCF_000005845.2 | *Escherichia coli* K-12 MG1655 |
| GCF_000009045.1 | *Bacillus subtilis* 168 |
| GCF_000006765.1 | *Pseudomonas aeruginosa* PAO1 |
| GCF_000021045.1 | *Azotobacter vinelandii* DJ |
| GCF_000009145.1 | *Nitrosomonas europaea* ATCC 19718 |
| GCF_000012725.1 | *Nitrobacter winogradskyi* Nb-255 |
| GCF_000195755.1 | *Desulfovibrio vulgaris* Hildenborough |
| GCF_000025485.1 | *Allochromatium vinosum* DSM 180 |
| GCF_000009725.1 | *Synechocystis* sp. PCC 6803 |
| GCF_000008325.1 | *Methylococcus capsulatus* Bath |

Provenance: the RefSeq protein FASTA of each assembly was annotated with eggNOG-mapper 2 in Galaxy
(October 2026) and the KO column was extracted from the Eggsplorer report. Only proteins with at least
one KO are kept; the RefSeq protein accessions identify the genes. All genomes are public.
