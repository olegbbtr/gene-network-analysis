# Data

| File | Tracked | Contents |
| --- | --- | --- |
| `graph.txt` | no (469 MB) | Directed gene network: `source<TAB>target<TAB>weight`, human Ensembl gene IDs. 34,033 genes and 10,100,294 edges; 5,562,481 have negative weight. Weights range from about −19 to 42. Only 10,155 genes have outgoing edges. |
| `tfs.txt` | yes | Reference list of 3,194 TF Ensembl gene IDs, one per line after an `Ensembl Gene ID` header. 1,639 of them have outgoing edges in the graph. Its source and selection criteria are not known. |
| `lambert-2018-tfs.txt` | yes | The 1,639 human TFs of Lambert et al., "The Human Transcription Factors", Cell 172(4):650–665 (2018), doi:10.1016/j.cell.2018.01.029: proteins with experimentally shown sequence-specific DNA binding, judged by an expert panel. Unmodified `TFs_Ensembl_v_1.01.txt` from <https://humantfs.ccbr.utoronto.ca/download/v_1.01/>, downloaded 2026-11-07; it equals the `Is TF? = Yes` rows of that release's database. Two entries (`DUX1_HUMAN`, `DUX3_HUMAN`) have no Ensembl ID. 741 have outgoing edges in the graph. |
| `hgnc-genes.tsv` | yes | The 42,640 HGNC genes with an Ensembl ID: `ensembl_gene_id`, `symbol`, `location`, `locus_group`, `gene_group` (`|`-separated). Used by `diagnostics` for protein-coding status, gene families and symbols. |
| `SHA256SUMS` | yes | Checksums of the data files, including the untracked `graph.txt` and `edges/`; `make check-data` verifies them. |

`hgnc-genes.tsv` comes from the HGNC complete set
(<https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt>)
with `gene-network-analysis import-genes HGNC OUTPUT`.

The two TF lists share 1,119 genes. `tfs.txt` has 2,075 genes that Lambert
2018 does not count as TFs; among the genes with outgoing edges, 497 are in
both lists.

`graph.txt` is not distributed with this repository; `SHA256SUMS` identifies
the copy the reports were made from.

Neither the graph's generating method nor the TF list's selection criteria are
known. A positive weight is read as activation and a negative one as
repression.

## Curated TF→target pairs (`edges/`)

Tab-separated `source_gene`, `target_gene`, `sign` (1 activation, −1
repression, 0 unknown or conflicting), in Ensembl gene IDs. The files are not in Git because the sources carry their own licences; rebuild
them from the downloads below with `gene-network-analysis import-edges`, which
maps them through the HGNC complete set (UniProt IDs for OmniPath, symbols for
TRRUST), then check them with `make check-data`:

```sh
uv run gene-network-analysis import-edges collectri collectri.tsv --hgnc hgnc_complete_set.txt data/edges/collectri.tsv
uv run gene-network-analysis import-edges dorothea-abc dorothea.tsv --hgnc hgnc_complete_set.txt data/edges/dorothea-abc.tsv
uv run gene-network-analysis import-edges trrust trrust_rawdata.human.tsv --hgnc hgnc_complete_set.txt data/edges/trrust.tsv
```

Downloads (OmniPath with `license=academic`, for research use):

- CollecTRI: <https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&fields=sources,references&license=academic>
- DoRothEA: <https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B,C,D,E&genesymbols=yes&fields=dorothea_level&license=academic>
- TRRUST v2: <https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv>

| File | Source | Pairs |
| --- | --- | ---: |
| `collectri.tsv` | CollecTRI via OmniPath (`datasets=collectri`), complexes split into member genes | 47,411 |
| `dorothea-abc.tsv` | DoRothEA via OmniPath, confidence levels A–C | 32,247 |
| `trrust.tsv` | TRRUST v2 human | 8,372 |

## DREAM5 network inference data (not in Git)

Used by `gene-network-analysis dream5` (see the main README) for networks with
a known gold standard: network 1 (in silico, 1,643 genes, 195 TFs, 4,012 true
edges, also signed) and network 3 (E. coli, 4,511 genes, 334 TFs, 2,066
RegulonDB edges), 805 expression samples each. Source: the DREAM5 challenge,
Marbach et al., "Wisdom of crowds for robust gene network inference", Nature
Methods 9:796–804 (2012), doi:10.1038/nmeth.2016; Synapse `syn2787209`, which
needs a free account and has its own terms of use, so the files are not
redistributed here.

The reports used copies fetched on 2026-11-08 from public mirrors, whose MD5s
match the Synapse manifest:

- expression and TF lists:
  <https://github.com/aertslab/arboreto/tree/master/resources/dream5>
- gold standards:
  <https://github.com/DuaaAlawad/AGRN/tree/main/Dataset/DREAM5/test%20data>
  (network 1's unsigned file is named `… - Copy.tsv` there)

```text
1a5c61fa117e002e65acc1c6c26a5599  net1/net1_expression_data.tsv
d4cfec1596a7dde09ce324d6476e2dbe  net1/net1_transcription_factors.tsv
a9f2da4f137a5827fedec5605acd401f  net1/DREAM5_NetworkInference_GoldStandard_Network1 - in silico.tsv
9bf6de800feb05c67fc04c5e30bac798  net1/DREAM5_NetworkInference_GoldStandardSigned_Network1 - in silico.tsv
d18b463871dd8e63532647a4cfc497fd  net3/net3_expression_data.tsv
489672e61c445c7d5d96841b13b36567  net3/net3_transcription_factors.tsv
5f7c7fbda344c38ca8b19057255bee17  net3/DREAM5_NetworkInference_GoldStandard_Network3 - E. coli.tsv
```

## Simulation

The simulated networks are not stored: `gene-network-analysis simulate` makes
them from a fixed seed (see the main README).
