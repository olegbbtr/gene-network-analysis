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
