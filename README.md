# Gene network analysis

Ranks the genes of a signed, weighted gene regulatory network by how likely
they are to be regulators, and checks each ranking against a list of known
transcription factors (TFs).

## Run

```sh
uv sync --group dev
make test
make lint
make reproduce        # checks data/SHA256SUMS, then runs every reported pipeline
uv run gene-network-analysis run --pipeline net-flow --output /tmp/one   # one pipeline
```

A run writes `results/report.json`, `results/report.pdf` and
`results/rankings.tsv` (every pipeline's score per gene, not tracked). The full
graph takes about 40 seconds.

## Data

See [data/README.md](data/README.md): `graph.txt` is the network;
`tfs.txt` and `lambert-2018-tfs.txt` are TF references. The report evaluates
every pipeline against each reference, each with its own random row and best
baseline; `--reference NAME=PATH` replaces them. `data/edges/` holds curated
TF→target pairs, rebuilt from their sources (see `data/README.md`); the report
checks every graph and DAG the pipelines build against them (`--edges
NAME=PATH` replaces them).

## Design

Steps are reusable functions of the graph, grouped by type:

- **Preprocessing** (`steps.py`, graph → graph): `resolve_reciprocals`,
  `drop_sinks`, `prune_triangles`, `orient_reciprocals`, `top_outgoing`,
  `keep_sign`, `absolute`.
- **De-cycling** (`steps.py`, graph → DAG and a level per gene, lower is
  further upstream): `greedy_layering`, `minimum_feedback_order`
  (Eades–Lin–Smyth ordering that keeps edges pointing forward).
- **Enrichment** (`pipelines.FEATURES`): per-gene scores computed on demand,
  such as degrees, net flow, PageRank and DAG level.
- **Selection** (`report.analyze`): every pipeline ranks the same genes, those
  with at least one outgoing edge in the input graph.
- **Metrics** (`evaluation.py`): TFs in the top k, ROC AUC, and p-values against
  2,000 random orderings.
- **Edge checks** (`edges.py`): for each graph and DAG, the share of curated
  TF→target pairs present as edges, how many positive and negative edges are
  curated pairs, and agreement with the curated sign, each against chance.

A pipeline (`pipelines.PIPELINES`) is one configured chain of steps that ends in
one score, so each pipeline gives one ranking. `baseline=True` marks the
pipelines the others are compared with; `reported=False` keeps a pipeline in
code but out of the default report. Pipelines that share leading steps reuse
their computed graphs.

| Pipeline | Steps |
| --- | --- |
| `out-degree`, `weighted-out-degree` | Baselines on the input graph. |
| `2016-pagerank`, `2016-layer` | Reciprocal subtraction, sink removal and same-sign triangle pruning, then lowest unweighted PageRank (damping 0.6) or earliest greedy layer. |
| `acyclic-order`, `acyclic-dag-pagerank` | Absolute weights, minimum-feedback ordering, then position or reverse PageRank on the DAG. |
| `net-flow` | Outgoing minus incoming absolute weight. |
| `oriented-order`, `top50-order`, `activation-order`, `repression-order` | One extra preprocessing step before the acyclic ordering. |

The report lists a random row (mean and 95% range of random orderings), every
pipeline's metrics, and its AUC difference from the best baseline.

## Pipelines that use known TFs

Marked ‡ in the report. Each scores a gene from the TF labels of other genes in
the reference being evaluated, never its own label, so they are computed once
per reference and are not label-free predictions (`neighbours.py`):

| Pipeline | Score |
| --- | --- |
| `tf-neighbour-share`, `tf-target-share`, `tf-sender-share` | Share of a gene's absolute edge weight (both directions, outgoing, incoming) that connects it to TFs. |
| `similar-targets` | Similarity-weighted TF share among the 20 genes whose 50 strongest targets are most alike by cosine. `similar-signed-targets` (not reported) compares target signs too. |

## Diagnostics

`gene-network-analysis diagnostics` (`diagnostics.py`) writes
`results/diagnostics.json`. For each reference it records:

- whether TFs are over- or under-represented among the protein-coding genes
  that send edges;
- the share of outgoing weight that TFs and other genes send to TFs;
- the top 100 genes of out-degree, repression-order and every reported ‡
  pipeline, with HGNC symbol, location and gene groups, and how many are in
  the "Zinc fingers C2H2-type" group;
- the same comparison and top genes without the HGNC "KRAB domain
  containing" genes, a TF family expressed together that dominates the ‡
  scores.

## Simulation

`gene-network-analysis simulate DIR --seed N --inference correlation|partial`
writes a network with known truth (`simulate.py`): 2,000 genes, 200 regulators
in 4 layers acting on later regulators and other genes, about 20 targets each,
30% repression, linear expression with noise, 500 samples. The inferred graph
z-scores each sender's correlation (or ridge partial correlation) with a target
against all senders for that target and keeps |z| > 2.35; senders are the
regulators plus 30% of other genes. `DIR` then holds `graph.txt`,
`regulators.txt`, `top-regulators.txt` and `true-edges.tsv` for `run`:

```sh
d=local/simulation/partial-seed-0
uv run gene-network-analysis simulate $d --seed 0 --inference partial
uv run gene-network-analysis run --graph $d/graph.txt \
  --reference regulators=$d/regulators.txt --reference top-layer=$d/top-regulators.txt \
  --edges true=$d/true-edges.tsv --output results/simulation/partial-seed-0
```

`results/simulation/` holds reports for seeds 0–2 of both inference methods.

## DREAM5

`gene-network-analysis dream5 N EXPRESSION TFS GOLD DIR [--signed SIGNED]
[--inference correlation|partial]` (`dream5.py`) infers a graph from a DREAM5
network inference dataset as `simulate` does, with every gene a candidate
sender, and writes `graph.txt`, `tfs.txt` (the DREAM5 TF list),
`regulators.txt` (TFs with a gold-standard edge) and `true-edges.tsv`. Genes
are renamed `ENSGD5N<N><gene>` because the report ranks only Ensembl-like IDs.
The input files come from the DREAM5 challenge (Marbach et al., Nature Methods
9:796–804, 2012; Synapse `syn2787209`). The reports used copies whose MD5s
match the Synapse manifest, from
<https://github.com/aertslab/arboreto/tree/master/resources/dream5>
(expression, TF lists) and
<https://github.com/DuaaAlawad/AGRN/tree/main/Dataset/DREAM5/test%20data>
(gold standards):

```sh
R=path/to/dream5/net1 d=local/dream5/net1-correlation
uv run gene-network-analysis dream5 1 $R/net1_expression_data.tsv \
  $R/net1_transcription_factors.tsv "$R/DREAM5_NetworkInference_GoldStandard_Network1 - in silico.tsv" $d \
  --signed "$R/DREAM5_NetworkInference_GoldStandardSigned_Network1 - in silico.tsv"
uv run gene-network-analysis run --graph $d/graph.txt --reference dream5-tfs=$d/tfs.txt \
  --reference true-regulators=$d/regulators.txt --edges gold=$d/true-edges.tsv \
  --output results/dream5/net1-correlation
```

`results/dream5/` holds reports for networks 1 (in silico) and 3 (E. coli)
with both inference methods.

## Limits

- The TF list marks known TFs, not master regulators; genes outside it are not
  confirmed non-regulators.
- How the network was generated is not known.
- `prune_triangles` decides every triangle on the same input snapshot, so the
  result does not depend on visiting order.
