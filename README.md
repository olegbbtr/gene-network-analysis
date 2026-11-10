# Gene network analysis

Studies what the shape of a gene network inferred from expression data says
about its genes. Networks are built from public data: 8 TCGA tumour types and
7 matching GTEx healthy tissues (`public.py`), plus simulations and DREAM5
with known truth. Every gene is a candidate sender, so no TF knowledge enters
a network.

## Research questions and answers

| # | Question | Answer | Result |
| --- | --- | --- | --- |
| 1 | Can network shape tell which genes regulate others (known TFs)? | **Negative** | Across all 16 public networks (Lambert 2018 TFs), out-degree gives AUC 0.39–0.52 and repression-order 0.42–0.50 against a random 95% range of 0.48–0.52; the best of the 19 label-free pipelines per network reaches only 0.48–0.56. Simulations and DREAM5 with known truth agree. Pipelines that use known TF labels reach 0.69–0.75 only because TFs are co-expressed with TFs. |
| 2 | What does network shape reflect instead? | Descriptive | Correlation, copy number and cell state. In all 16 networks 99.9–100% of signed triangles are balanced (random signs: 36–58%), while curated regulation (CollecTRI, DoRothEA) is 68–69% (random: 64–67%). Positive tumour edges join genes on the same chromosome arm 2.0–7.3× more often than shuffled locations (same band 3.1–11.9×), healthy tissue only 1.2–1.3× (1.5–2.0×); removing copy number in LUAD cuts 5.9× to 2.0×. In LUAD, 70 of the top 100 out-degree genes are cell-cycle genes or send mostly to them. |
| 3 | Do per-gene network descriptors predict essential genes beyond expression level? | **Positive** | In all 15 single-tissue networks, adding network descriptors to the gene's expression raises held-out AUC by +0.009 to +0.108 (every 95% interval above 0); with neighbours as well, AUC reaches 0.83–0.89. With expression from all 15 datasets (AUC 0.889), all networks still add +0.018 (95% interval +0.010 to +0.024). |
| 4 | Does it hold in healthy tissue as well as tumours? | **Positive** | Yes in all 7 healthy and 8 tumour networks. In healthy networks the neighbours of a gene carry most of it: neighbours alone give AUC 0.68–0.86, above expression alone (0.70–0.78) in lung, thyroid, prostate and liver. |
| 5 | Is the gain more than degree (number of links)? | **Positive, small** | Shape beyond degree adds +0.005 to +0.049 AUC per network (all intervals above 0), but not once all networks and all expression are combined (+0.003, interval −0.002 to +0.007); the 7 healthy networks together still give +0.026. |
| 6 | Is the gain robust to the strongest tumour artefact, copy number? | **Positive** | Removing each gene's copy-number effect (LUAD) leaves it unchanged: +0.068 against +0.070. |
| 7 | Do descriptors predict cancer genes (OncoKB) beyond expression? | **Weak** | Held-out AUC 0.56–0.69 at best; the network adds −0.002 to +0.050, often not significant. |
| 8 | Do a gene's neighbours predict its family and pathways? | **Positive** | Held-out AUC per HGNC gene group (124) and Reactome pathway (1,051): median 0.54–0.69 on single networks, 0.74–0.75 with all 15 combined; 90–94% of sets above 0.6. Ribosomal, proteasome, respiratory-chain, HOX and hnRNP groups reach 0.96–1.00. |
| 9 | Can that suggest roles for uncharacterised genes? | **Candidates** | 92 C\*orf\* genes without a Reactome pathway get a best-scoring pathway (`results/public/uncharacterised.tsv`). On known genes the top guess is right 20% of the time overall and 55% in the top score decile. |
| 10 | Does personalised PageRank after graph cleaning help? | **Positive in tumours** | For essential genes, keeping only two-way positive links and running PageRank restarted from known essential genes beats the plain neighbour share by +0.03 to +0.09 AUC in 7 of 9 tumour networks. In healthy networks the plain share is best. The 2016 cleaning drops AUC to 0.36–0.67. |
| 11 | Do cancer genes change their closest partners from healthy to tumour tissue? | **Negative** | Each gene's 50 closest partners change almost entirely (median 98–100%), most likely because GTEx and TCGA were processed differently; this adds nothing for OncoKB genes. Essential genes keep their partners more than other genes (AUC of rewiring 0.38–0.45). |
| 12 | Does a tissue's network best predict that tissue's cancer cell-line dependencies (DepMap)? | **Inconclusive** | Only 10–15 lineage-selective dependencies per lineage, testable for breast, colon, thyroid and prostate; the matched network is significantly best only for healthy colon. |
| 13 | Does a tissue's network best predict that tissue's GWAS genes? | **Negative** | The matched tissue ranks first with a significant gain only for forced expiratory volume in healthy lung (+0.049). GWAS genes score AUC 0.35–0.64 on the combined networks, like the tissue-free control, body height (0.58). |
| 14 | Do hubs reveal what the samples contain? | **Descriptive** | Marker modules (immune, stress, cell cycle, interferon) are 9–91× denser than average in every network except interferon in healthy colon (2×); the top hubs follow sample make-up: cell cycle in LUAD (70 of 100), immune genes in LUSC and healthy thyroid (36), fat-cell genes in healthy breast (PPARG, LPL), acid-secreting cell genes in stomach (KCNE2, KCNQ1). |

Questions 8–13 come from `gene-network-analysis applications` and question 10 from
`cleaning` (`applications.py`, `gene_sets.py`; `results/public/applications.json`,
`cleaning.json`); question 14 from `structure`. Questions 3–7 come from
`gene-network-analysis predict` (`prediction.py`,
`results/public/prediction.json`); question 1 from `run`, question 2 from
`structure` (`results/public/`). Essential genes are the 684 core essential
genes of Hart et al. 2017 (CEGv2), 606 of them among the 12,000 genes. Models
are logistic regressions scored on held-out chromosomes, so neighbouring genes
never sit on both sides of a split; gains are paired bootstraps over genes.

The earlier work on a 2016 graph of unknown origin is kept in the sections
below; that graph is not distributed and is not used for conclusions.

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

## Public data

`public.DATASETS` lists the inputs: UCSC Xena TCGA `HiSeqV2.gz` for LUAD,
LUSC, BRCA, COAD, THCA, PRAD, STAD and LIHC (primary tumours only; LUAD also
GISTIC2 copy number), and GTEx v8 TPMs for lung, breast, colon (transverse),
thyroid, prostate, stomach and liver. Sources and checksums are in
`~/workspace/data/bioinf/references/README.md`. The gene lists are CEGv2 core
essential genes (<https://github.com/hart-lab/bagel>) and the OncoKB cancer
gene list.
`applications` also reads Reactome pathways (`ReactomePathways.gmt.zip`), DepMap
24Q2 CRISPR gene effects with models and common essentials (figshare 25880521)
and the GWAS Catalog associations (`gwas-catalog-associations_ontology-annotated-full.zip`).
GWAS traits and DepMap lineages are mapped to tissues in `applications.py`,
fixed before any network was scored.

```sh
R=~/workspace/data/bioinf/references
uv run gene-network-analysis public $R local/public      # 16 graphs, about 1 minute
uv run gene-network-analysis predict local/public \
  --list essential=$R/gene-lists/CEGv2.txt \
  --list cancer=$R/gene-lists/oncokb-cancerGeneList.txt \
  --output results/public/prediction.json             # 10–20 minutes
uv run gene-network-analysis applications local/public $R --output results/public   # about 30 minutes
uv run gene-network-analysis cleaning local/public $R --output results/public       # about 20 minutes
uv run gene-network-analysis structure --graph tcga-luad=local/public/tcga-luad/graph.txt \
  --signed-edges collectri=data/edges/collectri.tsv --output results/public/structure.json
uv run gene-network-analysis run --graph local/public/tcga-luad/graph.txt \
  --reference lambert-2018=data/lambert-2018-tfs.txt --output results/public/ranking/tcga-luad
```

`uv run gene-network-analysis figures` draws the poster figures from these
results into `results/public/figures/`.

`public` keeps the 12,000 protein-coding genes that vary in every dataset and
have the highest mean expression quantile across them, infers each graph as
`simulate` does (z-scored correlations, |z| > 2.35) and writes
`expression.tsv` with each gene's mean and SD per dataset. For LUAD it also
writes a graph after removing each gene's linear dependence on its own copy
number. `results/public/structure.json` and `results/public/ranking/` were made
with every graph in `local/public/`.

`structure` (`structure.py`) measures, per graph: the share of edges joining
genes on the same chromosome arm or band against shuffled gene locations;
marker modules fixed in advance (sex, immune, immediate-early, cell cycle,
interferon) and how many of the top 100 out-degree genes belong to or send
mostly into one; and the share of balanced signed triangles against random
signs, also for curated and true networks.

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
