# Findings

Status on 2026-11-08. Numbers come from `report.json` (real graph),
`simulation/*/report.json` and `dream5/*/report.json`. AUC compares genes with outgoing edges: 10,155
genes, of which 741 are Lambert 2018 TFs and 1,639 are in `tfs.txt`.

## The graph

- 34,033 genes, 10,100,294 signed edges; only 10,155 genes send edges, 98% of
  them protein-coding. Many target-only genes are non-coding, pseudogenes or
  Ensembl IDs that no longer exist, so the graph predates current Ensembl.
- No weight lies between −2.35 and 2.35: scores were cut at a threshold.
- Weights are directional: no two-way pair has equal weights.
- Positive edges between senders are often two-way (46%); negative edges
  rarely are (21%).

## Edges against curated TF→target pairs

- Curated pairs (CollecTRI, DoRothEA A–C, TRRUST) are present 1.5–2.3× more
  often than chance.
- Positive edges are 2.6–4.5× enriched for curated pairs; negative edges are
  below chance (0.7×). Edge signs agree with curated signs only at chance level.
- The 2016 cleaning keeps mostly negative edges and drops nearly all positive
  ones, leaving curated-pair recall at or below chance.
- Each gene's 50 strongest edges are the most enriched (5–7×).

## Ranking regulators (Lambert 2018 reference)

Paired bootstrap over genes, 1,000 resamples; best baseline is out-degree,
AUC 0.571.

| Pipeline | AUC | vs best baseline (95%) |
| --- | ---: | --- |
| repression-order | 0.631 | +0.059 (+0.042, +0.077) |
| acyclic order | 0.611 | +0.040 (+0.028, +0.052) |
| clear-direction pairs, acyclic order | 0.611 | +0.039 (+0.027, +0.051) |
| net flow | 0.596 | +0.025 (+0.012, +0.038) |
| 2016 cleaning, lowest PageRank | 0.549 | −0.022 (−0.052, +0.009) |
| activation-only acyclic order | 0.507 | −0.065 (−0.088, −0.042) |

The same order holds against `tfs.txt` with smaller gains. PageRank, either
after the 2016 cleaning or on the acyclic DAG, is above random but no better
than out-degree.

PageRank variants fixed in advance (†) do not help either (Lambert AUC, 95%
interval against out-degree): reverse PageRank on the raw graph 0.554 (worse,
−0.026 to −0.008); PageRank balance, log reverse minus log forward weighted
PageRank at alpha 0.85, 0.579 (−0.002 to +0.018); the same on repressing
edges 0.570, with 2 TFs in the top 100. Alpha 0.6 and 0.95 give 0.570 and
0.582. PageRank divides each gene's score among its edges, so every gene
passes on the same total however much weight it sends; that removes the
out-against-in weight imbalance that net flow and repression-order use. The
log ratio also favours genes with almost no incoming edges.

Combining repression-order with PageRank (†, fixed in advance) lowers the AUC:
reverse PageRank on the repression-order DAG 0.540 (5 TFs in the top 100), and
the rank sum of order position and that PageRank 0.602 against 0.631 alone.
The rank sum does put more TFs at the very top (Lambert 20 in the top 100
against 13; `tfs.txt` 29 against 17), a post-hoc observation to confirm
elsewhere. The regulator
signal is in the direction of negative edges: genes that send more negative
weight than they receive are more often TFs (post hoc, AUC 0.599), although
negative edges are not curated regulatory pairs.

## Pipelines that use known TFs (‡, `report.json`; breakdown in `diagnostics.json`)

These pipelines use known TFs as input (never the scored gene's own label), so
they answer "does a gene look like the TFs around it?", not "is it a
regulator?".

| Score (Lambert 2018) | AUC | Top 100 | Without KRAB genes: AUC | Top 100 |
| --- | ---: | ---: | ---: | ---: |
| random | 0.500 | 7 | 0.500 | 6 |
| out-degree (best baseline) | 0.571 | 14 | 0.553 | 13 |
| repression-order | 0.631 | 13 | 0.615 | 13 |
| tf-neighbour-share ‡ | 0.622 | 48 | 0.573 | 18 |
| tf-target-share ‡ | 0.649 | 45 | 0.604 | 20 |
| similar-targets ‡ | 0.701 | 57 | 0.662 | 33 |

- Comparing target sets ranks TFs best, +0.13 AUC over out-degree
  (95% interval +0.10 to +0.16), and +0.11 without the KRAB family.
- The top of that ranking is two co-expressed groups, not TFs as such:
  46 of the top 100 are TFs in the HGNC "Zinc fingers C2H2-type" group and
  44 lie on chromosome 19, where these genes sit in clusters; without the
  KRAB family, the top 100 still holds 18 such TFs, plus immediate-early
  response genes (EGR1, FOS, JUN, ATF3, NR4A1) together with non-TF partners
  (DUSP1, IER2, ZFP36, SOCS3). The method finds TF-rich co-expression modules.
- The immediate-early group is a common signature of tissue handling and
  stress in bulk expression data, a clue to the graph's source data.
- TFs send 7.3% of their outgoing weight to TFs, other genes 5.2%; TFs are
  4.5% of all targets. The graph is only mildly TF-to-TF.
- Against `tfs.txt` all these scores are weaker (similar targets 0.583,
  neighbour shares near chance).

## How the sending genes were chosen

Among protein-coding genes in the graph, Lambert TFs are less likely to send
edges than other genes: 7.4% of senders against 10.4% of target-only genes
(odds ratio 0.69, Fisher p = 1e-11). The senders were not picked as candidate
regulators. For `tfs.txt` the opposite holds (16.4% against 12.6%, odds ratio
1.36), so that list may have been built with the graph's senders in view.

## Simulation with known truth

Three seeds each of two inference methods (`simulate.py`), settings fixed in
advance:

- Edges are recovered well: 71–90% of true edges, 23–60× chance, signs 99.9%
  right, direction mostly right.
- No pipeline ranks the true regulators reliably above the degree baseline.
  Under correlation, all methods are near chance and top-layer regulators rank
  below chance. Under partial correlation, regulators send fewer edges than
  other senders (37 against 64 on average) and rank below chance.

## DREAM5 networks with known truth

DREAM5 network 1 (in silico, 1,643 genes, 195 TFs, 4,012 true edges, signed)
and network 3 (E. coli, 4,511 genes, 334 TFs, 2,066 RegulonDB edges), 805
samples each (`results/dream5/`). Settings fixed before running: the graph is
inferred as in the simulation, |z| > 2.35, with every gene a candidate sender,
so no TF knowledge enters it. References: the DREAM5 TF list, and the TFs with
at least one true edge (178 and 141).

- Edges are recovered at 9–14× chance (recall 11–40%); on network 1 edge signs
  agree with the true signs 97% of the time.
- No label-free pipeline ranks TFs above the best degree baseline in any of
  the four graphs (network × inference), against either reference. On E. coli
  with correlation, out-degree ranks TFs below chance (AUC 0.40); the ordering
  pipelines beat it there but stay below chance (0.42–0.46).
- Repression-order, the best pipeline on the real graph, does not carry over:
  DREAM5-TF AUC 0.471, 0.491, 0.458, 0.621 against best baselines 0.544,
  0.564, 0.404, 0.632.
- The 2016 cleaning cuts curated-edge recall 2.4–9×, and on E. coli
  enrichment falls from 10.6× to 4.0× (correlation) and 8.9× to 2.9× (partial).
  As on the real graph it keeps mostly negative edges: on network 1 with
  correlation, regulator edges go from 4,571 positive and 3,823 negative to 670
  and 2,708, and true activation-edge recall from 38% to 4.5% while repression
  stays at 30%.
- The ‡ pipelines that use known TFs score high (AUC 0.63–0.86), as on the
  real graph: TFs sit near other TFs.

## What this means

1. Recovering edges well does not make regulators identifiable from graph
   structure: in the simulation a gene's out-degree reflects how its
   expression relates to the rest of the data, not its role. DREAM5 networks 1
   and 3 confirm this on a standard benchmark and on E. coli data.
2. On the real graph, ordering methods beat degree significantly but modestly
   (AUC up to 0.63), and the signal sits in edges that are not real regulatory
   pairs. Without knowing how the graph was made, this cannot be read as
   biology.
3. The 2016 method (cleaning, then low PageRank) is no better than degree, and
   its cleaning removes the edges that match curated regulation; on DREAM5 it
   removes true edges too, activation edges most.
4. Using known TFs, genes with TF-like target sets rank TFs best (AUC 0.70),
   but what it finds is TF-rich co-expression modules (zinc-finger clusters,
   immediate-early genes). It cannot separate a new TF from a co-expressed
   non-TF, so it does not discover TFs.

## Open

- Why direction in negative edges marks TFs on the real graph but not in the
  simulation; the simulation's inference may not match the real one.
- The inference method, expression data and sender selection of the real
  graph.
- Pipelines marked post hoc need confirming on independent data; on DREAM5
  none of them beats degree.
