# Findings

## Public data (2026-11-10)

The research questions and short answers are in the README. Networks: 8 TCGA
tumour types and 7 GTEx healthy tissues on one set of 12,000 genes
(`results/public/`). The 2016 graph below is of unknown origin and is not
used for conclusions.

### Positive: network descriptors predict essential genes beyond expression

606 CEGv2 core essential genes among 12,000. Held-out AUC (whole chromosomes
held out), logistic regression; "network gain" is expression + shape against
expression alone, paired bootstrap over genes.

| Networks | Expression | + degree | + shape | + neighbours | Network gain (95%) |
| --- | ---: | ---: | ---: | ---: | --- |
| tcga-luad | 0.775 | 0.822 | 0.844 | 0.870 | +0.070 (+0.055, +0.086) |
| tcga-luad-copy-number-removed | 0.775 | 0.826 | 0.842 | 0.878 | +0.068 (+0.054, +0.083) |
| tcga-lusc | 0.810 | 0.839 | 0.852 | 0.877 | +0.042 (+0.032, +0.051) |
| tcga-brca | 0.785 | 0.832 | 0.838 | 0.869 | +0.053 (+0.041, +0.065) |
| tcga-coad | 0.809 | 0.857 | 0.867 | 0.879 | +0.058 (+0.047, +0.069) |
| tcga-thca | 0.767 | 0.793 | 0.802 | 0.827 | +0.036 (+0.025, +0.046) |
| tcga-prad | 0.765 | 0.793 | 0.809 | 0.838 | +0.045 (+0.033, +0.057) |
| tcga-stad | 0.786 | 0.832 | 0.849 | 0.890 | +0.063 (+0.051, +0.074) |
| tcga-lihc | 0.796 | 0.833 | 0.839 | 0.853 | +0.043 (+0.032, +0.054) |
| gtex-lung | 0.716 | 0.733 | 0.761 | 0.861 | +0.045 (+0.031, +0.060) |
| gtex-breast | 0.780 | 0.779 | 0.794 | 0.826 | +0.015 (+0.008, +0.022) |
| gtex-colon | 0.776 | 0.776 | 0.786 | 0.845 | +0.009 (+0.002, +0.017) |
| gtex-thyroid | 0.775 | 0.786 | 0.794 | 0.848 | +0.019 (+0.010, +0.030) |
| gtex-prostate | 0.768 | 0.780 | 0.786 | 0.859 | +0.017 (+0.009, +0.026) |
| gtex-stomach | 0.758 | 0.757 | 0.793 | 0.849 | +0.035 (+0.022, +0.048) |
| gtex-liver | 0.705 | 0.763 | 0.813 | 0.868 | +0.108 (+0.087, +0.130) |
| all-tumour | 0.877 | 0.896 | 0.896 | 0.902 | +0.018 (+0.011, +0.025) |
| all-healthy | 0.822 | 0.844 | 0.870 | 0.892 | +0.047 (+0.036, +0.059) |
| all | 0.889 | 0.904 | 0.907 | 0.911 | +0.018 (+0.010, +0.024) |

- Every single network adds to expression, healthy and tumour alike.
- In healthy tissue the strongest feature is the share of a gene's neighbours
  that are essential (training labels only): essential genes are
  co-expressed with each other.
- Removing copy number (LUAD) does not change the gain.
- Shape beyond degree is small and disappears when all networks and all
  expression are combined.
- Cancer genes (OncoKB, 1,010 listed) are predicted weakly: best AUC 0.69
  with everything combined, against 0.675 from expression alone.

### Positive: neighbours predict gene function

`applications.json`, `function`: share of a gene's neighbours in each HGNC gene
group (124 sets of 20–500 genes) or Reactome pathway (1,051), members on
held-out chromosomes hidden. Single networks: median AUC 0.54–0.69 (healthy
stomach and colon lowest); mean over the 15 networks: 0.744 (HGNC) and 0.748
(Reactome), with 90% and 94% of sets above 0.6. Best: hnRNPs, respiratory
chain, ribosomal proteins, HOX genes, proteasome (0.96–1.00); worst near or
below 0.5: kinase-defined groups and signalling pathways.

`uncharacterised.tsv` lists the 92 C*orf* genes without a Reactome pathway with
their best-scoring pathway. On known genes the top guess is right 19.5% of the
time; precision by score decile rises from 6% to 55%, and each candidate
carries the precision of its decile. Candidates, not findings.

### Positive in tumours: personalised PageRank after cleaning

`cleaning.json`, essential genes, held-out AUC per network and cleaning. Two-way
positive links with PageRank restarted from training-fold essential genes beat
the plain neighbour share by +0.03 to +0.09 in LUAD, LUSC, BRCA, COAD, LIHC,
STAD and LUAD without copy number (95% intervals above 0), not in THCA or PRAD.
In healthy networks the plain neighbour share is best except breast, colon and
stomach, where PageRank on each gene's 50 strongest links helps. The 2016
cleaning gives 0.36–0.67, at or below chance for the neighbour share.

### Negative or inconclusive uses

- Rewiring (`rewiring`): each gene's 50 closest partners change by a median of
  98–100% between GTEx and TCGA networks of one tissue, too much to separate
  biology from processing differences. Adding it to expression and degree does
  not help predict OncoKB genes (gains −0.002 to +0.003). Essential genes change
  partners less (AUC of rewiring 0.38–0.45; gain +0.007 to +0.021).
- DepMap (`dependencies`): lineage-selective dependencies (mean effect ≤ −0.5,
  0.25 below other lineages, not common essential) number 10–15 for breast,
  colon, thyroid and prostate and fewer than 10 for lung, stomach and liver.
  The matched network is significantly best only for healthy colon.
- GWAS (`gwas`): the matched tissue ranks first with a significant gain only
  for forced expiratory volume in healthy lung (+0.049); most matched gains are
  near zero or negative. Combined-network AUCs (0.35–0.64) resemble the
  tissue-free control, body height (0.58).

### Sample make-up

From `structure.json`: the five marker modules are 9–91× denser than average
except interferon in healthy colon (2×). Top-100 out-degree hubs: LUAD 70
cell-cycle; LUSC and healthy thyroid 36 immune; healthy breast led by
adipocyte genes (PPP1R1A, PPARG, LPL, TIMP4); healthy stomach by parietal-cell
genes (KCNE2, KCNQ1); healthy lung by stress and inflammation genes. The cell
types are a reading of the gene lists, not a tested deconvolution.

### Negative: network shape does not find regulators

Across the 16 networks no label-free pipeline ranks Lambert TFs reliably above
random (`results/public/ranking/`); out-degree and repression-order are mostly
below it. See README question 1.

### What network shape reflects

`results/public/structure.json`: signed triangles 99.9–100% balanced in every
network (correlation structure); tumour edges concentrated on chromosome arms
and bands, explained largely by copy number; LUAD hubs dominated by the cell
cycle. See README question 2.

## 2016 graph (status 2026-11-08)

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
