# Findings

## Public data (2026-11-10)

The research questions and short answers are in the README. Networks: 8 TCGA
tumour types and 7 GTEx healthy tissues on one set of 12,000 genes
(`results/public/`). The 2016 graph below is of unknown origin and is not
used for conclusions.

### Harder targets for the same label-free features (2026-11-12)

`targets` (`results/public/targets.json`) refits the label-free groups, with the same nested
penalty selection and leave-one-out, on three targets fixed before the first run (definitions
in `src/gene_network_analysis/targets.py`).

| Target | Genes (positives) | AUC gene / local / global / biological / all four | Top-100 hits, all four (gene only; chance) |
| --- | --- | --- | --- |
| Essential genes (CEGv2) | 12,000 (606) | 0.890 / 0.879 / 0.888 / 0.802 / 0.914 | see above |
| Selective dependency (DepMap, >= 5% of lines, outside both lists) | 10,077 (970) | 0.804 / 0.776 / 0.781 / 0.718 / 0.820 | 65 (45; 9.6) |
| Knockdown impact (K562, non-essential, >= 50% knockdown, >= 289 genes changed) | 5,832 (178) | 0.653 / 0.626 / 0.636 / 0.634 / 0.688 | 12 (4; 3.1) |
| First focused experimental GO function after 2018 | 2,684 (406) | 0.619 / 0.610 / 0.606 / 0.559 / 0.634 | 31 (25; 15.1) |

Leave-one-out AUC loss (95% gene-bootstrap interval above 0 marked *):

| Removed | Essential | Selective | Knockdown | First function |
| --- | --- | --- | --- | --- |
| Gene | 0.007* | 0.008* | 0.009* | 0.009* |
| Graph (local + global) | 0.008* | 0.002 | -0.005 | 0.010 |
| Biological | 0.010* | 0.014* | 0.035* | 0.002 |

Graph features alone match gene features on every target, but add significantly to the
others only for essential genes. Biological features (pathway and chromosome-arm sharing)
carry the extra signal for selective dependency and knockdown impact.

### Label-free feature groups and leave-one-out (2026-11-12)

`label-free-groups` (`results/public/label-free-groups.json`) fits each group alone, all four,
and all four without gene, graph (local + global together) or biological features, with the
nested penalty selection of `feature-audit`. No feature is computed from essential-gene labels.

- Groups: gene (expression mean and SD); local graph (degree, signs, weights, flow,
  clustering, plus out- and in-link counts above the 25th, 50th and 75th percentile of each
  network's own |z|, recomputed without a cutoff by `cutoff-counts`: about 0.33, 0.69 and
  1.16 in every network); global graph (PageRank scores and positions, k-core depth,
  eigenvector centrality); biological (Reactome pathway sharing of neighbours and of partners
  linked in all 15 networks; chromosome-arm share of links).

  | | Gene | Local | Global | Biological | All four |
  | --- | --- | --- | --- | --- | --- |
  | CEGv2 AUC | 0.890 | 0.879 | 0.888 | 0.802 | 0.914 |
  | DepMap AUC outside CEGv2 | 0.840 | 0.839 | 0.844 | 0.749 | 0.871 |
  | DepMap essentials in the top 100 | 59 | 64 | 62 | 72 | 73 |

- Leave-one-out AUC loss with gene-bootstrap 95% intervals:

  | Removed | CEGv2 | DepMap outside CEGv2 |
  | --- | --- | --- |
  | Gene | 0.007 (0.005, 0.009) | 0.005 (0.002, 0.008) |
  | Graph (local + global) | 0.008 (0.001, 0.015) | 0.012 (0.006, 0.018) |
  | Biological | 0.010 (0.006, 0.014) | 0.007 (0.003, 0.011) |

  Removed alone, local or global costs little, because each carries much the same
  information as the other; they are therefore removed together.
- The percentile link counts raise local-graph-only top-100 hits from 50 to 61 (CEGv2 AUC
  0.874 to 0.878 before the biological regrouping).
- How hard the task is: mean expression alone gives CEGv2 AUC 0.731. Without the 904 genes
  of the large housekeeping complexes all four give 0.904; within deciles of mean expression
  the median AUC is 0.884. Average precision is 0.377 for all four, 0.284 for gene only
  (chance 0.051). The high AUC reflects how distinct CEGv2 genes are (steady expression in
  every tissue) and many easy negatives, not label leakage.
- Candidates outside CEGv2 and DepMap common essentials, screened in at least 1,000 lines
  and needed by at least 25%: SSBP1 (rank 49), MRPL20 (67), RPUSD4 (74), HNRNPA1, HNRNPA2B1,
  ARF1. Lineages most dependent (at least 20 lines): SSBP1 myeloid 56%, lymphoid 50%;
  MRPL20 lymphoid 52%, oesophagus/stomach 38%; RPUSD4 breast 45%, lymphoid 41%.

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

### Positive: the essential-gene result holds against harder baselines

`prediction.json` (list `depmap-common`) and `baselines.json`, held-out AUC as
above.

- Independent list: DepMap 24Q2 common essentials (1,654 genes; CEGv2 has 606).
  Network gain over expression +0.017 to +0.111 in all 15 single-tissue
  networks, every 95% interval above 0.
- Plain correlation: a gene's mean correlation with training-fold essential
  genes, no cutoff and no network. Expression + correlation reaches 0.81–0.90
  in tumours and 0.78–0.83 in healthy tissue (CEGv2). The network adds to it in
  all 15 networks on both lists (tumour +0.010 to +0.041, healthy +0.038 to
  +0.073). In healthy tissue correlation adds almost nothing beyond the network
  (+0.001 to +0.025, significant in 3 of 7 for CEGv2, 4 of 7 for DepMap): the
  cutoff and neighbourhood keep the useful part. With all 15 datasets combined
  the two are equivalent (network beyond correlation −0.001 CEGv2, +0.007
  DepMap).
- STRING v12, combined score ≥ 700 (150,845 links among 10,489 of the genes):
  expression + STRING neighbour share 0.87–0.90, stronger than any single
  network. Each network still adds +0.005 to +0.038 on top of STRING in all 15
  networks and both lists; all networks combined +0.005 (CEGv2, interval
  −0.001 to +0.010) and +0.006 (DepMap, +0.003 to +0.010). STRING includes
  co-expression and curated knowledge, so it is partly circular here.
- Fewer samples (`samples.json`): LUAD, BRCA, healthy lung and thyroid rebuilt
  from 25, 50, 100 and 200 random samples, two draws each. The gain over
  expression is significant at every size (+0.012 to +0.079) and at 50 samples
  is already close to the full-data gain (LUAD +0.070, BRCA +0.053, lung
  +0.045, thyroid +0.019). Gene-function median AUC (Reactome) needs more
  samples: 0.58–0.63 at 25, 0.64–0.67 with all.

### No clear gain: PageRank added to discovery outside the training list

`pagerank-discovery.json`, 2026-11-11, code revision `c72288f`. The same combined
15-network CEGv2 model and five chromosome folds are used as in `robustness`.
Only the 11,394 genes outside CEGv2 are scored against DepMap's 1,073 common
essentials among them. Two additions were fixed before this run: unweighted
forward/reverse PageRank (alpha 0.85), and personalised PageRank from
training-fold CEGv2 seeds on undirected links (alpha 0.85, 50 iterations).
Personalised training features use inner chromosome folds to hide each
gene's own label.

| Model | DepMap essentials in top 100 | AUC | AUC gain over original network model (95%) |
| --- | ---: | ---: | --- |
| Expression | 57 | 0.840283 | |
| Expression + network | 71 | 0.865009 | |
| + ordinary forward/reverse PageRank | 72 | 0.865614 | +0.000605 (−0.001408, +0.002741) |
| + personalised PageRank | 71 | 0.864998 | −0.000011 (−0.001223, +0.001358) |

The original AUCs and top-100 counts reproduce exactly. Ordinary PageRank
adds one hit, but neither AUC-gain interval excludes zero (1,000 paired
gene-bootstrap resamples). This supports reporting PageRank as an ablation;
it does not establish an improvement. The PageRank additions were tested
only on the combined 15-network model, not repeated across adjacent-normal
networks or the other robustness settings. Input hashes, settings and every
model's top 100 genes are saved with the result. No parameter search or
selection by DepMap performance was performed.

### Feature groups, PageRank positions and other graph descriptors

`feature-comparison.json`, 2026-11-11, revision `8ff0bb1`, clean source tree.
The three main models were fixed before this run. Five chromosome holdouts
give CEGv2 recovery AUC over all 12,000 genes; DepMap discovery is evaluated
only among the 11,394 genes outside CEGv2.

| Model | CEGv2 recovery AUC | DepMap discovery AUC | DepMap essentials in top 100 |
| --- | ---: | ---: | ---: |
| Expression | 0.889323 | 0.840283 | 57 |
| Network features | 0.897185 | 0.853863 | 67 |
| Expression + network | 0.910664 | 0.865009 | 71 |

Network-only removes direct expression mean/SD but retains graph descriptors
and training-fold essential-neighbour share. Original expression/combined
baselines reproduce. The main table was not chosen from the best top-100
outcome of the additional models.

Fixed additions to the combined model (1,000 paired gene-bootstrap resamples):

| Addition | Top-100 hits | DepMap AUC gain over original (95%) |
| --- | ---: | --- |
| PageRank scores | 72 | +0.000605 (−0.001408, +0.002741) |
| PageRank positions | 69 | +0.000036 (−0.001934, +0.002008) |
| k-core depth + eigenvector centrality | 72 | +0.001953 (−0.000353, +0.004265) |
| Core/eigenvector + PageRank positions | 65 | +0.001674 (−0.001418, +0.004792) |
| Edge-weight summaries | 62 | −0.000718 (−0.003996, +0.002602) |

PageRank positions average ties and are scaled to [0,1], with 0 the highest
score. They preserve score ordering, changing the scale given to the linear
model. Core depth and eigenvector centrality use undirected, unweighted
links. Edge summaries are incoming/outgoing nonzero-weight mean, SD, mean
absolute weight and positive share, without explicit edge counts.

Fixed feature replacements, top-100 DepMap hits:

| Removed block | No replacement | PageRank scores | PageRank positions |
| --- | ---: | ---: | ---: |
| Degree (log in/out count) | 72 | 68 | 70 |
| Other descriptors (sign/weight/flow/clustering/arm share) | 73 | 70 | 70 |
| Essential-neighbour share | 67 | 69 | 67 |

None of these replacements clearly beats the original combined model in
AUC. Removing the essential-neighbour share loses 0.0061 AUC; replacing it
with PageRank scores loses 0.0045 (95% interval −0.0076 to −0.0014), and with
positions 0.0045 (−0.0078 to −0.0015). PageRank does not recover the label-aware
neighbourhood information. Counts fluctuate without a clear AUC gain, so
the 73-hit feature removal is not promoted as a better model. Personalised
PageRank is not part of this extension or the revised poster. New descriptors
were not repeated across the other robustness settings.

### Five feature groups: gene, local graph and global graph

`feature-groups.json`, 2026-11-12, revision `13e4b5b`, clean source tree.
Five models were specified before this comparison, using the same 15
networks, 12,000 genes, chromosome holdouts and CEGv2 training labels.
Local graph features are degree, sign/weight/flow/arm descriptors, local
clustering and training-essential neighbour shares. Global features are
ordinary forward/reverse PageRank scores and normalised rank positions,
k-core depth and eigenvector centrality. Global only excludes expression
summaries and label-derived neighbour shares. No personalised PageRank is used.

| Feature set | CEGv2 recovery AUC | DepMap discovery AUC | Known DepMap essentials in top 100 outside CEGv2 |
| --- | ---: | ---: | ---: |
| Gene only | 0.889323 | 0.840283 | 57 |
| Local graph only | 0.897185 | 0.853863 | 67 |
| Global graph only | 0.887262 | 0.843408 | 63 |
| Gene + local graph | 0.910664 | 0.865009 | 71 |
| Gene + local graph + global graph | 0.907649 | 0.866309 | 67 |

The three existing baselines reproduce exactly. Gene + local has the
highest CEGv2 recovery AUC and top-100 count in these five models. Adding
global features gives a DepMap AUC gain of +0.001300 (95% paired
gene-bootstrap interval -0.001965 to +0.004534), providing no clear
improvement. The poster's specific gene examples retain the gene + local
ranking; they are not selected from whichever model ranks them highest.
The five-row poster figure uses this report, with the older three-row
comparison retained as historical evidence.

### Global-feature and regularisation audit

`feature-audit.json`, 2026-11-12, revision `05bb5ec`, clean source tree.
The fixed-C=1 global-only result reproduces exactly (CEGv2 AUC 0.887262,
DepMap AUC 0.843408, 63 top-100 hits), as does the original gene + local
baseline. The global block has 90 columns: ordinary forward/reverse PageRank
scores and positions, k-core depth and eigenvector centrality, across 15
networks. Essentiality labels do not enter these descriptors or graph
construction; the observed graph retains the held-out nodes, so evaluation
is transductive prediction of unseen labels, not unseen-node induction.

| Fixed-penalty inputs | CEGv2 AUC | DepMap AUC | Top-100 DepMap hits |
| --- | ---: | ---: | ---: |
| Degree | 0.815815 | 0.762801 | 33 |
| Degree + degree positions | 0.846236 | 0.784529 | 50 |
| PageRank scores | 0.830518 | 0.784205 | 47 |
| PageRank positions | 0.827426 | 0.770214 | 46 |
| Core depth | 0.776094 | 0.726854 | 14 |
| Eigenvector centrality | 0.821249 | 0.784524 | 44 |
| All global descriptors | 0.887262 | 0.843408 | 63 |
| Degree + all global | 0.888395 | 0.848871 | 62 |

Global descriptors overlap with degree but are not explained by the fitted
degree-only baseline. In LUAD, reverse PageRank correlates with outgoing
degree at Spearman rho 0.973; core depth with incoming degree at 0.706;
eigenvector centrality with outgoing degree at 0.577. All per-network
degree/expression correlations are saved in the audit. Ten CEGv2 label
permutations within chromosomes give global-model AUCs 0.470543–0.523956.
No label leakage was found in the checked path. Unit tests verify that
outer-held-out labels affect neither tuning nor predictions for that fold.

Nested L2 penalty selection uses five outer chromosome folds and four inner
folds. C = 0.001 / 0.01 / 0.1 / 1 / 10 is chosen separately per model and
outer fold by pooled inner CEGv2 AUC; DepMap labels never select C. Scaling
and neighbour shares are fitted/recomputed using each training fold only.

| Nested-penalty inputs | CEGv2 AUC | DepMap AUC | Top-100 DepMap hits |
| --- | ---: | ---: | ---: |
| Gene only | 0.889850 | 0.840362 | 59 |
| Local graph only | 0.904098 | 0.859170 | 69 |
| Global graph only | 0.887690 | 0.843799 | 62 |
| Gene + local graph | 0.908583 | 0.861823 | 73 |
| Gene + local graph + global graph | 0.912805 | 0.866444 | 73 |

The all-feature model selects C=0.01 in every outer fold; gene + local
selects 0.01 in four folds and 10 in one. Global-only selects 0.1 throughout,
local-only 0.01, and gene-only 10 (the upper grid boundary). Adding globals
to gene + local now gains +0.004222 CEGv2 AUC (paired gene-bootstrap 95%
interval +0.001365 to +0.007404) and +0.004621 DepMap AUC (+0.001468 to
+0.007704). The fixed penalty therefore affected the original comparison.
Neither L2 nor nested selection guarantees better held-out AUC from more
features; the augmented model only contains the smaller model as a feasible
training solution. These exploratory intervals resample genes, not the full
network/model-selection procedure. Robustness checks on the original model
have not been repeated for the nested models. The poster now shows nested
results, with the original gene examples explicitly attributed to C=1.

### PageRank scores and positions together

`pagerank-joint.json`, 2026-11-11, revision `218be58`, clean source tree.
The same chromosome splits reproduce the three original baselines and the
separate score/position additions. Adding both encodings to the combined
model gives 63 DepMap essentials in the top 100, versus 71 originally.
Discovery AUC is 0.864418 versus 0.865009, a gain of -0.000591 (95% paired
gene-bootstrap interval -0.003186 to +0.001871). CEGv2 recovery AUC is
0.909176 versus 0.910664. This provides no clear improvement. Personalised
PageRank is excluded from this comparison.

### Specific candidate examples outside both essential-gene lists

`candidate-evidence.json`, 2026-11-11, revision `2562d7b`, clean source tree.
The audit retains all 29 genes in the original combined model's top 100
outside both CEGv2 and DepMap 24Q2 common essentials. Ranks refer to the
11,394 genes outside CEGv2. Cell-line support is the fraction of finite
DepMap gene effects at or below -0.5; missing measurements are excluded.
The poster illustrates the three highest-ranked candidates screened in at
least 1,000 lines and supported in at least 25% of them:

| Gene | Rank outside CEGv2 | Dependent / screened lines | Fraction |
| --- | ---: | ---: | ---: |
| MRPL20 | 11 | 291 / 1,150 | 25.3% |
| ATIC | 19 | 374 / 1,150 | 32.5% |
| SSBP1 | 29 | 367 / 1,150 | 31.9% |

These are illustrative context-dependent candidates, not a new performance
test or newly established universal essentials. Absence from these two
lists does not establish absence from the literature: growth-related roles
for ATIC (PMID 34803509) and SSBP1 (PMCID PMC5479245) are already published.
MRPL20 encodes a mitochondrial ribosomal protein (NCBI Gene 55052).

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

`ontology.json`: 15 of the 92 have specific GO biological-process terms after
leaving out electronic, expression-pattern and Reactome-assigned annotations.
By GO similarity to each pathway's members, the guessed pathway ranks at a mean
percentile of 0.71 (95% interval 0.60–0.80; chance 0.5); 2 of 15 in the top
tenth. Weak support from few genes.

### Positive in tumours: personalised PageRank after cleaning

`cleaning.json`, essential genes, held-out AUC per network and cleaning. Two-way
positive links with PageRank restarted from training-fold essential genes beat
the plain neighbour share by +0.03 to +0.09 in LUAD, LUSC, BRCA, COAD, LIHC,
STAD and LUAD without copy number (95% intervals above 0), not in THCA or PRAD.
In healthy networks the plain neighbour share is best except breast, colon and
stomach, where PageRank on each gene's 50 strongest links helps. The 2016
cleaning gives 0.36–0.67, at or below chance for the neighbour share.

### Positive: the essential-gene model finds essential genes it never saw

`robustness.json`. Models: expression against expression + degree + shape +
neighbour share, held-out chromosomes, as above.

- Unseen list: trained on CEGv2 only (606 genes) and scored on the 11,394
  other genes against the 1,073 DepMap common essentials among them (9.4%),
  which count as negatives in training. Network gain +0.048 to +0.148 in all
  15 single networks (AUC 0.68–0.78 to 0.79–0.84), all intervals above 0; all
  15 together 0.840 to 0.865 (+0.025, interval +0.018 to +0.031).
- Top 100 new genes per model: 52–71 DepMap essentials with one network
  against 21–50 from expression alone; 71 against 57 with everything; 9.4
  expected by chance. Healthy liver: 64 against 21.
- The 29 of those 100 that DepMap does not list (all networks): 27 were
  screened; a cell line needs a gene when its effect is at most −0.5. They are
  needed by a median of 12% of lines and 22% of them by at least a quarter of
  lines; other genes outside both lists: 0.2% and 2%. The misses are mostly
  dependencies of some cell lines, such as DPH2, ADRM1, MRPL20 and EMG1.
- Housekeeping complexes removed (904 genes: Reactome Translation, Metabolism
  of RNA, DNA Replication; HGNC proteasome groups; 309 of the CEGv2 genes and
  625 of the DepMap ones): gain +0.059 to +0.180 (CEGv2, 297 left) and +0.044
  to +0.164 (DepMap, 1,029 left) in all 15 networks; with all networks +0.017
  (interval +0.005 to +0.029) and +0.027. The network does not only recognise
  the ribosome, spliceosome and proteasome.
- Transfer: a model fitted on one network and applied to another, each
  network's features standardised on its own training genes. Own AUC
  0.83–0.89; transferred 0.65–0.89; median loss 0.026, largest 0.217 (tumour
  to healthy liver); 157 of 210 pairs above the target's own expression-only
  model. Median loss: tumour to tumour 0.005, healthy to healthy 0.033,
  healthy to tumour 0.045, tumour to healthy 0.073. Healthy lung, thyroid and
  prostate behave like the tumours; healthy colon, stomach, liver and breast
  form a second group.
- Single descriptors, AUC for CEGv2 with no model (`descriptor_auc`): in
  tumours essential genes send more links (out-degree 0.58–0.71, negative
  out-degree 0.62–0.74) and are less clustered (0.41–0.49); in healthy tissue
  out-degree is 0.46–0.61, in-degree 0.27–0.56 and clustering 0.49–0.68. No
  descriptor points the same way everywhere, which is why transfer between
  the two kinds of network loses more. The section on adjacent normal tissue
  below separates tissue state from data source.
- Cutoff (`threshold-2.0/` and `threshold-3.09/robustness.json`): networks
  rebuilt with |z| > 2.0 and > 3.09 keep 5.1–7.4 and 0.15–0.49 million links
  instead of 1.6–3.8 million.

  | | 2.0 | 2.35 | 3.09 |
  | --- | --- | --- | --- |
  | Unseen list, gain per network | +0.046 to +0.152 | +0.048 to +0.148 | +0.026 to +0.097 |
  | Unseen list, all networks | +0.023 | +0.025 | +0.014 |
  | DepMap essentials in the top 100, per network | 48–70 | 52–71 | 50–73 |
  | The same, all networks | 65 | 71 | 75 |
  | Without housekeeping complexes, CEGv2, per network | +0.053 to +0.179 | +0.059 to +0.180 | +0.018 to +0.122 |
  | The same, all networks (95%) | +0.012 (−0.000, +0.022) | +0.017 (+0.005, +0.029) | +0.005 (−0.006, +0.017) |
  | Without housekeeping complexes, DepMap, per network | +0.047 to +0.169 | +0.044 to +0.164 | +0.030 to +0.130 |
  | The same, all networks | +0.029 | +0.027 | +0.021 |
  | Transfer, median loss | 0.031 | 0.026 | 0.021 |

  Every per-network interval is above 0 at all three cutoffs. The results do
  not depend on the 2.35 cutoff; the sparsest networks give smaller gains,
  most in healthy tissue, and the all-network CEGv2 gain without housekeeping
  complexes is significant only at 2.35.

### Positive: normal tissue next to the tumours, and what follows the data source

`normal/robustness.json`, `normal/constraint.json`. Healthy and tumour
networks differ in tissue state and in source (GTEx, TCGA). The TCGA cohorts
also hold normal tissue next to the tumour (sample type 11): 35 (STAD) to 114
(BRCA) samples, processed with the tumours. 8 networks on the same genes and
settings, 1.2–2.5 million links.

| | 8 tumour | 8 adjacent normal | 7 healthy (GTEx) |
| --- | --- | --- | --- |
| Unseen list, gain per network | +0.052 to +0.080 | +0.039 to +0.097 | +0.048 to +0.148 |
| DepMap essentials in the top 100 | 58–71 | 50–66 | 52–68 |
| The same from expression alone | 38–45 | 21–37 | 21–50 |
| Without housekeeping complexes, CEGv2 | +0.068 to +0.106 | +0.038 to +0.130 | +0.059 to +0.180 |
| Constrained genes, gain per network | −0.009 to +0.021 | +0.006 to +0.026 | +0.016 to +0.043 |
| Constrained genes, networks with interval above 0 | 3 of 8 | 8 of 8 | 7 of 7 |
| The same without cell-essential genes | 3 of 8 | 7 of 8 | 7 of 7 |
| Out-degree alone, AUC for CEGv2 | 0.58–0.71 | 0.54–0.66 | 0.46–0.61 |
| Negative out-degree alone | 0.62–0.74 | 0.53–0.69 | 0.34–0.50 |
| Clustering alone | 0.41–0.49 | 0.46–0.63 | 0.49–0.68 |
| Share of links on the gene's own arm, alone | 0.56–0.63 | 0.52–0.56 | 0.51–0.57 |

Median AUC loss when a model moves between kinds of network (source to
target): tumour to tumour 0.005, normal to normal 0.011, normal to tumour
0.017, tumour to normal 0.024, healthy to normal 0.030, healthy to healthy
0.033, healthy to tumour 0.045, normal to healthy 0.049, tumour to healthy
0.073.

- The essential-gene results hold in all 8 adjacent-normal networks, every
  interval above 0: 23 networks in all, and non-tumour tissue from two sources.
- Constraint follows the tissue: every non-tumour network gains, from either
  source, and most tumour networks do not. Per cohort the adjacent-normal
  gain is the larger in 6 of 8 (not in PRAD or STAD, whose tumour networks
  gain too).
- Transfer follows the source: adjacent-normal networks exchange models with
  tumours almost as well as tumours among themselves, and worse with GTEx.
  The single descriptors sit between the two, the copy-number one (own arm)
  with the healthy tissues.

### Positive, modest: neighbours predict functions published later

`prospective.json`. GO biological-process annotations (propagated through is_a
and part_of) are split by the PubMed ID of their reference at 30,000,000,
reached in about August 2018, after all expression samples were collected.
Scores use only annotations citing earlier papers: the share of a gene's
neighbours with the term, averaged over the 15 networks, for the 1,952 terms
with 20–500 genes. A positive is a gene that gained a term from a later
experiment (expression-pattern evidence left out) with no earlier or undated
support; genes with any other support for the term are left out; the rest are
negatives. 7,654 positive gene–term pairs, 1,126 genes, 1,359 terms.

| Genes | Network neighbours | Earlier annotations of the gene | Links of the gene |
| --- | --- | --- | --- |
| All (mean percentile, 95%) | 0.631 (0.615–0.650) | 0.611 (0.586–0.634) | 0.471 (0.447–0.495) |
| Studied later, 2,151 genes | 0.626 (0.609–0.644) | 0.524 (0.499–0.549) | 0.475 (0.451–0.499) |
| No earlier annotation, 226 genes, 1,473 pairs | 0.638 (0.604–0.673) | 0.500 | 0.484 (0.435–0.535) |
| Top 1% of candidates, all / studied / none earlier | 4.7× / 4.2× / 5.4× chance | 5.2× / 2.2× / 1.0× | 0.6× / 1.2× / 1.0× |

- Well-studied genes keep collecting annotations, so among all genes the count
  of earlier annotations ranks later ones almost as well as the network (gain
  of the network +0.021, interval −0.008 to +0.049). Among genes that received
  any later experimental annotation the network gains +0.102 (+0.073 to
  +0.130), and the number of links is below chance throughout.
- Single networks: 0.53–0.60 among all genes; the mean of 15 gives 0.63.
- Per term (198 terms with at least 10 later genes): median AUC 0.64, 62%
  above 0.6; among genes studied later 0.61 and 57%.
- Genes with no earlier annotation and their rank among all candidates:
  SCUBE3 1 of 11,941 (regulation of BMP signalling), RPUSD2 6 and DUS3L 7
  (tRNA processing and modification), HELLS 16 (double-strand break repair by
  homologous recombination), THUMPD2 17 (regulation of RNA splicing), MPHOSPH9
  22 (negative regulation of organelle assembly), TMEM223 59 (respiratory
  chain complex assembly), METTL4 94 (snRNA metabolism), ZMYND12 104 (cilium
  assembly).
- Positive links only (`positive_neighbours`, added after `consensus` showed
  that only positive links join related genes): mean percentile 0.655 among
  all genes (+0.024, interval +0.013 to +0.035) and 0.640 among genes studied
  later (+0.014, +0.003 to +0.025), but fewer later annotations in the top 1%
  of candidates (4.0× and 3.7× against 4.7× and 4.2×), because shares over
  fewer links are noisier at the top. Mixed; the main score stays as it was.
- Limits: intervals resample genes, not terms; a PubMed ID dates the paper, not
  the experiment's first report; 0.63 is far below the cross-validated 0.75
  for Reactome pathways above.

### Positive: links shared by more tissues are known relationships

`consensus.json`. Every pair of the 12,000 genes on different chromosomes
(68.1 million), by the number of the 15 networks linking it positively:

| Networks | Pairs | Same Reactome pathway | STRING link | Same HGNC group | Curated TF–target |
| ---: | ---: | ---: | ---: | ---: | ---: |
| all pairs | 68,129,518 | 1.26% | 0.21% | 0.33% | 0.043% |
| 0 | 60,542,927 | 1.13% | 0.14% | 0.28% | 0.038% |
| 1 | 4,768,613 | 1.70% | 0.32% | 0.50% | 0.065% |
| 5 | 178,148 | 3.97% | 1.60% | 1.19% | 0.158% |
| 10 | 21,623 | 8.51% | 6.69% | 2.71% | 0.416% |
| 15 | 7,439 | 27.3% | 42.1% | 12.5% | 0.928% |

- Every step from 0 to 15 raises the Reactome and STRING shares. Reactome
  pathways and HGNC groups of at most 200 genes; STRING score at least 700,
  partly built from co-expression, so not independent.
- Without the housekeeping complexes (58.2 million pairs): 13.2% in one
  pathway at 15 networks against 1.04% of all pairs (12.8×), STRING 30.9%
  against 0.15%.
- Pairs linked negatively are not enriched at any count: 0.8–1.3× for Reactome
  wherever at least 3,000 pairs exist. Functional information is in positive
  links only.
- Healthy and tumour networks both contribute: pairs linked in all 7 healthy
  and all 8 tumour networks are 13.2% in one pathway (without housekeeping
  complexes); in all 7 healthy and no tumour 5.5%; in all 8 tumour and no
  healthy 4.2%.

### Positive: the network helps most where STRING knows least

`understudied.json`. Genes split at the terciles of their number of STRING
links: 0–4 (4,026 genes), 5–22 (4,096), 23–660 (3,878). Models fitted on all
genes: expression + STRING (neighbour share and log number of links) against
the same plus the network; gains are measured inside each group.

| List, group | Listed | AUC with STRING | Gain per single network | Significant | All networks |
| --- | ---: | --- | --- | ---: | --- |
| DepMap, 0–4 links | 96 | 0.66–0.76 | +0.025 to +0.092 | 15 of 15 | +0.031 (+0.012, +0.051) |
| DepMap, 5–22 | 348 | 0.84–0.87 | −0.002 to +0.018 | 3 | +0.002 (−0.007, +0.013) |
| DepMap, 23–660 | 1,210 | 0.88–0.90 | +0.001 to +0.013 | 11 | +0.003 (−0.002, +0.008) |
| CEGv2, 0–4 links | 19 | 0.68–0.84 | −0.008 to +0.091 | 2 | +0.020 (−0.025, +0.060) |
| CEGv2, 5–22 | 84 | 0.87–0.90 | −0.003 to +0.038 | 3 | −0.004 (−0.023, +0.014) |
| CEGv2, 23–660 | 503 | 0.83–0.85 | +0.002 to +0.021 | 14 | +0.002 (−0.006, +0.012) |

A knowledge base cannot help with genes it has no links for; a network built
from one tissue's expression data does. CEGv2 has too few genes in the first
group to show it.

### Positive in healthy tissue: genes that people cannot lose

`constraint.json`. gnomAD v2.1.1: a gene is constrained when LOEUF is below
0.35; 9,665 of the 12,000 genes have at least 10 expected loss-of-function
variants, 2,367 of them constrained (37% of the cell-essential ones). The
baseline is expression plus the log expected variant count, so gene length
cannot pass as a network effect.

| Networks | Baseline AUC | Network gain | Significant | Without cell-essential genes |
| --- | --- | --- | ---: | --- |
| 7 healthy, one at a time | 0.70–0.73 | +0.016 to +0.043 | 7 of 7 | +0.015 to +0.043, 7 of 7 |
| 8 tumour, one at a time | 0.71–0.73 | −0.009 to +0.021 | 3 of 8 (THCA, PRAD, STAD) | −0.011 to +0.028, 3 of 8 |
| All 15 | 0.766 | +0.023 (+0.014, +0.031) | | +0.025 (+0.015, +0.033) |

Healthy-tissue networks carry information about constraint that tumour
networks mostly lack, and it is separate from cell essentiality.

### Known TFs and their curated targets

`regulons.json`. 388 TFs with at least 10 CollecTRI targets among the genes
(median 36). Per network, hypergeometric test of the overlap between a TF's
neighbours and its targets, Benjamini–Hochberg at 5%.

- Significant TFs per network: 30 (healthy stomach) to 80 (LUAD); median fold
  among them 3.8–6.3×. 189 TFs are significant somewhere, 45 in at least 8
  networks, 5 in all 15: E2F1, E2F2, FOXM1, SPI1, IRF1. In 13–14 networks:
  MYBL2, RUNX3, RELB, ATF3, FOS, JUN, JUNB, SREBF2, NFKB2, STAT1, EGR1, IRF9,
  IKZF1. These are the regulators of the cell cycle, interferon response,
  immune cells, immediate-early response and cholesterol synthesis, the
  programmes that vary most between samples. HNF4A is significant in 10
  networks, HNF1A in LUAD, COAD and STAD.
- A TF is seen when its own mRNA varies: its fold correlates with the spread of
  its expression in tumours (Spearman +0.28 to +0.37), less in healthy tissue
  (−0.09 to +0.23).
- Specificity is modest: a TF's fold for its own targets ranks at 0.54–0.65
  among its folds for the other 387 TFs' targets.
- Link signs match curated activation or repression at chance (0.65–0.88
  against 0.59–0.88 expected).
- Not tissue-specific: for 191 TFs at least 2-fold above their median of the
  other tissues in some GTEx tissue, target recovery in the tumour networks of
  those tissues is no higher than in the others (difference in log2 fold
  +0.05, interval −0.03 to +0.14; significant in 28% of own-tissue networks
  against 25%). PAX8 in thyroid and PPARA and NR1I3 in liver are among the
  exceptions.

So the network confirms regulators that are already known and whose targets
move with them; it still cannot tell a regulator from its targets
(question 1), their sign, or where a tissue's own regulators act.

### Tested by intervention: hubs, links and knockdown effects

`hubs.json`, `perturbation.json`, `partners.json`. Genome-scale Perturb-seq
(Replogle et al. 2022): CRISPRi knockdown of every expressed gene in K562
(8,033 knockdowns among the 12,000 genes, 6,739 of them measured), and of
essential genes in K562 (1,738) and RPE1 (1,809). One knockdown per gene, at
least 25 cells.

Hubs (WGCNA connectivity on each of the 15 datasets):

| Hubs against | AUC, unsigned | Among the top 100 hubs | Expected |
| --- | --- | --- | ---: |
| TFs (Lambert) | 0.37–0.62 | 0–18 | 9.6 |
| CEGv2 essential | 0.49–0.79 | 0–26 | 5.1 |
| DepMap common essential | 0.52–0.78 | 0–54 | 13.8 |
| Constrained (LOEUF < 0.35) | 0.51–0.65 | 9–71 | 24.5 |
| Broad-effect knockdown (K562) | 0.52–0.70 | 0–26 | 10.0 |

Which knockdowns change many genes (top 10% by differentially expressed
genes; K562 genome-wide, 289 genes or more):

- 38% of cell-essential genes against 3% of the others. TFs: 26 of 1,137, odds
  ratio 0.18 (p = 4 × 10⁻²⁷); without cell-essential genes 15 of 1,084, odds
  ratio 0.41 (p = 3 × 10⁻⁴). By this measure the genes that reshape a cell's
  transcriptome when removed are its housekeeping machinery, not its TFs.
- The network adds +0.039 to +0.121 AUC to expression plus knockdown strength,
  cell number and target expression in all 15 networks (all networks 0.784 to
  0.805), but without cell-essential genes only in 2 of 15 (all networks
  −0.019): it finds the essential genes again.
- The essential-gene screens select their targets, so they cannot repeat this
  (network gain significant in 0 and 1 of 15).
- Direction: among one-way edges where either knockdown affects the other gene,
  the edge points from the stronger one in 38–55% (K562 genome-wide).

Linked pairs, pooled over a screen (pairs on different chromosomes; affected =
top 1% of |z|):

| Networks linking the pair positively | K562 genome-wide | K562 essential | RPE1 essential |
| ---: | --- | --- | --- |
| 0 | 0.95% affected, 21% fall | 0.92%, 25% | 0.94%, 21% |
| 5 | 1.9%, 36% | 1.8%, 37% | 2.0%, 40% |
| 10 | 3.9%, 47% | 3.9%, 51% | 3.9%, 40% |
| 15 | 15.6%, 92% | 16.8%, 94% | 10.4%, 59% |

Negative links: 0.8–1.5% affected at every count. Taken alone this table says
that links shared by all tissues are causal and signed. The next one shows
why that reading is wrong.

The same question with the knockdown as the unit (`partners.json`): genes are
ranked within each knockdown, a gene's partners are those linked to it
positively in at least 8 networks (or all 15), and specificity asks whether
the partners respond more to this knockdown than to the others in the screen.
"Rank" is the median over knockdowns of the partners' mean rank (0.5 = like
any gene); "specific" is the share of knockdowns among the 5% most specific.

| K562 genome-wide | Knockdowns | Rank, falling | Rank, changing | Specific, falling / changing |
| --- | ---: | ---: | ---: | --- |
| Partners in ≥ 8 networks, all | 2,410 | 0.505 | 0.517 | 9.2% / 11.4% |
| Without cell-essential targets | 1,528 | 0.499 | 0.508 | 5.8% / 6.5% |
| Without housekeeping complexes | 1,746 | 0.498 | 0.508 | 5.8% / 7.6% |
| TF knockdowns | 225 | 0.502 | 0.506 | 7.6% / 6.2% |
| Partners in all 15, all | 188 | 0.557 | 0.573 | 24.5% / 18.6% |
| Without cell-essential targets | 66 | 0.501 | 0.516 | 4.5% / 10.6% |
| Without housekeeping complexes | 110 | 0.503 | 0.502 | 2.7% / 9.1% |

- The pooled enrichment comes from knockdowns of essential genes inside the
  housekeeping modules (ribosome, mitochondrial ribosome, respiratory chain,
  proteasome), whose members do fall specifically when one of them is knocked
  down: RPS19, RPS29, RPS9, MRPL27, MRPL43, NDUFB4, COX17, PSMG3 are among the
  25 most specific.
- Outside them a typical gene's partners do not respond to its knockdown.
  Exceptions: SREBF2 (18 partners, rank 0.97, the most specific of all 8,033),
  ZEB1, EGR3, ZKSCAN1, TSC2, DCUN1D1 and NEDD8.
- RPE1: partners in all 15 networks fall under most essential-gene knockdowns
  (rank 0.74; 0.88 without housekeeping complexes, 42 knockdowns), but only
  8% and 5% are among the most specific, and negatively linked genes rise in
  69% of knockdowns just as unspecifically. In these cells, which arrest when
  an essential gene is lost, the network's signs follow a proliferation
  programme: positive partners of growth genes fall and negative ones rise
  whatever the knockdown.
- Limits: two cell lines against tissue networks; a knockdown's effect is read
  days later, so indirect effects dominate; a partner that needs more than a
  CRISPRi knockdown to respond would be missed.

### Negative or inconclusive uses

- PageRank as model features (`prediction.json`): plain and reverse PageRank
  add +0.000 to +0.017 AUC beyond expression and shape on CEGv2 and DepMap
  common essentials (significant in 10 of 15 single networks), −0.004 to
  +0.017 on cancer genes, and +0.001 to +0.003 with all networks combined.
  Personalised PageRank from training-fold listed genes adds −0.003 to +0.011
  beyond the neighbour share (combined −0.002 to 0.000). It has to be
  cross-fitted: computed from the whole training fold, each seed's own restart
  mass enters its score and models lose 0.04–0.25 AUC on held-out genes
  (healthy liver 0.773 against 0.867 cross-fitted).
- Acyclic walks (`acyclic.json`): the Eades–Lin–Smyth order keeps 73–93% of
  edges. Personalised PageRank from essential genes along the kept edges gives
  AUC 0.53–0.71 forward and 0.62–0.80 backward, below the walk on all links in
  13 of 16 networks each; ignoring direction on the kept edges is within ±0.03.
  Collapsing strongly connected groups gives chance: one group holds 99–100% of
  genes in most networks, and 56–100% even among each gene's 50 strongest
  partners.
- Copy-number direction (`copy_number.json`, LUAD): for genes A and B on
  different chromosomes whose copy number moves their own expression (r ≥ 0.3,
  6,465 genes), the direction score is |corr(CN_A, E_B)| − |corr(CN_B, E_A)|.
  It points from TF to target in 49.7% of CollecTRI pairs (4,723), 49.7% of
  DoRothEA and 50.3% of TRRUST pairs; the same TFs with random non-targets give
  47.5%, and arm-level copy number gives the same picture (49.1% against
  47.2%), so the small difference is regional, not gene-specific. On the 212
  CollecTRI pairs that are network links: 56.6% (interval 50.0–63.2%). The
  share of a gene's links pointing away does not find Lambert TFs (AUC 0.515,
  random 95% range 0.476–0.523).
- DepMap levels (`dependency_levels.json`): neighbour means of each lineage's
  mean gene effect correlate with held-out genes' effects (Spearman 0.11–0.42,
  i.e. essentiality again), but the healthy prostate network is best for 6 of 7
  lineages and the matched network is never significantly better than the
  others. Lineage selectivity is barely predicted (Spearman ≤ 0.14).

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
