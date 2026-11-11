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
| 15 | Does the essential-gene result hold on an independent list? | **Positive** | DepMap 24Q2 common essentials (1,654 genes, CRISPR screens across cell lines): network descriptors add +0.017 to +0.111 AUC to expression in all 15 single-tissue networks, every 95% interval above 0. |
| 16 | Does the network add more than plain correlation with known essential genes? | **Positive** | Expression + mean correlation with training-fold essential genes (no cutoff, no network) against the same plus the network: the network adds +0.010 to +0.041 in all 8 tumour and +0.038 to +0.073 in all 7 healthy networks, on both lists. In healthy tissue correlation adds little beyond the network (+0.001 to +0.025, significant in 3–4 of 7). With all 15 datasets combined the network adds −0.001 (CEGv2) and +0.007 (DepMap). |
| 17 | Does it add to an existing interaction database (STRING)? | **Positive** | STRING v12 (combined score ≥ 700) is stronger alone (AUC 0.87–0.90 with expression), yet each single network adds +0.005 to +0.038 on top of it in all 15 networks and both lists; all networks combined add +0.005 (CEGv2, interval −0.001 to +0.010) and +0.006 (DepMap). STRING partly includes co-expression and curated knowledge. |
| 18 | How many samples does a network need? | **Positive** | Networks rebuilt from 25–200 random samples (LUAD, BRCA, healthy lung and thyroid, two draws each): the essential-gene gain is significant from 25 samples on (+0.012 to +0.079) and close to the full-data gain from 50; gene-function median AUC rises from 0.58–0.63 at 25 samples to 0.64–0.67 with all. |
| 19 | Do PageRank scores add to the essential-gene model? | **Negative (too small to matter)** | Plain and reverse PageRank added to expression and shape: +0.000 to +0.017 AUC on CEGv2 and DepMap (significant in 10 of 15 single networks each), −0.004 to +0.017 on cancer genes; all networks combined +0.001 to +0.003. Personalised PageRank from training-fold listed genes (cross-fitted so a gene's own label never enters its score) on top of the neighbour share: −0.003 to +0.011, combined −0.002 to 0.000. The neighbour share already carries what the walk finds. |
| 20 | Does PageRank on an acyclic version of the network help? | **Negative** | Eades–Lin–Smyth order keeps 73–93% of edges. Personalised PageRank from essential genes along the kept edges: forward AUC 0.53–0.71, backward 0.62–0.80, both below the walk on all links in 13 of 16 networks; ignoring direction on the kept edges changes ±0.03. Collapsing strongly connected groups is useless: one group holds 99–100% of genes in most networks (56–100% for each gene's 50 strongest partners). |
| 21 | Can copy number in tumours give links a direction? | **Negative** | LUAD, 512 tumours, 6,465 genes whose copy number moves their own expression (r ≥ 0.3). For curated TF→target pairs on different chromosomes the copy-number direction points from TF to target in 49.7% (CollecTRI, 4,723 pairs; DoRothEA 49.7%, TRRUST 50.3%), at chance; 56.6% on the 212 pairs that are network links (interval 50.0–63.2%), and 50.9% with arm-level copy number. The share of a gene's links pointing away from it does not find TFs (AUC 0.515, random 0.476–0.523). |
| 22 | Does a tissue's network predict how strongly that tissue's cell lines depend on each gene? | **Positive in general, negative per tissue** | Neighbour means of DepMap effects correlate with the held-out genes' mean effect per lineage (Spearman 0.11–0.42), but the healthy prostate network is best for 6 of 7 lineages; the matched network is never significantly better than the others. Selectivity (lineage minus others) is barely predicted (Spearman ≤ 0.14). |
| 23 | Do Gene Ontology annotations support the uncharacterised-gene guesses? | **Weak positive** | Only 15 of the 92 genes have specific non-electronic GO process terms not assigned by Reactome. Their guessed pathway ranks at a mean percentile of 0.71 (95% interval 0.60–0.80; chance 0.5) by GO similarity among all Reactome pathways with annotated members. |
| 24 | Does a model trained on one essential-gene list find the essential genes of a list it never saw? | **Positive** | Models trained on CEGv2 alone are scored on the 11,394 genes outside it against the 1,073 DepMap common essentials among them (9.4%). The network raises AUC by +0.048 to +0.148 in all 15 single-tissue networks (every 95% interval above 0) and from 0.840 to 0.865 with all 15. Of the 100 top-scoring new genes, 52–71 are DepMap essentials with one network (21–50 from expression alone) and 71 with all networks (57; 9 expected by chance). Of the 29 that are not, the 27 screened are needed by a median of 12% of DepMap cell lines, against 0.2% for other unlisted genes. |
| 25 | Does the essential-gene gain survive without the large housekeeping complexes? | **Positive** | Without the 904 genes of Reactome Translation, Metabolism of RNA and DNA Replication and the HGNC proteasome groups (309 of 606 CEGv2 and 625 of 1,654 DepMap essentials), the network still adds +0.059 to +0.180 (CEGv2) and +0.044 to +0.164 (DepMap) in all 15 single-tissue networks, and +0.017 and +0.027 with all 15; every interval above 0. |
| 26 | Does a model fitted on one tissue's network work on another's? | **Positive** | Each network's fitted model (expression, shape, neighbours) applied to every other network on held-out chromosomes: AUC 0.65–0.89 against 0.83–0.89 on its own network, median loss 0.026; 157 of 210 pairs still beat the target's own expression-only model. Median loss is 0.005 between tumour networks, 0.033 between healthy ones, 0.045 from healthy to tumour and 0.073 from tumour to healthy. The descriptors differ: essential genes have more links in tumours (out-degree alone AUC 0.58–0.71; healthy 0.46–0.61) and sit in tighter clusters in healthy tissue (clustering 0.49–0.68; tumours 0.41–0.49). Much of this difference follows the data source, not the tissue (question 33). |
| 27 | Do a gene's neighbours predict functions that were published later? | **Positive, modest** | GO biological-process annotations are split by PubMed ID at 30,000,000 (about August 2018). Scores use only the earlier papers; the test is the 7,654 gene–term pairs (1,126 genes, 1,359 terms) first supported by later experiments. Such a pair ranks at a mean percentile of 0.63 among never-annotated genes (95% interval 0.61–0.65; chance 0.5; single networks 0.53–0.60). A gene's count of earlier annotations does nearly as well among all genes (0.61), so the fair comparison is among genes studied later: 0.63 against 0.52 (+0.10, interval +0.07 to +0.13), with 4.2× the chance number of later annotations in the top 1% of candidates. For the 226 genes with no earlier annotation: 0.64 (0.60–0.67) and 5.4×. Examples: SCUBE3 ranked 1st of 11,941 for regulation of BMP signalling, HELLS 16th for repair by homologous recombination, TMEM223 59th for respiratory-chain assembly. |
| 28 | Are links found in more tissues more often known relationships? | **Positive** | Gene pairs on different chromosomes, by the number of networks that link them positively: the share in one Reactome pathway (at most 200 genes) rises at every step from 1.3% of all pairs to 4.0% at 5 networks, 8.5% at 10 and 27% at all 15 (7,439 pairs, 22×); STRING links from 0.2% to 42% (203×), same HGNC group from 0.3% to 12%, curated TF–target pairs from 0.04% to 0.9% (21×). Without the housekeeping complexes: 13% in one pathway (13×) and 31% in STRING at all 15. Negative links are not enriched (0.8–1.3× wherever at least 3,000 pairs exist). |
| 29 | Does the network help most where STRING knows least? | **Positive on DepMap** | Genes are split at the terciles of their number of STRING links (0–4, 5–22, 23–660). Among those with 0–4 links (96 DepMap common essentials in 4,026 genes) each single network adds +0.025 to +0.092 AUC to expression + STRING in all 15 networks (every interval above 0) and all networks together +0.031. With 5–22 links the gain is −0.002 to +0.018 (significant in 3 networks) and with more +0.001 to +0.013 (11). CEGv2 has only 19 genes in the first group: +0.02 with all networks (interval −0.02 to +0.06), significant in 2 single networks. |
| 30 | Do network descriptors predict genes that people cannot lose? | **Positive outside tumours** | Constrained genes: gnomAD v2.1.1 LOEUF below 0.35, 2,367 of the 9,665 genes with at least 10 expected loss-of-function variants. Against expression plus the expected variant count (AUC 0.70–0.73) the network adds +0.016 to +0.043 in all 7 healthy networks but only in 3 of 8 tumour networks (−0.009 to +0.021); all 15 together +0.023 (interval +0.014 to +0.031). Without cell-essential genes the gains are the same (healthy +0.015 to +0.043; together +0.025), so this is not the essential-gene result again. Normal tissue next to the tumours, processed by TCGA like the tumours, gains in all 8 networks (question 33), so the difference is between tumour and non-tumour tissue. |
| 31 | For a known TF, do its network neighbours recover its curated targets? | **Positive for general programmes, not tissue-specific** | 388 TFs have at least 10 CollecTRI targets among the genes. A TF's targets are over-represented among its neighbours (hypergeometric, 5% FDR per network) for 30–80 TFs per network, median 3.8–6.3× where significant; for 189 TFs in at least one network, 45 in at least 8, and E2F1, E2F2, FOXM1, SPI1 and IRF1 in all 15. The recurring TFs drive the programmes that vary most between samples: cell cycle (E2F, FOXM1, MYBL2), interferon and immune cells (IRF1, STAT1, IRF9, SPI1), immediate-early response (FOS, JUN, ATF3, EGR1), cholesterol synthesis (SREBF2). Limits: a TF's own targets rank at only 0.54–0.65 among all TFs' target sets (0.5 = no closer than to others'), link signs match curated signs at chance, and TFs elevated in a tissue (GTEx) do not recover their targets better in that tissue's tumour networks (difference +0.05, interval −0.03 to +0.14). |
| 32 | Do the essential-gene checks depend on the link cutoff? | **Positive (they hold)** | Networks rebuilt with cutoffs of 2.0 and 3.09 instead of 2.35 keep 5.1–7.4 and 0.15–0.49 million links instead of 1.6–3.8 million. Question 24 holds in all 15 networks at both cutoffs (gain +0.046 to +0.152 and +0.026 to +0.097, every interval above 0; 48–70 and 50–73 of the top 100 new genes are DepMap essentials, 65 and 75 with all networks). Question 25 holds in all 15 single networks at both; with all networks combined it holds on DepMap (+0.029 and +0.021) but not clearly on CEGv2 (+0.012, interval −0.000 to +0.022, and +0.005, −0.006 to +0.017). Transfer (question 26) loses a median of 0.031 and 0.021. The sparser networks give smaller gains, mostly in healthy tissue (lung +0.054 at 3.09 against +0.108). |
| 33 | Do the results hold in normal tissue next to the tumours, and do healthy–tumour differences follow the tissue or the data source? | **Positive; mixed for the differences** | 8 more networks from the adjacent normal samples of the TCGA cohorts (35–114 samples each; same genes and settings; same processing as the tumours, unlike GTEx). Question 24 holds in all 8 (gain +0.039 to +0.097, every interval above 0; 50–66 of the top 100 new genes are DepMap essentials against 21–37 from expression), and so does question 25 (+0.038 to +0.130 CEGv2, +0.041 to +0.107 DepMap). Constrained genes (question 30) gain in all 8 (+0.006 to +0.026; the tumours of the same cohorts in 3 of 8), so that result follows the tissue. Model transfer follows the data source instead: median loss 0.017–0.024 between adjacent-normal and tumour networks against 0.030–0.049 between adjacent-normal and GTEx healthy ones. Single descriptors lie in between (out-degree alone 0.54–0.66, clustering 0.46–0.63). |
| 34 | Are hubs regulators under the standard (WGCNA) definition, and what are they instead? | **Negative for regulators** | WGCNA connectivity (unsigned, power 6; signed, power 12) on the 15 datasets ranks Lambert TFs at AUC 0.37–0.62 (link counts 0.37–0.57); the 100 strongest hubs hold 0–18 TFs (9.6 expected). Hubs are more often cell-essential (CEGv2 AUC 0.49–0.79) and constrained in people (9–71 of the top 100, 24.5 expected; above that in 12 of 15 datasets). |
| 35 | Tested by intervention: which knockdowns change many genes, and does the network find them? | **TFs do not; the network finds essential genes again** | Genome-wide Perturb-seq in K562 (Replogle 2022; 8,033 knockdowns among the genes). The 10% of knockdowns with the most differentially expressed genes (289 or more) are 38% of the cell-essential genes and 3% of the others; TFs are depleted among them (26 of 1,137; odds ratio 0.18, p = 4 × 10⁻²⁷; 0.41 without essential genes, p = 3 × 10⁻⁴). The network adds +0.039 to +0.121 AUC to expression and experiment covariates in all 15 networks, but without cell-essential genes only in 2 of 15. One-way edges agree with the stronger of the two knockdown effects in 38–55% of pairs (chance 50%). |
| 36 | Do linked genes respond to each other's knockdown? | **Pooled yes; per knockdown mostly no** | Pooled over a screen, pairs linked positively in all 15 networks are affected 16× more often than unlinked pairs in K562 (15.6% against 0.95%; RPE1 11×) and the partner falls in 92% (RPE1 59%; 22–28% of all affected pairs); negative links show nothing. With the knockdown as the unit and genes ranked within each knockdown (`partners`), this mostly goes away: a gene's partners (linked in at least 8 networks; 2,410 knockdowns) have a median rank of 0.505 for falling and 0.517 for changing (0.5 = like any gene), and its own knockdown is among the 5% most specific for them in 9–11% of cases. Without cell-essential targets, or without the housekeeping complexes, the partners do not fall (0.499 and 0.498) and the specific share is 6–8%. In RPE1 the partners fall under most essential knockdowns (rank 0.74 for partners in all 15 networks) but under other knockdowns as well (most specific: 8%), and negative partners rise in 69% of knockdowns, equally unspecific: a proliferation programme, not a pair effect. Specific exceptions exist: SREBF2 (its 18 partners rank 0.97, the most specific of all 8,033 knockdowns), ribosomal, mitochondrial-ribosome, respiratory-chain and proteasome-assembly genes, TSC2. |
| 37 | Does adding PageRank improve discovery of essentials outside the training list? | **No clear gain** | Combined 15-network model, trained on CEGv2 and evaluated on other genes against DepMap: original 71 essentials in the top 100, ordinary forward/reverse PageRank added 72, cross-fitted personalised PageRank added 71. AUC 0.86501, 0.86561 and 0.86500; gain intervals over the original include zero for both additions. |
| 38 | How do expression-only, network-only and combined models compare, and do PageRank positions or more graph descriptors help? | **Combined best; no clear gain from extra descriptors** | CEGv2 recovery AUC 0.889 / 0.897 / 0.911; DepMap essentials among top 100 outside CEGv2 57 / 67 / 71. PageRank scores add one hit, positions give 69, and k-core/eigenvector additions 72, with no clear AUC gain from these additions to the combined model. Replacing essential-neighbour share with scores/positions loses AUC. |

The five-model comparison uses gene only, local graph only, global graph
only, gene + local graph, and gene + local graph + global graph. The poster
now uses nested penalty selection (`feature-audit.json`): CEGv2 recovery
AUCs are 0.890 / 0.904 / 0.888 / 0.909 / 0.913; known DepMap essentials in
the top 100 outside CEGv2 are 59 / 69 / 62 / 73 / 73. The original fixed-C=1
comparison remains in `feature-groups.json` (AUCs 0.889 / 0.897 / 0.887 /
0.911 / 0.908, top-100 counts 57 / 67 / 63 / 71 / 67).
The global group contains ordinary PageRank scores and positions, k-core
depth and eigenvector centrality. With nested penalty selection, adding
globals changes DepMap AUC from 0.861823 to 0.866444 (gain +0.004621;
paired gene-bootstrap 95% interval +0.001468 to +0.007704). Gene bootstraps
do not include all network/model-selection uncertainty.

```sh
uv run gene-network-analysis feature-groups local/public ../references --output results/public
uv run gene-network-analysis feature-audit local/public ../references --output results/public
uv run gene-network-analysis figures --poster-only
```

`feature-audit` checks the global-only result against degree, degree positions
and each centrality separately, plus ten within-chromosome label permutations.
It also repeats all five models with nested penalty selection: five outer
chromosome folds, four inner chromosome folds, and C in
0.001 / 0.01 / 0.1 / 1 / 10. Scaling and label-derived neighbour shares are
recomputed using each inner/outer training fold only. C is selected using
inner CEGv2 AUC; DepMap labels are excluded from selection. The original
fixed comparison uses C = 1. The global-only score reproduces (0.887262);
ten within-chromosome label shuffles give 0.470543–0.523956. Degree-only
scores 0.815815, and degree with rank positions 0.846236. Checked feature
paths show no essentiality-label leakage. Global features comprise 90
columns across 15 graphs, and remain high after nested tuning (0.887690).

```sh
uv run gene-network-analysis feature-audit local/public ../references --output results/public
```

The later `pagerank-joint` comparison (`results/public/pagerank-joint.json`)
tests scores and positions together: 63 top-100 DepMap hits and AUC 0.864418,
with no clear gain over the original combined model. The `candidates` audit
(`results/public/candidate-evidence.json`) lists all 29 original top-100
genes outside both CEGv2 and DepMap 24Q2 common essentials, with screened
cell-line support. MRPL20, ATIC and SSBP1 are illustrative context-dependent
candidates; absence from these lists does not establish biological novelty.

Questions 8–13 come from `gene-network-analysis applications` and question 10 from
`cleaning` (`applications.py`, `gene_sets.py`; `results/public/applications.json`,
`cleaning.json`); question 14 from `structure`. Questions 15 and 19 come from `predict`,
16–17 from `baselines` (`baselines.py`), 18 from `samples` (`samples.py`), 20 from
`acyclic`, 21 from `copy-number` (`copy_number.py`), 22 from `dependencies` and 23 from
`ontology` (`ontology.py`); each writes the JSON of the same name to `results/public/`.
Questions 24–26 come from `robustness`, 27 from `prospective`, 28 from `consensus`, 29 from
`understudied`, 30 from `constraint` and 31 from `regulons`, each in the module and
`results/public/` JSON of the same name; question 32 is `robustness` on networks built with
other cutoffs (`results/public/threshold-2.0/`, `threshold-3.09/`) and question 33 is
`robustness` and `constraint` on the networks that `normals` adds (`results/public/normal/`).
Question 34 comes from `hubs`, 35 from `perturbation` and 36 from `perturbation` and `partners`.
Question 37 comes from `pagerank-discovery` (`pagerank_discovery.py`,
`results/public/pagerank-discovery.json`); question 38 from `feature-comparison`
(`feature_comparison.py`, `results/public/feature-comparison.json`). Questions 3–7 come from
`gene-network-analysis predict` (`prediction.py`,
`results/public/prediction.json`); question 1 from `run`, question 2 from
`structure` (`results/public/`). Essential genes are the 684 core essential
genes of Hart et al. 2017 (CEGv2), 606 of them among the 12,000 genes. Models
are logistic regressions scored on held-out chromosomes, so neighbouring genes
never sit on both sides of a split; gains are paired bootstraps over genes.

The earlier work on a 2016 graph of unknown origin is kept in the sections
below; that graph is not distributed and is not used for conclusions.

## Run

### PageRank addition to essential-gene discovery

The `pagerank-discovery` check compares the original combined 15-network
CEGv2-to-DepMap discovery model with two fixed additions: unweighted forward
and reverse PageRank (alpha 0.85), and personalised PageRank restarted from
training-fold CEGv2 essentials (alpha 0.85, 50 iterations). Personalised scores
for training genes use inner chromosome folds to hide their own labels.
These variants and settings were fixed before running this comparison on
2026-11-11; no DepMap-driven parameter search is performed.

The original expression and expression-plus-network models are rerun with
the same five chromosome folds. Only genes outside CEGv2 are scored against
DepMap. Report top-100 hits and AUC for every variant, with paired gene-bootstrap
intervals for the PageRank additions relative to the original network model.
This additional analysis was requested during poster preparation. Neither
addition gives a clear gain: ordinary PageRank adds one top-100 hit (72 versus
71), with an AUC-gain interval spanning zero; personalised PageRank gives 71.
The original result stays in `results/public/robustness.json`.

```sh
uv run gene-network-analysis pagerank-discovery local/public ../references --output results/public
```

The report, `results/public/pagerank-discovery.json`, includes input hashes,
the code revision, model definitions and each model's top 100 candidates.

### Main commands

`candidates` audits every gene in the original combined model's top 100 that
is absent from both CEGv2 and DepMap common essentials. It checks reference
membership again, keeps the original rank, and reports the fraction of
screened cell lines with DepMap gene effect at most -0.5, excluding missing
measurements from the denominator. Absence from the two lists is not a claim
of novelty or absence of published essentiality evidence.

```sh
uv run gene-network-analysis candidates results/public/feature-comparison.json ../references --output results/public
```

The `pagerank-joint` check tests scores and rank positions in the same model,
next to the original baseline and the two individual encodings. Network-only
and combined joint models are fixed before this run; no personalised variant
or parameter search is included. Rank and score convey the same ordering,
but their combination can give the linear classifier a nonlinear dependence
on centrality. This does not guarantee a gain on the biological data.

```sh
uv run gene-network-analysis pagerank-joint local/public ../references --output results/public
```

The poster's three-row comparison uses expression-only, network-only and
expression-plus-network models. `feature-comparison` reproduces the existing
expression/combined baselines, adds the network-only row, and evaluates known
CEGv2 reconstruction separately from DepMap discovery outside CEGv2. Network
features include the training-only essential-neighbour share.

Before running the 2026-11-11 extension, the following comparisons were fixed:
ordinary PageRank scores versus top-first percentile positions; replacing
degree, other structural descriptors or neighbour share with each; adding
k-core depth and eigenvector centrality; adding nonzero edge-weight summaries.
The original three models remain the main table, so it is not selected by
the best DepMap outcome. Additional models are exploratory comparisons.
Personalised PageRank is excluded from this extension and from the poster.

```sh
uv run gene-network-analysis feature-comparison local/public ../references --output results/public
uv run gene-network-analysis figures --poster-only
```

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
GISTIC2 copy number; `normals` builds 8 further graphs from the adjacent
normal samples in the same files), and GTEx v8 TPMs for lung, breast, colon (transverse),
thyroid, prostate, stomach and liver. Sources and checksums are in
`~/workspace/data/bioinf/references/README.md`. The gene lists are CEGv2 core
essential genes (<https://github.com/hart-lab/bagel>) and the OncoKB cancer
gene list.
`applications` also reads Reactome pathways (`ReactomePathways.gmt.zip`), DepMap
24Q2 CRISPR gene effects with models and common essentials (figshare 25880521)
and the GWAS Catalog associations (`gwas-catalog-associations_ontology-annotated-full.zip`).
`baselines` reads STRING v12 human links and protein names; `ontology` reads the
GO human annotations (`goa_human.gaf.gz`) and `go-basic.obo`.
`prospective` reads the same GO files, `consensus` and `understudied` STRING, `constraint`
the gnomAD v2.1.1 gene constraint table (`gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz`), and
`consensus` and `regulons` the CollecTRI pairs in `data/edges/`.
`perturbation`, `partners` and `hubs` read the pseudobulk files of the genome-scale
Perturb-seq screens of Replogle et al. 2022 (`perturb-seq-replogle-2022/`).
GWAS traits and DepMap lineages are mapped to tissues in `applications.py`,
fixed before any network was scored.

```sh
R=~/workspace/data/bioinf/references
uv run gene-network-analysis public $R local/public      # 16 graphs, about 1 minute
uv run gene-network-analysis predict local/public \
  --list essential=$R/gene-lists/CEGv2.txt \
  --list depmap-common=$R/depmap-24q2/CRISPRInferredCommonEssentials.csv \
  --list cancer=$R/gene-lists/oncokb-cancerGeneList.txt \
  --output results/public/prediction.json             # about 1 hour
uv run gene-network-analysis applications local/public $R --output results/public   # about 30 minutes
uv run gene-network-analysis cleaning local/public $R --output results/public       # about 20 minutes
uv run gene-network-analysis acyclic local/public $R --output results/public        # about 1 hour
uv run gene-network-analysis baselines local/public $R --output results/public \
  --list essential=$R/gene-lists/CEGv2.txt \
  --list depmap-common=$R/depmap-24q2/CRISPRInferredCommonEssentials.csv
uv run gene-network-analysis samples local/public $R --output results/public        # about 1 hour
uv run gene-network-analysis copy-number local/public $R --output results/public    # about 30 seconds
uv run gene-network-analysis dependencies local/public $R --output results/public
uv run gene-network-analysis ontology $R --output results/public   # reads uncharacterised.tsv
for check in robustness prospective consensus understudied constraint; do
  uv run gene-network-analysis $check local/public $R --output results/public   # 5–10 minutes each
done
uv run gene-network-analysis regulons local/public --output results/public     # about 1 minute
for check in perturbation partners hubs; do     # Perturb-seq and WGCNA hubs; 1–8 minutes each
  uv run gene-network-analysis $check local/public $R --output results/public
done
uv run gene-network-analysis normals $R local/public local/public-normal   # 8 adjacent-normal graphs
for check in robustness constraint; do
  uv run gene-network-analysis $check local/public-normal $R --output results/public/normal
done
for z in 2.0 3.09; do   # other cutoffs; 2.0 needs 4 GB of disk and 15 GB of memory
  uv run gene-network-analysis public $R local/public-z$z --threshold $z
  uv run gene-network-analysis robustness local/public-z$z $R --output results/public/threshold-$z
done
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
