"""Given a known TF, do its network neighbours recover its curated targets?

Network shape does not tell which genes are regulators. This asks the easier question:
for a TF that is already known, with its targets curated from the literature (CollecTRI),
are those targets over-represented among the genes linked to it?

For every TF with at least MIN_TARGETS curated targets among the network's genes, and in
every network, the overlap between its neighbours (linked in either direction) and its
targets is tested against drawing the neighbours at random (hypergeometric), with
Benjamini-Hochberg correction over TFs inside a network.

Two checks on what a significant overlap means:

- own against others: targets of many TFs are popular, well-connected genes, so each
  TF's enrichment for its own targets is ranked among its enrichments for every other
  TF's targets (0.5 if a TF's neighbours are no closer to its own targets);
- visibility: the rank correlation between a TF's enrichment and the spread of its own
  expression in that dataset, since a TF whose mRNA does not vary cannot co-vary with
  anything.

MIN_TARGETS, the test and both checks were fixed before the first run. A third check was
added after it and fixed before it was computed:

- tissue match: a TF is elevated in a tissue when its mean expression there (GTEx,
  healthy) is at least ELEVATED log2 TPM above its median in the other six tissues. If
  the network sees a regulator where it is at work, such a TF should recover its targets
  better in the tumour networks of its tissues than in the other tumour networks. Tumour
  data play no part in calling a TF elevated. The same comparison inside the healthy
  networks is given for reference. The measure is log2((shared + 1) / (expected + 1)).
"""

import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.stats import false_discovery_control, hypergeom, spearmanr

from .applications import load
from .references import read_edges

MIN_TARGETS = 10
FDR = 0.05
LISTED = 15
ELEVATED = 1.0


def regulons(
    pairs: list[tuple[str, str, int]], genes: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """TF indices, TFs x genes target membership and TFs x genes curated sign."""
    index = {g: i for i, g in enumerate(genes)}
    targets: dict[int, dict[int, int]] = {}
    for source, target, sign in pairs:
        if source in index and target in index and source != target:
            targets.setdefault(index[source], {})[index[target]] = sign
    chosen = sorted(tf for tf, found in targets.items() if len(found) >= MIN_TARGETS)
    member = np.zeros((len(chosen), len(genes)), dtype=bool)
    signs = np.zeros((len(chosen), len(genes)), dtype=np.int8)
    for row, tf in enumerate(chosen):
        member[row, list(targets[tf])] = True
        signs[row, list(targets[tf])] = list(targets[tf].values())
    return np.array(chosen), member, signs


def enrichment(links: csr_array, tfs: np.ndarray, member: np.ndarray) -> dict[str, np.ndarray]:
    """Per TF: neighbours, targets among them, fold over chance, p-value and the rank of
    its own fold among its folds for the other TFs' targets."""
    neighbours = (links[tfs] != 0).astype(np.float64)
    shared = np.asarray(neighbours @ member.T.astype(np.float64))
    degree = np.asarray(neighbours.sum(axis=1)).ravel()
    sizes = member.sum(axis=1)
    others = links.shape[0] - 1
    expected = degree[:, None] * sizes[None, :] / others
    fold = np.divide(shared, expected, out=np.zeros_like(shared), where=expected > 0)
    own = np.diag(shared)
    below = (fold < np.diag(fold)[:, None]).sum(axis=1)
    ties = (fold == np.diag(fold)[:, None]).sum(axis=1) - 1
    return {
        "neighbours": degree,
        "shared": own,
        "expected": np.diag(expected),
        "fold": np.diag(fold),
        "p": hypergeom.sf(own - 1, others, sizes, degree),
        "own_rank": (below + 0.5 * ties) / (len(tfs) - 1),
    }


def sign_agreement(adjacency: csr_array, tfs: np.ndarray, signs: np.ndarray) -> dict:
    """Among linked curated pairs with a known sign: does the link's sign match?"""
    both = (adjacency + adjacency.T).tocsr()[tfs].toarray()
    judged = (both != 0) & (signs != 0)
    if not judged.any():
        return {"pairs": 0}
    link, curated = np.sign(both[judged]), signs[judged]
    positive, activation = (link > 0).mean(), (curated > 0).mean()
    return {
        "pairs": int(judged.sum()),
        "agreement": float((link == curated).mean()),
        "chance": float(positive * activation + (1 - positive) * (1 - activation)),
    }


def elevated(means: np.ndarray) -> np.ndarray:
    """TFs x tissues: the TF's mean is at least ELEVATED above its median elsewhere."""
    return np.column_stack(
        [
            means[:, j] - np.median(np.delete(means, j, axis=1), axis=1) >= ELEVATED
            for j in range(means.shape[1])
        ]
    )


def tissue_match(
    strength: np.ndarray,
    significant: np.ndarray,
    tissue: np.ndarray,
    high: np.ndarray,
    runs: int = 2000,
    seed: int = 0,
) -> dict:
    """Networks x TFs `strength` and `significant`, each network's tissue index and the
    TFs x tissues `high`: every TF's networks in its own tissues against its others."""
    own = high[:, tissue].T
    judged = own.any(axis=0) & ~own.all(axis=0)
    own, strength, significant = own[:, judged], strength[:, judged], significant[:, judged]

    def means(values: np.ndarray, chosen: np.ndarray) -> np.ndarray:
        return (values * chosen).sum(axis=0) / chosen.sum(axis=0)

    gain = means(strength, own) - means(strength, ~own)
    rng = np.random.default_rng(seed)
    draws = gain[rng.integers(0, len(gain), (runs, len(gain)))].mean(axis=1)
    return {
        "tfs": int(judged.sum()),
        "mean_difference": float(gain.mean()),
        "ci95": [float(q) for q in np.quantile(draws, [0.025, 0.975])],
        "tfs_higher_in_own_tissues": int((gain > 0).sum()),
        "significant_in_own_tissues": float(means(significant, own).mean()),
        "significant_elsewhere": float(means(significant, ~own).mean()),
        "judged": judged,
        "gain": gain,
    }


def run(root: Path, gene_table: Path, data: Path, output: Path) -> dict:
    context = load(root, gene_table)
    genes = context.genes
    tfs, member, signs = regulons(read_edges(data / "edges/collectri.tsv"), genes)
    symbols = [context.table[genes[i]]["symbol"] for i in tfs]
    report: dict = {
        "tfs": len(tfs),
        "min_targets": MIN_TARGETS,
        "median_targets": float(np.median(member.sum(axis=1))),
        "networks": {},
    }
    significant = np.zeros((len(context.plain), len(tfs)), dtype=bool)
    strength = np.zeros(significant.shape)
    for row, name in enumerate(context.plain):
        found = enrichment(context.links[name], tfs, member)
        strength[row] = np.log2((found["shared"] + 1) / (found["expected"] + 1))
        testable = found["neighbours"] > 0
        q = np.ones(len(tfs))
        q[testable] = false_discovery_control(found["p"][testable])
        significant[row] = q <= FDR
        spread = context.expression[f"{name}_sd"][tfs]
        order = np.argsort(found["p"], kind="stable")[:LISTED]
        chosen = significant[row]
        report["networks"][name] = {
            "tfs_with_links": int(testable.sum()),
            "significant": int(chosen.sum()),
            "median_fold": float(np.median(found["fold"][testable])),
            "median_fold_significant": float(np.median(found["fold"][chosen]))
            if chosen.any()
            else None,
            "mean_own_rank": float(found["own_rank"][testable].mean()),
            "fold_against_expression_spread": float(
                spearmanr(found["fold"][testable], spread[testable]).statistic
            ),
            "sign_agreement": sign_agreement(context.adjacency[name], tfs[chosen], signs[chosen]),
            "top": [
                {
                    "tf": symbols[i],
                    "targets": int(member[i].sum()),
                    "neighbours": int(found["neighbours"][i]),
                    "shared": int(found["shared"][i]),
                    "fold": float(found["fold"][i]),
                    "q": float(q[i]),
                }
                for i in order
            ],
        }
        print(
            name,
            json.dumps({k: v for k, v in report["networks"][name].items() if k != "top"}),
            " ".join(t["tf"] for t in report["networks"][name]["top"]),
            flush=True,
        )
    count = significant.sum(axis=0)
    report["across_networks"] = {
        "significant_in_any": int((count > 0).sum()),
        "significant_in_half": int((count >= (len(context.plain) + 1) // 2).sum()),
        "significant_in_all": int((count == len(context.plain)).sum()),
        "most_often": [
            {"tf": symbols[i], "networks": int(count[i])}
            for i in np.argsort(-count, kind="stable")[:40]
            if count[i]
        ],
    }
    print(json.dumps(report["across_networks"]), flush=True)
    state = np.array([context.stats[n]["state"] for n in context.plain])
    tissues = sorted({context.stats[n]["tissue"] for n in context.plain})
    tissue = np.array([tissues.index(context.stats[n]["tissue"]) for n in context.plain])
    high = elevated(np.column_stack([context.expression[f"gtex-{t}_mean"][tfs] for t in tissues]))
    report["tissue_match"] = {
        "tissues": tissues,
        "elevated_tfs_per_tissue": {t: int(high[:, j].sum()) for j, t in enumerate(tissues)},
    }
    for name in ("tumour", "healthy"):
        chosen = state == name
        found = tissue_match(strength[chosen], significant[chosen], tissue[chosen], high)
        judged, gain = found.pop("judged"), found.pop("gain")
        order = np.flatnonzero(judged)[np.argsort(-gain, kind="stable")]
        found["largest"] = [
            {
                "tf": symbols[i],
                "tissues": [t for j, t in enumerate(tissues) if high[i, j]],
                "difference": float(d),
            }
            for i, d in zip(order[:LISTED], np.sort(gain)[::-1][:LISTED], strict=True)
        ]
        report["tissue_match"][name] = found
        print("tissue match", name, json.dumps(found), flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output / "regulons.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
