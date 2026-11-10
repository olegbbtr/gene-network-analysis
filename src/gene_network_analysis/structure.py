"""What a graph's structure reflects: genome position, sample make-up, correlation.

- `locality`: share of edges joining genes on the same chromosome arm or band, against
  the same share after shuffling gene locations.
- `modules`: marker gene sets fixed in advance, and how many top out-degree genes are
  members of a module or send most of their edges into one.
- `balance`: share of signed triangles that are balanced (an even number of negative
  edges), against the share expected if the same signs were placed at random.
"""

import re

import numpy as np
from scipy.sparse import csr_array

from .graph import Graph
from .references import read_edges

ARM = re.compile(r"^(\d+|X|Y)([pq])")
BAND = re.compile(r"^(\d+|X|Y)([pq]\d+)")

MODULES = {
    "sex": ["XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"],
    "immune": ["PTPRC", "CD2", "CD3E", "CD247", "LCK", "CD48", "CCL5", "CD53", "LCP2", "CORO1A"],
    "immediate_early": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "NR4A1", "IER2", "DUSP1"],
    "cell_cycle": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "CENPF", "KIF11", "TPX2"],
    "interferon": ["ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OAS1", "RSAD2", "IFI44L", "STAT1"],
}
MODULE_ENRICHMENT = 5.0
MODULE_MIN_EDGES = 3


def location_codes(nodes: tuple[str, ...], genes: dict[str, dict], pattern: re.Pattern):
    """An integer per node for its arm or band; -1 where HGNC gives no usable location."""
    labels = []
    for node in nodes:
        match = pattern.match(genes.get(node, {}).get("location", ""))
        labels.append("".join(match.groups()) if match else "")
    names = sorted(set(labels) - {""})
    index = {name: i for i, name in enumerate(names)}
    return np.array([index.get(label, -1) for label in labels])


def locality(
    graph: Graph, genes: dict[str, dict], permutations: int = 100, seed: int = 0, top_k: int = 100
) -> dict:
    edges = graph.adjacency.tocoo()
    rng = np.random.default_rng(seed)
    result = {}
    for level, pattern in (("arm", ARM), ("band", BAND)):
        codes = location_codes(graph.nodes, genes, pattern)
        known = np.flatnonzero(codes >= 0)
        located = (codes[edges.row] >= 0) & (codes[edges.col] >= 0)
        rows, columns, weights = edges.row[located], edges.col[located], edges.data[located]
        subsets = {"all": slice(None), "positive": weights > 0, "negative": weights < 0}
        shuffled = []
        for _ in range(permutations):
            permuted = codes.copy()
            permuted[known] = rng.permutation(codes[known])
            shuffled.append(permuted)
        level_result = {}
        for name, subset in subsets.items():
            r, c = rows[subset], columns[subset]
            observed = float((codes[r] == codes[c]).mean())
            null = np.array([(p[r] == p[c]).mean() for p in shuffled])
            level_result[name] = {
                "edges": int(len(r)),
                "same_share": observed,
                "shuffled_mean": float(null.mean()),
                "shuffled_95": [float(np.quantile(null, q)) for q in (0.025, 0.975)],
                "fold": observed / float(null.mean()) if null.mean() else None,
                "p_value": float((1 + (null >= observed).sum()) / (1 + permutations)),
            }
        same = (codes[rows] == codes[columns]).astype(float)
        degree = np.bincount(rows, minlength=len(codes))
        local = np.bincount(rows, weights=same, minlength=len(codes))
        senders = np.flatnonzero(degree > 0)
        share = local[senders] / degree[senders]
        top = senders[np.argsort(-degree[senders], kind="stable")[:top_k]]
        level_result["sender_local_share_mean"] = float(share.mean())
        level_result[f"top_{top_k}_out_degree_local_share_mean"] = float(
            (local[top] / degree[top]).mean()
        )
        result[level] = level_result
    return result


def modules(graph: Graph, genes: dict[str, dict], top_k: int = 100) -> dict:
    symbols = [genes.get(node, {}).get("symbol", node) for node in graph.nodes]
    index = {symbol: i for i, symbol in enumerate(symbols)}
    adjacency = graph.adjacency.tocsr()
    present = (adjacency != 0).astype(np.float64)
    nodes = len(graph.nodes)
    degree = np.asarray(present.sum(axis=1)).ravel()
    density = present.nnz / (nodes * (nodes - 1))
    members = {
        name: np.array(sorted(index[g] for g in markers if g in index), dtype=int)
        for name, markers in MODULES.items()
    }
    enrichment = np.zeros((nodes, len(MODULES)))
    counts = np.zeros((nodes, len(MODULES)))
    result: dict = {"modules": {}}
    for column, (name, member) in enumerate(members.items()):
        if len(member) < 2:
            result["modules"][name] = {"present": [symbols[i] for i in member]}
            continue
        counts[:, column] = np.asarray(present[:, member].sum(axis=1)).ravel()
        share = np.divide(counts[:, column], degree, out=np.zeros(nodes), where=degree > 0)
        enrichment[:, column] = share / (len(member) / nodes)
        block = adjacency[member][:, member]
        result["modules"][name] = {
            "present": [symbols[i] for i in member],
            "within_density_fold": float(block.nnz / (len(member) * (len(member) - 1)) / density),
            "within_positive_share": float((block.data > 0).mean()) if block.nnz else None,
        }
    assigned = np.full(nodes, "", dtype=object)
    names = list(MODULES)
    strongest = enrichment.argmax(axis=1)
    for node in range(nodes):
        column = strongest[node]
        if (
            enrichment[node, column] >= MODULE_ENRICHMENT
            and counts[node, column] >= MODULE_MIN_EDGES
        ):
            assigned[node] = names[column]
    for name, member in members.items():
        assigned[member] = name
    top = np.argsort(-degree, kind="stable")[:top_k]
    result[f"top_{top_k}_out_degree"] = {
        name: int((assigned[top] == name).sum()) for name in names
    } | {"none": int((assigned[top] == "").sum())}
    result["top_20"] = [
        {"symbol": symbols[i], "out_degree": int(degree[i]), "module": assigned[i] or None}
        for i in top[:20]
    ]
    sex = members["sex"]
    xist = index.get("XIST")
    if xist is not None and len(sex) > 1:
        y_genes = sex[sex != xist]
        weights = np.concatenate(
            [adjacency[[xist]][:, y_genes].data, adjacency[y_genes][:, [xist]].data]
        )
        result["xist_y_edges"] = {
            "edges": int(len(weights)),
            "possible": int(2 * len(y_genes)),
            "negative": int((weights < 0).sum()),
        }
    result["settings"] = {"enrichment": MODULE_ENRICHMENT, "min_edges": MODULE_MIN_EDGES}
    return result


def undirected_signs(adjacency: csr_array) -> csr_array:
    """Each linked pair once: the sign of the summed weights in both directions."""
    total = (adjacency + adjacency.T).tocsr()
    total.setdiag(0)
    total.eliminate_zeros()
    total.data = np.sign(total.data)
    return total


def balance(signs: csr_array) -> dict:
    dense = signs.toarray().astype(np.float32)
    magnitude = np.abs(dense)
    signed = float(((dense @ dense) * dense).sum(dtype=np.float64))
    total = float(((magnitude @ magnitude) * magnitude).sum(dtype=np.float64))
    positive = float((signs.data > 0).mean())
    negative = 1 - positive
    return {
        "nodes": int(signs.shape[0]),
        "edges": int(signs.nnz // 2),
        "positive_share": positive,
        "triangles": int(round(total / 6)),
        "balanced_share": (total + signed) / (2 * total) if total else None,
        "random_sign_balanced_share": positive**3 + 3 * positive * negative**2,
    }


def signed_pairs(path) -> csr_array:
    """Curated pairs with a known sign as a symmetric sign matrix."""
    pairs = [(s, t, sign) for s, t, sign in read_edges(path) if sign != 0 and s != t]
    names = sorted({g for s, t, _ in pairs for g in (s, t)})
    index = {g: i for i, g in enumerate(names)}
    rows = [index[s] for s, _, _ in pairs]
    columns = [index[t] for _, t, _ in pairs]
    matrix = csr_array(
        (np.array([sign for _, _, sign in pairs], dtype=float), (rows, columns)),
        shape=(len(names), len(names)),
    )
    return undirected_signs(matrix)
