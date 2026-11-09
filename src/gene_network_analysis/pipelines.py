from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import csr_array
from scipy.stats import rankdata

from .graph import component_summary, pagerank
from .neighbours import neighbour_vote, similarity_vote
from .steps import (
    absolute,
    drop_sinks,
    greedy_layering,
    keep_sign,
    minimum_feedback_order,
    orient_reciprocals,
    prune_triangles,
    resolve_reciprocals,
    top_outgoing,
)

ALPHA = 0.6


@dataclass
class State:
    graph: csr_array
    dag: csr_array | None = None
    level: np.ndarray | None = None
    features: dict[str, np.ndarray] = field(default_factory=dict)


def _negative_net_flow(graph: csr_array) -> np.ndarray:
    negative = graph.multiply(graph < 0).tocsr()
    weight = abs(negative)
    return np.asarray(weight.sum(axis=1)).ravel() - np.asarray(weight.sum(axis=0)).ravel()


def _pagerank_balance(graph: csr_array, alpha: float) -> np.ndarray:
    """log(reverse PageRank) - log(PageRank) on absolute weights: PageRank's net flow."""
    outgoing = pagerank(graph, alpha=alpha, reverse=True, weighted=True)[0]
    incoming = pagerank(graph, alpha=alpha, weighted=True)[0]
    return np.log(outgoing) - np.log(incoming)


# Enrichment: per-node metadata, computed on demand.
FEATURES: dict[str, Callable[[State], np.ndarray]] = {
    "out_degree": lambda s: np.diff(s.graph.indptr).astype(float),
    "weighted_out_degree": lambda s: np.asarray(abs(s.graph).sum(axis=1)).ravel(),
    "net_flow": lambda s: (
        np.asarray(abs(s.graph).sum(axis=1)).ravel() - np.asarray(abs(s.graph).sum(axis=0)).ravel()
    ),
    "positive_out_degree": lambda s: np.asarray((s.graph > 0).sum(axis=1)).ravel().astype(float),
    "negative_net_flow": lambda s: _negative_net_flow(s.graph),
    "sign_aware": lambda s: (
        rankdata(feature(s, "positive_out_degree")) + rankdata(feature(s, "negative_net_flow"))
    ),
    "pagerank": lambda s: pagerank(s.graph, alpha=ALPHA)[0],
    "dag_reverse_pagerank": lambda s: pagerank(s.dag, alpha=ALPHA, reverse=True)[0],
    "reverse_pagerank": lambda s: pagerank(s.graph, alpha=ALPHA, reverse=True)[0],
    "order_pagerank_rank_sum": lambda s: (
        rankdata(-s.level) + rankdata(feature(s, "dag_reverse_pagerank"))
    ),
    "pagerank_balance": lambda s: _pagerank_balance(s.graph, 0.85),
    "pagerank_balance_alpha_0.6": lambda s: _pagerank_balance(s.graph, 0.6),
    "pagerank_balance_alpha_0.95": lambda s: _pagerank_balance(s.graph, 0.95),
    "level": lambda s: s.level.astype(float),
}


def _similar_targets(state: State, labels: np.ndarray, universe: np.ndarray, signed: bool):
    values = np.zeros(len(labels))
    values[universe] = similarity_vote(
        state.graph, universe, labels[universe], targets=50, neighbours=20, signed=signed
    )
    return values


# Enrichment from the reference: labels of every node in, one value per node out.
LABEL_FEATURES: dict[str, Callable[[State, np.ndarray, np.ndarray], np.ndarray]] = {
    "tf_neighbour_share": lambda s, labels, _: neighbour_vote(s.graph, labels, "both"),
    "tf_target_share": lambda s, labels, _: neighbour_vote(s.graph, labels, "out"),
    "tf_sender_share": lambda s, labels, _: neighbour_vote(s.graph, labels, "in"),
    "similar_targets": lambda s, labels, universe: _similar_targets(s, labels, universe, False),
    "similar_signed_targets": lambda s, labels, universe: _similar_targets(
        s, labels, universe, True
    ),
}


def feature(state: State, name: str) -> np.ndarray:
    if name not in state.features:
        state.features[name] = FEATURES[name](state)
    return state.features[name]


@dataclass(frozen=True)
class Step:
    """A reusable graph step with its settings."""

    function: Callable
    params: tuple[tuple[str, object], ...] = ()

    def __call__(self, adjacency: csr_array):
        return self.function(adjacency, **dict(self.params))

    def describe(self) -> dict:
        return {"step": self.function.__name__, **dict(self.params)}


def step(function: Callable, **params) -> Step:
    return Step(function, tuple(sorted(params.items())))


@dataclass(frozen=True)
class Pipeline:
    """Configured steps ending in one score; higher_first says which end of it ranks first."""

    name: str
    title: str
    preprocess: tuple[Step, ...]
    decycle: Step | None
    score: str
    higher_first: bool
    baseline: bool = False
    reported: bool = True
    post_hoc: bool = False  # designed after seeing these results; needs independent data
    uses_labels: bool = False  # scored from other genes' labels in each reference


CLEAN_2016 = (step(resolve_reciprocals), step(drop_sinks), step(prune_triangles))
ABSOLUTE = step(absolute)
ACYCLIC = step(minimum_feedback_order)


def acyclic_variant(name: str, title: str, *preprocess: Step) -> Pipeline:
    return Pipeline(name, title, (*preprocess, ABSOLUTE), ACYCLIC, "level", False)


PIPELINES = (
    Pipeline("out-degree", "Most outgoing edges (raw graph)", (), None, "out_degree", True, True),
    Pipeline(
        "weighted-out-degree",
        "Largest total outgoing weight (raw graph)",
        (),
        None,
        "weighted_out_degree",
        True,
        True,
    ),
    Pipeline(
        "2016-pagerank", "2016 cleaning, lowest PageRank", CLEAN_2016, None, "pagerank", False
    ),
    Pipeline(
        "2016-layer",
        "2016 cleaning, earliest greedy layer",
        CLEAN_2016,
        step(greedy_layering),
        "level",
        False,
    ),
    acyclic_variant("acyclic-order", "Earliest in acyclic order"),
    Pipeline("net-flow", "Outgoing minus incoming weight", (ABSOLUTE,), None, "net_flow", True),
    Pipeline(
        "acyclic-dag-pagerank",
        "Reverse PageRank on the acyclic DAG",
        (ABSOLUTE,),
        ACYCLIC,
        "dag_reverse_pagerank",
        True,
    ),
    acyclic_variant(
        "oriented-order",
        "Pairs with a 2x stronger side only, earliest in acyclic order",
        step(orient_reciprocals, ratio=2.0),
    ),
    acyclic_variant(
        "top50-order",
        "Top 50 outgoing edges per gene, earliest in acyclic order",
        step(top_outgoing, k=50),
    ),
    acyclic_variant(
        "activation-order",
        "Activating edges only, earliest in acyclic order",
        step(keep_sign, sign=1),
    ),
    acyclic_variant(
        "repression-order",
        "Repressing edges only, earliest in acyclic order",
        step(keep_sign, sign=-1),
    ),
    Pipeline(
        "positive-out-degree",
        "Most positive outgoing edges",
        (),
        None,
        "positive_out_degree",
        True,
        post_hoc=True,
    ),
    Pipeline(
        "negative-net-flow",
        "Negative outgoing minus incoming weight",
        (),
        None,
        "negative_net_flow",
        True,
        post_hoc=True,
    ),
    Pipeline(
        "sign-aware",
        "Rank sum: positive out-degree and negative net flow",
        (),
        None,
        "sign_aware",
        True,
        post_hoc=True,
    ),
    # PageRank variants, fixed before running: reverse PageRank, and the balance of
    # reverse against forward PageRank (alpha 0.85; 0.6 and 0.95 as unreported checks).
    Pipeline(
        "reverse-pagerank",
        "Highest reverse PageRank (raw graph)",
        (),
        None,
        "reverse_pagerank",
        True,
        post_hoc=True,
    ),
    Pipeline(
        "pagerank-balance",
        "Reverse minus forward log PageRank, weighted, alpha 0.85",
        (),
        None,
        "pagerank_balance",
        True,
        post_hoc=True,
    ),
    Pipeline(
        "repression-pagerank-balance",
        "As pagerank-balance on repressing edges only",
        (step(keep_sign, sign=-1),),
        None,
        "pagerank_balance",
        True,
        post_hoc=True,
    ),
    # Repression order combined with PageRank, fixed before running.
    Pipeline(
        "repression-dag-pagerank",
        "Reverse PageRank on the repression-order DAG",
        (step(keep_sign, sign=-1), ABSOLUTE),
        ACYCLIC,
        "dag_reverse_pagerank",
        True,
        post_hoc=True,
    ),
    Pipeline(
        "repression-order-pagerank",
        "Rank sum: repression-order position and reverse PageRank on its DAG",
        (step(keep_sign, sign=-1), ABSOLUTE),
        ACYCLIC,
        "order_pagerank_rank_sum",
        True,
        post_hoc=True,
    ),
    *(
        Pipeline(
            f"pagerank-balance-alpha-{alpha}",
            f"As pagerank-balance with alpha {alpha}",
            (),
            None,
            f"pagerank_balance_alpha_{alpha}",
            True,
            reported=False,
            post_hoc=True,
        )
        for alpha in ("0.6", "0.95")
    ),
    *(
        Pipeline(name, title, (), None, score, True, uses_labels=True, reported=reported)
        for name, title, score, reported in (
            (
                "tf-neighbour-share",
                "Share of edge weight to or from TFs",
                "tf_neighbour_share",
                True,
            ),
            ("tf-target-share", "Share of outgoing weight to TFs", "tf_target_share", True),
            ("tf-sender-share", "Share of incoming weight from TFs", "tf_sender_share", True),
            (
                "similar-targets",
                "TF share among the 20 genes with the most similar top-50 targets",
                "similar_targets",
                True,
            ),
            (
                "similar-signed-targets",
                "As similar-targets, comparing target signs too",
                "similar_signed_targets",
                False,
            ),
        )
    ),
)


def graph_stats(matrix: csr_array) -> dict:
    components = component_summary(matrix)
    return {
        "edges": int(matrix.nnz),
        "negative_edges": int((matrix.data < 0).sum()),
        "nodes_with_edges": components["active_nodes"],
        "largest_strong_component": components["largest_strong_component"],
    }


class Runner:
    """Runs pipelines, reusing graphs and features of shared leading steps."""

    def __init__(self, raw: csr_array):
        self.states: dict[tuple, State] = {(): State(raw)}

    def state(self, steps: tuple[Step, ...], decycle: Step | None = None) -> State:
        key = steps + ((decycle,) if decycle else ())
        if key not in self.states:
            if decycle:
                parent = self.state(steps)
                dag, level = decycle(parent.graph)
                self.states[key] = State(parent.graph, dag, level)
            else:
                self.states[key] = State(steps[-1](self.state(steps[:-1]).graph))
        return self.states[key]

    def record(self, pipeline: Pipeline) -> dict:
        state = self.state(pipeline.preprocess, pipeline.decycle)
        return {
            "name": pipeline.name,
            "title": pipeline.title,
            "baseline": pipeline.baseline,
            "post_hoc": pipeline.post_hoc,
            "uses_labels": pipeline.uses_labels,
            "steps": {
                "preprocess": [s.describe() for s in pipeline.preprocess],
                "decycle": pipeline.decycle.describe() if pipeline.decycle else None,
                "score": pipeline.score,
                "order": "descending" if pipeline.higher_first else "ascending",
            },
            "graph": graph_stats(state.graph),
            **({"dag": graph_stats(state.dag)} if state.dag is not None else {}),
        }

    def score(
        self, pipeline: Pipeline, universe: np.ndarray, labels: np.ndarray | None = None
    ) -> np.ndarray:
        """Scores of the universe genes, highest first; labels cover every node."""
        state = self.state(pipeline.preprocess, pipeline.decycle)
        if pipeline.uses_labels:
            if labels is None:
                raise ValueError(f"{pipeline.name} needs reference labels")
            values = LABEL_FEATURES[pipeline.score](state, labels, universe)
        else:
            values = feature(state, pipeline.score)
        # Selection: the shared universe, ordered so the first-ranked gene scores highest.
        return (values if pipeline.higher_first else -values)[universe]
