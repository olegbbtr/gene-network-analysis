import csv
import datetime
import hashlib
import importlib.metadata
import json
import subprocess
from pathlib import Path

import numpy as np

from .edges import edge_metrics
from .evaluation import (
    bootstrap_aucs,
    normalize_gene_id,
    p_value,
    random_rankings,
    ranking_metrics,
    read_reference,
    resolve_labels,
)
from .graph import read_graph
from .pipelines import PIPELINES, Runner
from .references import read_edges


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def revision() -> dict:
    root = Path(__file__).resolve().parents[2]
    head = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True)
    diff = subprocess.run(
        ["git", "-C", root, "status", "--porcelain", "--", "src", "pyproject.toml", "uv.lock"],
        capture_output=True,
        text=True,
    )
    return {
        "git_commit": head.stdout.strip() or None,
        "source_modified": bool(diff.stdout) if diff.returncode == 0 else None,
    }


def compare(
    records: list[dict], scores: dict, reference: str, labels: np.ndarray, top_k: int
) -> list[dict]:
    """Random row, then every pipeline's metrics against random and the best baseline.

    The AUC difference from the best baseline gets a 95% interval and a two-sided p-value
    from a paired bootstrap over genes.
    """
    null = random_rankings(labels, top_k)
    resampled = bootstrap_aucs(scores, labels)
    rows = [
        {
            "pipeline": "random",
            "title": f"Random order (mean of {null['runs']:,}, seed {null['seed']})",
            "baseline": True,
            "top_k_reference_genes": float(null["top_k_reference_genes"].mean()),
            "top_k_95_range": np.percentile(null["top_k_reference_genes"], [2.5, 97.5]).tolist(),
            "roc_auc": float(null["roc_auc"].mean()),
            "roc_auc_95_range": np.percentile(null["roc_auc"], [2.5, 97.5]).tolist(),
        }
    ]
    best = max(
        (r for r in records if r["baseline"]),
        key=lambda record: record["metrics"][reference]["roc_auc"],
        default=None,
    )
    for r in records:
        metrics = r["metrics"][reference]
        difference = resampled[r["name"]] - resampled[best["name"]] if best else None
        rows.append(
            {
                "pipeline": r["name"],
                "title": r["title"],
                "baseline": r["baseline"],
                "post_hoc": r["post_hoc"],
                "uses_labels": r["uses_labels"],
                "top_k_reference_genes": metrics["top_k_reference_genes"],
                "roc_auc": metrics["roc_auc"],
                "p_top_k_vs_random": p_value(
                    null["top_k_reference_genes"], metrics["top_k_reference_genes"]
                ),
                "p_auc_vs_random": p_value(null["roc_auc"], metrics["roc_auc"]),
                "auc_minus_best_baseline": (
                    metrics["roc_auc"] - best["metrics"][reference]["roc_auc"] if best else None
                ),
                "auc_difference_95_interval": (
                    np.percentile(difference, [2.5, 97.5]).tolist() if best else None
                ),
                "p_auc_vs_best_baseline": (
                    float(min(1.0, 2 * min((difference <= 0).mean(), (difference >= 0).mean())))
                    if best and r["name"] != best["name"]
                    else None
                ),
            }
        )
    return rows


def chain_text(steps: tuple, decycle) -> str:
    parts = [step.describe() for step in steps] + ([decycle.describe()] if decycle else [])
    return (
        " → ".join(
            entry["step"] + "".join(f" {k}={v}" for k, v in entry.items() if k != "step")
            for entry in parts
        )
        or "raw graph"
    )


def check_edges(runner: Runner, pipelines: list, nodes: tuple, references: dict) -> list[dict]:
    """Edge metrics for every distinct graph and DAG the pipelines build."""
    index = {normalize_gene_id(node): i for i, node in enumerate(nodes)}
    pairs = {}
    for name, path in references.items():
        known = [
            (index[normalize_gene_id(s)], index[normalize_gene_id(t)], sign)
            for s, t, sign in read_edges(path)
            if normalize_gene_id(s) in index and normalize_gene_id(t) in index
        ]
        pairs[name] = tuple(np.array(column, dtype=np.int64) for column in zip(*known, strict=True))
    checks: dict[tuple, dict] = {}
    for pipeline in pipelines:
        for decycle in (None, pipeline.decycle) if pipeline.decycle else (None,):
            key = (pipeline.preprocess, decycle)
            if key not in checks:
                state = runner.state(*key)
                matrix = state.dag if decycle else state.graph
                checks[key] = {
                    "graph": chain_text(*key),
                    "pipelines": [],
                    "references": {name: edge_metrics(matrix, *p) for name, p in pairs.items()},
                }
            checks[key]["pipelines"].append(pipeline.name)
    return list(checks.values())


def analyze(
    graph_path: Path,
    references: dict[str, Path],
    output: Path,
    *,
    edge_references: dict[str, Path] | None = None,
    pipelines: list[str] | None = None,
    top_k: int = 100,
) -> dict:
    graph = read_graph(graph_path)
    raw = graph.adjacency
    nodes = list(graph.nodes)
    _, known = resolve_labels(graph.nodes, set())
    # Result set: genes with an outgoing edge in the raw graph and an Ensembl identity.
    universe = np.flatnonzero((np.diff(raw.indptr) > 0) & known)
    universe_nodes = [nodes[i] for i in universe]
    selected = [p for p in PIPELINES if (p.name in pipelines if pipelines else p.reported)]
    runner = Runner(raw)
    records = [dict(runner.record(p), metrics={}) for p in selected]
    shared = {p.name: runner.score(p, universe) for p in selected if not p.uses_labels}
    edge_checks = check_edges(runner, selected, graph.nodes, edge_references or {})
    evaluations, labels, scores = {}, {}, {}
    for name, path in references.items():
        every_node = resolve_labels(graph.nodes, read_reference(path))[0]
        labels[name] = every_node[universe]
        scores[name] = {
            p.name: shared[p.name] if p.name in shared else runner.score(p, universe, every_node)
            for p in selected
        }
        for record in records:
            record["metrics"][name] = ranking_metrics(
                scores[name][record["name"]], labels[name], universe_nodes, top_k
            )
        evaluations[name] = {
            "path": path.name,
            "sha256": sha256(path),
            "reference_genes": int(labels[name].sum()),
            "comparison": compare(records, scores[name], name, labels[name], top_k),
        }
    report = {
        "date": datetime.date.today().isoformat(),
        "implementation": revision(),
        "versions": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "matplotlib")},
        "input": {
            "graph": {"path": graph_path.name, "sha256": sha256(graph_path)},
            "nodes": len(nodes),
            "edges": int(raw.nnz),
        },
        "selection": {
            "universe": "genes with at least one outgoing edge in the raw graph",
            "genes": len(universe),
            "top_k": top_k,
        },
        "references": evaluations,
        "pipelines": records,
        "edge_references": {
            name: {"path": path.name, "sha256": sha256(path)}
            for name, path in (edge_references or {}).items()
        },
        "edge_checks": edge_checks,
        "limitations": [
            "A reference marks known transcription factors; absence is not a negative label.",
            "The graph's experimental or simulation origin is unverified.",
            "Pipelines that use labels score each gene from other genes' labels in the same "
            "reference, never its own; they are not label-free predictions.",
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    columns = {p.name: shared[p.name] for p in selected if p.name in shared}
    for name in references:
        for p in selected:
            if p.uses_labels:
                columns[f"{p.name}[{name}]"] = scores[name][p.name]
    with (output / "rankings.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["gene_id", *(f"tf_{name}" for name in labels), *columns])
        for row, gene in enumerate(universe_nodes):
            writer.writerow(
                [
                    gene,
                    *(int(values[row]) for values in labels.values()),
                    *(float(values[row]) for values in columns.values()),
                ]
            )
    write_pdf(report, output / "report.pdf")
    return report


def p_text(value: float) -> str:
    return "< 0.001" if value < 0.001 else f"{value:.3f}"


def describe_steps(steps: dict) -> str:
    def one(entry: dict) -> str:
        params = ", ".join(f"{k}={v}" for k, v in entry.items() if k != "step")
        return f"{entry['step']}({params})"

    chain = [one(s) for s in steps["preprocess"]]
    if steps["decycle"]:
        chain.append(one(steps["decycle"]))
    chain.append(f"score {steps['score']}, {steps['order']}")
    return " → ".join(chain)


def comparison_page(plt, report: dict, name: str):
    evaluation = report["references"][name]
    rows = evaluation["comparison"]
    random_row = rows[0]
    top_k = report["selection"]["top_k"]
    figure = plt.figure(figsize=(8.27, 11.69))
    figure.text(0.06, 0.95, f"Regulator ranking against reference '{name}'", size=15)
    figure.text(
        0.06,
        0.925,
        f"{report['date']} · {report['selection']['genes']:,} genes ranked, "
        f"{evaluation['reference_genes']:,} in {evaluation['path']} · top k = {top_k}",
        size=9,
    )
    cells = [
        [
            "random",
            "{:.1f} ({:.0f}–{:.0f})".format(
                random_row["top_k_reference_genes"], *random_row["top_k_95_range"]
            ),
            "{:.3f} ({:.3f}–{:.3f})".format(random_row["roc_auc"], *random_row["roc_auc_95_range"]),
            "",
            "",
            "",
            "",
        ]
    ]
    for r in rows[1:]:
        delta = r["auc_minus_best_baseline"]
        cells.append(
            [
                r["pipeline"]
                + (" †" if r["post_hoc"] else "")
                + (" ‡" if r["uses_labels"] else ""),
                str(r["top_k_reference_genes"]),
                f"{r['roc_auc']:.3f}",
                p_text(r["p_top_k_vs_random"]),
                p_text(r["p_auc_vs_random"]),
                ""
                if delta is None
                else "{:+.3f} ({:+.3f}, {:+.3f})".format(delta, *r["auc_difference_95_interval"]),
                "" if r["p_auc_vs_best_baseline"] is None else p_text(r["p_auc_vs_best_baseline"]),
            ]
        )
    table = figure.add_axes((0.03, 0.56, 0.94, 0.35))
    table.axis("off")
    grid = table.table(
        cellText=cells,
        colLabels=[
            "Pipeline",
            f"Reference genes in top {top_k}",
            "AUC",
            "p top-k",
            "p AUC",
            "AUC vs best baseline (95%)",
            "p vs baseline",
        ],
        loc="upper center",
        cellLoc="left",
    )
    grid.auto_set_font_size(False)
    grid.set_fontsize(6.5)
    grid.auto_set_column_width(list(range(7)))
    # Rows shrink as pipelines are added so the table stays above the footnote.
    grid.scale(1, min(1.15, 24 / len(cells)))
    figure.text(
        0.06,
        0.455,
        "Random: mean and 95% range over random orderings. p: share of random orderings "
        "scoring at least as well. Grey bars: random and baselines.\n"
        "AUC vs best baseline: paired bootstrap over genes (1,000 resamples). "
        "†: designed after seeing these results; needs confirming on other data.\n"
        "‡: uses known TFs: scored from other genes' labels in this reference (lighter bars).",
        size=7,
    )
    chart = figure.add_axes((0.3, 0.05, 0.65, 0.36))
    chart.barh(
        [r["pipeline"] for r in rows][::-1],
        [r["roc_auc"] for r in rows][::-1],
        color=[
            "#8a8a8a" if r["baseline"] else "#8cc3cf" if r.get("uses_labels") else "#20788c"
            for r in rows
        ][::-1],
    )
    chart.axvspan(*random_row["roc_auc_95_range"], color="#b06b29", alpha=0.15, label="random 95%")
    chart.set(xlabel="ROC AUC", xlim=(0.4, 0.75))
    chart.tick_params(axis="y", labelsize=7)
    chart.legend(loc="lower right")
    return figure


def steps_page(plt, report: dict):
    figure = plt.figure(figsize=(8.27, 11.69))
    figure.text(0.06, 0.95, "Pipeline steps", size=15)
    y = 0.92
    for record in report["pipelines"]:
        graph = record["graph"]
        lines = [
            describe_steps(record["steps"]),
            f"Graph: {graph['edges']:,} edges ({graph['negative_edges']:,} negative), "
            f"largest strong component {graph['largest_strong_component']:,}"
            + (f"; DAG keeps {record['dag']['edges']:,} edges." if "dag" in record else "."),
        ]
        note = " (uses known TFs)" if record["uses_labels"] else ""
        figure.text(0.06, y, f"{record['name']}: {record['title']}{note}", size=9, weight="bold")
        for line in lines:
            y -= 0.018
            figure.text(0.08, y, line, size=7)
        y -= 0.03
    return figure


def write_pdf(report: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        pages = [comparison_page(plt, report, name) for name in report["references"]]
        if report["edge_checks"]:
            pages.append(edges_page(plt, report))
        for figure in [*pages, steps_page(plt, report)]:
            pdf.savefig(figure)
            plt.close(figure)


def edges_page(plt, report: dict):
    figure = plt.figure(figsize=(8.27, 11.69))
    figure.text(0.06, 0.95, "Edge checks against curated TF→target pairs", size=15)
    figure.text(
        0.06,
        0.915,
        "Recall ×: share of curated pairs present as edges, over chance.\n"
        "Positive/negative ×: share of edges of that sign that are curated pairs, "
        "over chance. Sign: agreement with the curated sign (chance).",
        size=7,
    )
    names = list(report["edge_references"])
    height = 0.85 / max(len(names), 1)
    for position, name in enumerate(names):
        cells = []
        for check in report["edge_checks"]:
            m = check["references"][name]
            if not m.get("reference_pairs"):
                continue
            by_sign = m["precision_by_graph_sign"] or {}
            agreement = m["sign_agreement"]
            cells.append(
                [
                    check["graph"],
                    f"{m['reference_pairs']:,}",
                    f"{m['recall']:.3f} ({m['recall_enrichment']:.2f}×)",
                    *(
                        f"{by_sign[k]['enrichment']:.2f}×"
                        if by_sign.get(k) and by_sign[k]["edges"]
                        else ""
                        for k in ("positive", "negative")
                    ),
                    f"{agreement['agreement']:.2f} ({agreement['chance']:.2f})"
                    if agreement
                    else "",
                ]
            )
        top = 0.9 - position * height
        figure.text(0.06, top, name, size=10, weight="bold")
        axes = figure.add_axes((0.03, top - height + 0.02, 0.94, height - 0.04))
        axes.axis("off")
        grid = axes.table(
            cellText=cells,
            colLabels=["Graph", "Pairs", "Recall (×)", "Positive ×", "Negative ×", "Sign (chance)"],
            loc="upper center",
            cellLoc="left",
        )
        grid.auto_set_font_size(False)
        grid.set_fontsize(5.5)
        grid.auto_set_column_width(list(range(6)))
    return figure
