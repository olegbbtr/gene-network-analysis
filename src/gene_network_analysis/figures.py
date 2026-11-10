"""Poster figures drawn from the committed results in results/public/."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, BAND = "#0b0b0b", "#52514e", "#dddcd7"

# Tumour and healthy networks of one tissue sit together.
NETWORKS = [
    ("tcga-luad", "Lung adenocarcinoma"),
    ("tcga-luad-copy-number-removed", "Lung adenocarcinoma, CN removed"),
    ("tcga-lusc", "Lung squamous carcinoma"),
    ("gtex-lung", "Lung, healthy"),
    ("tcga-brca", "Breast carcinoma"),
    ("gtex-breast", "Breast, healthy"),
    ("tcga-coad", "Colon adenocarcinoma"),
    ("gtex-colon", "Colon, healthy"),
    ("tcga-thca", "Thyroid carcinoma"),
    ("gtex-thyroid", "Thyroid, healthy"),
    ("tcga-prad", "Prostate adenocarcinoma"),
    ("gtex-prostate", "Prostate, healthy"),
    ("tcga-stad", "Stomach adenocarcinoma"),
    ("gtex-stomach", "Stomach, healthy"),
    ("tcga-lihc", "Liver carcinoma"),
    ("gtex-liver", "Liver, healthy"),
]


def ranking_rows(ranking: Path, reference: str = "lambert-2018") -> list[dict]:
    """Per network: random 95% range, out-degree, repression-order, best label-free."""
    rows = []
    for name, label in NETWORKS:
        report = json.loads((ranking / name / "report.json").read_text())
        comparison = report["references"][reference]["comparison"]
        by_name = {row["pipeline"]: row for row in comparison}
        label_free = [
            row for row in comparison if row["pipeline"] != "random" and not row["uses_labels"]
        ]
        best = max(label_free, key=lambda row: row["roc_auc"])
        rows.append(
            {
                "network": name,
                "label": label,
                "random_95": by_name["random"]["roc_auc_95_range"],
                "out_degree": by_name["out-degree"]["roc_auc"],
                "repression_order": by_name["repression-order"]["roc_auc"],
                "best_label_free": best["roc_auc"],
                "best_pipeline": best["pipeline"],
                "label_free_pipelines": len(label_free),
            }
        )
    return rows


def prediction_rows(prediction: Path, gene_list: str = "essential") -> list[dict]:
    """Per single-tissue network: held-out AUC with and without network descriptors."""
    results = json.loads(prediction.read_text())["lists"][gene_list]["results"]
    rows = []
    for name, label in NETWORKS:
        if name not in results:
            continue
        auc = results[name]["auc"]
        gain = results[name]["gains"]["network_beyond_expression"]
        rows.append(
            {
                "network": name,
                "label": label,
                "expression": auc["expression"],
                "with_shape": auc["expression+shape"],
                "with_neighbours": auc["expression+shape+neighbours"],
                "gain": gain["gain"],
                "gain_95": gain["ci95"],
            }
        )
    return rows


def locality_rows(structure: Path) -> list[dict]:
    """Per network: fold enrichment of positive edges within one chromosome arm."""
    report = json.loads(structure.read_text())
    return [
        {
            "network": name,
            "label": label,
            "arm_fold": report[name]["locality"]["arm"]["positive"]["fold"],
            "band_fold": report[name]["locality"]["band"]["positive"]["fold"],
            "balanced_share": report[name]["balance"]["balanced_share"],
        }
        for name, label in NETWORKS
    ]


def _style(size: float) -> None:
    plt.rcParams.update(
        {
            "font.family": ["Helvetica Neue", "DejaVu Sans"],
            "font.size": size,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": MUTED,
            "xtick.color": MUTED,
            "ytick.color": INK,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )


def _network_axis(ax, labels: list[str]) -> list[int]:
    positions = list(range(len(labels)))[::-1]
    ax.set_yticks(positions, labels)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=BAND, linewidth=1)
    ax.set_axisbelow(True)
    return positions


def draw_ranking(rows: list[dict], output: Path, size: tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    y = _network_axis(ax, [row["label"] for row in rows])
    for position, row in zip(y, rows, strict=True):
        low, high = row["random_95"]
        ax.fill_betweenx([position - 0.38, position + 0.38], low, high, color=BAND, linewidth=0)
    series = [
        ("out_degree", "Out-degree", BLUE, "o"),
        ("repression_order", "Repression order", ORANGE, "s"),
        ("best_label_free", f"Best of {rows[0]['label_free_pipelines']} pipelines", AQUA, "^"),
    ]
    for key, label, colour, marker in series:
        ax.scatter(
            [row[key] for row in rows],
            y,
            s=220,
            marker=marker,
            color=colour,
            edgecolor="white",
            linewidth=2,
            label=label,
            zorder=3,
        )
    ax.axvline(0.5, color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)))
    ax.set_xlim(0.35, 0.65)
    ax.set_xlabel("ROC AUC, known TFs (Lambert 2018)")
    handles, _ = ax.get_legend_handles_labels()
    handles.append(Patch(color=BAND, label="Random order, 95% range"))
    fig.legend(handles=handles, loc="outside upper center", ncol=2, frameon=False)
    fig.savefig(output)
    plt.close(fig)


def draw_prediction(rows: list[dict], output: Path, size: tuple[float, float]) -> None:
    fig, (left, right) = plt.subplots(
        1, 2, figsize=size, sharey=True, layout="constrained", width_ratios=[3, 2]
    )
    y = _network_axis(left, [row["label"] for row in rows])
    _network_axis(right, [row["label"] for row in rows])
    for position, row in zip(y, rows, strict=True):
        left.plot(
            [row["expression"], row["with_neighbours"]],
            [position, position],
            color=MUTED,
            linewidth=1.5,
            zorder=2,
        )
    series = [
        ("expression", "Expression", BLUE, "o"),
        ("with_shape", "+ shape", ORANGE, "s"),
        ("with_neighbours", "+ neighbours", AQUA, "^"),
    ]
    for key, label, colour, marker in series:
        left.scatter(
            [row[key] for row in rows],
            y,
            s=220,
            marker=marker,
            color=colour,
            edgecolor="white",
            linewidth=2,
            label=label,
            zorder=3,
        )
    left.set_xlabel("Held-out AUC")
    fig.legend(loc="outside upper center", ncol=3, frameon=False)
    for position, row in zip(y, rows, strict=True):
        low, high = row["gain_95"]
        right.plot([low, high], [position, position], color=ORANGE, linewidth=4)
        right.scatter(row["gain"], position, s=120, color=ORANGE, edgecolor="white", zorder=3)
    right.axvline(0, color=MUTED, linewidth=1.5)
    right.set_xlabel("Gain, 95% CI")
    fig.savefig(output)
    plt.close(fig)


def draw_locality(rows: list[dict], output: Path, size: tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    y = _network_axis(ax, [row["label"] for row in rows])
    for position, row in zip(y, rows, strict=True):
        healthy = row["network"].startswith("gtex")
        ax.barh(
            position,
            row["arm_fold"],
            height=0.7,
            color=AQUA if healthy else BLUE,
            edgecolor="white",
            linewidth=2,
        )
        ax.text(row["arm_fold"] + 0.1, position, f"{row['arm_fold']:.1f}×", va="center")
    ax.axvline(1, color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)))
    ax.set_xlabel("Same-arm positive edges, fold over shuffled")
    handles = [Patch(color=BLUE, label="Tumour (TCGA)"), Patch(color=AQUA, label="Healthy (GTEx)")]
    fig.legend(handles=handles, loc="outside upper center", ncol=2, frameon=False)
    fig.savefig(output)
    plt.close(fig)


def write_figures(results: Path, output: Path, font_size: float = 24) -> dict:
    _style(font_size)
    output.mkdir(parents=True, exist_ok=True)
    ranking = ranking_rows(results / "ranking")
    prediction = prediction_rows(results / "prediction.json")
    locality = locality_rows(results / "structure.json")
    draw_ranking(ranking, output / "tf-ranking.pdf", (12.5, 9))
    draw_prediction(prediction, output / "essential-genes.pdf", (12.5, 9.5))
    draw_locality(locality, output / "chromosome-locality.pdf", (12.5, 8.5))
    return {"ranking": ranking, "prediction": prediction, "locality": locality}
