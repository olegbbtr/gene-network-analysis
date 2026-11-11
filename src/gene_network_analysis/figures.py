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


def discovery_rows(robustness: Path) -> list[dict]:
    """Per network and for all together: DepMap essentials among the top new predictions
    of a model trained on CEGv2, and the AUC gain from the network."""
    report = json.loads(robustness.read_text())["discovery"]
    chosen = [(name, label) for name, label in NETWORKS if name in report]
    rows = []
    for name, label in [*chosen, ("all", f"All {len(chosen)} networks")]:
        entry = report[name]
        rows.append(
            {
                "network": name,
                "label": label,
                "chance": 100 * entry["base_rate"],
                "expression": entry["top_hits"]["expression"],
                "with_network": entry["top_hits"]["expression+network"],
                "gain": entry["gain"],
            }
        )
    return rows


def consensus_rows(consensus: Path, kind: str = "positive") -> dict[str, dict]:
    """Per kind of known pair: its share (%) among gene pairs linked in 1..N networks, for
    all pairs on different chromosomes and for those outside the housekeeping complexes."""
    links = json.loads(consensus.read_text())["links"][kind]
    series = {
        "Pairs on different chromosomes": links["different_chromosomes"],
        "Without housekeeping complexes": links["different_chromosomes_without_core"],
    }
    return {
        known: {
            label: {
                "base": 100 * table["base"][known],
                "share": [100 * row[known] for row in table["by_networks"][1:]],
            }
            for label, table in series.items()
        }
        for known in ("reactome", "string")
    }


def understudied_rows(understudied: Path, gene_list: str = "depmap") -> tuple[list[dict], dict]:
    """Per network: the network's AUC gain over expression + STRING in each STRING group;
    and each group's range of STRING links per gene."""
    report = json.loads(understudied.read_text())
    results = report["lists"][gene_list]
    rows = [
        {"network": name, "label": label}
        | {group: results[name][group]["gain"] for group in report["string_links_per_group"]}
        for name, label in NETWORKS
        if name in results
    ]
    return rows, report["string_links_per_group"]


def prospective_rows(prospective: Path) -> list[dict]:
    """Per gene universe: mean percentile of later annotations under each ranking."""
    universes = json.loads(prospective.read_text())["universes"]
    labels = {
        "all": "All genes",
        "studied_after": "Genes studied afterwards",
        "first_annotated_after": "Not annotated before",
    }
    return [
        {"universe": name, "label": label, "pairs": universes[name]["pairs"]}
        | {
            score: (universes[name][score]["mean_percentile"], universes[name][score]["ci95"])
            for score in ("neighbours", "terms_before", "links")
        }
        for name, label in labels.items()
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


def _markers(ax, values: list[float], y: list[int], label: str, colour: str, marker: str) -> None:
    ax.scatter(
        values, y, s=220, marker=marker, color=colour, edgecolor="white", linewidth=2,
        label=label, zorder=3,
    )  # fmt: skip


def _intervals(ax, rows: list[dict], y: list[int], key: str, colour: str, shift: float = 0) -> None:
    """Each row's `key` is a paired gain: {"gain": value, "ci95": [low, high]}."""
    for position, row in zip(y, rows, strict=True):
        ax.plot(row[key]["ci95"], [position + shift] * 2, color=colour, linewidth=4)
        ax.scatter(
            row[key]["gain"], position + shift, s=120, color=colour, edgecolor="white", zorder=3
        )


def draw_discovery(rows: list[dict], output: Path, size: tuple[float, float]) -> None:
    fig, (left, right) = plt.subplots(
        1, 2, figsize=size, sharey=True, layout="constrained", width_ratios=[3, 2]
    )
    y = _network_axis(left, [row["label"] for row in rows])
    _network_axis(right, [row["label"] for row in rows])
    for position, row in zip(y, rows, strict=True):
        left.plot(
            [row["expression"], row["with_network"]], [position] * 2, color=MUTED, linewidth=1.5
        )
    _markers(left, [row["expression"] for row in rows], y, "Expression", BLUE, "o")
    _markers(left, [row["with_network"] for row in rows], y, "+ network", AQUA, "^")
    left.axvline(rows[0]["chance"], color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)))
    left.set_xlim(0, 100)
    left.set_xlabel("Confirmed among top 100 new genes")
    fig.legend(loc="outside upper center", ncol=2, frameon=False)
    _intervals(right, rows, y, "gain", ORANGE)
    right.axvline(0, color=MUTED, linewidth=1.5)
    right.set_xlabel("AUC gain, 95% CI")
    fig.savefig(output)
    plt.close(fig)


def draw_consensus(rows: dict[str, dict], output: Path, size: tuple[float, float]) -> None:
    titles = {"reactome": "Same Reactome pathway", "string": "Linked in STRING"}
    fig, axes = plt.subplots(1, 2, figsize=size, layout="constrained")
    for ax, (known, series) in zip(axes, rows.items(), strict=True):
        for (label, values), colour, marker in zip(
            series.items(), (BLUE, AQUA), ("o", "^"), strict=True
        ):
            x = range(1, len(values["share"]) + 1)
            ax.plot(x, values["share"], color=colour, linewidth=3, marker=marker, markersize=11,
                    markeredgecolor="white", label=label)  # fmt: skip
            ax.axhline(values["base"], color=colour, linewidth=1.5, linestyle=(0, (4, 3)))
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(lambda value, _: f"{value:g}%")
        ax.grid(axis="y", color=BAND, linewidth=1)
        ax.set_axisbelow(True)
        ax.set_xticks([1, 5, 10, 15])
        ax.set_xlabel("Networks linking the pair")
        ax.set_title(titles[known], fontsize=plt.rcParams["font.size"], color=INK)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=2, frameon=False)
    fig.savefig(output)
    plt.close(fig)


def draw_understudied(
    rows: list[dict], groups: dict, output: Path, size: tuple[float, float]
) -> None:
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    y = _network_axis(ax, [row["label"] for row in rows])
    for (key, (low, high)), colour, shift in zip(
        groups.items(), (ORANGE, BLUE, AQUA), (0.25, 0.0, -0.25), strict=True
    ):
        _intervals(ax, rows, y, key, colour, shift)
        ax.plot([], [], color=colour, linewidth=4, label=f"{low}–{high}")
    ax.axvline(0, color=MUTED, linewidth=1.5)
    ax.set_xlabel("AUC gain over expression + STRING, 95% CI")
    fig.legend(loc="outside upper center", ncol=3, frameon=False, title="STRING links of the gene")
    fig.savefig(output)
    plt.close(fig)


def draw_prospective(rows: list[dict], output: Path, size: tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    y = _network_axis(ax, [f"{row['label']}\n({row['pairs']:,} pairs)" for row in rows])
    series = [
        ("neighbours", "Network neighbours", AQUA, 0.25),
        ("terms_before", "Earlier annotations", BLUE, 0.0),
        ("links", "Links", MUTED, -0.25),
    ]
    for key, label, colour, shift in series:
        for position, row in zip(y, rows, strict=True):
            value, (low, high) = row[key]
            ax.plot([low, high], [position + shift] * 2, color=colour, linewidth=4)
            ax.scatter(value, position + shift, s=160, color=colour, edgecolor="white", zorder=3)
        ax.plot([], [], color=colour, linewidth=4, label=label)
    ax.axvline(0.5, color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)))
    ax.set_xlabel("Mean percentile of later functions")
    fig.legend(loc="outside upper center", ncol=3, frameon=False)
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
    discovery = discovery_rows(results / "robustness.json")
    consensus = consensus_rows(results / "consensus.json")
    understudied, groups = understudied_rows(results / "understudied.json")
    prospective = prospective_rows(results / "prospective.json")
    draw_discovery(discovery, output / "essential-discovery.pdf", (12.5, 9.5))
    draw_consensus(consensus, output / "consensus-links.pdf", (12.5, 7))
    draw_understudied(understudied, groups, output / "understudied-genes.pdf", (12.5, 9.5))
    draw_prospective(prospective, output / "later-functions.pdf", (12.5, 6))
    return {
        "ranking": ranking,
        "prediction": prediction,
        "locality": locality,
        "discovery": discovery,
        "consensus": consensus,
        "understudied": understudied,
        "prospective": prospective,
    }
