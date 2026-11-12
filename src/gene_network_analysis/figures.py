"""Poster figures drawn from the committed results in results/public/."""

import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, BAND = "#0b0b0b", "#52514e", "#dddcd7"
IMPACT = "shapley"

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


def feature_comparison_rows(report: Path, audit: Path | None = None) -> tuple[list[dict], float]:
    data = json.loads(report.read_text())
    rows = data["main_comparison"]
    if audit is not None:
        tuned = json.loads(audit.read_text())
        if not data.get("inputs") or data["inputs"] != tuned.get("inputs"):
            raise ValueError("Nested and fixed comparisons must use identical input hashes")
        rows = [row | tuned["nested"][row["key"]] for row in rows]
    return rows, 100 * data["discovery"]["base_rate"]


def draw_impact(ax, rows, positions, colours, impact: dict) -> None:
    """Leave-one-out AUC loss per block on CEGv2 (solid) and DepMap outside CEGv2 (light),
    with 95% intervals; the graph block spans the local and global rows."""
    y = {row["key"]: pos for row, pos in zip(rows, positions, strict=True)}
    colour = {row["key"]: c for row, c in zip(rows, colours, strict=True)}
    centres = {
        "gene": (y["gene"], colour["gene"]),
        "graph": ((y["local"] + y["global"]) / 2, colour["local"]),
        "biological": (y["biological"], colour["biological"]),
    }
    values = [v for lists in impact.values() for v in lists.values()]
    top = max(v["auc_loss_ci95"][1] for v in values)
    for label, shift, alpha in (("cegv2", 0.17, 1.0), ("depmap", -0.17, 0.45)):
        for block, (centre, face) in centres.items():
            value = impact[label][block]
            at = centre + shift
            ax.barh(at, value["auc_loss"], height=0.3, color=face, alpha=alpha)
            low, high = value["auc_loss_ci95"]
            ax.plot([low, high], [at, at], color=INK, linewidth=2)
            ax.text(high + 0.04 * top, at, f"{value['auc_loss']:.3f}", va="center", fontsize=19)
    ax.text(0.03 * top, min(y.values()), "Graph bars: local and global\nremoved together."
            "\nDark: CEGv2; light: DepMap.", va="center", fontsize=17, color=MUTED)
    ax.axvline(0, color=MUTED, linewidth=2)
    low = min(v["auc_loss_ci95"][0] for v in values)
    ax.set_xlim(min(0, low) - 0.05 * top, 1.45 * top)
    ax.set_xlabel("Leave-one-out impact\nAUC lost when removed")


def write_poster_figures(results: Path, output: Path) -> None:
    _style(24)
    output.mkdir(parents=True, exist_ok=True)
    label_free = results / "label-free-groups.json"
    impact = {}
    if label_free.exists():
        data = json.loads(label_free.read_text())
        rows, chance = data["main_comparison"], 100 * data["discovery_base_rate"]
        impact = {k: v for k, v in data.get("leave_one_out", {}).items() if k != "method"}
    else:
        comparison = results / "feature-groups.json"
        if not comparison.exists():
            comparison = results / "feature-comparison.json"
        audit = results / "feature-audit.json"
        rows, chance = feature_comparison_rows(comparison, audit if audit.exists() else None)
    ratios = [1, 1.1, 0.9] if impact else [1, 1.1]
    fig, axes = plt.subplots(
        1,
        len(ratios),
        figsize=(13.5 * sum(ratios) / 2.1, 1.12 * len(rows)),
        layout="constrained",
        sharey=True,
        width_ratios=ratios,
    )
    left, right = axes[:2]
    positions = list(reversed(range(len(rows))))
    colours = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#000000", "#E69F00", "#56B4E9"][
        : len(rows)
    ]
    markers = ["o", "s", "D", "^", "P", "X", "v"][: len(rows)]
    labels = [
        textwrap.fill(row.get("label", row.get("key", "Model")), width=21, break_long_words=False)
        for row in rows
    ]
    left.set_yticks(positions, labels)
    left.tick_params(axis="y", length=0)
    for row, y, colour, marker in zip(rows, positions, colours, markers, strict=True):
        left.plot([0.5, row["known_auc"]], [y, y], color=colour, linewidth=4)
        left.scatter(row["known_auc"], y, color=colour, marker=marker, s=220, zorder=3)
        left.annotate(
            f"{row['known_auc']:.3f}",
            (row["known_auc"], y),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=21,
        )
        right.barh(y, row["top_hits"], color=colour, height=0.32)
        right.text(row["top_hits"] + 2, y, str(row["top_hits"]), va="center", weight="bold")
    if impact:
        draw_impact(axes[2], rows, positions, colours, impact)
    left.set_xlim(0.5, 1.0)
    left.set_xticks([0.5, 0.7, 0.9, 1.0])
    left.set_xlabel("Known essentials\nROC AUC on held-out genes")
    right.set_xlim(0, 100)
    right.set_xticks([0, 25, 50, 75, 100])
    right.set_xlabel("Known essentials\nin the top 100")
    right.axvline(chance, color=MUTED, linestyle="--", linewidth=2)
    for ax in axes:
        ax.set_ylim(-0.6, len(rows) - 0.3)
        ax.grid(axis="x", color=BAND, linewidth=1)
        ax.set_axisbelow(True)
    fig.savefig(output / "feature-comparison.pdf")
    plt.close(fig)
    if (results / "targets.json").exists():
        draw_targets(target_rows(results), output)

    series = consensus_rows(results / "consensus.json")["reactome"]
    fig, ax = plt.subplots(figsize=(13.5, 4.2), layout="constrained")
    for (_, values), label, colour, marker in zip(
        series.items(),
        ["All eligible gene pairs", "Housekeeping excluded"],
        ["#0072B2", "#D55E00"],
        ["o", "^"],
        strict=True,
    ):
        ax.plot(
            range(1, 16),
            values["share"],
            color=colour,
            linewidth=3,
            marker=marker,
            markersize=9,
            label=label,
        )
        ax.axhline(values["base"], color=colour, linestyle="--", linewidth=1.5)
        ax.annotate(
            f"{values['share'][-1]:.1f}%",
            (15, values["share"][-1]),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            color=colour,
        )
    ax.set_xlim(0.5, 17)
    ax.set_ylim(0, 31)
    ax.set_xticks([1, 5, 10, 15])
    ax.set_xlabel("Networks sharing a positive link")
    ax.set_ylabel("Pairs in the same\nReactome pathway (%)")
    ax.legend(loc="upper left", frameon=False, fontsize=20)
    ax.grid(axis="y", color=BAND, linewidth=1)
    ax.set_axisbelow(True)
    fig.savefig(output / "functional-consensus.pdf")
    plt.close(fig)


TARGET_PANELS = [
    ("essentiality", "Essential genes\n(CEGv2)"),
    ("selective_dependency", "Selective cancer\ndependency (DepMap)"),
    ("knockdown_impact", "Knockdown impact\n(Perturb-seq)"),
    ("first_function", "First experimental\nfunction (GO)"),
]
GROUP_COLOURS = {
    "gene": "#0072B2",
    "local": "#D55E00",
    "global": "#009E73",
    "biological": "#CC79A7",
    "gene+local+global+biological": "#000000",
}


def target_rows(results: Path) -> list[dict]:
    """Per target: AUC per model, leave-one-out losses, genes, positives, top-100 hits."""
    essential = json.loads((results / "label-free-groups.json").read_text())
    rows = [
        {
            "key": "essentiality",
            "auc": {r["key"]: r["known_auc"] for r in essential["main_comparison"]},
            "loss": essential["leave_one_out"]["cegv2"],
            "genes": 12000,
            "positives": 606,
            "base": 606 / 12000,
        }
    ]
    targets = json.loads((results / "targets.json").read_text())
    for key, _ in TARGET_PANELS[1:]:
        t = targets[key]
        rows.append(
            {
                "key": key,
                "auc": {k: v["auc"] for k, v in t["models"].items()},
                "loss": t["leave_one_out"],
                "genes": t["genes"],
                "positives": t["positives"],
                "base": t["base_rate"],
            }
        )
    return rows


def draw_targets(rows: list[dict], output: Path) -> None:
    """One AUC panel per target, drawn at poster size (36.6 cm wide)."""
    _style(19)
    models = list(GROUP_COLOURS)
    labels = ["Gene", "Local graph", "Global graph", "Biological", "All four"]
    fig, axes = plt.subplots(1, len(rows), figsize=(14.4, 4.3), layout="constrained", sharey=True)
    y = list(reversed(range(len(models))))
    for column, (ax, row, (_, title)) in enumerate(zip(axes, rows, TARGET_PANELS, strict=True)):
        for yy, m in zip(y, models, strict=True):
            v = row["auc"][m]
            ax.plot([0.5, v], [yy, yy], color=GROUP_COLOURS[m], linewidth=3.5)
            ax.scatter(v, yy, color=GROUP_COLOURS[m], s=110, zorder=3)
            ax.text(v + 0.015, yy, f"{v:.2f}", va="center", fontsize=16)
        ax.set_xlim(0.5, 1.0)
        ax.set_xticks([0.6, 0.8, 1.0])
        ax.set_title(f"{title}\n{row['positives']:,} of {row['genes']:,} genes", fontsize=17)
        if column == 0:
            ax.set_yticks(y, labels)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=BAND, linewidth=1)
        ax.set_axisbelow(True)
    fig.supxlabel("ROC AUC on held-out genes", fontsize=18)
    fig.savefig(output / "target-comparison.pdf")
    plt.close(fig)
