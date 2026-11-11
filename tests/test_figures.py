import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from matplotlib.figure import Figure

from gene_network_analysis.figures import (
    NETWORKS,
    consensus_rows,
    discovery_rows,
    prediction_rows,
    ranking_rows,
    write_poster_figures,
)


def rows_from(function, report: dict):
    with tempfile.NamedTemporaryFile("w", suffix=".json") as handle:
        handle.write(json.dumps(report))
        handle.flush()
        return function(Path(handle.name))


def comparison(auc: dict[str, float], labels: dict[str, float]) -> list[dict]:
    rows = [{"pipeline": "random", "roc_auc": 0.5, "roc_auc_95_range": [0.48, 0.52]}]
    rows += [{"pipeline": n, "roc_auc": a, "uses_labels": False} for n, a in auc.items()]
    rows += [{"pipeline": n, "roc_auc": a, "uses_labels": True} for n, a in labels.items()]
    return rows


class RankingRowsTests(unittest.TestCase):
    def test_best_pipeline_ignores_label_using_pipelines(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, _ in NETWORKS:
                (root / name).mkdir()
                rows = comparison(
                    {"out-degree": 0.45, "repression-order": 0.42, "net-flow": 0.53},
                    {"similar-targets": 0.75},
                )
                report = {"references": {"lambert-2018": {"comparison": rows}}}
                (root / name / "report.json").write_text(json.dumps(report))
            row = ranking_rows(root)[0]
        self.assertEqual(row["best_pipeline"], "net-flow")
        self.assertEqual(row["label_free_pipelines"], 3)
        self.assertEqual(row["random_95"], [0.48, 0.52])


class PredictionRowsTests(unittest.TestCase):
    def test_skips_networks_without_results_and_keeps_order(self):
        result = {
            "auc": {
                "expression": 0.7,
                "expression+shape": 0.75,
                "expression+shape+neighbours": 0.8,
            },
            "gains": {"network_beyond_expression": {"gain": 0.05, "ci95": [0.03, 0.07]}},
        }
        names = ["gtex-liver", "tcga-luad"]
        report = {"lists": {"essential": {"results": {n: result for n in names}}}}
        with tempfile.NamedTemporaryFile("w", suffix=".json") as handle:
            handle.write(json.dumps(report))
            handle.flush()
            rows = prediction_rows(Path(handle.name))
        self.assertEqual([row["network"] for row in rows], ["tcga-luad", "gtex-liver"])
        self.assertEqual(rows[0]["gain_95"], [0.03, 0.07])


class NewResultRowsTests(unittest.TestCase):
    def test_nested_metrics_require_the_same_inputs_as_the_fixed_comparison(self):
        from gene_network_analysis.figures import feature_comparison_rows

        fixed = {
            "inputs": {"graph": "abc"},
            "discovery": {"base_rate": 0.1},
            "main_comparison": [
                {"key": "global", "label": "Global graph only", "known_auc": 0.8, "top_hits": 60}
            ],
        }
        nested = {
            "inputs": {"graph": "abc"},
            "nested": {"global": {"known_auc": 0.9, "top_hits": 70, "discovery_auc": 0.85}},
        }
        with tempfile.TemporaryDirectory() as directory:
            report, audit = Path(directory) / "fixed.json", Path(directory) / "audit.json"
            report.write_text(json.dumps(fixed))
            audit.write_text(json.dumps(nested))
            rows, chance = feature_comparison_rows(report, audit)
            self.assertEqual(rows[0]["known_auc"], 0.9)
            self.assertEqual(rows[0]["top_hits"], 70)
            self.assertEqual(rows[0]["label"], "Global graph only")
            self.assertEqual(chance, 10)
            nested["inputs"]["graph"] = "changed"
            audit.write_text(json.dumps(nested))
            with self.assertRaisesRegex(ValueError, "identical input hashes"):
                feature_comparison_rows(report, audit)

    def test_main_plot_keeps_known_label_auc_separate_from_discovery_auc(self):
        rows = [
            {"known_auc": a, "discovery_auc": b, "top_hits": hits}
            for a, b, hits in [(0.9, 0.6, 57), (0.8, 0.7, 70), (0.95, 0.65, 71)]
        ]
        observed = {}

        def inspect(figure, path):
            if path.name == "feature-comparison.pdf":
                left, right = figure.axes
                observed["auc"] = [float(c.get_offsets()[0, 0]) for c in left.collections]
                observed["hits"] = [p.get_width() for p in right.patches]
                observed["labels"] = [tick.get_text() for tick in left.get_yticklabels()]

        series = {k: {"base": 1, "share": list(range(1, 16))} for k in ("all", "excluded")}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "feature-comparison.json").write_text(
                json.dumps({"main_comparison": rows, "discovery": {"base_rate": 0.1}})
            )
            with (
                patch(
                    "gene_network_analysis.figures.consensus_rows",
                    return_value={"reactome": series},
                ),
                patch.object(Figure, "savefig", autospec=True, side_effect=inspect),
            ):
                write_poster_figures(root, root)
        self.assertEqual(observed["auc"], [0.9, 0.8, 0.95])
        self.assertEqual(observed["hits"], [57, 70, 71])

    def test_poster_prefers_five_group_report_and_preserves_row_labels(self):
        labels = [
            "Gene only",
            "Local graph only",
            "Global graph only",
            "Gene + local graph",
            "Gene + local graph + global graph",
        ]
        rows = [
            {"label": label, "known_auc": 0.6 + i / 20, "top_hits": i + 40}
            for i, label in enumerate(labels)
        ]
        observed = {}

        def inspect(figure, path):
            if path.name == "feature-comparison.pdf":
                left, right = figure.axes
                observed["labels"] = [
                    tick.get_text().replace("\n", " ") for tick in left.get_yticklabels()
                ]
                observed["hits"] = [p.get_width() for p in right.patches]

        series = {k: {"base": 1, "share": list(range(1, 16))} for k in ("all", "excluded")}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "feature-comparison.json").write_text("{}")
            (root / "feature-groups.json").write_text(
                json.dumps({"main_comparison": rows, "discovery": {"base_rate": 0.1}})
            )
            with (
                patch(
                    "gene_network_analysis.figures.consensus_rows",
                    return_value={"reactome": series},
                ),
                patch.object(Figure, "savefig", autospec=True, side_effect=inspect),
            ):
                write_poster_figures(root, root)
        self.assertEqual(observed["labels"], labels)
        self.assertEqual(observed["hits"], list(range(40, 45)))

    def test_discovery_rows_end_with_all_networks_and_give_chance_per_hundred(self):
        entry = {
            "base_rate": 0.094,
            "top_hits": {"expression": 45, "expression+network": 67},
            "gain": {"gain": 0.08, "ci95": [0.06, 0.09]},
        }
        report = {"discovery": {n: entry for n in ("all", "gtex-liver", "tcga-luad")}}
        rows = rows_from(discovery_rows, report)
        self.assertEqual([row["network"] for row in rows], ["tcga-luad", "gtex-liver", "all"])
        self.assertEqual(rows[-1]["label"], "All 2 networks")
        self.assertAlmostEqual(rows[0]["chance"], 9.4)
        self.assertEqual((rows[0]["expression"], rows[0]["with_network"]), (45, 67))

    def test_consensus_rows_drop_unlinked_pairs_and_give_percentages(self):
        table = {
            "base": {"reactome": 0.01, "string": 0.002},
            "by_networks": [
                {"networks": k, "reactome": 0.01 * (k + 1), "string": 0.1 * k} for k in range(3)
            ],
        }
        report = {
            "links": {
                "positive": {
                    "different_chromosomes": table,
                    "different_chromosomes_without_core": table,
                }
            }
        }
        rows = rows_from(consensus_rows, report)
        series = rows["reactome"]["Pairs on different chromosomes"]
        self.assertAlmostEqual(series["base"], 1.0)
        self.assertEqual(len(series["share"]), 2)
        self.assertAlmostEqual(series["share"][0], 2.0)
        self.assertAlmostEqual(rows["string"]["Without housekeeping complexes"]["share"][1], 20.0)


if __name__ == "__main__":
    unittest.main()
