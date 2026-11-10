import json
import tempfile
import unittest
from pathlib import Path

from gene_network_analysis.figures import NETWORKS, prediction_rows, ranking_rows


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


if __name__ == "__main__":
    unittest.main()
