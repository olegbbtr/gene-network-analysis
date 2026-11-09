import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.evaluation import (
    bootstrap_aucs,
    p_value,
    random_rankings,
    ranking_metrics,
    read_reference,
    resolve_labels,
)
from gene_network_analysis.graph import component_summary, pagerank, read_graph
from gene_network_analysis.pipelines import PIPELINES
from gene_network_analysis.report import analyze
from gene_network_analysis.steps import (
    drop_sinks,
    greedy_layering,
    keep_sign,
    minimum_feedback_order,
    orient_reciprocals,
    prune_triangles,
    resolve_reciprocals,
    top_outgoing,
)


def matrix(edges, size):
    rows, columns, weights = zip(*edges, strict=True)
    return csr_array((np.asarray(weights, dtype=float), (rows, columns)), shape=(size, size))


class GraphTests(unittest.TestCase):
    def test_reader_preserves_identity_and_sign(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "graph.txt"
            path.write_text("A10 A2 -2\nA2 A3 1\n")
            graph = read_graph(path)
            self.assertEqual(graph.nodes, ("A10", "A2", "A3"))
            self.assertEqual(graph.adjacency[0, 1], -2)
            path.write_text("A10 A2 1\nA10 A2 2\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                read_graph(path)
            for contents in ("a b nan\n", "a b inf\n", "a b 0\n", "a b\n", ""):
                path.write_text(contents)
                with self.assertRaises(ValueError):
                    read_graph(path)

    def test_pagerank_direction_and_dangling_mass(self):
        graph = matrix([(0, 1, 1), (1, 2, 1)], 3)
        forward, _ = pagerank(graph)
        reverse, _ = pagerank(graph, reverse=True)
        self.assertGreater(forward[2], forward[0])
        self.assertGreater(reverse[0], reverse[2])
        self.assertAlmostEqual(float(forward.sum()), 1)
        self.assertAlmostEqual(float(reverse.sum()), 1)
        np.testing.assert_allclose(forward, reverse[::-1])

    def test_pagerank_balance_puts_sources_above_sinks(self):
        from gene_network_analysis.pipelines import _pagerank_balance

        graph = matrix([(0, 1, 2.0), (0, 2, -1.0), (1, 2, 1.0)], 3)
        balance = _pagerank_balance(graph, 0.85)
        self.assertGreater(balance[0], balance[1])
        self.assertGreater(balance[1], balance[2])

    def test_order_pagerank_rank_sum_favours_early_reaching_genes(self):
        from gene_network_analysis.pipelines import PIPELINES, Runner

        graph = matrix([(0, 1, -2.0), (1, 2, -1.0), (0, 3, -1.0), (2, 0, -0.5)], 4)
        pipeline = next(p for p in PIPELINES if p.name == "repression-order-pagerank")
        scores = Runner(graph).score(pipeline, np.arange(4))
        self.assertEqual(int(np.argmax(scores)), 0)

    def test_pagerank_stationary_equation(self):
        graph = matrix([(0, 1, -2), (0, 2, 1), (1, 0, 1)], 4)
        ranks, _ = pagerank(graph, weighted=True, alpha=0.85)
        dense = abs(graph.toarray())
        degree = dense.sum(axis=1)
        for node in range(4):
            if degree[node]:
                dense[node] /= degree[node]
            else:
                dense[node] = 0.25
        np.testing.assert_allclose(ranks, 0.85 * dense.T @ ranks + 0.15 / 4, atol=1e-11)
        with self.assertRaises(RuntimeError):
            pagerank(graph, max_iterations=1)

    def test_components_use_correct_finishing_order(self):
        graph = matrix([(0, 1, 1), (1, 2, 1)], 3)
        self.assertEqual(component_summary(graph)["strong_components"], 3)
        graph = matrix([(0, 1, 1), (1, 0, 1), (1, 2, 0.1)], 4)
        summary = component_summary(graph, 0.5)
        self.assertEqual(summary["strong_components"], 3)
        self.assertEqual(summary["largest_strong_component"], 2)
        self.assertEqual(summary["isolated_nodes"], 2)

    def test_layering_keeps_heavy_edges_and_is_acyclic(self):
        graph = matrix([(0, 1, 10), (1, 2, 9), (2, 0, 1), (3, 4, 3), (4, 4, 20)], 5)
        dag, layers = greedy_layering(graph)
        edges = dag.tocoo()
        selected = set(zip(edges.row, edges.col, strict=True))
        self.assertIn((0, 1), selected)
        self.assertIn((1, 2), selected)
        self.assertNotIn((2, 0), selected)
        self.assertNotIn((4, 4), selected)
        self.assertTrue(np.all(layers[edges.row] < layers[edges.col]))

    def test_reciprocal_subtraction_handles_signs_and_ties(self):
        graph = matrix([(0, 1, -5), (1, 0, -2), (1, 2, 3), (2, 1, -4)], 3)
        cleaned = resolve_reciprocals(graph)
        self.assertEqual(cleaned[0, 1], -3)
        self.assertEqual(cleaned[1, 0], 0)
        self.assertEqual(cleaned[1, 2], 3)
        self.assertEqual(cleaned[2, 1], -4)
        tied = resolve_reciprocals(matrix([(0, 1, 2), (1, 0, 2)], 2))
        self.assertEqual(tied.nnz, 0)
        self.assertEqual(graph.nnz, 4)

    def test_sink_filter_is_one_pass(self):
        graph = matrix([(0, 1, 1), (1, 2, 1)], 3)
        cleaned = drop_sinks(graph)
        self.assertEqual(cleaned.nnz, 1)
        self.assertEqual(cleaned[0, 1], 1)

    def test_triangles_need_one_sign_on_all_three_edges(self):
        same = matrix([(0, 1, -2), (1, 2, -3), (0, 2, -1)], 3)
        pruned = prune_triangles(same, block=1)
        self.assertEqual(pruned[0, 2], 0)
        self.assertEqual(pruned.nnz, 2)
        mixed = matrix([(0, 1, -2), (1, 2, -3), (0, 2, 1)], 3)
        self.assertEqual(prune_triangles(mixed).nnz, 3)

    def test_orientation_keeps_only_a_clearly_stronger_side(self):
        graph = matrix([(0, 1, 6), (1, 0, -3), (1, 2, 4), (2, 1, 3), (2, 3, 1)], 4)
        oriented = orient_reciprocals(graph, ratio=2.0)
        self.assertEqual(set(zip(*oriented.nonzero(), strict=True)), {(0, 1), (2, 3)})
        self.assertEqual(oriented[0, 1], 6)

    def test_top_outgoing_keeps_the_strongest_edges_per_gene(self):
        graph = matrix([(0, 1, 1), (0, 2, -5), (0, 3, 3), (1, 0, 2)], 4)
        kept = top_outgoing(graph, k=2)
        self.assertEqual(set(zip(*kept.nonzero(), strict=True)), {(0, 2), (0, 3), (1, 0)})

    def test_sign_filter(self):
        graph = matrix([(0, 1, 1), (1, 2, -2)], 3)
        self.assertEqual(keep_sign(graph, 1).nnz, 1)
        self.assertEqual(keep_sign(graph, -1)[1, 2], -2)

    def test_feedback_order_is_topological_on_a_dag(self):
        graph = matrix([(3, 1, 1), (1, 0, 2), (3, 0, 1), (0, 2, 5)], 4)
        dag, position = minimum_feedback_order(graph)
        self.assertEqual(dag.nnz, graph.nnz)
        self.assertEqual(sorted(position.tolist()), [0, 1, 2, 3])
        self.assertEqual(int(np.argmin(position)), 3)

    def test_feedback_order_drops_the_lighter_edge_of_a_cycle(self):
        graph = matrix([(0, 1, 5), (1, 2, -4), (2, 0, 1), (2, 3, 1)], 4)
        dag, position = minimum_feedback_order(graph)
        kept = set(zip(*dag.nonzero(), strict=True))
        self.assertEqual(kept, {(0, 1), (1, 2), (2, 3)})
        self.assertEqual(dag[1, 2], -4)
        self.assertEqual(component_summary(dag)["strong_components"], 4)


class EvaluationTests(unittest.TestCase):
    def test_reference_and_graph_identities_share_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / "tf.txt"
            reference.write_text("Ensembl Gene ID\nENSG0001.3\nENSG0001\n")
            known_tfs = read_reference(reference)
            labels, known = resolve_labels(("ENSG0001.7", "ENSG0002", "X1"), known_tfs)
            self.assertEqual(known_tfs, {"ENSG0001"})
            self.assertEqual(labels.tolist(), [True, False, False])
            self.assertEqual(known.tolist(), [True, True, False])

    def test_version_aliases_cannot_double_count_a_gene(self):
        with self.assertRaisesRegex(ValueError, "same gene"):
            resolve_labels(("ENSG1.2", "ENSG1.3"), {"ENSG1"})

    def test_perfect_and_reverse_auc_and_chance(self):
        scores = np.array([4.0, 3.0, 2.0, 1.0])
        labels = np.array([True, True, False, False])
        result = ranking_metrics(scores, labels, ["a", "b", "c", "d"], k=2)
        self.assertEqual(result["roc_auc"], 1)
        self.assertEqual(result["top_k_reference_genes"], 2)
        self.assertEqual(result["expected_by_chance"], 1)
        self.assertEqual(result["fold_enrichment"], 2)
        self.assertEqual(ranking_metrics(-scores, labels, list("abcd"), 2)["roc_auc"], 0)

    def test_ties_do_not_inflate_auc(self):
        result = ranking_metrics(np.ones(4), np.array([True, False, True, False]), list("abcd"), 1)
        self.assertEqual(result["roc_auc"], 0.5)
        self.assertEqual(result["tie_averaged_top_k_reference_genes"], 0.5)
        self.assertEqual(result["boundary_tie_size"], 4)


class RandomBaselineTests(unittest.TestCase):
    def test_random_orderings_centre_on_chance(self):
        labels = np.zeros(1000, dtype=bool)
        labels[:100] = True
        null = random_rankings(labels, k=100, runs=500)
        self.assertAlmostEqual(float(null["roc_auc"].mean()), 0.5, delta=0.01)
        self.assertAlmostEqual(float(null["top_k_reference_genes"].mean()), 10, delta=0.5)
        self.assertEqual(p_value(null["roc_auc"], 1.0), 1 / 501)
        self.assertEqual(p_value(null["roc_auc"], 0.0), 1.0)


class BootstrapTests(unittest.TestCase):
    def test_bootstrap_matches_exact_auc_on_average(self):
        rng = np.random.default_rng(3)
        labels = rng.random(400) < 0.3
        score = labels + rng.normal(size=400)
        exact = ranking_metrics(score, labels, [str(i) for i in range(400)])["roc_auc"]
        resampled = bootstrap_aucs({"s": score, "t": np.round(score)}, labels, runs=300)
        self.assertAlmostEqual(float(resampled["s"].mean()), exact, delta=0.01)
        tied = ranking_metrics(np.round(score), labels, [str(i) for i in range(400)])["roc_auc"]
        self.assertAlmostEqual(float(resampled["t"].mean()), tied, delta=0.01)


class PipelineTests(unittest.TestCase):
    def test_pipeline_evaluates_signed_graph_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            graph, first, second = base / "graph.txt", base / "a.txt", base / "b.txt"
            graph.write_text("ENSG1 ENSG2 9\nENSG2 ENSG3 -4\nENSG1 ENSG3 2\n")
            first.write_text("Ensembl Gene ID\nENSG1\n")
            second.write_text("ENSG2\n")
            references = {"a": first, "b": second}
            summary = analyze(graph, references, base / "out", top_k=1)
            names = [p["name"] for p in summary["pipelines"]]
            self.assertEqual(names, [p.name for p in PIPELINES if p.reported])
            self.assertEqual(summary["selection"]["genes"], 2)
            comparison = summary["references"]["a"]["comparison"]
            rows = {r["pipeline"]: r for r in comparison}
            self.assertEqual(comparison[0]["pipeline"], "random")
            self.assertEqual(rows["out-degree"]["top_k_reference_genes"], 1)
            self.assertEqual(rows["weighted-out-degree"]["auc_minus_best_baseline"], 0)
            self.assertIn("p_auc_vs_random", rows["net-flow"])
            other = {r["pipeline"]: r for r in summary["references"]["b"]["comparison"]}
            self.assertEqual(other["out-degree"]["top_k_reference_genes"], 0)
            steps = summary["pipelines"][names.index("top50-order")]["steps"]
            self.assertEqual(steps["preprocess"][0], {"step": "top_outgoing", "k": 50})
            only = analyze(graph, references, base / "one", pipelines=["net-flow"], top_k=1)
            self.assertEqual([p["name"] for p in only["pipelines"]], ["net-flow"])
            self.assertTrue((base / "out/report.pdf").stat().st_size > 0)
            self.assertEqual(json.loads((base / "out/report.json").read_text()), summary)


if __name__ == "__main__":
    unittest.main()
