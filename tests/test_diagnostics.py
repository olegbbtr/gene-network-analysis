import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.diagnostics import diagnose, sender_selection, top_genes
from gene_network_analysis.neighbours import neighbour_vote, similarity_vote
from gene_network_analysis.pipelines import PIPELINES, Runner
from gene_network_analysis.references import read_gene_table, write_gene_table


def matrix(edges, size):
    rows, columns, weights = zip(*edges, strict=True)
    return csr_array((np.asarray(weights, dtype=float), (rows, columns)), shape=(size, size))


class NeighbourVoteTests(unittest.TestCase):
    def test_uses_absolute_weight_and_ignores_own_label(self):
        graph = matrix([(0, 1, 3.0), (0, 2, -1.0), (0, 0, 5.0), (3, 0, 2.0)], 4)
        labels = np.array([True, True, False, False])
        out = neighbour_vote(graph, labels, "out")
        self.assertAlmostEqual(out[0], 0.75)
        self.assertEqual(out[1], 0.0)
        incoming = neighbour_vote(graph, labels, "in")
        self.assertEqual(incoming[0], 0.0)
        self.assertEqual(incoming[1], 1.0)
        both = neighbour_vote(graph, labels, "both")
        self.assertAlmostEqual(both[0], 3 / 6)


class SimilarityVoteTests(unittest.TestCase):
    def test_votes_with_genes_sharing_targets(self):
        # Genes 0 and 1 share targets 4 and 5; genes 2 and 3 share targets 6 and 7.
        graph = matrix(
            [(0, 4, 3), (0, 5, 3), (1, 4, 3), (1, 5, -3), (2, 6, 3), (2, 7, 3), (3, 6, 3)], 8
        )
        rows = np.arange(4)
        labels = np.array([True, False, False, False])
        scores = similarity_vote(graph, rows, labels, targets=2, neighbours=1, block=3)
        np.testing.assert_allclose(scores, [0.0, 1.0, 0.0, 0.0])

    def test_signed_profiles_separate_opposite_targets(self):
        graph = matrix([(0, 2, 3), (1, 2, -3)], 3)
        labels = np.array([True, False])
        unsigned = similarity_vote(graph, np.arange(2), labels, neighbours=1)
        signed = similarity_vote(graph, np.arange(2), labels, neighbours=1, signed=True)
        self.assertEqual(unsigned[1], 1.0)
        self.assertEqual(signed[1], 0.0)


class SenderSelectionTests(unittest.TestCase):
    def test_counts_protein_coding_genes_only(self):
        senders = np.array([True, True, True, False, False, False])
        labels = np.array([True, True, False, False, True, True])
        coding = np.array([True, True, True, True, True, False])
        result = sender_selection(senders, labels, coding)
        self.assertEqual(result["table"]["counts"], [[2, 1], [1, 1]])
        self.assertAlmostEqual(result["tf_share_among_senders"], 2 / 3)
        self.assertAlmostEqual(result["odds_ratio"], 2.0)


class LabelPipelineTests(unittest.TestCase):
    def test_label_pipelines_need_labels_and_never_see_their_own(self):
        graph = matrix([(0, 2, 3), (1, 2, 3), (2, 0, -1)], 3)
        runner = Runner(graph)
        pipeline = next(p for p in PIPELINES if p.name == "tf-target-share")
        universe = np.arange(3)
        with self.assertRaises(ValueError):
            runner.score(pipeline, universe)
        labels = np.array([True, False, True])
        np.testing.assert_allclose(runner.score(pipeline, universe, labels), [1, 1, 1])
        labels[2] = False
        np.testing.assert_allclose(runner.score(pipeline, universe, labels), [0, 0, 1])


class TopGenesTests(unittest.TestCase):
    def test_lists_first_ranked_genes_with_family_and_chromosome(self):
        genes = {
            "ENSG1": {
                "symbol": "ZNF1",
                "location": "19q13.4",
                "gene_group": ["Zinc fingers C2H2-type"],
            },
            "ENSG2": {"symbol": "FOS", "location": "14q24.3", "gene_group": []},
        }
        result = top_genes(
            np.array([1.0, 2.0, 2.0]),
            np.array([True, True, False]),
            ["ENSG1", "ENSG3", "ENSG2"],
            genes,
            2,
        )
        self.assertEqual([g["gene"] for g in result["genes"]], ["ENSG2", "ENSG3"])
        self.assertEqual(result["reference_genes"], 1)
        self.assertEqual(result["family_genes"], 0)
        self.assertEqual(result["chromosomes"], {"14": 1, "unknown": 1})


class EndToEndTests(unittest.TestCase):
    def test_diagnostics_and_gene_table(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            genes = [f"ENSG{i:011d}" for i in range(8)]
            edges = [(0, 4), (0, 5), (1, 4), (1, 5), (2, 6), (2, 7), (3, 6), (3, 7), (4, 0)]
            (base / "graph.txt").write_text(
                "".join(f"{genes[a]}\t{genes[b]}\t{3 if a % 2 else -3}\n" for a, b in edges)
            )
            (base / "tfs.txt").write_text(f"{genes[0]}\n{genes[1]}\n")
            (base / "hgnc.tsv").write_text(
                "hgnc_id\tensembl_gene_id\tsymbol\tlocation\tlocus_group\tgene_group\n"
                + "".join(
                    f"H{i}\t{gene}\tG{i}\t19q13\tprotein-coding gene\t"
                    f"{'KRAB domain containing|Zinc fingers C2H2-type' if i == 2 else ''}\n"
                    for i, gene in enumerate(genes)
                )
                + "H9\t\tNOID\t1p\tprotein-coding gene\t\n"
            )
            table = base / "genes.tsv"
            self.assertEqual(write_gene_table(base / "hgnc.tsv", table), 8)
            self.assertEqual(
                read_gene_table(table)[genes[2]]["gene_group"],
                ["KRAB domain containing", "Zinc fingers C2H2-type"],
            )
            report = diagnose(base / "graph.txt", {"t": base / "tfs.txt"}, table, base / "out", 1)
            result = report["references"]["t"]
            self.assertEqual(report["selection"]["genes"], 5)
            self.assertEqual(result["excluded_genes"], 1)
            rows = {r["pipeline"]: r for r in result["all_senders"]["comparison"]}
            self.assertEqual(rows["similar-targets"]["roc_auc"], 1.0)
            self.assertTrue(rows["similar-targets"]["uses_labels"])
            kept = result["without_excluded"]["top_genes"]["out-degree"]["genes"]
            self.assertNotIn(genes[2], [g["gene"] for g in kept])
            self.assertTrue((base / "out/diagnostics.json").exists())


if __name__ == "__main__":
    unittest.main()
