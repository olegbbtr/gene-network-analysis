import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.edges import edge_metrics
from gene_network_analysis.graph import edge_weights
from gene_network_analysis.references import (
    GeneMap,
    combine,
    omnipath_pairs,
    read_edges,
    trrust_pairs,
    write_edges,
)

HGNC = (
    "symbol\tprev_symbol\talias_symbol\tensembl_gene_id\tuniprot_ids\n"
    "JUN\t\tAP1\tENSG1\tP1\n"
    "FOS\tOLDFOS\tAP1\tENSG2\tP2\n"
    "TERT\t\t\tENSG3\tP3\n"
    "NOENS\t\t\t\tP4\n"
)


def matrix(edges, size):
    rows, columns, weights = zip(*edges, strict=True)
    return csr_array((np.asarray(weights, dtype=float), (rows, columns)), shape=(size, size))


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.base = Path(self.directory.name)
        (self.base / "hgnc.txt").write_text(HGNC)
        self.genes = GeneMap(self.base / "hgnc.txt")

    def tearDown(self):
        self.directory.cleanup()

    def test_symbols_resolve_through_previous_names_but_not_ambiguous_aliases(self):
        self.assertEqual(self.genes.by_symbol("fos"), "ENSG2")
        self.assertEqual(self.genes.by_symbol("OLDFOS"), "ENSG2")
        self.assertIsNone(self.genes.by_symbol("AP1"))
        self.assertIsNone(self.genes.by_uniprot("P4"))

    def test_complexes_split_into_members_and_signs_combine(self):
        path = self.base / "omnipath.tsv"
        header = "source\ttarget\tconsensus_stimulation\tconsensus_inhibition\tdorothea_level\n"
        path.write_text(
            header
            + "COMPLEX:P1_P2\tP3\tTrue\tFalse\tA\n"
            + "P1\tP3\tFalse\tTrue\tD\n"
            + "P4\tP3\tTrue\tFalse\tA\n"
        )
        pairs, stats = combine(omnipath_pairs(path, self.genes))
        self.assertEqual(pairs, {("ENSG1", "ENSG3"): 0, ("ENSG2", "ENSG3"): 1})
        self.assertEqual(stats, {"input_pairs": 4, "dropped": 1, "pairs": 2})
        pairs, _ = combine(omnipath_pairs(path, self.genes, {"A"}))
        self.assertEqual(pairs[("ENSG1", "ENSG3")], 1)

    def test_trrust_round_trip(self):
        path = self.base / "trrust.tsv"
        path.write_text("JUN\tTERT\tRepression\t1\nJUN\tTERT\tUnknown\t2\nFOS\tJUN\tUnknown\t3\n")
        pairs, _ = combine(trrust_pairs(path, self.genes))
        write_edges(self.base / "out.tsv", pairs)
        self.assertEqual(
            read_edges(self.base / "out.tsv"), [("ENSG1", "ENSG3", -1), ("ENSG2", "ENSG1", 0)]
        )


class EdgeMetricTests(unittest.TestCase):
    def test_lookup_finds_present_and_absent_edges(self):
        graph = matrix([(0, 1, 2.0), (2, 0, -1.0)], 3)
        weights = edge_weights(graph, np.array([0, 2, 1]), np.array([1, 0, 0]))
        np.testing.assert_array_equal(weights, [2.0, -1.0, 0.0])

    def test_recall_precision_and_sign_agreement(self):
        graph = matrix([(0, 1, 1.0), (0, 2, -1.0), (0, 3, 1.0), (4, 0, 1.0)], 5)
        sources, targets = np.array([0, 0, 0, 3]), np.array([1, 2, 4, 0])
        result = edge_metrics(graph, sources, targets, np.array([1, 1, -1, 1]))
        self.assertEqual(result["reference_pairs"], 3)
        self.assertAlmostEqual(result["recall"], 2 / 3)
        self.assertAlmostEqual(result["chance_recall"], 3 / 4)
        self.assertEqual(result["direction"]["forward_only"], 2)
        self.assertEqual(result["sign_agreement"]["agreement"], 0.5)
        by_sign = result["precision_by_graph_sign"]
        self.assertEqual((by_sign["positive"]["edges"], by_sign["positive"]["curated"]), (2, 1))
        self.assertEqual((by_sign["negative"]["edges"], by_sign["negative"]["curated"]), (1, 1))
        self.assertIsNone(edge_metrics(abs(graph), sources, targets, np.ones(4))["sign_agreement"])


if __name__ == "__main__":
    unittest.main()
