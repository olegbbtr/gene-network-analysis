import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array

from gene_network_analysis.prospective import (
    CUTOFF,
    examples,
    gene_bootstrap,
    kind,
    percentiles,
    read_dated,
    score_pairs,
    split,
    top_hits,
)

OLD, NEW = f"PMID:{CUTOFF - 1}", f"PMID:{CUTOFF}"


def gaf(symbol, term, reference, evidence, qualifier="involved_in", aspect="P"):
    fields = ["UniProtKB", "P1", symbol, qualifier, term, reference, evidence, "", aspect]
    return "\t".join(fields + ["name", "", "protein", "taxon:9606", "20200101", "UniProt"])


class ProspectiveTests(unittest.TestCase):
    def test_rows_are_dated_by_their_pubmed_id(self):
        self.assertEqual(kind(OLD, "IEA"), "before")
        self.assertEqual(kind(NEW, "IDA"), "after")
        self.assertEqual(kind(NEW, "IEP"), "later")  # expression pattern is not counted
        self.assertEqual(kind(NEW, "TAS"), "later")
        self.assertEqual(kind("GO_REF:0000033", "IBA"), "undated")
        self.assertEqual(kind(f"{NEW}|{OLD}", "IDA"), "before")

    def test_annotations_are_propagated_and_negations_skipped(self):
        parents = {"GO:1": set(), "GO:2": {"GO:1"}}
        rows = [
            "!gaf-version: 2.2",
            gaf("A", "GO:2", OLD, "IDA"),
            gaf("A", "GO:1", NEW, "IDA", qualifier="NOT|involved_in"),
            gaf("B", "GO:2", NEW, "IMP"),
            gaf("B", "GO:9", NEW, "IMP"),
            gaf("B", "GO:1", NEW, "IDA", aspect="F"),
            gaf('C "quoted', "GO:1", "GO_REF:0000002", "IEA"),
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "goa.gaf.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("\n".join(rows) + "\n")
            dated = read_dated(path, {"A": "a", "B": "b", 'C "QUOTED': "c"}.get, parents)
        self.assertEqual(dated["before"], {"a": {"GO:1", "GO:2"}})
        self.assertEqual(dated["after"], {"b": {"GO:1", "GO:2"}})
        self.assertEqual(dated["undated"], {})  # the lookup is given the symbol as written

    def test_only_unsupported_new_terms_are_positives(self):
        genes = [f"g{i}" for i in range(26)]
        members = {g: {"T"} for g in genes[:20]}
        dated = {
            "before": members,
            "after": {"g20": {"T"}, "g21": {"T"}, "g0": {"T"}},
            "later": {"g22": {"T"}},
            "undated": {"g21": {"T"}, "g23": {"T"}},
        }
        data = split(dated, genes)
        self.assertEqual(data.terms, ["T"])
        self.assertEqual(list(np.flatnonzero(data.positive[:, 0])), [20])
        self.assertEqual(list(np.flatnonzero(data.negative[:, 0])), [24, 25])
        self.assertEqual(list(np.flatnonzero(data.studied)), [0, 20, 21])
        self.assertTrue(data.unannotated[20] and not data.unannotated[0])

    def test_percentile_counts_ties_as_half(self):
        scores = np.array([0.0, 0.0, 1.0, 2.0, 1.0, 0.0])
        positive = np.array([False, False, False, False, True, True])
        negative = np.array([True, True, True, True, False, False])
        np.testing.assert_allclose(percentiles(scores, positive, negative), [0.625, 0.25])

    def test_top_hits_share_tied_places(self):
        scores = np.array([3.0, 1.0, 1.0, 1.0, 0.0])
        positive = np.array([True, True, False, False, True])
        self.assertAlmostEqual(top_hits(scores, positive, 2), 1 + 1 / 3)

    def test_a_gene_linked_to_annotated_genes_ranks_first(self):
        genes = [f"g{i}" for i in range(40)]
        dated = {
            "before": {g: {"T"} for g in genes[:20]},
            "after": {"g20": {"T"}, "g21": {"T"}},
            "later": {},
            "undated": {},
        }
        data = split(dated, genes)
        dense = np.zeros((40, 40))
        dense[20, :5] = dense[:5, 20] = 1  # g20 is linked to annotated genes, g21 to none
        dense[21, 30:35] = dense[30:35, 21] = 1
        scores = csr_array(dense) @ data.before.astype(float)
        scored = score_pairs(scores, data, np.ones(40, dtype=bool))
        np.testing.assert_allclose(scored.percentile, [1.0, 0.5])
        self.assertEqual(list(scored.gene), [20, 21])
        self.assertAlmostEqual(scored.hits, 1.0)  # top 1% of 20 candidates is one gene
        chosen = np.arange(40) != 20
        rows = examples(scored, data, chosen, genes, {"T": "a term"})
        self.assertEqual(len(rows), 1)  # g20 is not chosen; g21 ties with all 18 negatives
        self.assertEqual(
            (rows[0]["symbol"], rows[0]["rank"], rows[0]["candidates"]), ("g21", 10, 19)
        )
        self.assertEqual((rows[0]["name"], rows[0]["annotated_before"]), ("a term", 20))

    def test_bootstrap_resamples_whole_genes(self):
        low, high = gene_bootstrap(np.array([1, 1, 1, 2]), np.array([1.0, 1.0, 1.0, 0.0]), 500)
        self.assertEqual((low, high), (0.0, 1.0))


if __name__ == "__main__":
    unittest.main()
