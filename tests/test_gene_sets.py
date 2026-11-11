import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from gene_network_analysis.applications import closest, rewiring_scores, walks
from gene_network_analysis.gene_sets import (
    condensed_pagerank,
    depmap_selective,
    gwas_traits,
    held_out_scores,
    membership,
    neighbour_scores,
    personalised_pagerank,
    read_gmt,
    restrict,
    set_aucs,
)


def chain(size: int) -> csr_array:
    links = csr_array(np.diag(np.ones(size - 1), 1))
    return (links + links.T).tocsr()


class EngineTests(unittest.TestCase):
    def test_neighbour_scores_are_member_shares_per_set(self):
        labels = np.array([[True, False], [False, True], [True, True]])
        scores = neighbour_scores(chain(3), labels)
        np.testing.assert_allclose(scores[1], [1.0, 0.5])
        np.testing.assert_allclose(scores[0], [0.0, 1.0])

    def test_held_out_scores_never_use_the_genes_own_fold(self):
        rng = np.random.default_rng(0)
        labels = rng.random((300, 3)) < 0.3
        groups = np.arange(300) % 5
        links = csr_array(np.eye(300))  # each gene is its own only neighbour
        scores = held_out_scores(links, labels, groups)
        np.testing.assert_array_equal(scores, 0)

    def test_a_connected_set_scores_high(self):
        labels = np.zeros((40, 1), dtype=bool)
        labels[:20] = True
        block = np.zeros((40, 40))
        block[:20, :20] = 1
        np.fill_diagonal(block, 0)
        scores = held_out_scores(csr_array(block), labels, np.arange(40) % 5)
        self.assertEqual(set_aucs(scores, labels)[0], 1.0)

    def test_sets_are_cut_to_the_genes_and_filtered_by_size(self):
        sets = {"small": {"a"}, "kept": {"a", "b", "z"}}
        chosen = restrict(sets, ["a", "b", "c"], smallest=2, largest=5)
        self.assertEqual(chosen, {"kept": {"a", "b"}})
        np.testing.assert_array_equal(
            membership(["a", "b", "c"], chosen)[:, 0], [True, True, False]
        )


class PagerankTests(unittest.TestCase):
    def test_personalised_pagerank_reaches_two_steps_but_favours_one(self):
        labels = np.array([[True], [False], [False], [False]])
        scores = personalised_pagerank(chain(4), labels)[:, 0]
        self.assertGreater(scores[1], scores[2])
        self.assertGreater(scores[2], scores[3])
        self.assertGreater(scores[2], 0)

    def test_held_out_pagerank_hides_the_test_fold(self):
        labels = np.zeros((10, 1), dtype=bool)
        labels[0] = True
        groups = np.array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4])
        isolated = csr_array(np.eye(10))
        scores = held_out_scores(isolated, labels, groups, scorer=personalised_pagerank)
        self.assertEqual(scores[0, 0], 0)

    def test_condensed_pagerank_shares_a_group_score_and_flows_downstream(self):
        # 4 -> {0 <-> 1} -> 2 -> 3, walk seeded at gene 0
        edges = [(4, 0), (0, 1), (1, 0), (1, 2), (2, 3)]
        links = csr_array((np.ones(5), tuple(np.array(edges).T)), shape=(5, 5))
        labels = np.zeros((5, 1), dtype=bool)
        labels[0] = True
        scores = condensed_pagerank(links, labels)[:, 0]
        self.assertEqual(scores[0], scores[1])
        self.assertGreater(scores[2], scores[3])
        self.assertGreater(scores[3], 0)
        self.assertEqual(scores[4], 0)

    def test_acyclic_walks_have_no_cycles(self):
        rng = np.random.default_rng(0)
        dense = rng.normal(size=(30, 30)) * (rng.random((30, 30)) < 0.3)
        np.fill_diagonal(dense, 0)
        variants = walks(csr_array(dense))
        for name in ("acyclic-forward/pagerank", "acyclic-backward/pagerank"):
            count, _ = connected_components(variants[name][0], directed=True, connection="strong")
            self.assertEqual(count, 30)


class ReaderTests(unittest.TestCase):
    def test_gmt_in_a_zip_maps_symbols(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sets.gmt.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("sets.gmt", "Pathway A\tR-1\tAAA\tBBB\tNOPE\n")
            self.assertEqual(
                read_gmt(path, {"AAA": "E1", "BBB": "E2"}), {"Pathway A": {"E1", "E2"}}
            )

    def test_depmap_selects_lineage_selective_non_common_dependencies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            models = [f"ACH-{i}" for i in range(12)]
            lineage = ["Lung"] * 6 + ["Skin"] * 6
            (root / "Model.csv").write_text(
                "ModelID,OncotreeLineage\n"
                + "".join(f"{m},{x}\n" for m, x in zip(models, lineage, strict=True))
            )
            (root / "common.csv").write_text("Essentials\nCOMMON (3)\n")
            rows = ["," + ",".join(["SELECTIVE (1)", "EVERYWHERE (2)", "COMMON (3)"])]
            for x in lineage:
                selective = -1.0 if x == "Lung" else 0.0
                rows.append(f"m,{selective},-1.0,{selective}")
            rows = [rows[0]] + [r.replace("m,", f"{m},", 1) for m, r in zip(models, rows[1:])]
            (root / "effects.csv").write_text("\n".join(rows) + "\n")
            symbols = {"SELECTIVE": "E1", "EVERYWHERE": "E2", "COMMON": "E3"}
            sets, counts = depmap_selective(
                root / "effects.csv",
                root / "Model.csv",
                root / "common.csv",
                symbols,
                {"Lung": "lung", "Skin": "skin"},
            )
            self.assertEqual(sets["lung"], {"E1"})
            self.assertEqual(sets["skin"], set())
            self.assertEqual(counts["lung"], 6)

    def test_gwas_keeps_significant_mapped_genes_of_chosen_traits(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gwas.tsv"
            path.write_text(
                "MAPPED_GENE\tP-VALUE\tMAPPED_TRAIT\n"
                "AAA - BBB\t1e-9\tasthma, body height\n"
                "CCC\t1e-5\tasthma\n"
                "DDD\t1e-10\tother trait\n"
            )
            sets = gwas_traits(
                path, {s: s.lower() for s in ("AAA", "BBB", "CCC", "DDD")}, {"asthma"}
            )
            self.assertEqual(sets, {"asthma": {"aaa", "bbb"}})


class RewiringTests(unittest.TestCase):
    def test_closest_partners_and_rewiring(self):
        healthy = csr_array(np.array([[0, 5.0, 1.0], [0, 0, 0], [0, 0, 0]]))
        tumour = csr_array(np.array([[0, 0, -4.0], [0, 0, 0], [0, 0, 0]]))
        self.assertEqual(closest(healthy, k=1)[0], {1})
        change = rewiring_scores(healthy, tumour)
        self.assertEqual(change[0], 0.5)  # {1, 2} -> {2}
        self.assertEqual(change[1], 1.0)  # {0} -> nothing in common


if __name__ == "__main__":
    unittest.main()
