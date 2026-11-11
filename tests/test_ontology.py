import tempfile
import unittest
from pathlib import Path

from gene_network_analysis.ontology import ancestors, percentiles, read_ontology

OBO = """[Term]
id: GO:0000001
name: process
namespace: biological_process

[Term]
id: GO:0000002
name: child
namespace: biological_process
is_a: GO:0000001 ! process

[Term]
id: GO:0000003
name: part
namespace: biological_process
relationship: part_of GO:0000002 ! child

[Term]
id: GO:0000004
name: component
namespace: cellular_component
"""


class OntologyTests(unittest.TestCase):
    def test_ancestors_follow_is_a_and_part_of_within_processes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "go.obo"
            path.write_text(OBO)
            parents = read_ontology(path)
        self.assertNotIn("GO:0000004", parents)
        found = ancestors("GO:0000003", parents, {})
        self.assertEqual(found, {"GO:0000001", "GO:0000002", "GO:0000003"})

    def test_a_guess_sharing_specific_terms_ranks_first(self):
        # 100 background genes with a common term; pathway "match" shares gene x's rare term
        terms = {f"b{i}": {"common"} for i in range(100)}
        terms |= {"m1": {"rare", "common"}, "m2": {"rare"}, "x": {"rare", "common"}}
        pathways = {
            "match": {"m1", "m2"},
            "other": {"b1", "b2"},
            "third": {"b3", "b4"},
        }
        ranks = percentiles({"x": "match", "y": "match"}, terms, pathways)
        self.assertEqual(ranks, {"x": 1.0})  # y has no terms; common is not specific


if __name__ == "__main__":
    unittest.main()
