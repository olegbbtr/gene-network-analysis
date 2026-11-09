import argparse
import json
from pathlib import Path

from .diagnostics import diagnose
from .dream5 import prepare
from .pipelines import PIPELINES
from .references import (
    GeneMap,
    combine,
    omnipath_pairs,
    trrust_pairs,
    write_edges,
    write_gene_table,
)
from .report import analyze
from .simulate import Settings, simulate

REFERENCES = {"tfs": Path("data/tfs.txt"), "lambert-2018": Path("data/lambert-2018-tfs.txt")}
GENE_TABLE = Path("data/hgnc-genes.tsv")
EDGE_REFERENCES = {
    name: Path(f"data/edges/{name}.tsv") for name in ("collectri", "dorothea-abc", "trrust")
}


def named_path(value: str) -> tuple[str, Path]:
    name, separator, path = value.partition("=")
    if not separator or not name or not path:
        raise argparse.ArgumentTypeError("expected NAME=PATH")
    return name, Path(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Rank regulators in a gene regulatory network")
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run pipelines and write the report")
    run.add_argument("--graph", type=Path, default=Path("data/graph.txt"))
    run.add_argument(
        "--reference",
        type=named_path,
        action="append",
        help="NAME=PATH of a TF list (repeatable); default: the lists in data/",
    )
    run.add_argument(
        "--edges",
        type=named_path,
        action="append",
        help="NAME=PATH of curated TF->target pairs (repeatable); default: data/edges/",
    )
    run.add_argument("--output", type=Path, default=Path("results"))
    run.add_argument(
        "--pipeline",
        action="append",
        choices=[p.name for p in PIPELINES],
        help="run this pipeline (repeatable); default: every reported pipeline",
    )
    run.add_argument("--top-k", type=int, default=100)
    convert = commands.add_parser(
        "import-edges", help="convert a downloaded TF->target reference to Ensembl IDs"
    )
    convert.add_argument("format", choices=["collectri", "dorothea-abc", "trrust"])
    convert.add_argument("source", type=Path)
    convert.add_argument("--hgnc", type=Path, required=True)
    convert.add_argument("output", type=Path)
    genes = commands.add_parser(
        "import-genes", help="reduce the HGNC complete set to the gene table diagnostics read"
    )
    genes.add_argument("hgnc", type=Path)
    genes.add_argument("output", type=Path)
    diagnostics = commands.add_parser(
        "diagnostics",
        help="sender selection, TF partners, top genes, and rankings without the KRAB family",
    )
    diagnostics.add_argument("--graph", type=Path, default=Path("data/graph.txt"))
    diagnostics.add_argument("--genes", type=Path, default=GENE_TABLE)
    diagnostics.add_argument("--output", type=Path, default=Path("results"))
    sim = commands.add_parser("simulate", help="write a simulated network with known truth")
    sim.add_argument("output", type=Path)
    sim.add_argument("--seed", type=int, default=0)
    sim.add_argument("--inference", choices=["correlation", "partial"], default="correlation")
    dream = commands.add_parser(
        "dream5", help="infer a graph from a DREAM5 dataset and write its truth files"
    )
    dream.add_argument("network", help="DREAM5 network number, used in gene names")
    dream.add_argument("expression", type=Path)
    dream.add_argument("tfs", type=Path)
    dream.add_argument("gold", type=Path)
    dream.add_argument("output", type=Path)
    dream.add_argument("--signed", type=Path, help="signed gold standard, if there is one")
    dream.add_argument("--inference", choices=["correlation", "partial"], default="correlation")
    args = parser.parse_args(argv)
    try:
        if args.command == "dream5":
            stats = prepare(
                args.network,
                args.expression,
                args.tfs,
                args.gold,
                args.output,
                signed=args.signed,
                inference=args.inference,
            )
            print(json.dumps(stats))
            return 0
        if args.command == "simulate":
            print(
                json.dumps(
                    simulate(Settings(seed=args.seed, inference=args.inference), args.output)
                )
            )
            return 0
        if args.command == "import-genes":
            print(json.dumps({"genes": write_gene_table(args.hgnc, args.output)}))
            return 0
        if args.command == "diagnostics":
            report = diagnose(args.graph, REFERENCES, args.genes, args.output)
            summary = {
                name: {
                    "sender_odds_ratio": result["sender_selection"]["odds_ratio"],
                    **{
                        f"auc_{subset}": {
                            row["pipeline"]: round(row["roc_auc"], 3)
                            for row in result[subset]["comparison"]
                        }
                        for subset in ("all_senders", "without_excluded")
                    },
                }
                for name, result in report["references"].items()
            }
            print(json.dumps(summary, indent=2))
            return 0
        if args.command == "import-edges":
            genes = GeneMap(args.hgnc)
            pairs = {
                "collectri": lambda: omnipath_pairs(args.source, genes),
                "dorothea-abc": lambda: omnipath_pairs(args.source, genes, {"A", "B", "C"}),
                "trrust": lambda: trrust_pairs(args.source, genes),
            }[args.format]()
            edges, stats = combine(pairs)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            write_edges(args.output, edges)
            print(json.dumps(stats))
            return 0
        edges = (
            dict(args.edges)
            if args.edges
            else {name: path for name, path in EDGE_REFERENCES.items() if path.exists()}
        )
        report = analyze(
            args.graph,
            dict(args.reference) if args.reference else REFERENCES,
            args.output,
            edge_references=edges,
            pipelines=args.pipeline,
            top_k=args.top_k,
        )
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        parser.exit(1, f"error: {error}\n")
    summary = {
        name: [{k: row[k] for k in ("pipeline", "roc_auc")} for row in evaluation["comparison"]]
        for name, evaluation in report["references"].items()
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
