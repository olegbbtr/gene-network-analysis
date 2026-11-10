import argparse
import json
from pathlib import Path

from .diagnostics import diagnose
from .dream5 import prepare
from .figures import write_figures
from .graph import read_graph
from .pipelines import PIPELINES
from .prediction import predict
from .public import build
from .references import (
    GeneMap,
    combine,
    omnipath_pairs,
    read_gene_table,
    trrust_pairs,
    write_edges,
    write_gene_table,
)
from .report import analyze
from .simulate import Settings, simulate
from .structure import balance, locality, modules, signed_pairs, undirected_signs

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
    public = commands.add_parser(
        "public", help="infer one graph per TCGA and GTEx dataset in public.DATASETS"
    )
    public.add_argument("references", type=Path, help="folder holding the downloaded datasets")
    public.add_argument("output", type=Path)
    public.add_argument("--genes", type=Path, default=GENE_TABLE)
    public.add_argument("--size", type=int, default=12000)
    predict_parser = commands.add_parser(
        "predict", help="predict listed genes from expression and network descriptors"
    )
    predict_parser.add_argument("networks", type=Path, help="output folder of `public`")
    predict_parser.add_argument(
        "--list", type=named_path, action="append", required=True, help="NAME=PATH gene list"
    )
    predict_parser.add_argument("--genes", type=Path, default=GENE_TABLE)
    predict_parser.add_argument("--output", type=Path, required=True)
    structure = commands.add_parser(
        "structure", help="genome locality, marker modules and signed triangle balance"
    )
    structure.add_argument("--graph", type=named_path, action="append", default=[])
    structure.add_argument(
        "--signed-edges",
        type=named_path,
        action="append",
        default=[],
        help="NAME=PATH of curated or true signed pairs; triangle balance only",
    )
    structure.add_argument("--genes", type=Path, default=GENE_TABLE)
    structure.add_argument("--output", type=Path, required=True)
    figures = commands.add_parser("figures", help="poster figures from the public-data results")
    figures.add_argument("--results", type=Path, default=Path("results/public"))
    figures.add_argument("--output", type=Path, default=Path("results/public/figures"))
    args = parser.parse_args(argv)
    try:
        if args.command == "figures":
            write_figures(args.results, args.output)
            print(json.dumps(sorted(p.name for p in args.output.glob("*.pdf"))))
            return 0
        if args.command == "public":
            stats = build(args.references, args.genes, args.output, size=args.size)
            print(json.dumps(stats, indent=2))
            return 0
        if args.command == "predict":
            report = predict(args.networks, args.genes, dict(args.list))
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
            print(
                json.dumps(
                    {
                        name: {s: r["auc"] for s, r in result["results"].items()}
                        for name, result in report["lists"].items()
                    },
                    indent=2,
                )
            )
            return 0
        if args.command == "structure":
            genes = read_gene_table(args.genes)
            report: dict = {}
            for name, path in args.graph:
                graph = read_graph(path)
                report[name] = {
                    "locality": locality(graph, genes),
                    "modules": modules(graph, genes),
                    "balance": balance(undirected_signs(graph.adjacency)),
                }
            for name, path in args.signed_edges:
                report[name] = {"balance": balance(signed_pairs(path))}
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps(report, indent=2))
            return 0
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
