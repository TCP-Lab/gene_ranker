from gene_ranker.ranker import run_method
from gene_ranker.ranking_methods import RANKING_METHODS
from pathlib import Path
import sys

def bin(args = None):
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--list-methods",
        help="List all available methods and a brief description.",
        action="store_true",
    )

    parser.add_argument(
        "case_matrix", help="Expression Matrix with case samples.", type=Path
    )
    parser.add_argument(
        "control_matrix", help="Expression Matrix with control samples.", type=Path
    )
    parser.add_argument(
        "method",
        help="Ranking method",
        choices=list(RANKING_METHODS.keys()),
    )

    parser.add_argument(
        "--output-file", help="Output file path", type=Path, default=None
    )

    args = parser.parse_args(args)

    if args.list_methods:
        res = "Available methods:\n"
        for key, method in RANKING_METHODS.items():
            res += f"[{key}] - {method.name}: {method.desc}\n"
        sys.stdout.write(res)
        return

    result = run_method(
        case_matrix=args.case_matrix,
        control_matrix=args.control_matrix,
        method=RANKING_METHODS[args.method],
    )

    if args.output_file:
        with args.output_file.open("w+") as stream:
            result.to_csv(stream, index=False)
        return

    result.to_csv(sys.stdout)

