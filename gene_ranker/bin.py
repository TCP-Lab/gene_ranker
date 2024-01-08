from gene_ranker.ranker import run_method
from gene_ranker.ranking_methods import RANKING_METHODS
from pathlib import Path
import sys
import argparse
import logging

log = logging.getLogger("archon")

class ListMethodsAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        res = "Available methods:\n"
        for key, method in RANKING_METHODS.items():
            res += f"\t'{key}' - {method.name}: {method.desc}\n"
        sys.stdout.write(res)
        parser.exit()
        return

def bin(args = None):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--list-methods",
        help="List all available methods and a brief description.",
        nargs=0,
        action=ListMethodsAction,
    )

    parser.add_argument(
        "case_matrix", help="Expression Matrix with log2 expression of case samples.", type=Path
    )
    parser.add_argument(
        "control_matrix", help="Expression Matrix with log2 expression of control samples.", type=Path
    )

    # TODO: Edit this to add the custom subparsers from each ranking method.
    parser.add_argument(
        "method",
        help="Ranking method",
        choices=list(RANKING_METHODS.keys()),
    )

    parser.add_argument(
        "--output-file", help="Output file path", type=Path, default=None
    )

    parser.add_argument(
        "--id-col", help="Name of shared ID comlumn between files", type=str, default="gene_id"
    )
    args = parser.parse_args(args)

    result = run_method(
        case_matrix=args.case_matrix,
        control_matrix=args.control_matrix,
        method=RANKING_METHODS[args.method],
        shared_col=args.id_col
    )
    
    log.info("Writing output to {}".format(args.output_file if args.output_file else "stdout"))
    
    out_stream = args.output_file.open("w+") if args.output_file else sys.stdout
    result.to_csv(out_stream, index=False)

