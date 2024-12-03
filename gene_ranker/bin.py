import argparse
import logging
import sys
from pathlib import Path

from gene_ranker import __version__
from gene_ranker.methods import RANKING_METHODS
from gene_ranker.ranker import run_method

log = logging.getLogger(__name__)


class ListMethodsAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        res = "Available methods:\n"
        for key, method in RANKING_METHODS.items():
            res += f"\t'{key}' - {method.name}: {method.desc}\n"
        sys.stdout.write(res)
        parser.exit()
        return


class PrintVersionAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(__version__)
        parser.exit()
        return


def bin(args=None):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--list-methods",
        help="List all available methods and a brief description.",
        nargs=0,
        action=ListMethodsAction,
    )

    parser.add_argument(
        "--version",
        "-v",
        help="Print version and exit",
        nargs=0,
        action=PrintVersionAction,
    )

    parser.add_argument(
        "case_matrix",
        help="Expression Matrix with log2 expression of case samples.",
        type=Path,
    )
    parser.add_argument(
        "control_matrix",
        help="Expression Matrix with log2 expression of control samples.",
        type=Path,
    )

    parser.add_argument(
        "--output-file", help="Output file path", type=Path, default=None
    )

    parser.add_argument(
        "--id-col",
        help="Name of shared ID comlumn between files",
        type=str,
        default="gene_id",
    )

    general_args = [x.dest for x in parser._actions] + ["method"]

    # Add the individual parsers
    subparsers = parser.add_subparsers(dest="method")
    for key, values in RANKING_METHODS.items():
        if not values.parser:
            log.warn(f"Parser {key} has no parser. Skipping...")
            continue
        subparsers.add_parser(key, parents=[values.parser], add_help=False)

    args = parser.parse_args(args)
    extra_args = {k: v for k, v in vars(args).items() if k not in general_args}

    result = run_method(
        case_matrix=args.case_matrix,
        control_matrix=args.control_matrix,
        method=RANKING_METHODS[args.method],
        shared_col=args.id_col,
        extra_args=extra_args,
    )

    log.info(
        "Writing output to {}".format(
            args.output_file if args.output_file else "stdout"
        )
    )

    out_stream = args.output_file.open("w+") if args.output_file else sys.stdout
    result.to_csv(out_stream, index=False)
