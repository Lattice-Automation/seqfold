"""Fold sequences through seqfold's CLI"""

import argparse
import sys
from typing import List

from . import __version__
from .fold import calc_dg, calc_dg_verbose


def run():
    """Entry point for console_scripts.

    Simply fold the DNA and report the minimum free energy
    """

    args = parse_args(sys.argv[1:])

    if args.verbose:
        dg, desc = calc_dg_verbose(args.seq, temp=args.temp)
        print(desc)
        print(dg)
    else:
        print(calc_dg(args.seq, temp=args.temp))


def parse_args(args: List[str]) -> argparse.Namespace:
    """Parse command line parameters

    Created and based on an example from pyscaffold:
    https://github.com/pyscaffold/pyscaffold/blob/master/src/pyscaffold/templates/skeleton.template

    Args:
        args ([str]): List of parameters as strings

    Returns:
        `argparse.Namespace`: command line parameters namespace
    """

    parser = argparse.ArgumentParser(
        description="Predict the minimum free energy (kcal/mol) of a nucleic acid sequence"
    )

    parser.add_argument(
        "seq", type=str, metavar="SEQ", help="nucleic acid sequence to fold",
    )
    parser.add_argument(
        "-t",
        dest="temp",
        type=float,
        metavar="FLOAT",
        default=37.0,
        help="temperature in Celcius",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="log a 2D folding description",
    )
    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__),
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
