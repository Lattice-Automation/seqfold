"""Fold sequences through seqfold's CLI"""

import argparse
import sys
from typing import List

from . import __version__
from .fold import calc_dg, fold, Struct


def run():
    """Entry point for console_scripts.
    """

    args = parse_args(sys.argv[1:])
    structs = fold(args.seq, temp=args.temp)

    if args.verbose:
        # log structure with dot-bracket notation
        desc = ["."] * len(args.seq)
        for s in structs:
            if len(s.ij) == 1:
                i, j = s.ij[0]
                desc[i] = "("
                desc[j] = ")"
        print(args.seq)
        print("".join(desc))

    if args.log:
        # log each structure
        print(Struct.fmt.format("i", "j", "dg", "description"))
        for s in structs:
            print(s)

    # log free energy
    print(round(sum(s.e for s in structs), 2))


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
        "-l", "--log", action="store_true", help="log each structure",
    )
    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__),
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
