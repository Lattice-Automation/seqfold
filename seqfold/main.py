"""Fold sequences through seqfold's CLI"""

import argparse
import sys
from typing import List

from . import __version__
from .fold import fold, dot_bracket, Struct


def run():
    """Entry point for console_scripts."""

    args = parse_args(sys.argv[1:])
    structs = fold(args.seq, temp=args.temp)

    if args.verbose or args.dot_bracket:
        # log structure with dot-bracket notation
        print(args.seq)
        print(dot_bracket(structs))

    if args.log or args.sub_structures:
        # log each structure
        print(Struct.fmt.format("i", "j", "ddg", "description"))
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
        "seq", type=str, metavar="SEQ", help="nucleic acid sequence to fold"
    )
    parser.add_argument(
        "-t",
        "--celcius",
        dest="temp",
        type=float,
        metavar="FLOAT",
        default=37.0,
        help="temperature in Celsius",
    )

    # I'm hiding these flags because 1. they're an embarassing testament
    # to me not really knowing at the time what the --verbose flag is for
    # 2. both point out the glaring lack of logging in this tool
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-d",
        "--dot-bracket",
        action="store_true",
        help="write a dot-bracket of the MFE folding to stdout",
    )
    parser.add_argument(
        "-r",
        "--sub-structures",
        action="store_true",
        help="write each substructure of the MFE folding to stdout",
    )

    parser.add_argument(
        "--version", action="version", version="seqfold {ver}".format(ver=__version__)
    )

    return parser.parse_args(args)


if __name__ == "__main__":
    run()
