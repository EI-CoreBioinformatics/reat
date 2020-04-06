#!/usr/bin/env python3
import pyfaidx
import argparse


__doc__ = """Quick script to determine whether to use long or short gmap."""


def main():

    limit = 2 ** 32 - 1
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-s", "--gsnap", action="store_true", default=False)
    parser.add_argument("genome")
    args = parser.parse_args()

    fa = pyfaidx.Fasta(args.genome)
    if args.gsnap:
        prefix = "gsnap"
    else:
        prefix = "gmap"

    if sum(len(fa.records[_]) for _ in fa.records) < limit:
        print(prefix)
    else:
        print(prefix + "l")

    return


if __name__ == "__main__":
    main()
