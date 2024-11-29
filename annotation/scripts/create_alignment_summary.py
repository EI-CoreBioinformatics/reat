#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create alignment summary from multi_genome_compare --detailed tab file
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"

# import libraries
import argparse
from pathlib import Path
from collections import defaultdict

class CreateAlignmentSummary:
    bins = {
        "0-10": [],
        "11-20": [],
        "21-30": [],
        "31-40": [],
        "41-50": [],
        "51-60": [],
        "61-70": [],
        "71-80": [],
        "81-90": [],
        "91-100": [],
    }

    def __init__(self, args):
        self.tab_files = args.tab_file
        self.header = args.header
        self.header_index = None
        self.bins = {key: [] for key in CreateAlignmentSummary.bins}

    def check_and_update_bin(self, value):
        for bin_range in self.bins:
            start, end = map(float, bin_range.split("-"))
            if start <= value <= end:
                self.bins[bin_range].append(value)
                break

    def run(self, tab_file):
        with open(tab_file) as f:
            for line in f:
                if line.startswith("Group"):
                    if self.header in line:
                        header = line.strip().split("\t")
                        self.header_index = header.index(self.header)
                        continue
                    else:
                        print(
                            f"Error: {self.header} not found in {tab_file}. Please provide a valid input file (see help)."
                        )
                        break
                if self.header_index:
                    line = line.strip().split("\t")
                    value = float(line[self.header_index])
                    self.check_and_update_bin(value)
        return self.bins


def main():

    parser = argparse.ArgumentParser(
        description="Script to create alignment summary from multi_genome_compare --detailed tab file"
    )
    parser.add_argument(
        "tab_file",
        nargs="+",
        help="Provide multi_genome_compare --detailed tab file(s)",
    )
    parser.add_argument(
        "--header",
        type=str,
        default="Precision(Junction)",
        help="Provide the header to plot (default: %(default)s)",
    )
    args = parser.parse_args()

    # store each input file bin values in a dictionary
    collective_bin_values = defaultdict(dict)
    for tab_file in args.tab_file:
        bin_values = CreateAlignmentSummary(args).run(tab_file)
        if tab_file not in collective_bin_values:
            collective_bin_values[tab_file] = bin_values

    # print collective_bin_values into a table
    header = f"Bin {args.header}"
    for tab_file in args.tab_file:
        header += f"\t{Path(tab_file).name}"
    print(header)
    for bin_range in CreateAlignmentSummary.bins:
        row = bin_range
        for tab_file in args.tab_file:
            row += f"\t{len(collective_bin_values[tab_file][bin_range])}"
        print(row)



if __name__ == "__main__":
    main()
