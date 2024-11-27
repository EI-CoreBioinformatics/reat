#!/usr/bin/env python3

# This script reads the GFF file from the alignments and creates a groups.txt file. At this point, the GFF file should have the mRNA lines with the protein_alias attribute and the ID attribute and that is all that is needed to create the groups.txt file for the multi genome compare step. This script can also be used generically.

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"


from collections import defaultdict
import argparse
import re
import os
import sys

script = os.path.basename(sys.argv[0])


def main():
    parser = argparse.ArgumentParser(
        prog=script,
        description="Create groups.txt file for multi genome compare from the alignments GFF file",
    )
    parser.add_argument(
        "--alignments_gff",
        type=argparse.FileType("r"),
        help="Input GFF from the alignments",
    )
    parser.add_argument(
        "--feature",
        type=str,
        default="mRNA",
        help="Provide the feature type in the GFF3 file which contains the '--protein_alias'. For example, if the 'mRNA' type contains the original protein name in 'protein_alias' attribute, then provide 'mRNA' for this option (default: %(default)s).",
    )
    parser.add_argument(
        "--protein_alias",
        type=str,
        default="protein_alias",
        help="Provide the GFF3 attribute field in the '--feature' type which can be used to link the actual input protein_fasta to the alignment gff_file. For example, if the 'mRNA' attribute has the original protein name in 'protein_alias=ABCD_trans_v1;' attribute, then provide 'protein_alias' for this option (default: %(default)s).",
    )
    # get the output file name groups.txt
    parser.add_argument(
        "--output",
        type=str,
        default="groups.txt",
        help="Output file name to write the groups.txt file (default: %(default)s).",
    )

    args = parser.parse_args()
    alignments_filename = args.alignments_gff.name
    output_filename = args.output
    feature_type = args.feature
    protein_alias = args.protein_alias
    pattern = protein_alias + "=([^;]+)"
    model_to_group = defaultdict(list)
    with open(alignments_filename) as gff:
        for line in gff:
            line = line.rstrip("\n")
            if re.match(r"^\s*$", line) or line.startswith("#"):
                continue
            x = line.split("\t")
            if len(x) != 9:
                continue
            feature = x[2]
            if feature.lower() != feature_type.lower():
                continue
            attribute = x[8]
            align_id = re.search(r"ID=([^;]+)", attribute).group(1)
            protein_alias_id = re.search(pattern, attribute).group(1)

            if not align_id or not protein_alias_id:
                print(
                    "Error: align id or protein alias not found in the mRNA line: ",
                    line,
                )
                exit(1)
            model_to_group[protein_alias_id].append(align_id)

    groupid = 0
    with open(output_filename, "w") as gf:
        for k, v in model_to_group.items():
            for a in v:
                groupid += 1
                print(k, groupid, file=gf)
                print(a, groupid, file=gf)


if __name__ == "__main__":
    main()
