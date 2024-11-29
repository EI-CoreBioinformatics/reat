#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to reformat miniprot GFF from mRNA->CDS to gene->mRNA>[exon|CDS] OR match->match_part format

# Friday, 30 September 2022, 10:28
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"


# import libraries
import argparse
from itertools import count
import os
import re
import logging
import sys
from collections import defaultdict
import pyfaidx

# get script name
script = os.path.basename(sys.argv[0])

# get the GFF3 attributes
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTE = range(9)

logging.basicConfig(format="%(message)s", datefmt="%d-%b-%y %H:%M:%S")

# process file


class Miniprot2GFF3:
    @staticmethod
    def get_id(file, line, attribute, field):
        """
        Extract query from GFF3 attribute column
        """
        pattern = field + "=([^;]+)"
        id_field = str()
        try:
            # Check for GFF3 file
            id_search = re.search(pattern, attribute)
            if id_search:
                id_field = id_search.group(1)
            else:
                logging.error(
                    f"Error: Cannot extract attribute type '{field}=' from the file '{file}'. Please check if the attribute is present in line:'\n{line}\n"
                )
        except AttributeError as err:
            logging.error(
                f"Error: {err}. '{field}' cannot be extracted from the file '{file}' line below\n{line}"
            )

        if not id_field:
            logging.error(
                f"Error: Cannot find {field} from the file '{file}', exiting.."
            )
            sys.exit(1)

        return id_field

    def __init__(self, args):
        self.args = args
        logging.warning(f"# INFO: Loading FAI index for '{args.protein_fasta}'")
        self.fai = pyfaidx.Fasta(args.protein_fasta)
        logging.warning(f"# INFO: Finished loading FAI index")
        self.gff_file = args.gff_file
        self.source = args.source
        self.identity = args.identity
        self.coverage = args.coverage
        self.print_match = args.print_match

    def process_gff(self):
        rtinfo = defaultdict()  # real transcripts
        atinfo = defaultdict()  # aligned transcripts
        frtinfo = defaultdict()  # filtered real transcripts
        fatinfo = defaultdict()  # filtered aligned transcripts
        print("##gff-version 3")
        counter = 1
        add_directive = False
        with open(self.gff_file, "r") as filehandle:
            for line in filehandle:
                line = line.rstrip("\n")
                if re.match(r"^\s*$", line) or line.startswith("#"):
                    continue
                else:
                    x = line.split("\t")
                    if len(x) != 9:
                        logging.error(
                            f"Error: Not a standard GFF3 9 column line, see below:\n{line}"
                        )

                    # change source if requested
                    if self.source:
                        x[SOURCE] = self.source

                    # start processing main types
                    if x[TYPE] in ["mRNA"]:

                        mrna_id = Miniprot2GFF3.get_id(
                            self.gff_file, line, x[ATTRIBUTE], "ID"
                        )
                        mIdentity = (
                            float(
                                Miniprot2GFF3.get_id(
                                    self.gff_file, line, x[ATTRIBUTE], "Identity"
                                )
                            )
                            * 100
                        )
                        mIdentity = f"{mIdentity:.2f}"
                        mPositive = (
                            float(
                                Miniprot2GFF3.get_id(
                                    self.gff_file, line, x[ATTRIBUTE], "Positive"
                                )
                            )
                            * 100
                        )
                        mPositive = f"{mPositive:.2f}"
                        mStopCodon = (
                            Miniprot2GFF3.get_id(
                                self.gff_file, line, x[ATTRIBUTE], "StopCodon"
                            )
                            if "StopCodon" in line
                            else ""
                        )
                        mFrameshift = (
                            Miniprot2GFF3.get_id(
                                self.gff_file, line, x[ATTRIBUTE], "Frameshift"
                            )
                            if "Frameshift" in line
                            else ""
                        )
                        mTarget = Miniprot2GFF3.get_id(
                            self.gff_file, line, x[ATTRIBUTE], "Target"
                        )
                        tid, tstart, tend = mTarget.split()

                        if tid not in self.fai:
                            logging.error(
                                f"Error: Transcript '{tid}' is not in the input fasat file, exiting.."
                            )
                            sys.exit(1)

                        # store aligned transcripts
                        if mrna_id not in atinfo:
                            atinfo[mrna_id] = 1
                        # store real transcripts
                        if tid not in rtinfo:
                            rtinfo[tid] = 1

                        aln = f"{tstart}-{tend}"
                        alnLen = int(tend) - int(tstart) + 1
                        tcov = (int(tend) - int(tstart) + 1) / len(self.fai[tid]) * 100
                        tcov = f"{tcov:.2f}"
                        # print(tid, len(self.fai[tid]), tstart, tend, tcov)

                        # new ID
                        new_mrna_id = mrna_id + "_" + tid
                        new_gene_id = new_mrna_id + "-gene"

                        # mRNA/match note
                        mnote = f"aln:{aln}|alnLen:{alnLen}|len:{len(self.fai[tid])}|cov:{tcov}|id:{mIdentity}|Positive:{mPositive}"
                        mnote += f"|StopCodon:{mStopCodon}" if mStopCodon else ""
                        mnote += f"|Frameshift:{mFrameshift}" if mFrameshift else ""
                        mnote += f"|Target:{mTarget}"

                        x[ATTRIBUTE] = x[ATTRIBUTE].replace("ID=", "Prev_ID=").replace("Identity=", "mp_Identity=")
                        gattr = f"ID={new_gene_id};Name={new_gene_id};Note={tid}"
                        mattr = f"ID={new_mrna_id};Parent={new_gene_id};Name={new_mrna_id};Note={mnote};protein_alias={tid};cov={tcov};identity={mIdentity};{x[ATTRIBUTE]}"
                        mmattr = f"ID={new_mrna_id};Name={new_mrna_id};Note={mnote};protein_alias={tid};cov={tcov};identity={mIdentity};{x[ATTRIBUTE]}"
                        # print(gattr) # gene attribute
                        # print(mattr) # mRNA attribute
                        # print(mmattr) # match attribute

                        # apply filter
                        if float(mIdentity) >= float(self.identity) and float(
                            tcov
                        ) >= float(self.coverage):

                            if add_directive:
                                print("###")
                            add_directive = True
                            # store filtered aligned transcripts
                            if mrna_id not in fatinfo:
                                fatinfo[mrna_id] = new_mrna_id
                            # store real filtered aligned transcripts
                            if tid not in frtinfo:
                                frtinfo[tid] = 1

                            if self.print_match:
                                x[TYPE] = "match"
                                new_line = x[:8]
                                new_line.append(mmattr)
                                print("\t".join(new_line))
                            else:
                                x[TYPE] = "gene"
                                new_line = x[:8]
                                new_line.append(gattr)
                                print("\t".join(new_line))
                                x[TYPE] = "mRNA"
                                new_line = x[:8]
                                new_line.append(mattr)
                                print("\t".join(new_line))

                            counter = 1

                    # all other features
                    elif x[TYPE] in ["CDS"]:
                        cParent = Miniprot2GFF3.get_id(
                            self.gff_file, line, x[ATTRIBUTE], "Parent"
                        )
                        if cParent in fatinfo:
                            new_cParent = fatinfo[cParent]
                            x[ATTRIBUTE] = x[ATTRIBUTE].replace(
                                "Parent=", "Prev_Parent="
                            )
                            eattr = f"ID={new_cParent}.exon{counter};Parent={new_cParent};{x[ATTRIBUTE]}"
                            cattr = f"ID={new_cParent}.CDS{counter};Parent={new_cParent};{x[ATTRIBUTE]}"
                            if self.print_match:
                                x[TYPE] = "match_part"
                                new_line = x[:8]
                                new_line.append(eattr)
                                print("\t".join(new_line))
                            else:
                                x[TYPE] = "exon"
                                new_line = x[:8]
                                new_line.append(eattr)
                                print("\t".join(new_line))
                                x[TYPE] = "CDS"
                                new_line = x[:8]
                                new_line.append(cattr)
                                print("\t".join(new_line))
                            counter += 1
        print("###")

        # Summarise
        logging.warning("# INFO: Summary stats")
        logging.warning(f"# INFO: Total transcripts:\t{len(self.fai.keys())}")
        logging.warning(
            f"# INFO: Total transcripts aligned:\t{len(rtinfo.keys())} [{(100 * len(rtinfo.keys()) / len(self.fai.keys())):.2f}%]"
        )
        logging.warning(f"# INFO: Total transcript alignments:\t{len(atinfo.keys())}")
        logging.warning(
            f"# INFO: Summary stats pass filter [coverage >= {self.coverage} and identity >= {self.identity}]"
        )
        logging.warning(
            f"# INFO: Total transcripts aligned (pass filter):\t{len(frtinfo.keys())} [{(100 * len(frtinfo.keys()) / len(self.fai.keys())):.2f}%]"
        )
        logging.warning(
            f"# INFO: Total transcript alignments (pass filter):\t{len(fatinfo.keys())}"
        )
        logging.warning(f"# INFO: Analysis complete")

    def run(self):
        logging.warning(f"# INFO: Processing input file '{self.gff_file}'")
        self.process_gff()


def main():

    parser = argparse.ArgumentParser(
        description="Script to reformat miniprot GFF from mRNA->CDS to gene->mRNA>[exon|CDS] OR match->match_part format"
    )
    parser.add_argument("protein_fasta", help="Provide query protein fasta file")
    parser.add_argument("gff_file", help="Provide GFF file")
    parser.add_argument(
        "--source",
        type=str,
        default=None,
        help="Provide new source for GFF3 (default: %(default)s)",
    )
    parser.add_argument(
        "--identity",
        type=int,
        default=0,
        help="identify cutoff [INT 0 - 100] (default: %(default)s)",
    )
    parser.add_argument(
        "--coverage",
        type=int,
        default=0,
        help="coverage cutoff [INT 0 - 100] (default: %(default)s)",
    )
    parser.add_argument(
        "--print_match",
        action="store_true",
        help="print output GFF3 file as match->match_part for loading to the browser (default: %(default)s)",
    )
    args = parser.parse_args()
    Miniprot2GFF3(args).run()


if __name__ == "__main__":
    main()
