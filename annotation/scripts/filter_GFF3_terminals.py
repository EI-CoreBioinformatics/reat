#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to filter GFF3 terminals. Script expects GFF3 input in
gene->mRNA>[exon|CDS] format. CDS coordinates are used to filter
the terminals and intron coordinates are computed of the CDS coordinates

# Friday, 30 September 2022, 10:28
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"


# import libraries
import argparse
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

logging.basicConfig(
    format="%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.DEBUG,
)

# CONTANT FOR MINIMUM
MIN_LENGTH = 1_000_000_000_000


class FilterGFF3Terminals:
    @staticmethod
    def get_id(input_file, line, attribute, field, required=True):
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
                if required:
                    raise ValueError(
                        f"Error: Cannot extract attribute type '{field}=' from the file '{input_file.name}'. Please check if the attribute is present in line:'\n{line}\n"
                    )
                else:
                    id_field = str()
        except AttributeError as err:
            if required:
                raise ValueError(
                    f"Error: {err}. '{field}' cannot be extracted from the file '{input_file.name}' line below\n{line}"
                )
            else:
                id_field = str()
        if required:
            if not id_field:
                raise ValueError(
                    f"Error: Cannot find {field} from the file '{input_file.name}', exiting.."
                )
        return id_field

    @staticmethod
    def compute_coords_sum(coords):
        # logging.info("Compute coordinate sum ...")
        sum_size = 0
        for start, end, info in coords:
            sum_size += end - start + 1
        return sum_size

    @staticmethod
    def get_start_end(coords):
        """
        Compute start and end for given coords
        """
        # logging.info("Compute start and end ...")
        starts = list()
        ends = list()
        for start, end, info in coords:
            starts.append(start)
            ends.append(end)
        # print(f"start:{min(starts)}, end:{max(ends)}")
        return min(starts), max(ends)

    @staticmethod
    def get_terminal_len_multi(mrna_id, strand, coords_info):
        """
        Get terminal length for multi-CDS, multi-intron models.
        It will error if coords_info files are multi-exonic
        [IndexError: list index out of range]. We are not planning
        to filter the models on exons yet.
        """
        term5, term3 = 0, 0
        if mrna_id in coords_info:
            # DEBUG
            # logging.info("*MULTI*")
            # logging.info(coords_info[mrna_id])

            # get five prime terminal
            x, y, _ = (
                coords_info[mrna_id][0]
                if strand in ["+", "."]
                else coords_info[mrna_id][-1]
            )
            term5 = (y - x) + 1
            # DEBUG
            # logging.info("*MULTI*term5*")
            # logging.info(f"{term5} = ({y} - {x}) + 1 -- {(y - x) + 1}")

            # get three prime terminal
            x, y, _ = (
                coords_info[mrna_id][-1]
                if strand in ["+", "."]
                else coords_info[mrna_id][0]
            )
            term3 = (y - x) + 1
            # DEBUG
            # logging.info("*MULTI*term3*")
            # logging.info(f"{term3} = ({y} - {x}) + 1 -- {(y - x) + 1}")
        else:
            raise ValueError(
                f"Error: Transcript '{mrna_id}' not present in child attribute '{coords_info}'"
            )
        return term5, term3

    def __init__(self, args):
        self.args = args
        logging.info(f"Loading FAI index for '{args.protein_fasta}'")
        self.fai = pyfaidx.Fasta(args.protein_fasta)
        logging.info(f"Finished loading FAI index")
        self.gff_file = args.gff_file
        self.protein_alias = args.protein_alias
        self.source = args.source
        self.term5c_len = args.term5c_len
        self.term3c_len = args.term3c_len
        self.term5i_len = args.term5i_len
        self.term3i_len = args.term3i_len
        self.clip = args.clip
        # self.identity = args.identity
        # self.coverage = args.coverage
        # self.print_match = args.print_match
        self.gene_info = defaultdict(dict)
        self.mrna_info = defaultdict(dict)
        self.exon_info = defaultdict(list)
        self.cds_info = defaultdict(list)
        self.intron_info = defaultdict(list)

    def process_gff(self):
        for line in self.gff_file:
            line = line.rstrip("\n")
            if re.match(r"^\s*$", line) or line.startswith("#"):
                continue
            x = line.split("\t")
            if len(x) != 9:
                raise ValueError(
                    f"Error: Not a standard GFF3 9 column line, see below:\n{line}"
                )

            start, end = (
                (int(x[END]), int(x[START]))
                if int(x[START]) > int(x[END])
                else (int(x[START]), int(x[END]))
            )

            if x[TYPE].lower() in ["gene"]:
                gene_id = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], "ID"
                )
                if gene_id in self.gene_info:
                    raise ValueError(
                        f"Error: Duplicate gene encountered '{gene_id}'. Please make sure that the GFF3 file is sorted using genometools. Exiting..."
                    )
                else:
                    self.gene_info[gene_id]["start"] = start
                    self.gene_info[gene_id]["end"] = end
                    self.gene_info[gene_id]["line"] = line

            elif x[TYPE].lower() in ["mrna"]:
                mrna_id = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], "ID"
                )
                parent_id = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], "Parent"
                )
                protein_alias = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], self.protein_alias
                )
                if mrna_id in self.mrna_info:
                    raise ValueError(
                        f"Error: Duplicate mRNA encountered '{mrna_id}'. Please make sure that the GFF3 file is sorted using genometools. Exiting..."
                    )
                else:
                    self.mrna_info[mrna_id]["parent_id"] = parent_id
                    self.mrna_info[mrna_id]["protein_alias"] = protein_alias
                    self.mrna_info[mrna_id]["start"] = start
                    self.mrna_info[mrna_id]["end"] = end
                    self.mrna_info[mrna_id]["strand"] = x[STRAND]
                    self.mrna_info[mrna_id]["line"] = line

            elif x[TYPE].lower() in ["cds"]:
                parent_id = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], "Parent"
                )
                # store line information to get the score and phase of CDS
                self.cds_info[parent_id].append([start, end, line])

            elif x[TYPE].lower() in ["exon"]:
                parent_id = FilterGFF3Terminals.get_id(
                    self.gff_file, line, x[ATTRIBUTE], "Parent"
                )
                # store line information to get the score of exon, not really necessary here as we focus on CDS coords
                self.exon_info[parent_id].append([start, end, line])

            else:
                continue

            # # change source if requested
            # if self.source:
            #     x[SOURCE] = self.source

    def get_intron_coords(self):
        """
        Get intron coordinates from CDS coordinates and NOT from exon
        coordinates. Because we need to filter short CDS exons and long
        intron models that spans the terminal CDS exons. In addition,
        sort the CDS exon coordinates.
        """
        logging.info("Compute intron coordinates ...")
        for transcript, coords in self.cds_info.items():
            # VERY IMPORTANT
            # sanity - sort coordinates
            coords.sort(key=lambda interval: interval[0])
            # print(transcript, coords)
            if len(coords) == 1:  # For mono-cds models
                self.intron_info[transcript] = []
                # update mrna info with mono-cds models
                self.mrna_info[transcript]["mono_multi"] = "mono"
                # print(coords)
                self.mrna_info[transcript]["len"] = (coords[0][1] - coords[0][0]) + 1
            else:
                # update mrna info with multi-cds models
                self.mrna_info[transcript]["mono_multi"] = "multi"
                for i, _ in enumerate(coords[:-1]):
                    # adding dummy 'transcript' to the intron_info to match the three value list as we have in cds_info or exon_info
                    self.intron_info[transcript].append(
                        [coords[i][1] + 1, coords[i + 1][0] - 1, transcript]
                    )
                    # print(
                    #     f"i:{i} {coords[i][1] + 1 } {coords[i+1][0] - 1}, {transcript}"
                    # )

    def compute_terminal_sizes(self):
        """
        Get:
        term5c_len - terminal five prime CDS length
        term3c_len - terminal three prime CDS length
        term5i_len - terminal five prime intron length
        term3i_len - terminal three prime intron length
        """
        logging.info("Compute terminal sizes [CDS, intron] ...")

        for mrna_id, feat in self.mrna_info.items():
            # print(mrna_id, feat)
            # compute terminal size only for multi-cds models
            if feat["mono_multi"] == "multi":
                (
                    self.mrna_info[mrna_id]["term5c_len"],
                    self.mrna_info[mrna_id]["term3c_len"],
                ) = FilterGFF3Terminals.get_terminal_len_multi(
                    mrna_id, feat["strand"], self.cds_info
                )

                (
                    self.mrna_info[mrna_id]["term5i_len"],
                    self.mrna_info[mrna_id]["term3i_len"],
                ) = FilterGFF3Terminals.get_terminal_len_multi(
                    mrna_id, feat["strand"], self.intron_info
                )
            else:
                (
                    self.mrna_info[mrna_id]["term5c_len"],
                    self.mrna_info[mrna_id]["term3c_len"],
                ) = (0, 0)
                (
                    self.mrna_info[mrna_id]["term5i_len"],
                    self.mrna_info[mrna_id]["term3i_len"],
                ) = (0, 0)

            # print(mrna_id, self.mrna_info[mrna_id])

    def update_mrna_info(self):
        """
        Bring all together
        """
        logging.info("Update mRNA records with CDS and intron coordinates ...")
        for mrna_id, feat in self.mrna_info.items():
            if mrna_id in self.cds_info:
                self.mrna_info[mrna_id]["cds_coords"] = self.cds_info[mrna_id]
            else:
                raise ValueError(
                    f"Error: Transcript '{mrna_id}' does not have CDS coordinates"
                )
            if mrna_id in self.intron_info:
                self.mrna_info[mrna_id]["intron_coords"] = self.intron_info[mrna_id]
            else:
                raise ValueError(
                    f"Error: Transcript '{mrna_id}' does not have intron coordinates"
                )
            # add gene line to mRNA
            if feat["parent_id"] in self.gene_info:
                self.mrna_info[mrna_id]["gline"] = self.gene_info[feat["parent_id"]][
                    "line"
                ]
            else:
                raise ValueError(f"Error: Transcript '{mrna_id}' Parent not available")

    def apply_filter(self):
        """
        apply the filter
        """
        logging.info("Tagging/Clipping mRNA records ...")
        # print GFF3 header
        print("##gff-version 3")
        for mrna_id, feat in self.mrna_info.items():

            # print(mrna_id, feat)
            gx = feat["gline"].split("\t")
            mx = feat["line"].split("\t")
            mstrand = feat["strand"]

            # apply source
            if self.source:
                gx[SOURCE] = self.source
                mx[SOURCE] = self.source

            # mstart, mend = (
            #     (int(mx[END]), int(mx[START]))
            #     if int(mx[START]) > int(mx[END])
            #     else (int(mx[START]), int(mx[END]))
            # )

            # print mono exonic models as they are
            if feat["mono_multi"] == "mono":

                # add defaults for mono exonic models
                mx[
                    ATTRIBUTE
                ] += ";Q5EXON=0;Q5INTN=0;Q5TERM=0;Q3EXON=0;Q3INTN=0;Q3TERM=0;CLIP5TERM=0;CLIP3TERM=0"

            # apply filter on multi-exonic models
            elif feat["mono_multi"] == "multi":
                # print("**Processing multi**")
                # print(mrna_id, feat)

                q5exon = False
                q5intn = False
                q5term = False
                q3exon = False
                q3intn = False
                q3term = False
                if feat["term5c_len"] <= self.term5c_len:
                    q5exon = True
                    mx[ATTRIBUTE] += ";Q5EXON=1"
                else:
                    q5exon = False
                    mx[ATTRIBUTE] += ";Q5EXON=0"
                if feat["term5i_len"] >= self.term5i_len:
                    q5intn = True
                    mx[ATTRIBUTE] += ";Q5INTN=1"
                else:
                    q5intn = False
                    mx[ATTRIBUTE] += ";Q5INTN=0"
                if q5exon and q5intn:
                    q5term = True
                    mx[ATTRIBUTE] += ";Q5TERM=1"
                else:
                    q5term = False
                    mx[ATTRIBUTE] += ";Q5TERM=0"

                if feat["term3c_len"] <= self.term3c_len:
                    q3exon = True
                    mx[ATTRIBUTE] += ";Q3EXON=1"
                else:
                    q3exon = False
                    mx[ATTRIBUTE] += ";Q3EXON=0"
                if feat["term3i_len"] >= self.term3i_len:
                    q3intn = True
                    mx[ATTRIBUTE] += ";Q3INTN=1"
                else:
                    q3intn = False
                    mx[ATTRIBUTE] += ";Q3INTN=0"
                if q3exon and q3intn:
                    q3term = True
                    mx[ATTRIBUTE] += ";Q3TERM=1"
                else:
                    q3term = False
                    mx[ATTRIBUTE] += ";Q3TERM=0"

                # # USE FOR DEBUGGING
                # print(
                #     f"{mrna_id, feat}\n"
                #     f'feat["term5c_len"] <= self.term5c_len:{feat["term5c_len"]} <= {self.term5c_len}:{feat["term5c_len"] <= self.term5c_len};\n'
                #     f'feat["term5i_len"] >= self.term5i_len:{feat["term5i_len"]} >= {self.term5i_len}:{feat["term5i_len"] >= self.term5i_len};\n'
                #     f'feat["term3c_len"] <= self.term3c_len:{feat["term3c_len"]} <= {self.term3c_len}:{feat["term3c_len"] <= self.term3c_len};\n'
                #     f'feat["term3i_len"] >= self.term3i_len:{feat["term3i_len"]} >= {self.term3i_len}:{feat["term3i_len"] >= self.term3i_len};\n'
                # )
                # print(
                #     f"q5exon:{q5exon};q5intn:{q5intn};q3exon:{q3exon};q3intn:{q3intn};"
                # )

                # apply clipping to multi exonic models

                # for --clip clip_term_exon option we need to meet q5exon or q3exon
                if self.clip in ["clip_term_exon"]:
                    # clip 5'
                    if q5exon:
                        # remove first coord set from 5' for +ve
                        if mstrand in ["+", "."]:
                            logging.warning(
                                f"Clipping 5prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][0]}'"
                            )
                            del self.cds_info[mrna_id][0]
                        # remove last coord set from 5' for -ve
                        else:
                            logging.warning(
                                f"Clipping 5prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][-1]}'"
                            )
                            del self.cds_info[mrna_id][-1]
                        mx[ATTRIBUTE] += ";CLIP5TERM=1"
                    else:
                        mx[ATTRIBUTE] += ";CLIP5TERM=0"

                    # clip 3'
                    if q3exon:
                        # remove last coord set from 3' for +ve
                        if mstrand in ["+", "."]:
                            logging.warning(
                                f"Clipping 3prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][-1]}'"
                            )
                            del self.cds_info[mrna_id][-1]
                        # remove first coord set from 3' for -ve
                        else:
                            logging.warning(
                                f"Clipping 3prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][0]}'"
                            )
                            del self.cds_info[mrna_id][0]
                        mx[ATTRIBUTE] += ";CLIP3TERM=1"
                    else:
                        mx[ATTRIBUTE] += ";CLIP3TERM=0"

                # for --clip clip_term_intron-exon option we need to meet q5term or q3term
                if self.clip in ["clip_term_intron-exon"]:
                    # clip 5'
                    if q5term:
                        # remove first coord set from 5' for +ve
                        if mstrand in ["+", "."]:
                            logging.warning(
                                f"Clipping 5prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][0]}'"
                            )
                            del self.cds_info[mrna_id][0]
                        # remove last coord set from 5' for -ve
                        else:
                            logging.warning(
                                f"Clipping 5prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][-1]}'"
                            )
                            del self.cds_info[mrna_id][-1]
                        mx[ATTRIBUTE] += ";CLIP5TERM=1"
                    else:
                        mx[ATTRIBUTE] += ";CLIP5TERM=0"

                    # clip 3'
                    if q3term:
                        # remove last coord set from 3' for +ve
                        if mstrand in ["+", "."]:
                            logging.warning(
                                f"Clipping 3prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][-1]}'"
                            )
                            del self.cds_info[mrna_id][-1]
                        # remove first coord set from 3' for -ve
                        else:
                            logging.warning(
                                f"Clipping 3prime '{mrna_id}' '{mstrand}' '{self.cds_info[mrna_id][0]}'"
                            )
                            del self.cds_info[mrna_id][0]
                        mx[ATTRIBUTE] += ";CLIP3TERM=1"
                    else:
                        mx[ATTRIBUTE] += ";CLIP3TERM=0"

                if self.clip in ["no_clip"]:
                    mx[ATTRIBUTE] += ";CLIP5TERM=0;CLIP3TERM=0"

            # extract existing cov and identity based on SPALN
            # For example, current spaln2gff Note attribute format is:
            # 'Note=AL1G11410.t1.v2.1|cov:100.00|id:85.20'
            # note_id = AL1G11410.t1.v2.1
            # note_cov = cov:100.00
            # note_iden = id:85.20
            note_id, note_cov, note_iden = None, None, None
            note_search = re.search("Note=([^;]+)", mx[ATTRIBUTE])
            if note_search:
                id_field = note_search.group(1)
                # AL1G11360.t1.v2.1|cov:100.00|id:89.40
                try:
                    note_id, note_cov, note_iden = id_field.split("|")
                # rather than exiting with error 'ValueError: not enough values to unpack (expected 3, got 2)', do not use the information if we do not follow the spaln2gff format
                # AL1G11360.t1.v2.1|cov:100.00|id:89.40
                except ValueError:
                    note_id, note_cov, note_iden = None, None, None

            if note_cov and note_cov.startswith("cov:"):
                note_cov = note_cov.replace("cov:", "cov=")
                mx[ATTRIBUTE] += f";{note_cov}"
            if note_iden and note_iden.startswith("id:"):
                note_iden = note_iden.replace("id:", "identity=")
                mx[ATTRIBUTE] += f";{note_iden}"
            # print(note_id, note_cov, note_iden)

            # compute actual coverage
            # coverage is based on cds coords, NOT exons
            sum_size = FilterGFF3Terminals.compute_coords_sum(self.cds_info[mrna_id])
            # print(f"mrna_id:'{mrna_id}';sum_size:'{sum_size}'")
            mcov = 100 * ((sum_size / 3) / len(self.fai[feat["protein_alias"]]))
            mcov = f"{mcov:.2f}"

            # DEBUG
            # print(f"Align_len (genomic bps):{sum_size}")
            # print(f"Align_len (aa bps):{sum_size / 3}")
            # print(f"Query_len (aa bps):{len(self.fai[feat['protein_alias']])}")
            # print(f"Coverage:{mcov}")

            # check if there a cov=*; information in the GFF and use it
            mrna_cov = FilterGFF3Terminals.get_id(
                self.gff_file, feat["line"], mx[ATTRIBUTE], "cov", required=False
            )
            if float(mrna_cov):
                if float(mcov) >= float(mrna_cov):
                    mcov = mrna_cov

            mx[ATTRIBUTE] += f";recal_cov={mcov}"

            # recompute mRNA start and end
            rc_mstart, rc_mend = FilterGFF3Terminals.get_start_end(
                self.cds_info[mrna_id]
            )

            # print the GFF3
            print(*gx[:3], rc_mstart, rc_mend, *gx[5:9], sep="\t")
            print(*mx[:3], rc_mstart, rc_mend, *mx[5:9], sep="\t")

            # PRINT CDS AND EXON COORDS BASED ON CDS COORDS
            for count, coords in enumerate(self.cds_info[mrna_id], start=1):
                # print(count, coords)
                cx = coords[2].split("\t")  # extract the GFF line
                cds_attrib = (
                    f"ID={mrna_id}.cds{count};{cx[8].replace('ID=', 'prev_ID=')}"
                )
                # print(*mx[:2], "CDS", *cx[3:9], sep="\t")
                print(*mx[:2], "CDS", *cx[3:8], cds_attrib, sep="\t")
                exon_attrib = f"ID={mrna_id}.exon{count};Parent={mrna_id}"
                print(*mx[:2], "exon", *cx[3:8], exon_attrib, sep="\t")
                # print(*mx[:2], "exon", *cx[3:9], sep="\t")

            # print GFF3 directive
            print("###")

    def run(self):
        logging.info(f"Processing input file '{self.gff_file}'")
        self.process_gff()
        self.get_intron_coords()
        self.compute_terminal_sizes()
        self.update_mrna_info()
        self.apply_filter()
        logging.info("Analysis complete")


class HelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


def main():

    parser = argparse.ArgumentParser(
        prog=script,
        # formatter_class=lambda prog: argparse.HelpFormatter(prog, width=80),
        formatter_class=HelpFormatter,
        description="""
        Script to filter GFF3 terminals. Script expects GFF3 input in
        gene->mRNA>[exon|CDS] format. CDS coordinates are used to filter
        the terminals and intron coordinates are computed of the CDS coordinates
        """,
        epilog=f"Contact: {__author__} ({__email__})",
    )
    parser.add_argument("protein_fasta", help="Provide query protein fasta file")
    parser.add_argument(
        "gff_file",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Provide GFF file",
    )
    parser.add_argument(
        "--protein_alias",
        type=str,
        required=True,
        default=None,
        help="Provide the GFF3 attribute field in the 'mRNA' type which can be used to link the actual\ninput protein_fasta to the alignment gff_file. For example, if the 'mRNA' attribute has\nthe original protein name in 'protein_alias=ABCD_trans_v1;' attribute, then provide\n'protein_alias' for this option. This option is required because protein aligners adds\nadditional prefix/suffix to the original protein fasta in the output GFF and makes it\ndifficult to recalculate the coverage.",
    )
    parser.add_argument(
        "--source", type=str, default=None, help="Provide new source for GFF3",
    )
    parser.add_argument(
        "--term5c_len",
        type=int,
        default=10,
        help="Minimum five prime terminal CDS exon length (bps).\nAny five prime terminal CDS exon less than the value (default: %(default)s)\nwill be flagged/removed in the GFF3 mRNA attribute (Q5TERM=1) (Q5TERM=0, if not met),\nprovided it also meets the --term5i_len parameter",
    )
    parser.add_argument(
        "--term5i_len",
        type=int,
        default=200000,
        help="Maximum five prime terminal intron length (bps).\nAny five prime terminal intron more than the value (default: %(default)s)\nwill be flagged/removed in the GFF3 mRNA attribute (Q5TERM=1) (Q5TERM=0, if not met),\nprovided it also meets the --term5c_len parameter",
    )
    parser.add_argument(
        "--term3c_len",
        type=int,
        default=10,
        help="Minimum three prime terminal CDS exon length (bps).\nAny three prime terminal CDS exon less than the value (default: %(default)s)\nwill be flagged/removed in the GFF3 mRNA attribute (Q3TERM=1) (Q3TERM=0, if not met),\nprovided it also meets the --term3i_len parameter",
    )
    parser.add_argument(
        "--term3i_len",
        type=int,
        default=200000,
        help="Maximum three prime terminal intron length (bps).\nAny three prime terminal intron more than the value (default: %(default)s)\nwill be flagged/removed in the GFF3 mRNA attribute (Q3TERM=1) (Q3TERM=0, if not met),\nprovided it also meets the --term3c_len parameter",
    )
    parser.add_argument(
        "--clip",
        type=str,
        choices=("no_clip", "clip_term_exon", "clip_term_intron-exon",),
        default="no_clip",
        help="Filter questionable alignment results by the filter types specified\n"
        "Select from one of the choices:\n"
        "'no_clip': Do not clip terminal features.\n"
        "'clip_term_exon': Clip short terminal exon if either --term5c_len or --term3c_len are\n   met (preferred option for plants like annotation).\n"
        "'clip_term_intron-exon': Clip short terminal exon/long terminal intron if both\n   --term5c_len,--term5i_len or --term3c_len,--term3i_len are met (preferred\n   option for human like annotation).\n",
    )
    # TODO: Disabiling for now below options
    # parser.add_argument(
    #     "--identity", type=float, default=0.0, help="identify cutoff [INT 0.00 - 1.00]",
    # )
    # parser.add_argument(
    #     "--coverage", type=float, default=0.0, help="coverage cutoff [INT 0.00 - 1.00]",
    # )
    # parser.add_argument(
    #     "--print_match",
    #     action="store_true",
    #     help="print output GFF3 file as match->match_part for loading to the browser",
    # )
    args = parser.parse_args()
    FilterGFF3Terminals(args).run()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
