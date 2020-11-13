#!/usr/bin/env python3

import Bio.SeqIO
import Bio.Seq
import sys
import argparse
from collections import Counter
import logging
import re


"""Simple utility to discard any non-description information out of the FASTA file, so to reduce
the chances of invalid characters creeping in. It will also check the consistency of
the identifiers and remove duplicated ones."""


formatter = logging.Formatter(
        "{asctime} - {name} - {filename}:{lineno} - {levelname} - {funcName} \
- {processName} - {message}",
        style="{"
        )


def create_default_logger(name, level="WARN"):
    """Default logger
    :param name: string used to give a name to the logger.
    :type name: str

    :param level: level of the logger. Default: WARN
    """

    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def main():

    logger = create_default_logger("sanitizer")

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("-cstop", "--correct-stops", default=False, action="store_true",
                        help="Transform plain '.' characters into '*' in protein sequences.")
    parser.add_argument("fasta", nargs="+", type=argparse.FileType("rt"))
    args = parser.parse_args()

    found_ids = Counter()

    starter = 96

    for fasta in args.fasta:
        if len(args.fasta) > 1:
            starter += 1
            prefix = "{}_".format(chr(starter))
        else:
            prefix = ""
        for record in Bio.SeqIO.parse(fasta, "fasta"):
            if record.id in found_ids:
                logger.warning("ID found other {} time{} in the input files!".format(
                    found_ids[record.id], "s" if found_ids[record.id] > 1 else ""))
            record.id = "{}{}".format(prefix, re.sub("\|", "_", record.id))
            record.description = ""
            if '.' in record.seq:
                record.seq = Bio.Seq.Seq(re.sub('\.', '*', str(record.seq)))

            Bio.SeqIO.write(record, args.out, "fasta")

    args.out.close()

main()