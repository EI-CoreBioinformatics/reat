# coding: utf-8
"""
Module to parse BED12 objects. Much more basic than what PyBedtools could offer,
but at the same time more pythonic.
"""
import sys
from pathlib import Path

from ..models.Transcript import Exon, Gene, Transcript
from ..parsers import get_handle


class BEDReader:

    def __init__(self, filepath):
        self.file_name = filepath
        self.handle = get_handle(filepath)
        self.done = False
        self.num_records = 1
        self.file_size = Path(filepath).stat().st_size
        self.check_format()
        self.t = self.parse_line()  # Current transcript, helps accumulate transcripts in genes when file is gene sorted

    def parse_line(self):
        # The line should now contain 3 to 12 or 13 fields
        #  0) chrom
        #  1) start
        #  2) end
        #  3) name
        #  4) score
        #  5) strand
        #  6) thickStart
        #  7) thickEnd
        #  8) itemRgb
        #  9) blockCount
        # 10) blockSizes
        # 12) blockStarts
        # We will consider thickStart and thickEnd to be coding sequence coordinates

        cur = self.handle.tell()
        line = self.handle.readline()
        if line == '\n' and self.handle.tell() < self.file_size:
            raise ValueError("Unexpected empty line in " + self.handle.name + ":" + str(cur))
        if line == '\n' or line == '':
            self.done = True
            return None

        raw = line.strip().split()

        num_fields = len(raw)
        if num_fields < 3:
            raise ValueError("Found " +
                             str(num_fields) +
                             " fields when expecting at least 3 for file\n" +
                             self.handle.name + ":" + str(cur) +
                             "\nLine:\n" +
                             line)
        t = self.make_transcript_from_bed(raw, self.file_name + "_" + str(self.num_records))
        self.num_records += 1
        return t

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handle.close()

    def __iter__(self):
        return self

    def __next__(self):
        if self.done:
            raise StopIteration
        return self.read()

    def check_format(self):
        line = self.handle.readline().strip()
        if len(line.split()) < 3:
            raise ValueError("Invalid format found in " + self.handle.name)
        self.handle.seek(0)

    def read(self):
        t = self.parse_line()
        ts = dict()
        ts[self.t.uid] = self.t
        gene_uid = self.t.attr["geneID"]
        gene_start = self.t.start
        gene_end = self.t.end
        gene_chrom = self.t.chrom
        gene_strand = self.t.strand
        while t and t.attr["geneID"] == self.t.attr["geneID"]:
            ts[t.uid] = t
            self.t = t
            gene_start = min(gene_start, t.start)
            gene_end = max(gene_end, t.end)
            t = self.parse_line()
        self.t = t
        return Gene(gene_uid, self.file_name, gene_chrom,
                    gene_start, gene_end, gene_strand, 0, ts, [], dict())

    @classmethod
    def parse_bed_attr(cls, blob):
        d = dict()
        attr_order = []
        for pair in blob.split(';'):
            name, value = pair.split('=', 1)
            d[sys.intern(name)] = value
            attr_order.append(name)
        return d

    @classmethod
    def make_transcript_from_bed(cls, raw, uid):
        num_fields = len(raw)
        chrom, cstart, cend = raw[0], int(raw[1]), int(raw[2])
        name = ""
        score = 0
        strand = '.'
        attr = dict()
        # block_count = 0
        block_sizes = []
        block_starts = []
        if num_fields >= 4:
            name = raw[3]
        if num_fields >= 5:
            score = raw[4]
        if num_fields >= 6:
            strand = raw[5]
        if num_fields >= 12:
            # Parse comma separated values to list of ints (the if x handles trailing commas)
            block_sizes = [int(x) for x in raw[10].split(',') if x]
            block_starts = [int(x) for x in raw[11].split(',') if x]
        if num_fields == 13:  # Process potential ad-hoc comments
            attr = cls.parse_bed_attr(raw[12])
            attr['colour'] = raw[8]
        if num_fields >= 8:
            attr['colour'] = raw[8]
        # Transform the block information into exons
        exons = [
            Exon(uid, name + "_" + str(i), cstart + p[0] + 1, cstart + p[0] + p[1], 0, 0, strand, None)
            for i, p in enumerate(zip(block_starts, block_sizes), 1)
        ]

        fp_utr = None
        tp_utr = None
        # Extract from the attr CDS;CDSPhase;geneID, calculate the 5', 3' and if possible coding exons coordinates
        cds_coords = attr.get('CDS', None)
        if cds_coords is not None:
            cds_coords = cds_coords.split(':')
            cds_coords = (int(cds_coords[0]), int(cds_coords[1]))
            fp_utr = cds_coords[0]
            tp_utr = cds_coords[1]
            if strand == '-':
                fp_utr, tp_utr = tp_utr, fp_utr
                tp_utr += 1
            else:
                fp_utr += 1

        return Transcript(name, "bed", "mRNA", chrom, cstart + 1, cend, strand, score, exons, [], fp_utr, tp_utr, attr)
