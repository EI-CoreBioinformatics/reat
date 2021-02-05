# coding: utf_8
"""
Module to serialize GFF files.
"""
import sys

from ..models.Transcript import Exon, Gene, Transcript
from ..parsers import get_handle


class EmptyFileError(Exception):
    pass


class GFFReader(object):
    """
    A GFFReader provides a bit more than just parsing a GFF file, it also collates attributes
    and unique identifiers for the features in the file.

    Whilst working with the elements of a GFFReader generated objects the GFFReader should
    remain in scope.
    """

    def __init__(self, filename):
        self.filename = filename
        self.fh = get_handle(filename)
        self.line = ""
        self.raw = [None, None, None]
        self.records_count = 0
        self.uids = set()  # Collection of all unique identifiers defined in a file
        self.attr_orders_to_id = dict()  # Attribute attr_order to id
        self.order_id_to_attr_order = dict()  # Reverse order_id to attr_order tuple
        self.num_orders = 0
        self.attr_file_pos = dict()  # Associates the position in file of an attr line to a uid
        self.header = dict()
        # Skip header
        while True:
            cur = self.fh.tell()
            line = self.fh.readline()
            if not line.startswith('#') or self.fh.tell() == cur:
                self.fh.seek(cur)
                break
            else:
                if line[:2] == '##' and line[2] != '#':
                    header_kv = line.split(maxsplit=1)
                    self.header[header_kv[0][2:]] = header_kv[1]

        self.format_check()

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self.fh.close()

    def __iter__(self):
        return self

    def __next__(self):
        return self.read()

    def format_check(self):
        """
        GFF3 formats are supposed to be tab-delimited and 9 required fields
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        """
        cur = self.fh.tell()
        raw = self.fh.readline().strip().split('\t')
        if raw[0] == '':
            raise EmptyFileError(f"{self.filename} appears to be empty or contain a leading empty line")
        if len(raw) != 9:
            raise IOError(
                "ERROR reading {0}, expected 9 tab-delimited fields but saw {1}".format(self.filename, len(raw))
            )
        self.fh.seek(cur)

    def read(self):
        """
        GFF files
        (0) chr
        (1) annotation source
        (2) type: gene|transcript|CDS|exon|UTR|match_part|cDNA_match
        (3) 1-based start
        (4) 1-based end
        (5) score
        (6) strand: +|-|.
        (7) phase
        (8) attributes

        Entries need to be delimited by '###' or closed when encountering a feature without a parent attribute
        """
        self.skip_comments()
        assert 'Parent' not in self.raw[8]
        gene_chr = sys.intern(self.raw[0])
        gene_source = sys.intern(self.raw[1])
        gene_start = int(self.raw[3])
        gene_end = int(self.raw[4])
        if gene_start < 0 or gene_end < 0:
            print(f"WARNING: Negative coordinates ({self.raw[3]}, {self.raw[4]}) in file {self.filename}\n"
                  f"Line:\n{self.line}", file=sys.stderr)
            gene_start = abs(gene_start)
            gene_end = abs(gene_end)
        try:
            gene_score = float(self.raw[5])
        except ValueError:
            gene_score = '.'
        gene_strand = self.raw[6]
        gene_attr = dict()
        gene_uid = self.parse_attributes(gene_attr, "gene")
        gene_mrnas = {}
        gene_features = []
        while True:  # Read features until the feature we are reading has no parent
            self.line = self.fh.readline().strip()
            cur = self.fh.tell()
            if self.line == '' and cur == self.fh.tell() or self.line.startswith('###'):
                return Gene(gene_uid, gene_source, gene_chr, gene_start, gene_end, gene_strand,
                            gene_score, gene_mrnas, gene_features, gene_attr)
            if self.line.startswith('#'):
                continue

            self.raw = self.line.split('\t')
            assert len(self.raw) == 9
            entry_type = self.raw[2].lower()
            if 'region' == entry_type:
                continue

            entry_attr = dict()
            entry_uid = self.parse_attributes(entry_attr, entry_type)
            if 'Parent' not in entry_attr:
                return Gene(gene_uid, gene_source, gene_chr, gene_start, gene_end, gene_strand,
                            gene_score, gene_mrnas, gene_features, gene_attr)

            # Here we are passing in the 'raw' entry_type to store it in the Transcript object for printing
            self.collect_gene_entries(cur, entry_attr, self.raw[2], entry_uid, gene_features, gene_mrnas)

    def collect_gene_entries(self, cur, entry_attr, entry_type, entry_uid, gene_features, gene_mrnas):
        # Get the lowercase version of the entry_type for checking against
        entry_type_lc = entry_type.lower()
        entry_chr = self.raw[0]
        entry_source = sys.intern(self.raw[1])
        entry_start = int(self.raw[3])
        entry_end = int(self.raw[4])
        if entry_start < 0 or entry_end < 0:
            print(f"WARNING: Negative coordinates ({self.raw[3]}, {self.raw[4]}) in file {self.filename}\n"
                  f"Line:\n{self.line}", file=sys.stderr)
            entry_start = abs(entry_start)
            entry_end = abs(entry_end)
        entry_strand = self.raw[6]
        entry_score = 0 if self.raw[5] == '.' else float(self.raw[5])
        entry_phase = self.raw[7]
        if entry_type_lc == 'exon':
            for parent in entry_attr['Parent']:
                if parent not in gene_mrnas:
                    gene_mrnas[parent] = Transcript(entry_uid, entry_type, entry_source, entry_chr, entry_start, entry_end,
                                                    entry_strand, entry_score, [], [], None, None,
                                                    entry_attr)
                gene_mrnas[parent].exons.append(
                    Exon(entry_uid, entry_chr, entry_start, entry_end,
                         entry_score, entry_phase, entry_strand, entry_attr)
                )
        elif entry_type_lc == 'cds':
            for parent in entry_attr['Parent']:
                if parent not in gene_mrnas:
                    gene_mrnas[parent] = Transcript(entry_uid, entry_type, entry_source, entry_chr, entry_start, entry_end,
                                                    entry_strand, entry_score, [], [], None, None,
                                                    entry_attr)
                gene_mrnas[parent].cds_exons.append(
                    Exon(entry_uid, entry_chr, entry_start, entry_end,
                         entry_score, entry_phase, entry_strand, entry_attr)
                )
        elif entry_type_lc == 'five_prime_utr' or entry_type_lc == '5\'utr':
            pass
        elif entry_type_lc == 'three_prime_utr' or entry_type_lc == '3\'utr':
            pass
        elif entry_type_lc == 'mrna' or entry_type_lc.endswith("_gene_segment"):
            score = 0 if self.raw[5] == '.' else float(self.raw[5])
            gene_mrnas[entry_uid] = Transcript(entry_uid, entry_type, entry_source, entry_chr, entry_start,
                                               entry_end, entry_strand, score, [], [], None, None,
                                               entry_attr)
        else:
            gene_features.append(self.line)
            print("WARNING: Unhandled type {0}, in file {3} at position {1}\nLine:\n{2}"
                  .format(entry_type, cur, self.line, self.filename), file=sys.stderr)

    def parse_attributes(self, attr, entry_type):
        # Parent, Alias, Note, Dbxref and Ontology_term attributes can have multiple values
        attr_order = []
        self.records_count += 1
        uid = self.records_count
        for blob in [x.strip() for x in self.raw[8].split(';')]:
            blob = blob.strip()
            if blob.startswith('Name='):
                attr['Name'] = blob[5:]
                attr_order.append('Name')
            elif blob.startswith('Parent='):
                attr['Parent'] = blob[7:].split(',')
                attr_order.append('Parent')
                for parent in attr['Parent']:
                    if parent not in self.uids:
                        raise Exception("Formatting error, unseen parent {}. In file {}, at position {}\nLine: {}"
                                        .format(parent, self.filename, self.fh.tell(), self.line))
            elif blob.startswith('Alias='):
                attr['Alias'] = [x.strip() for x in blob[6:].split(',')]
                attr_order.append('Alias')
            elif blob.startswith('Dbxref='):
                attr['Dbxref'] = [x.strip() for x in blob[7:].split(',')]
                attr_order.append('Dbxref')
            elif blob.startswith('Ontology_term='):
                attr['Ontology_term'] = [x.strip() for x in blob[14:].split(',')]
                attr_order.append('Ontology_term')
            elif blob.startswith('ID='):
                uid = blob[3:]
                attr_order.append('ID')
            elif blob.startswith('Note='):
                attr["Note"] = blob[5:]  # Keep the 'note' attribute intact
                attr_order.append('Note')
            else:
                parts = blob.split('=', 1)
                # Handle trailing semi-colons and attributes with no values
                parts_len = len(parts)
                if parts_len > 1:
                    attr[sys.intern(parts[0])] = parts[1]
                    attr_order.append(parts[0])

        if uid not in self.uids:
            self.uids.add(uid)
        else:
            if 'Parent' not in attr:
                pass
            elif entry_type != 'cds':
                raise Exception("Repeated ID {} in file {} at position {}\nCurrent line: {}\nuids: {}".format(
                    uid, self.filename, self.fh.tell(), self.line, self.uids))

        return uid

    def skip_comments(self):
        if self.raw[0] is None:  # Just for the first pass through
            self.line = self.fh.readline().strip()
        self.raw = self.line.strip().split('\t')
        while self.raw[0].startswith('#') or (len(self.raw) > 2 and self.raw[2] in ['region', 'repeat_region']):
            line = self.fh.readline().strip()
            self.raw = line.strip().split('\t')
        if len(self.raw) == 0 or self.raw[0] == '':
            raise StopIteration("EOF")
