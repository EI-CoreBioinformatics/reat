import sys
from copy import deepcopy
from functools import total_ordering


@total_ordering
class GenomicSegment:
    """
    A genomic segment defines a totally ordered hashable base type
    """

    def __init__(self, chrom, start, end, strand):
        self.chrom = sys.intern(chrom)  # Chromosome names are interned to save space and speedup comparisons
        self.start = start
        self.end = end
        self.strand = strand

    def __eq__(self, other):
        return (self.start, self.end, self.strand, self.chrom) == (other.start, other.end, other.strand, self.chrom)

    def __lt__(self, other):
        return (self.chrom, self.start, self.end, self.strand) < (self.chrom, other.start, other.end, other.strand)

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.strand))


class UniquelyIdentifiableSegment(GenomicSegment):
    """
    A UniquelyIdentifiableSegment inherits all behaviour from GenomicSegment adding the
    uid member which will be used for both hashing and equality, the uniqueness of the uid
    is not enforced by this base class
    """

    def __init__(self, uid, chrom, start, end, strand, attr):
        super().__init__(chrom, start, end, strand)
        self.uid = str(uid)
        self.attr = attr

    def __eq__(self, other):
        return self.uid == other.uid

    def __hash__(self):
        return hash(self.uid)

    def print_attributes(self, output_file=sys.stdout):
        if self.uid:
            print(f"ID={self.uid}", file=output_file, end=';')
        attrl = list(self.attr.items())
        if len(attrl) == 0:
            return
        for k, v in attrl[:-1]:
            print(f"{k}=", file=output_file, end='')
            if isinstance(v, list):
                print(f"{','.join(v)}", file=output_file, end=';')
            else:
                print(f"{v}", file=output_file, end=';')

        k, v = attrl[-1]
        print(f"{k}=", file=output_file, end='')
        if isinstance(v, list):
            print(f"{','.join(v)}", file=output_file)
        else:
            print(f"{v}", file=output_file)


class Exon(UniquelyIdentifiableSegment):
    def __init__(self, uid, chrom, gstart, gend, score, phase, strand, attr):
        super().__init__(uid, chrom, gstart, gend, strand, attr)
        self._score = score
        self._phase = phase

    def __eq__(self, other):
        return (self.start, self.end, self.strand) == (other.start, other.end, other.strand)

    def print_gff_line(self, source, ftype, file=sys.stdout):
        print(f"{self.chrom}\t{source}\t{ftype}\t{self.start}\t{self.end}\t"
              f"{self.score if self.score else '.'}\t"
              f"{self.strand if self.strand else '.'}\t"
              f"{self.phase}\t",
              file=file, end='')
        self.print_attributes(file)

    @property
    def phase(self):
        return '.' if self._phase is None else self._phase

    @phase.setter
    def phase(self, val):
        self._phase = None if val != '.' else int(val)

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, val):
        self._score = None if val != '.' else float(val)


class Transcript(UniquelyIdentifiableSegment):
    def __init__(self, uid, type, source, chrom, gstart, gend, strand, score, exons: list, cds_exons: list,
                 css, ces, attr: dict):
        super().__init__(uid, chrom, gstart, gend, strand, attr)

        # If exon and cds_exon reference the same list, make a full copy of exons onto cds_exons
        # See test_gene.py:test_translate_to:214, without the following it would fail
        if exons is cds_exons:
            cds_exons = deepcopy(exons)
        self.type = type
        self.score = score
        self.source = sys.intern(source)
        self.exons = exons
        self.exons.sort()
        self.cds_exons = cds_exons
        if strand == '+':
            self.exons.sort()
            self.cds_exons.sort()
        else:
            self.exons.sort(reverse=True)
            self.cds_exons.sort(reverse=True)
        self.css = css
        self.ces = ces

        if self.css is None and len(self.cds_exons) > 0:
            if self.strand == '+':
                self.css = self.cds_exons[0].start
            else:
                self.css = self.cds_exons[0].end

        if self.ces is None and len(self.cds_exons) > 0:
            if self.strand == '+':
                self.ces = self.cds_exons[-1].end
            else:
                self.ces = self.cds_exons[-1].start

        # Check cds_exons and css, ces

    def __eq__(self, altr):
        return (self.start, self.end, self.chrom, self.css, self.ces, self.exons) == \
               (altr.start, altr.end, altr.chrom, altr.css, altr.ces, altr.exons)

    def __hash__(self):
        return hash((self.start, self.end, self.chrom, self.css, self.ces))

    @property
    def coverage(self):
        return self.attr.get('coverage', 0)

    @property
    def identity(self):
        return self.attr.get('identity', 0)

    @property
    def tp_utr(self):
        if self.strand == '-':
            return self.css
        return self.ces

    @property
    def fp_utr(self):
        if self.strand == '-':
            return self.ces
        return self.css

    @property
    def thick_start(self):
        if self.strand == '-':
            return self.ces - 1 if self.ces else None
        return self.css - 1 if self.css else None

    @property
    def thick_end(self):
        if self.strand == '-':
            return self.css if self.css else None
        return self.ces if self.ces else None

    @property
    def colour(self):
        return self.attr.get('colour', None)

    @property
    def exon_starts(self):
        if self.strand == '-':
            return [e.start - self.start for e in reversed(self.exons)]
        return [e.start - self.start for e in self.exons]

    @property
    def exon_lens(self):
        if self.strand == '-':
            return [e.end - e.start + 1 for e in reversed(self.exons)]
        return [e.end - e.start + 1 for e in self.exons]

    @property
    def cds_starts(self):
        if self.strand == '-':
            return [e.start - self.ces for e in reversed(self.cds_exons)]
        return [e.start - self.css for e in self.cds_exons]

    @property
    def cds_lens(self):
        if self.strand == '-':
            return [e.end - e.start + 1 for e in reversed(self.cds_exons)]
        return [e.end - e.start + 1 for e in self.cds_exons]

    @property
    def bed_attr(self):
        return dict([(key, value) for key, value in self.attr.items() if key != 'colour'])

    @property
    def introns(self):
        """
        Returns:
            List of intron start, end pairs in 5' to 3' order
        """
        if self.strand == '+':
            return [(e0.end + 1, e1.start - 1) for e0, e1 in zip(self.exons, self.exons[1:])]
        else:
            return [(e1.end + 1, e0.start - 1) for e0, e1 in zip(self.exons, self.exons[1:])]

    @property
    def cds_introns(self):
        """
        Returns:
            List of cds_intron start, end pairs in 5' to 3' order
        """
        if self.strand == '+':
            return [(e0.end + 1, e1.start - 1) for e0, e1 in zip(self.cds_exons, self.cds_exons[1:])]
        else:
            return [(e0.end + 1, e1.start - 1) for e0, e1 in zip(self.cds_exons[::-1], self.cds_exons[-2::-1])]

    def print_gff(self, file=sys.stdout):
        self.print_gff_line(file)
        for e in self.exons:
            e.print_gff_line(self.source, 'exon', file=file)
        for c in self.cds_exons:
            c.print_gff_line(self.source, 'CDS', file=file)

    def print_gff_line(self, file=sys.stdout):
        print(f"{self.chrom}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t"
              f"{self.score if self.score else '.'}\t"
              f"{self.strand if self.strand else '.'}\t"
              f".\t",
              file=file, end='')
        self.print_attributes(file)

    def to_bed(self, cds_only=False):
        s = self.chrom
        if cds_only:
            s += '\t' + str(self.thick_start)
            s += '\t' + str(self.thick_end)
        else:
            s += '\t' + str(self.start - 1)
            s += '\t' + str(self.end)

        if self.uid is None:
            return s
        s += '\t' + self.uid
        if self.score is None:
            return s
        s += '\t' + str(self.score)
        if self.strand is None:
            return s
        s += '\t' + self.strand
        if self.thick_start is None:
            return s
        s += '\t' + str(self.thick_start)
        if self.thick_end is None:
            return s
        s += '\t' + str(self.thick_end)
        if self.colour is None:
            s += "\t0,0,0"
        else:
            s += '\t' + self.colour

        if cds_only:
            if len(self.cds_exons) == 0:
                return s
            s += '\t' + str(len(self.cds_exons))
            s += '\t' + ','.join([str(ln) for ln in self.cds_lens])
            s += '\t' + ','.join([str(st) for st in self.cds_starts])
        else:
            if len(self.exons) == 0:
                return s
            s += '\t' + str(len(self.exons))
            s += '\t' + ','.join([str(ln) for ln in self.exon_lens])
            s += '\t' + ','.join([str(st) for st in self.exon_starts])

        if self.bed_attr is None:
            return s

        s += '\t'
        if len(self.cds_exons) > 0:
            s += 'CDSphase=' + self.cds_exons[0].phase + ';'

        s += ';'.join(
            [key + '='
             +
             str(value) if not isinstance(value, list)
             else key + '=' + ','.join([v for v in value])
             for key, value in self.bed_attr.items()
             ]
        )
        return s

    @staticmethod
    def _merge_exons(exons_it, distance):
        merged = []
        for higher in exons_it:
            if not merged:
                merged.append(higher)
            else:
                lower = merged[-1]
                # test for intersection
                if higher.start - lower.end < distance:
                    upper_bound = max(lower.end, higher.end)
                    merged[-1].start, merged[-1].end = (lower.start, upper_bound)
                else:
                    merged.append(higher)
        return merged

    def merge_exons(self, distance=4):
        num_exons = len(self.exons)
        num_cds_exons = len(self.cds_exons)
        self.exons = self._merge_exons(reversed(self.exons) if self.strand == '-' else iter(self.exons),
                                       distance)
        self.cds_exons = self._merge_exons(reversed(self.cds_exons) if self.strand == '-' else iter(self.cds_exons),
                                           distance)
        if self.strand == '-':
            self.exons.reverse()
            self.cds_exons.reverse()
        return num_exons - len(self.exons), num_cds_exons - len(self.cds_exons)


class Gene(UniquelyIdentifiableSegment):
    # A gene has uid, source, chrom, start, end, strand, score
    # mrnas (aka Transcripts)
    # potentially multiple other GFF features
    # and attributes
    def __init__(self, uid, source, chrom, start, end, strand, score, mrnas, features, attr):
        super().__init__(uid, chrom, start, end, strand, attr)
        self.source = source
        self.score = score
        for k, v in mrnas.items():
            if len(v.exons) == 0 and len(v.cds_exons) > 0:
                v.exons = deepcopy(v.cds_exons)
                for e in v.exons:
                    e.uid = f"virt_exon-{e.uid}"
            if v.strand == '-':
                v.exons.sort(reverse=True)
                v.cds_exons.sort(reverse=True)
            else:
                v.exons.sort()
                v.cds_exons.sort()

            if v.css is None and len(v.cds_exons) > 0:
                if v.strand == '+':
                    v.css = v.cds_exons[0].start
                else:
                    v.css = v.cds_exons[0].end

            if v.ces is None and len(v.cds_exons) > 0:
                if v.strand == '+':
                    v.ces = v.cds_exons[-1].end
                else:
                    v.ces = v.cds_exons[-1].start

        self.mrnas = mrnas
        self.features = features
        self.attr = attr

    def __eq__(self, altr):
        if (self.start, self.end, self.chrom, self.strand) == (altr.start, altr.end, altr.chrom, altr.strand):
            if len(self.mrnas.values()) == len(altr.mrnas.values()):
                for s, a in zip(self.mrnas.values(), altr.mrnas.values()):
                    if s != a:
                        return False
                return True
        return False

    def __hash__(self):
        return hash((self.uid, self.start, self.end, self.chrom, self.strand))

    def extend_three_prime_end_cds(self, mrna_id, num_nts):
        mrna = self.mrnas[mrna_id]
        last_exon = mrna.exons[-1]
        last_cds_exon = mrna.cds_exons[-1]
        if mrna.strand == '-':
            last_cds_exon.start -= num_nts
            mrna.ces = last_cds_exon.start
            last_exon.start = min(last_exon.start, last_cds_exon.start)
            mrna.start = min(mrna.start, last_cds_exon.start)
            self.start = min(self.start, last_cds_exon.start)
        else:
            last_cds_exon.end += num_nts
            mrna.ces = last_cds_exon.end
            last_exon.end = max(last_exon.end, last_cds_exon.end)
            mrna.end = max(mrna.end, last_cds_exon.end)
            self.end = max(self.end, last_cds_exon.end)

    def extend_five_prime_end_cds(self, mrna_id, num_nts):
        mrna = self.mrnas[mrna_id]
        first_exon = mrna.exons[0]
        first_cds_exon = mrna.cds_exons[0]
        if mrna.strand == '-':
            first_cds_exon.end += num_nts
            mrna.css = first_cds_exon.end
            first_exon.end = max(first_exon.end, first_cds_exon.end)
            mrna.end = max(mrna.end, first_cds_exon.end)
            self.end = max(self.end, first_cds_exon.end)
        else:
            first_cds_exon.start -= num_nts
            mrna.css = first_cds_exon.start
            first_exon.start = min(first_exon.start, first_cds_exon.start)
            mrna.start = min(mrna.start, first_cds_exon.start)
            self.start = min(self.start, first_cds_exon.start)

    def translate_to(self, new_start, chrom=None):
        chrom = self.chrom if chrom is None else chrom
        for mrna in self.mrnas.values():
            for e in mrna.exons:
                e.chrom = chrom
                e.end -= mrna.start
                e.end += new_start
                e.start -= mrna.start
                e.start += new_start
            for c in mrna.cds_exons:
                c.chrom = chrom
                c.end -= mrna.start
                c.end += new_start
                c.start -= mrna.start
                c.start += new_start

            if mrna.strand == '-':
                mrna.ces = mrna.cds_exons[-1].start
                mrna.css = mrna.cds_exons[0].end
            else:
                mrna.ces = mrna.cds_exons[-1].end
                mrna.css = mrna.cds_exons[0].start
            mrna.chrom = chrom
            mrna.end -= mrna.start
            mrna.end += new_start
            mrna.start -= mrna.start
            mrna.start += new_start

        self.chrom = chrom
        self.end -= self.start
        self.end += new_start
        self.start = new_start

    def change_id(self, new_uid):
        for m in self.mrnas.values():
            m.attr["Parent"] = [p.replace(self.uid, str(new_uid)) for p in m.attr["Parent"]]
        self.uid = str(new_uid)

    def print_gff_line(self, file=sys.stdout):
        print(f"{self.chrom}\t{self.source}\tgene\t{self.start}\t{self.end}\t"
              f"{self.score if self.score else '.'}\t"
              f"{self.strand if self.strand else '.'}\t"
              f".\t",
              file=file, end='')
        self.print_attributes(file)

    def print_gff(self, file=sys.stdout):
        # Format the gene as gff
        self.print_gff_line(file=file)
        for mrna_id, mrna in self.mrnas.items():
            mrna.print_gff(file=file)
