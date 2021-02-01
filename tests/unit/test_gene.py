from annotation.lib.models.Transcript import Exon, Gene, GenomicSegment, Transcript
from annotation.lib.parsers.BED import BEDReader


def test_segment_hashable():
    s = GenomicSegment("chr1", 0, 10, 1)
    d = dict()
    d[s] = 0
    assert d[s] == 0


def test_identity():
    attr = {'identity': 20}
    g = Transcript("uid", "source", "chr1", 0, 10, '+', 0, [], [], None, None, attr)
    assert g.identity == 20


def test_coverage():
    attr = {'coverage': 20}
    g = Transcript("uid", "source", "chr1", 0, 10, '+', 0, [], [], None, None, attr)
    assert g.coverage == 20


def test_missing_identity():
    g = Transcript("uid", "source", "chr1", 0, 10, '+', 0, [], [], None, None, dict())
    assert g.identity == 0


def test_missing_coverge():
    g = Transcript("uid", "source", "chr1", 0, 10, '+', 0, [], [], None, None, dict())
    assert g.coverage == 0


def test_bed_representation():
    bed_line = "Chr1\t3759\t6272\tAL1G11530.t1.v2.1.m1\t100\t+\t3759\t6272\t0,0,0\t8\t154,281,120,390,153,96,7,101\t0,236,726,946,1414,1679,2326,2412\tCDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1"
    g = BEDReader.make_transcript_from_bed(bed_line.strip().split(), "test")
    s = g.to_bed()
    assert s == bed_line


def test_gene_hashable():
    g = Gene("uid", "source", "chr1", 0, 10, '+', 0, dict(), dict(), dict())
    d = dict()
    d[g] = "test"
    assert d[g] == "test"


def test_gene_lt():
    g = Gene("uid", "source", "chr1", 0, 10, '+', 0, dict(), dict(), dict())
    assert g < Gene("uid", "source", "chr1", 10, 20, '+', 0, dict(), dict(), dict())


def test_gene_equality():
    g = Gene("uid", "source", "chr1", 0, 10, '+', 0, dict(), dict(), dict())
    assert g == Gene("uid", "source", "chr1", 0, 10, '+', 0, dict(), dict(), dict())
    assert g != Gene("uid", "source", "chr1", 2, 10, '+', 0, dict(), dict(), dict())

    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, dict())

    assert g != Gene("uid", "source", "chr1", 0, 10, '+', 0, t, dict(), dict())

    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 82)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t2 = dict()
    t2["muid2"] = Transcript("muid2", "source", "chr1", 5, 82, '-', 0, exons, exons, None, None, dict())

    assert Gene("uid", "source", "chr1", 0, 10, '+', 0, t, dict(), dict()) != Gene("uid", "source", "chr1",
                                                                                   0, 10, '+', 0,
                                                                                   t2, dict(), dict())


def test_segment_equality():
    s = GenomicSegment("chr1", 0, 10, 1)
    assert s == GenomicSegment("chr1", 0, 10, 1)


def test_segment_lt():
    s = GenomicSegment("chr1", 10, 20, 1)
    assert GenomicSegment("chr1", 0, 10, 1) < s


def test_intron_plus_strand():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = Transcript("uid", "source", "chr1", 5, 80, '+', 0, exons, [], None, None, dict())

    assert t.introns == [(11, 19), (31, 49)]


def test_intron_minus_strand():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = Transcript("uid", "source", "chr1", 5, 80, '-', 0, exons, [], None, None, dict())

    assert t.introns == [(31, 49), (11, 19)]


def test_extending_three_prime_end():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert g.start == 5
    assert t["muid"].start == 5
    assert t["muid"].ces == 5
    assert exons[-1].start == 5

    g.extend_three_prime_end_cds("muid", 3)
    assert g.start == 2
    assert t["muid"].start == 2
    assert t["muid"].ces == 2
    assert exons[-1].start == 2

    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '+', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert g.end == 80
    assert t["muid"].end == 80
    assert t["muid"].ces == 80
    assert exons[-1].end == 80

    g.extend_three_prime_end_cds("muid", 3)
    assert g.end == 83
    assert t["muid"].end == 83
    assert t["muid"].ces == 83
    assert exons[-1].end == 83


def test_extending_five_prime_end():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert g.end == 80
    assert t["muid"].end == 80
    assert exons[0].end == 80
    assert t["muid"].css == 80

    g.extend_five_prime_end_cds("muid", 3)
    assert g.end == 83
    assert t["muid"].end == 83
    assert exons[0].end == 83
    assert t["muid"].css == 83

    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '+', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert g.start == 5
    assert t["muid"].start == 5
    assert exons[0].start == 5
    assert t["muid"].css == 5

    g.extend_five_prime_end_cds("muid", 3)
    assert g.start == 2
    assert t["muid"].start == 2
    assert exons[0].start == 2
    assert t["muid"].css == 2


def test_merge_exons():
    exons = []
    for i, pair in enumerate([(5, 10), (12, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    for e, p in zip(t['muid'].exons, reversed([(5, 10), (12, 30), (50, 80)])):
        assert (e.start, e.end) == p

    assert len(g.mrnas["muid"].exons) == 3
    g.mrnas["muid"].merge_exons()

    for e, p in zip(g.mrnas['muid'].exons, reversed([(5, 30), (50, 80)])):
        assert (e.start, e.end) == p
    assert len(g.mrnas["muid"].exons) == 2

    exons = []
    for i, pair in enumerate([(5, 10), (20, 48), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '+', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert len(g.mrnas["muid"].exons) == 3
    g.mrnas["muid"].merge_exons()
    merged_exons = [(5, 10), (20, 80)]
    for e, p in zip(g.mrnas['muid'].exons, merged_exons):
        assert (e.start, e.end) == p
    assert len(g.mrnas["muid"].exons) == 2

    exons = []
    for i, pair in enumerate([(5, 15), (16, 48), (50, 80), (83, 95)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '+', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert len(t['muid'].exons) == 4
    t['muid'].merge_exons()
    assert len(t['muid'].exons) == 1
    assert t['muid'].exons[0].start, t['muid'].exons[0].end == (5, 95)


def test_translate_to():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    new_start = 1000
    assert g.chrom == "chr1"
    assert t["muid"].chrom == "chr1"
    assert t["muid"].exons[0].chrom == "chr1"
    g.translate_to(new_start)
    assert g.chrom == "chr1"
    assert t["muid"].chrom == "chr1"
    assert t["muid"].exons[0].chrom == "chr1"

    g.translate_to(new_start, "other")
    assert g.chrom == "other"
    assert g.start == new_start
    assert g.end == 1075
    assert t["muid"].start == new_start
    assert t["muid"].end == 1075
    assert t["muid"].chrom == "other"

    assert t["muid"].exons[0].start == 1045
    assert t["muid"].exons[0].end == 1075
    assert t["muid"].exons[0].chrom == "other"

    assert t["muid"].css == 1075
    assert t["muid"].ces == 1000

    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '+', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '+', 0, exons, exons, None, None, dict())
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    new_start = 1000
    g.translate_to(new_start, "other2")
    assert g.chrom == "other2"
    assert g.start == new_start
    assert g.end == 1075
    assert t["muid"].start == new_start
    assert t["muid"].end == 1075
    assert t["muid"].chrom == "other2"

    assert t["muid"].exons[0].start == 1000
    assert t["muid"].exons[0].end == 1005
    assert t["muid"].exons[0].chrom == "other2"

    assert t["muid"].css == 1000
    assert t["muid"].ces == 1075


def test_change_id():
    exons = []
    for i, pair in enumerate([(5, 10), (20, 30), (50, 80)]):
        exons.append(Exon(f"e{i}", "chr1", pair[0], pair[1], 0, 0, '-', None))
    t = dict()
    t["muid"] = Transcript("muid", "source", "chr1", 5, 80, '-', 0, exons, exons, None, None, {"Parent": ["guid"]})
    g = Gene("guid", "source", "chr1", t["muid"].start, t["muid"].end, t["muid"].strand, t["muid"].score, t, [], dict())

    assert g.uid == "guid"
    assert g.mrnas["muid"].uid == "muid"
    assert len(g.mrnas["muid"].attr["Parent"]) == 1
    assert "guid" in g.mrnas["muid"].attr["Parent"]

    g.change_id("guid2")
    assert g.uid == "guid2"
    assert g.mrnas["muid"].uid == "muid"
    assert len(g.mrnas["muid"].attr["Parent"]) == 1
    assert "guid2" in g.mrnas["muid"].attr["Parent"]

# TODO: Check constructor, it should catch any mismatching coordinates on the exon->mrna->gene hierarchy
