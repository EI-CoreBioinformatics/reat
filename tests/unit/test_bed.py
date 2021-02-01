import tempfile
import pytest

from annotation.lib.parsers.BED import BEDReader


def test_invalid_format():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print("Chr1    6917    8666", file=tmp_bed)
        print("a\tb", file=tmp_bed)
        tmp_bed.flush()
        with pytest.raises(ValueError) as ex:
            for _ in BEDReader(tmp_bed.name):
                pass
        assert ex.match(r".* fields when expecting at least 3 for *.")


def test_empty_file():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        tmp_bed.flush()
        with pytest.raises(ValueError) as ex:
            for _ in BEDReader(tmp_bed.name):
                pass
        assert ex.match("Invalid format .*")


def test_empty_line():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_empty_line_in_contents, file=tmp_bed)
        tmp_bed.flush()
        with pytest.raises(ValueError) as ex:
            for _ in BEDReader(tmp_bed.name):
                pass

        assert ex.match("Unexpected empty line .*")


def test_bed_with_construct():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_contents, file=tmp_bed)
        tmp_bed.flush()
        cont = 0
        with BEDReader(tmp_bed.name) as bed_parser:
            for _ in bed_parser:
                if cont == 3:
                    result = _
                cont += 1

        assert result.chrom == "Chr1"
        assert result.start == 11867
        assert result.end == 12940
        assert result.mrnas["AL181U10040.t1.v2.1.m1"].css == 12940
        assert result.mrnas['AL181U10040.t1.v2.1.m1'].ces == 11867
        assert len(result.mrnas['AL181U10040.t1.v2.1.m1'].exons) == 1
        assert result.mrnas['AL181U10040.t1.v2.1.m1'].exons[0].end == result.start + 1073

        assert _.chrom == "Chr1"
        assert _.start == 38982
        assert _.end == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].css == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].ces == 38982
        assert len(_.mrnas['AL1G22550.t1.v2.1.m3'].exons) == 4
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].exons[-1].end == _.start + 72


def test_bed_parser():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_contents, file=tmp_bed)
        tmp_bed.flush()
        cont = 0
        for _ in BEDReader(tmp_bed.name):
            if cont == 3:
                result = _
            cont += 1

        assert result.chrom == "Chr1"
        assert result.start == 11867
        assert result.end == 12940
        assert result.mrnas['AL181U10040.t1.v2.1.m1'].css == 12940
        assert result.mrnas['AL181U10040.t1.v2.1.m1'].ces == 11867
        assert len(result.mrnas['AL181U10040.t1.v2.1.m1'].exons) == 1
        assert result.mrnas['AL181U10040.t1.v2.1.m1'].exons[0].end == result.start + 1073

        assert _.chrom == "Chr1"
        assert _.start == 38982
        assert _.end == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].css == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].ces == 38982
        assert len(_.mrnas['AL1G22550.t1.v2.1.m3'].exons) == 4
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].exons[-1].end == _.start + 72


def test_no_newline_at_end():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_no_newline_at_end, file=tmp_bed)
        tmp_bed.flush()
        for _ in BEDReader(tmp_bed.name):
            pass

        assert _.chrom == "Chr1"
        assert _.start == 38982
        assert _.end == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].css == 40877
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].ces == 38982
        assert len(_.mrnas['AL1G22550.t1.v2.1.m3'].exons) == 4
        assert _.mrnas['AL1G22550.t1.v2.1.m3'].exons[-1].end == _.start + 72


def test_bed_parser_groups_by_gene_id():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_contents_collapse, file=tmp_bed)
        tmp_bed.flush()
        cont = 0
        for _ in BEDReader(tmp_bed.name):
            if cont == 1:
                result = _
            cont += 1

        assert len(result.mrnas) == 2


def test_bed_parser_groups_by_gene_id_at_end():
    with tempfile.NamedTemporaryFile("wt") as tmp_bed:
        print(gffread_bed_contents_collapse_at_end, file=tmp_bed)
        tmp_bed.flush()
        for _ in BEDReader(tmp_bed.name):
            pass

        # Just check that the last transcripts were collapsed onto the _ gene
        assert len(_.mrnas) == 2


gffread_bed_contents = """Chr1    3759    6272    AL1G11530.t1.v2.1.m1    100     +       3759    6272    0,0,0   8       154,281,120,390,153,96,7,101,   0,236,726,946,1414,1679,2326,2412,      CDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t1.v2.1.m1    100     -       6917    8666    0,0,0   9       152,76,67,86,74,46,90,48,96,    0,239,466,646,844,1024,1318,1499,1653,  CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t2.v2.1.m1    100     -       6917    8666    0,0,0   8       152,76,67,86,117,90,48,96,      0,239,466,646,844,1318,1499,1653,       CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t2.v2.1.g1
Chr1    11866   12940   AL181U10040.t1.v2.1.m1  100     -       11866   12940   0,0,0   1       1074,   0,      CDS=11866:12940;CDSphase=0;geneID=AL181U10040.t1.v2.1.g1
Chr1    38981   40877   AL1G22550.t1.v2.1.m3    100     -       38981   40877   0,0,0   4       73,152,194,232, 0,154,427,1664, CDS=38981:40877;CDSphase=0;geneID=AL1G22550.t1.v2.1.g3
"""

gffread_bed_contents_collapse = """Chr1    3759    6272    AL1G11530.t1.v2.1.m1    100     +       3759    6272    0,0,0   8       154,281,120,390,153,96,7,101,   0,236,726,946,1414,1679,2326,2412,      CDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t1.v2.1.m1    100     -       6917    8666    0,0,0   9       152,76,67,86,74,46,90,48,96,    0,239,466,646,844,1024,1318,1499,1653,  CDS=6917:8666;CDSphase=0;geneID=AL1G11510.v2.1
Chr1    6917    8666    AL1G11510.t2.v2.1.m1    100     -       6917    8666    0,0,0   8       152,76,67,86,117,90,48,96,      0,239,466,646,844,1318,1499,1653,       CDS=6917:8666;CDSphase=0;geneID=AL1G11510.v2.1
Chr1    11866   12940   AL181U10040.t1.v2.1.m1  100     -       11866   12940   0,0,0   1       1074,   0,      CDS=11866:12940;CDSphase=0;geneID=AL181U10040.t1.v2.1.g1
Chr1    38981   40877   AL1G22550.t1.v2.1.m3    100     -       38981   40877   0,0,0   4       73,152,194,232, 0,154,427,1664, CDS=38981:40877;CDSphase=0;geneID=AL1G22550.t1.v2.1.g3
"""

gffread_bed_contents_collapse_at_end = """Chr1    3759    6272    AL1G11530.t1.v2.1.m1    100     +       3759    6272    0,0,0   8       154,281,120,390,153,96,7,101,   0,236,726,946,1414,1679,2326,2412,      CDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1
Chr1    11866   12940   AL181U10040.t1.v2.1.m1  100     -       11866   12940   0,0,0   1       1074,   0,      CDS=11866:12940;CDSphase=0;geneID=AL181U10040.t1.v2.1.g1
Chr1    38981   40877   AL1G22550.t1.v2.1.m3    100     -       38981   40877   0,0,0   4       73,152,194,232, 0,154,427,1664, CDS=38981:40877;CDSphase=0;geneID=AL1G22550.t1.v2.1.g3
Chr1    6917    8666    AL1G11510.t1.v2.1.m1    100     -       6917    8666    0,0,0   9       152,76,67,86,74,46,90,48,96,    0,239,466,646,844,1024,1318,1499,1653,  CDS=6917:8666;CDSphase=0;geneID=AL1G11510.v2.1
Chr1    6917    8666    AL1G11510.t2.v2.1.m1    100     -       6917    8666    0,0,0   8       152,76,67,86,117,90,48,96,      0,239,466,646,844,1318,1499,1653,       CDS=6917:8666;CDSphase=0;geneID=AL1G11510.v2.1
"""

gffread_bed_empty_line_in_contents = """Chr1    3759    6272    AL1G11530.t1.v2.1.m1    100     +       3759    6272    0,0,0   8       154,281,120,390,153,96,7,101,   0,236,726,946,1414,1679,2326,2412,      CDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t1.v2.1.m1    100     -       6917    8666    0,0,0   9       152,76,67,86,74,46,90,48,96,    0,239,466,646,844,1024,1318,1499,1653,  CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t1.v2.1.g1

Chr1    6917    8666    AL1G11510.t2.v2.1.m1    100     -       6917    8666    0,0,0   8       152,76,67,86,117,90,48,96,      0,239,466,646,844,1318,1499,1653,       CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t2.v2.1.g1
Chr1    11866   12940   AL181U10040.t1.v2.1.m1  100     -       11866   12940   0,0,0   1       1074,   0,      CDS=11866:12940;CDSphase=0;geneID=AL181U10040.t1.v2.1.g1
Chr1    38981   40877   AL1G22550.t1.v2.1.m3    100     -       38981   40877   0,0,0   4       73,152,194,232, 0,154,427,1664, CDS=38981:40877;CDSphase=0;geneID=AL1G22550.t1.v2.1.g3
"""

gffread_bed_no_newline_at_end = """Chr1    3759    6272    AL1G11530.t1.v2.1.m1    100     +       3759    6272    0,0,0   8       154,281,120,390,153,96,7,101,   0,236,726,946,1414,1679,2326,2412,      CDS=3759:6272;CDSphase=0;geneID=AL1G11530.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t1.v2.1.m1    100     -       6917    8666    0,0,0   9       152,76,67,86,74,46,90,48,96,    0,239,466,646,844,1024,1318,1499,1653,  CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t1.v2.1.g1
Chr1    6917    8666    AL1G11510.t2.v2.1.m1    100     -       6917    8666    0,0,0   8       152,76,67,86,117,90,48,96,      0,239,466,646,844,1318,1499,1653,       CDS=6917:8666;CDSphase=0;geneID=AL1G11510.t2.v2.1.g1
Chr1    11866   12940   AL181U10040.t1.v2.1.m1  100     -       11866   12940   0,0,0   1       1074,   0,      CDS=11866:12940;CDSphase=0;geneID=AL181U10040.t1.v2.1.g1
Chr1    38981   40877   AL1G22550.t1.v2.1.m3    100     -       38981   40877   0,0,0   4       73,152,194,232, 0,154,427,1664, CDS=38981:40877;CDSphase=0;geneID=AL1G22550.t1.v2.1.g3"""
