import os
import tempfile

import pytest

from annotation.lib.parsers.GFF import EmptyFileError, GFFReader


def test_invalid_file():
    with pytest.raises(FileNotFoundError):
        _ = GFFReader("foo")


def test_insufficient_permissions_file():
    with tempfile.NamedTemporaryFile("wt") as tmp:
        os.chmod(tmp.name, 0o000)
        with pytest.raises(PermissionError):
            _ = GFFReader(tmp.name)


def test_empty_file():
    with tempfile.NamedTemporaryFile() as tmp:
        with pytest.raises(EmptyFileError):
            for _ in GFFReader(tmp.name):
                pass


def test_invalid_mid_file():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_not_valid_file_content, file=gff_temp)
        gff_temp.flush()
        with pytest.raises(AssertionError):
            for _ in GFFReader(gff_temp.name):
                pass


def test_with_construct():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_content, file=gff_temp)
        gff_temp.flush()
        with GFFReader(gff_temp.name) as reader:
            result = next(reader)
    assert result.start == 2


def test_for_construct():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_content, file=gff_temp)
        gff_temp.flush()
        count = 0
        with GFFReader(gff_temp.name) as parser:
            for _ in parser:
                if count == 1:
                    result = _
                count += 1
        assert count == 3
        assert result.start == 3306
        assert _.uid == "AL1G10030.v2.1"
        assert _.mrnas["AL1G10030.t1.v2.1"].attr["Note"] == "cov:100|id:97.65"
        assert _.mrnas["AL1G10030.t1.v2.1"].fp_utr == 7449
        assert _.mrnas["AL1G10030.t1.v2.1"].tp_utr == 7610


def test_with_directives():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_content_with_directives, file=gff_temp)
        gff_temp.flush()
        count = 0
        for _ in GFFReader(gff_temp.name):
            if count == 1:
                result = _
            count += 1
    assert result.start == 3306
    assert _.uid == "AL1G10030.v2.1"


def test_ends_at_eof():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_ends_at_eof, file=gff_temp)
        gff_temp.flush()
        for _ in GFFReader(gff_temp.name):
            result = _
        assert result.start == 2


def test_ends_mid_comment():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_ends_at_comment, file=gff_temp)
        gff_temp.flush()
        for _ in GFFReader(gff_temp.name):
            result = _
    assert result.start == 2


def test_ends_at_directive():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_ends_at_directive, file=gff_temp)
        gff_temp.flush()
        for _ in GFFReader(gff_temp.name):
            result = _
    assert result.start == 2
    assert result.mrnas[next(iter(result.mrnas.keys()))].attr['Alias'] == ['gene1', 'gene2', 'gnee3']
    assert result.mrnas[next(iter(result.mrnas.keys()))].attr['Ontology_term'] == ['t1', 't2']
    assert result.attr['Dbxref'] == ['db1']


def test_gene_parented():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_gene_parented_CDS, file=gff_temp)
        gff_temp.flush()
        for _ in GFFReader(gff_temp.name):
            pass
        assert _


def test_eden():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_eden, file=gff_temp)
        gff_temp.flush()
        gene_count = 0
        for _ in GFFReader(gff_temp.name):
            gene_count += 1
    assert gene_count == 1
    assert _.uid == "gene00001"
    assert len(_.mrnas) == 3
    assert _.mrnas[next(iter(_.mrnas.keys()))].cds_exons[0].uid == "cds00001"


def test_eden_integer_chr():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_eden_integerableChr, file=gff_temp)
        gff_temp.flush()
        gene_count = 0
        for _ in GFFReader(gff_temp.name):
            gene_count += 1
    assert gene_count == 1
    assert _.uid == "gene00001"
    assert len(_.mrnas) == 3
    assert _.mrnas[next(iter(_.mrnas.keys()))].cds_exons[0].uid == "cds00001"
    _.print_gff()
    for mrna in _.mrnas.values():
        print(mrna.to_bed(cds_only=True))


def test_eden_noattr():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_eden_noattr, file=gff_temp)
        gff_temp.flush()
        gene_count = 0
        for _ in GFFReader(gff_temp.name):
            gene_count += 1
    assert gene_count == 1
    assert _.uid == "gene00001"
    assert len(_.mrnas) == 3
    assert _.mrnas[next(iter(_.mrnas.keys()))].cds_exons[0].uid == "cds00001"
    _.print_gff()


def test_eden_noattr_trailingsemicolon():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_file_eden_noattr_trailingsemicolon, file=gff_temp)
        gff_temp.flush()
        gene_count = 0
        for _ in GFFReader(gff_temp.name):
            gene_count += 1
    assert gene_count == 1
    assert _.uid == "gene00001"
    assert len(_.mrnas) == 3
    assert _.mrnas[next(iter(_.mrnas.keys()))].cds_exons[0].uid == "cds00001"
    _.print_gff()


def test_wrong_number_of_fields():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_wrong_number_of_fields, file=gff_temp)
        gff_temp.flush()
        with pytest.raises(IOError):
            for _ in GFFReader(gff_temp.name):
                pass


def test_dicts():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_ends_at_directive, file=gff_temp)
        gff_temp.flush()
        r1 = GFFReader(gff_temp.name)
        r2 = GFFReader(gff_temp.name)

        r1.uids.add("Test")

        assert r1.uids != r2.uids


def test_ids_are_unique():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as gff_temp:
        print(gff_non_unique_ids, file=gff_temp)
        gff_temp.flush()
        with pytest.raises(Exception):
            for _ in GFFReader(gff_temp.name):
                pass


def test_gff_gnomon():
    num_genes = 0
    total_transcripts = 0
    total_exons = 0
    monoexonic_cds = 0
    monoexonic = 0
    for gene in GFFReader("tests/unit/data/ncbi_annot.gff"):
        num_transcripts = 0
        for mrna_id, mrna in gene.mrnas.items():
            num_coding_exons = len(mrna.cds_exons)
            if num_coding_exons > 0:
                total_exons += len(mrna.exons)
                num_transcripts += 1
                if len(mrna.cds_exons) == 1:
                    monoexonic_cds += 1
                if len(mrna.exons) == 1:
                    monoexonic += 1
        if num_transcripts > 0:
            num_genes += 1
            total_transcripts += num_transcripts
    assert num_genes == 3
    assert total_transcripts == 3
    assert monoexonic_cds == 2
    assert monoexonic == 2
    assert total_exons == 11

    s = gene.mrnas["rna-XM_003689506.3"].to_bed(cds_only=True)
    bed_line = "NW_003789112.1\t4846\t25552\trna-XM_003689506.3\t0\t-\t4846\t25552\t0,0,0\t9\t38,131,245,228,185,186,198,115,315\t0,548,1486,3275,5756,6363,7036,16397,20391"
    assert s.split('\t')[:12] == bed_line.split('\t')[:12]


def test_mitochondrial_gff():
    genes = {}
    for gene in GFFReader("tests/unit/data/mitochondrial.gff"):
        genes[gene.uid] = gene

    assert genes['gene-CYTB'] is not None
    assert len(genes['gene-CYTB'].mrnas) > 0
    assert len(genes['gene-CYTB'].mrnas['gene-CYTB'].exons) > 0
    assert genes['gene-CYTB'].mrnas['gene-CYTB'].exons[0].uid[:5] == 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].cds_exons[0].uid[:5] != 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].exons[0].uid[:5] == 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].cds_exons[0].uid[:5] != 'virt_'


def test_int_exon_mitochondrial_gff():
    genes = {}
    for gene in GFFReader("tests/unit/data/int_mitochondrial.gff"):
        genes[gene.uid] = gene

    assert genes['gene-CYTB'] is not None
    assert len(genes['gene-CYTB'].mrnas) > 0
    assert len(genes['gene-CYTB'].mrnas['gene-CYTB'].exons) > 0
    assert genes['gene-CYTB'].mrnas['gene-CYTB'].exons[0].uid[:5] == 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].cds_exons[0].uid[:5] != 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].exons[0].uid[:5] == 'virt_'
    assert genes['gene-ND1'].mrnas['gene-ND1'].cds_exons[0].uid[:5] != 'virt_'


def test_cds_w_no_phase():
    with tempfile.NamedTemporaryFile("wt", suffix=".gff") as tmp:
        print(gff_CDS_w_no_phase, file=tmp)
        tmp.flush()
        for _ in GFFReader(tmp.name):
            pass

        for mrna in _.mrnas.values():
            for cds in mrna.cds_exons:
                assert cds.phase != '.'


gff_header_only = "##gff-version 3"

gff_wrong_number_of_fields = \
    """scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1\tName=AL1G10010\tancestorIdentifier=918720.v1
"""

gff_not_valid_file_content = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\terror@Q_$(G)QPV\t2444\t25\tt41\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.invalid.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
"""

gff_ends_at_comment = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
###
"""

gff_ends_at_eof = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745"""

gff_ends_at_directive = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1;Dbxref=db1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1;Alias=gene1, gene2, gnee3;Ontology_term=t1,t2
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
"""

gff_non_unique_ids = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
"""

gff_ends_at_directive_no_header = \
    """scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
"""

gff_file_content = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2444\t2503\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tfive_prime_UTR\t2504\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.five_prime_UTR.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t2124\t2347\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.2;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2124\t2347\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.2;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t1803\t2035\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.3;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t1803\t2035\t.\t-\t1\tID=AL1G10010.t1.v2.1.CDS.3;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t1423\t1642\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.4;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t1423\t1642\t.\t-\t2\tID=AL1G10010.t1.v2.1.CDS.4;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t407\t782\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.5;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t407\t782\t.\t-\t1\tID=AL1G10010.t1.v2.1.CDS.5;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t145\t252\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.6;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t145\t252\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.6;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t2\t61\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.7;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2\t61\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.7;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tgene\t3306\t6161\t.\t-\t.\tID=AL1G10020.v2.1;Name=AL1G10020;ancestorIdentifier=470048.v1
scaffold_1\tphytozomev11\tmRNA\t3306\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1;Name=AL1G10020.t1;pacid=35931978;longest=1;ancestorIdentifier=470048.v1.107;Parent=AL1G10020.v2.1
scaffold_1\tphytozomev11\texon\t5631\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t5631\t6123\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tfive_prime_UTR\t6124\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1.five_prime_UTR.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t5319\t5538\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.2;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t5319\t5538\t.\t-\t2\tID=AL1G10020.t1.v2.1.CDS.2;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t4851\t5228\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.3;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t4851\t5228\t.\t-\t1\tID=AL1G10020.t1.v2.1.CDS.3;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t4164\t4362\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.4;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t4164\t4362\t.\t-\t1\tID=AL1G10020.t1.v2.1.CDS.4;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t3810\t4082\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.5;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t3810\t4082\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.5;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t3306\t3708\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.6;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tthree_prime_UTR\t3306\t3528\t.\t-\t.\tID=AL1G10020.t1.v2.1.three_prime_UTR.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t3529\t3708\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.6;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tgene\t7429\t7630\t.\t+\t.\tID=AL1G10030.v2.1;Name=AL1G10030;ancestorIdentifier=918722.v1
scaffold_1\tphytozomev11\ttranscript\t7429\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1;Name=AL1G10030.t1;pacid=35935000;longest=1;ancestorIdentifier=918722.v1.107;Parent=AL1G10030.v2.1;Note=cov:100|id:97.65
scaffold_1\tphytozomev11\texon\t7429\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1.exon.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tfive_prime_UTR\t7429\t7448\t.\t+\t.\tID=AL1G10030.t1.v2.1.five_prime_UTR.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tCDS\t7449\t7610\t.\t+\t0\tID=AL1G10030.t1.v2.1.CDS.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tthree_prime_UTR\t7611\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1.three_prime_UTR.1;Parent=AL1G10030.t1.v2.1;pacid=35935000"""

gff_file_content_with_directives = \
    """##gff-version 3
##annot-version	v2.1
##species Arabidopsis lyrata
scaffold_1\tphytozomev11\tgene\t2\t2541\t.\t-\t.\tID=AL1G10010.v2.1;Name=AL1G10010;ancestorIdentifier=918720.v1
scaffold_1\tphytozomev11\tmRNA\t2\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1;Name=AL1G10010.t1;pacid=35930745;longest=1;ancestorIdentifier=918720.v1.107;Parent=AL1G10010.v2.1
scaffold_1\tphytozomev11\texon\t2444\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2444\t2503\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tfive_prime_UTR\t2504\t2541\t.\t-\t.\tID=AL1G10010.t1.v2.1.five_prime_UTR.1;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t2124\t2347\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.2;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2124\t2347\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.2;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t1803\t2035\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.3;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t1803\t2035\t.\t-\t1\tID=AL1G10010.t1.v2.1.CDS.3;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t1423\t1642\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.4;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t1423\t1642\t.\t-\t2\tID=AL1G10010.t1.v2.1.CDS.4;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t407\t782\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.5;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t407\t782\t.\t-\t1\tID=AL1G10010.t1.v2.1.CDS.5;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t145\t252\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.6;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t145\t252\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.6;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\texon\t2\t61\t.\t-\t.\tID=AL1G10010.t1.v2.1.exon.7;Parent=AL1G10010.t1.v2.1;pacid=35930745
scaffold_1\tphytozomev11\tCDS\t2\t61\t.\t-\t0\tID=AL1G10010.t1.v2.1.CDS.7;Parent=AL1G10010.t1.v2.1;pacid=35930745
###
scaffold_1\tphytozomev11\tgene\t3306\t6161\t.\t-\t.\tID=AL1G10020.v2.1;Name=AL1G10020;ancestorIdentifier=470048.v1
scaffold_1\tphytozomev11\tmRNA\t3306\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1;Name=AL1G10020.t1;pacid=35931978;longest=1;ancestorIdentifier=470048.v1.107;Parent=AL1G10020.v2.1
scaffold_1\tphytozomev11\texon\t5631\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t5631\t6123\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tfive_prime_UTR\t6124\t6161\t.\t-\t.\tID=AL1G10020.t1.v2.1.five_prime_UTR.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t5319\t5538\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.2;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t5319\t5538\t.\t-\t2\tID=AL1G10020.t1.v2.1.CDS.2;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t4851\t5228\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.3;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t4851\t5228\t.\t-\t1\tID=AL1G10020.t1.v2.1.CDS.3;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t4164\t4362\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.4;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t4164\t4362\t.\t-\t1\tID=AL1G10020.t1.v2.1.CDS.4;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t3810\t4082\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.5;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t3810\t4082\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.5;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\texon\t3306\t3708\t.\t-\t.\tID=AL1G10020.t1.v2.1.exon.6;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tthree_prime_UTR\t3306\t3528\t.\t-\t.\tID=AL1G10020.t1.v2.1.three_prime_UTR.1;Parent=AL1G10020.t1.v2.1;pacid=35931978
scaffold_1\tphytozomev11\tCDS\t3529\t3708\t.\t-\t0\tID=AL1G10020.t1.v2.1.CDS.6;Parent=AL1G10020.t1.v2.1;pacid=35931978
###
scaffold_1\tphytozomev11\tgene\t7429\t7630\t.\t+\t.\tID=AL1G10030.v2.1;Name=AL1G10030;ancestorIdentifier=918722.v1
scaffold_1\tphytozomev11\tmRNA\t7429\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1;Name=AL1G10030.t1;pacid=35935000;longest=1;ancestorIdentifier=918722.v1.107;Parent=AL1G10030.v2.1
scaffold_1\tphytozomev11\texon\t7429\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1.exon.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tfive_prime_UTR\t7429\t7448\t.\t+\t.\tID=AL1G10030.t1.v2.1.five_prime_UTR.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tCDS\t7449\t7610\t.\t+\t0\tID=AL1G10030.t1.v2.1.CDS.1;Parent=AL1G10030.t1.v2.1;pacid=35935000
scaffold_1\tphytozomev11\tthree_prime_UTR\t7611\t7630\t.\t+\t.\tID=AL1G10030.t1.v2.1.three_prime_UTR.1;Parent=AL1G10030.t1.v2.1;pacid=35935000"""

gff_file_eden = """##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN
ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3
ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003
ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002
ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003
ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
"""


gff_file_eden_integerableChr = """##gff-version 3.1.26
##sequence-region ctg123 1 1497228
1\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN
1\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001
1\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1
1\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2
1\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3
1\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003
1\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002
1\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003
1\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
1\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
1\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
1\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
1\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
1\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1
1\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
1\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
1\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2
1\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
1\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
1\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3
1\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
1\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
1\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4
"""


gff_file_eden_noattr = """##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001
ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001
ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001
ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003
ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002
ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003
ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001
ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002
ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003
ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003
"""


gff_file_eden_noattr_trailingsemicolon = """##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;
ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;
ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;
ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;
ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003
ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002
ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003
ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;
ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001;
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001;
ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;
ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;
ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002;
ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003;
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003;
ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003;
ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003;
ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003;
"""


gff_gene_parented_CDS = """##gff-version 3
NC_000014.9\tCurated Genomic\tgene\t105741473\t105743070\t.\t-\t.\tID=gene28689;geneID=gene28689;gene_name=IGHG1;Dbxref=GeneID:3500,HGNC:HGNC:5525,IMGT/GENE-DB:IGHG1,MIM:147100;Name=IGHG1;description=immuno
NC_000014.9\tCurated Genomic\tCDS\t105741473\t105741795\t.\t-\t2\tParent=gene28689;Dbxref=GeneID:3500,HGNC:HGNC:5525,IMGT/GENE-DB:IGHG1,MIM:147100
NC_000014.9\tCurated Genomic\tCDS\t105741893\t105742222\t.\t-\t2\tParent=gene28689;Dbxref=GeneID:3500,HGNC:HGNC:5525,IMGT/GENE-DB:IGHG1,MIM:147100
NC_000014.9\tCurated Genomic\tCDS\t105742341\t105742385\t.\t-\t2\tParent=gene28689;Dbxref=GeneID:3500,HGNC:HGNC:5525,IMGT/GENE-DB:IGHG1,MIM:147100
NC_000014.9\tCurated Genomic\tCDS\t105742777\t105743070\t.\t-\t2\tParent=gene28689;Dbxref=GeneID:3500,HGNC:HGNC:5525,IMGT/GENE-DB:IGHG1,MIM:147100
"""


gff_CDS_w_no_phase = """##gff-version 3
NC_000014.9\tCuratedGenomic\tgene\t105765914\t105771405\t.\t-\t.\tID=gene28690;geneID=gene28690;gene_name=IGHG3
NC_000014.9\tCuratedGenomic\tCDS\t105765914\t105765997\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105767813\t105767943\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105769237\t105769559\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105769657\t105769986\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105770105\t105770149\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105770293\t105770337\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105770481\t105770525\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105770669\t105770719\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105771112\t105771405\t.\t-\t.\tParent=gene28690
NC_000014.9\tCuratedGenomic\tCDS\t105765914\t105765997\t.\t-\t0\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105767813\t105767943\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105769237\t105769559\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105769245\t105769559\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105769657\t105769986\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105769657\t105769986\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770105\t105770149\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770105\t105770149\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770293\t105770337\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770293\t105770337\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770481\t105770525\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770481\t105770525\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770669\t105770719\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105770669\t105770719\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105771112\t105771405\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
NC_000014.9\tCuratedGenomic\tCDS\t105771112\t105771405\t.\t-\t2\tParent=gene28690;Dbxref=GeneID:3502,HGNC:HGNC
"""