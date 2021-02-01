import os
from unittest import TestCase

import pyfaidx

from annotation.lib.sequence.RefSeq import check_splicing_sites, get_spliced_cds_seq, stop_codons
from annotation.lib.parsers.GFF import GFFReader


class TestRefSeq(TestCase):

    def test_check_splicing_sites(self):
        genes = [g for g in GFFReader("tests/unit/data/refseq/pstrand_c.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/pstrand_c.fa")
        long_intron, short_intron, nc_splicing, intron_len = check_splicing_sites(genome,
                                                                                  genes[0].mrnas["rna-XM_625167.6.m1"],
                                                                                  200000, 59)
        os.remove("tests/unit/data/refseq/pstrand_c.fa.fai")
        assert nc_splicing == 0
        assert short_intron == 1

        genes = [g for g in GFFReader("tests/unit/data/refseq/mstrand_c.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/mstrand_c.fa")
        long_intron, short_intron, nc_splicing, intron_len = check_splicing_sites(genome,
                                                                                  genes[0].mrnas["rna-XM_392640.7.m1"],
                                                                                  200000)
        os.remove("tests/unit/data/refseq/mstrand_c.fa.fai")
        assert nc_splicing == 0
        assert short_intron == 0

        genes = [g for g in GFFReader("tests/unit/data/refseq/pstrand_nc.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/pstrand_nc.fa")
        long_intron, short_intron, nc_splicing, intron_len = check_splicing_sites(genome,
                                                                                  genes[0].mrnas[
                                                                                      "rna-XM_006567179.3.m1"],
                                                                                  200000)
        os.remove("tests/unit/data/refseq/pstrand_nc.fa.fai")
        assert nc_splicing == 1

        genes = [g for g in GFFReader("tests/unit/data/refseq/mstrand_nc.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/mstrand_nc.fa")
        long_intron, short_intron, nc_splicing, intron_len = check_splicing_sites(genome,
                                                                                  genes[0].mrnas[
                                                                                      "rna-XM_026443867.1.m1"],
                                                                                  200000)
        os.remove("tests/unit/data/refseq/mstrand_nc.fa.fai")
        assert nc_splicing == 1

    def test_get_spliced_cds_seq(self):
        genes = [g for g in GFFReader("tests/unit/data/refseq/mstrand_e.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/mstrand_e.fa")
        seq, seq_len, intron_len, last_codon = get_spliced_cds_seq(genome,
                                                                   genes[0].mrnas["rna-XM_033456960.1.m1"],
                                                                   True)
        os.remove("tests/unit/data/refseq/mstrand_e.fa.fai")
        assert last_codon in stop_codons

        genes = [g for g in GFFReader("tests/unit/data/refseq/pstrand_e.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/pstrand_e.fa")
        seq, seq_len, intron_len, last_codon = get_spliced_cds_seq(genome,
                                                                   genes[0].mrnas["rna-XM_033441273.1.m1"],
                                                                   True)
        os.remove("tests/unit/data/refseq/pstrand_e.fa.fai")
        assert last_codon in stop_codons

        genes = [g for g in GFFReader("tests/unit/data/refseq/mstrand_e2.gff")]
        genome = pyfaidx.Fasta("tests/unit/data/refseq/mstrand_e2.fa")
        seq, seq_len, intron_len, last_codon = get_spliced_cds_seq(genome,
                                                                   genes[0].mrnas["rna-XM_392669.6.m1"],
                                                                   True)
        os.remove("tests/unit/data/refseq/mstrand_e2.fa.fai")
        assert last_codon in stop_codons
