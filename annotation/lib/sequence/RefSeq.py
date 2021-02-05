import pyfaidx

start_codons = {"ATG"}  # {"ATG", "CTG", "TTG"}
stop_codons = {"TAA", "TGA", "TAG"}


def check_splicing_sites(genome: pyfaidx.Fasta, mrna, max_intron=200000, min_intron=20):
    noncanonical_splicing = 0
    long_intron = 0
    short_intron = 0
    intron_length = 0
    for i in mrna.introns:
        ilen = i[1] - i[0]
        intron_length += ilen
        if ilen > max_intron:
            long_intron += 1
        if ilen < min_intron:
            short_intron += 1

        if mrna.strand == '+':
            donor = genome.get_seq(mrna.chrom, i[0], i[0] + 1).seq
            try:
                acceptor = genome.get_seq(mrna.chrom, i[1] - 1, i[1]).seq
            except pyfaidx.FetchError:
                acceptor = "NN"
        else:
            acceptor = (-genome.get_seq(mrna.chrom, i[0], i[0] + 1)).seq
            try:
                donor = (-genome.get_seq(mrna.chrom, i[1] - 1, i[1])).seq
            except pyfaidx.FetchError:
                donor = "NN"

        # canonical (GT AG), (GC AG), (AT AC)
        if (donor.upper(), acceptor.upper()) not in (("GT", "AG"), ("GC", "AG"), ("AT", "AC")):
            noncanonical_splicing += 1
    return long_intron, short_intron, noncanonical_splicing, intron_length


def get_spliced_cds_seq(genome: pyfaidx.Fasta, mrna, is_alignment=False):
    aln_seq = []
    cds_intron_len = 0
    last_codon = ""
    for e in mrna.cds_exons:
        if mrna.strand == '+':
            aln_seq.append(genome.get_seq(e.chrom, e.start, e.end).seq)
        else:
            aln_seq.append(genome.get_seq(e.chrom, e.start, e.end).reverse.complement.seq)

    if is_alignment:
        last_cds_exon = mrna.cds_exons[-1]
        try:
            if mrna.strand == '-':
                last_codon = (-genome.get_seq(
                    last_cds_exon.chrom, mrna.ces - 3, mrna.ces - 1)).seq.upper()
            else:
                last_codon = genome.get_seq(
                    last_cds_exon.chrom, mrna.ces + 1, mrna.ces + 3).seq.upper()
        except pyfaidx.FetchError:
            last_codon = "NNN"

    for ei in mrna.cds_introns:
        cds_intron_len += ei[1] - ei[0]

    aln_cds = ''.join(aln_seq)
    return aln_cds.upper(), len(aln_cds), cds_intron_len, last_codon
