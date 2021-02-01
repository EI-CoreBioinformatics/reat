from annotation.lib.parsers.BED import BEDReader
from annotation.lib.parsers.GFF import GFFReader


def test_gff_bv():
    genes_bed = []
    genes_gff = []
    for gene_bed in BEDReader("tests/unit/data/bv.bed"):
        genes_bed.append(gene_bed)

    for gene_gff in GFFReader("tests/unit/data/bv.gff"):
        genes_gff.append(gene_gff)

    assert genes_bed == genes_gff

    for gene_bed, gene_gff in zip(genes_bed, genes_gff):
        for mrna_bed, mrna_gff in zip(gene_bed.mrnas.values(), gene_gff.mrnas.values()):
            bed_string_repr = mrna_bed.to_bed().split('\t')
            gff_string_repr = mrna_gff.to_bed().split('\t')
            bed_string_repr = bed_string_repr[0:4] + bed_string_repr[5:12]
            gff_string_repr = gff_string_repr[0:4] + gff_string_repr[5:12]

            assert bed_string_repr == gff_string_repr


def test_gff_crubella():
    genes_bed = []
    genes_gff = []
    for gene_bed in BEDReader("tests/unit/data/crubella.bed"):
        genes_bed.append(gene_bed)

    for gene_gff in GFFReader("tests/unit/data/crubella.gff"):
        genes_gff.append(gene_gff)

    # assert genes_bed == genes_gff

    for gene_bed, gene_gff in zip(genes_bed, genes_gff):
        for mrna_bed, mrna_gff in zip(gene_bed.mrnas.values(), gene_gff.mrnas.values()):
            bed_string_repr = mrna_bed.to_bed().split('\t')
            gff_string_repr = mrna_gff.to_bed(cds_only=True).split('\t')
            bed_string_repr = bed_string_repr[0:4] + bed_string_repr[5:12]
            gff_string_repr = gff_string_repr[0:4] + gff_string_repr[5:12]

            assert bed_string_repr == gff_string_repr


def test_gff_ncbi():
    genes_bed = []
    genes_gff = []
    for gene_bed in BEDReader("tests/unit/data/ncbi_annot.bed"):
        genes_bed.append(gene_bed)

    for gene_gff in GFFReader("tests/unit/data/ncbi_annot.gff"):
        genes_gff.append(gene_gff)

    assert genes_bed == genes_gff

    for gene_bed, gene_gff in zip(genes_bed, genes_gff):
        for mrna_bed, mrna_gff in zip(gene_bed.mrnas.values(), gene_gff.mrnas.values()):
            bed_string_repr = mrna_bed.to_bed().split('\t')
            gff_string_repr = mrna_gff.to_bed().split('\t')
            bed_string_repr = bed_string_repr[0:4] + bed_string_repr[5:12]
            gff_string_repr = gff_string_repr[0:4] + gff_string_repr[5:12]

            assert bed_string_repr == gff_string_repr

    # assert bed[:12] == gff[:12]
    # bed_fields = bed.split()[:12]
    # original_fields = "NW_003789112.1	4548	25766	rna-XM_003689506.3	100	-	4846	25552	0,0,0	9	" \
    #                   "336,131,245,228,185,186,198,115,529,	0,846,1784,3573,6054,6661,7334,16695,20689,".split()[:12]
    # assert bed_fields == original_fields


def test_bter():
    genes_bed = []
    genes_gff = []
    for gene_bed in BEDReader("tests/unit/data/bter_top.bed"):
        genes_bed.append(gene_bed)

    for gene_gff in GFFReader("tests/unit/data/bter_top.gff"):
        genes_gff.append(gene_gff)

    assert genes_bed == genes_gff

    for gene_bed, gene_gff in zip(genes_bed, genes_gff):
        for mrna_bed, mrna_gff in zip(gene_bed.mrnas.values(), gene_gff.mrnas.values()):
            bed_string_repr = mrna_bed.to_bed().split('\t')
            gff_string_repr = mrna_gff.to_bed().split('\t')
            bed_string_repr = bed_string_repr[0:4] + bed_string_repr[5:12]
            gff_string_repr = gff_string_repr[0:4] + gff_string_repr[5:12]
            assert bed_string_repr == gff_string_repr


def test_bter_CDS():
    genes_bed = []
    genes_gff = []
    for gene_bed in BEDReader("tests/unit/data/bter_top.coding.bed"):
        genes_bed.append(gene_bed)

    for gene_gff in GFFReader("tests/unit/data/bter_top.gff"):
        genes_gff.append(gene_gff)

    for gene_bed, gene_gff in zip(genes_bed, genes_gff):
        for mrna_bed, mrna_gff in zip(gene_bed.mrnas.values(), gene_gff.mrnas.values()):
            bed_string_repr = mrna_bed.to_bed().split('\t')

            coding_gff_string_repr = mrna_gff.to_bed(cds_only=True).split('\t')
            bed_string_repr = bed_string_repr[0:4] + bed_string_repr[5:12]
            coding_gff_string_repr = coding_gff_string_repr[0:4] + coding_gff_string_repr[5:12]
            assert bed_string_repr == coding_gff_string_repr


if __name__ == "__main__":
    test_gff_crubella()
    test_gff_bv()
    test_gff_ncbi()
