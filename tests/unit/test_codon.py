from annotation.lib.sequence.Codon import translate


def test_dna():
    dna = "ATGATTCTTTGA"
    assert translate(dna) == "MIL*"


def test_ambiguous_dna():
    dna = "ATGATNCTTTGA"
    assert translate(dna) == "MXL*"


def test_x():
    dna = "ATGAT-CTTTGA"
    assert translate(dna) == "MXL*"


def test_GYG():
    dna = "ATGAGYGTTTGA"
    assert translate(dna) == "MXV*"
