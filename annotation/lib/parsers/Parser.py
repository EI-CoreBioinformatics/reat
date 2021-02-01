import io

from ..parsers import get_handle


class Parser:
    """Generic parser iterator. Base parser class."""

    def __init__(self, handle):
        self._handle = get_handle(handle)

    def __iter__(self):
        return self

    def __next__(self):
        cur = self._handle.tell()
        line = self._handle.readline()
        if line == "" and self._handle.tell() == cur:
            raise StopIteration
        return line

    def __enter__(self):
        if self.closed is True:
            raise ValueError('I/O operation on closed file.')
        return self

    def __exit__(self, *args):
        _ = args
        self._handle.close()

    def close(self):
        """
        Alias for __exit__
        """
        self.__exit__()

    @property
    def name(self):
        """
        Return the filename.
        """
        return self._handle.name

    @property
    def closed(self):
        """
        Boolean flag. If True, the file has been closed already.
        """
        return self._handle.closed


def get_parser(handle, input_format=None):
    """
    Function to recognize the input file type (GFF, GTF, BAM or BED).

    Args:
        handle: (str|io.TextIOWrapper|io.BytesIO|io.BufferedReader|IO)
        :rtype: (Mikado.parsers.GTF.GTF | Mikado.parsers.GFF.GFF3)
        input_format: Optional annotation of the expected file format
    """

    if isinstance(handle, io.TextIOWrapper):
        fname = "-"
    elif isinstance(handle, (io.BytesIO, io.BufferedReader)):
        fname = "-"
        handle = io.TextIOWrapper(handle)
    else:
        fname = handle

    if input_format == "bam" or fname.endswith(".bam"):
        pass
        # return bam_parser.BamParser(handle)
    if input_format == "gtf" or ".gtf" in fname:
        from annotation.lib.parsers import GTF
        return GTF.GTF(handle)
    elif input_format == "gff3" or ".gff" in fname or ".gff3" in fname:
        from annotation.lib.parsers import GFF
        return GFF.GFF3(handle)
    elif input_format == "bed12" or ".bed12" in fname or ".bed" in fname:
        from annotation.lib.parsers import BED
        return BED.Bed12Parser(handle)
    else:
        raise ValueError('Unrecognized format for {}'.format(fname))
