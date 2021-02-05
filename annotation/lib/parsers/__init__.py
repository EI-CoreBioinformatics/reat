# coding: utf_8

"""
    This module defines the iterators that will parse BED12, GTF, GFF files.
"""
import bz2
import gzip
import io
import sys
from functools import partial

from ..utilities.file_type import filetype


class HeaderError(Exception):
    """
    Mock exception which is raised when a header/comment line (e.g. starting with "#") is found.
    """
    pass


def get_handle(handle, position=None):
    if not isinstance(handle, io.IOBase):
        file_type = filetype(handle)
        if handle.endswith(".gz") or file_type == b"application/gzip":
            opener = gzip.open
        elif handle.endswith(".bz2") or file_type == b"application/x-bzip2":
            opener = bz2.open
        else:
            # TODO: Check how effective is this buffering policy
            opener = partial(open, **{"buffering": 1})
        try:
            handle = opener(handle, "rt")
        except IOError as e:
            print("Cannot open file ", handle.name, file=sys.stderr)
            print("I/O error({0}): {1}".format(e.errno, e.strerror), file=sys.stderr)

    if position is not None:
        handle.seek(position)
    return handle
