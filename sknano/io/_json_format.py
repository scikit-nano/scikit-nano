# -*- coding: utf-8 -*-
"""
======================================================
JSON format (:mod:`sknano.io._json_format`)
======================================================

.. currentmodule:: sknano.io._json_format

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

#from sknano.core import get_fpath

from ._base import StructureIO  # , default_comment_line

__all__ = ['JSONDATA', 'JSONReader', 'JSONWriter']


class JSONReader(StructureIO):
    """Class for reading json file format.

    Parameters
    ----------
    fpath : str
        json file

    """
    def __init__(self, fpath):
        super(JSONReader, self).__init__(fpath=fpath)

        if fpath is not None:
            self.read()

    def read(self):
        pass


class JSONWriter:
    """Class for writing json chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, atoms=None, comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str
        outpath : str, optional
            Output path for structure data file.
        atoms : :class:`~sknano.io.atoms.Atoms`
            An :class:`~sknano.io.atoms.Atoms` instance.
        comment_line : str, optional
            A string written to the first line of `json` file. If `None`,
            then it is set to the full path of the output `json` file.

        """
        pass


class JSONDATA(JSONReader):
    """Class for reading and writing structure data in JSON data format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None):
        super(JSONDATA, self).__init__(fpath=fpath)

    def write(self, jsonfile=None):
        """Write json file.

        Parameters
        ----------
        jsonfile : {None, str}, optional

        """
        pass
