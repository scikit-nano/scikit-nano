# -*- coding: utf-8 -*-
"""
============================================================
ZMATRIX format (:mod:`sknano.io._zmatrix_format`)
============================================================

.. currentmodule:: sknano.io._zmatrix_format

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import os

from sknano.core import get_fpath
from ._base import Atom, StructureIO, StructureConverter, \
    StructureFormat, StructureIOError, default_comment_line

__all__ = ['ZMATRIXDATA', 'ZMATRIXReader', 'ZMATRIXWriter',
           'ZMATRIX2DATAConverter']


class ZMATRIXReader(StructureIO):
    """Class for reading zmatrix chemical file format.

    Parameters
    ----------
    zmatrixfile : str
        zmatrix structure file

    """
    def __init__(self, fpath):
        super(ZMATRIXReader, self).__init__(fpath=fpath)

        if fpath is not None:
            self.read()

    def read(self):
        pass


class ZMATRIXWriter:
    """Class for writing zmatrix chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, atoms=None,
              comment_line=None):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : :class:`~sknano.io.atoms.Atoms`
            An :class:`~sknano.io.atoms.Atoms` instance.
        comment_line : str, optional
            A string written to the first line of `zmatrix` file. If `None`,
            then it is set to the full path of the output `zmatrix` file.

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='zmatrix', outpath=outpath,
                              overwrite=True, add_fnum=False)
            print('fpath: {}'.format(fpath))

        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero_coords()


class ZMATRIXDATA(ZMATRIXReader):
    """Class for reading and writing structure data in ZMATRIX data format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None):
        super(ZMATRIXDATA, self).__init__(fpath=fpath)

    def write(self, zmatrixfile=None):
        """Write zmatrix file.

        Parameters
        ----------
        zmatrixfile : {None, str}, optional

        """
        try:
            if (zmatrixfile is None or zmatrixfile == '') and \
                    (self.fpath is None or self.fpath == ''):
                error_msg = '`zmatrixfile` must be a string at least 1 ' + \
                    'character long.'
                if zmatrixfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                zmatrixfile = self.fpath
            ZMATRIXWriter.write(fname=zmatrixfile,
                                atoms=self.atoms,
                                comment_line=self.comment_line)
        except (TypeError, ValueError) as e:
            print(e)


class ZMATRIX2DATAConverter(StructureConverter):
    """
    Class for converting structure data from `zmatrix` to LAMMPS `data` format.

    Parameters
    ----------
    zmatrixfile : str

    """
    def __init__(self, zmatrixfile, boxbounds=None, pad_box=True,
                 xpad=10., ypad=10., zpad=10.):

        self._zmatrixfile = zmatrixfile
        self._datafile = os.path.splitext(self._zmatrixfile)[0] + '.data'

        super(ZMATRIX2DATAConverter, self).__init__(
            infile=self._zmatrixfile, outfile=self._datafile)

        self._boxbounds = boxbounds
        self._pad_box = pad_box
        self._xpad = xpad
        self._ypad = ypad
        self._zpad = zpad

        self._new_atoms = []
        self._add_new_atoms = False

        self._new_types = []
        self._add_new_types = False

    @property
    def zmatrixfile(self):
        return self.infile

    @property
    def datafile(self):
        """LAMMPS data file name."""
        return self.outfile

    def add_atom(self, atom=None):
        """Add new atom to atoms.

        Parameters
        ----------
        atom : :class:`~sknano.io.atoms.Atom`

        """
        self._new_atoms.append(atom)
        if not self._add_new_atoms:
            self._add_new_atoms = True

    def add_type(self, atom=None):
        """Add new atom type to atom type dictionary.

        Parameters
        ----------
        atom : :class:`~sknano.io.atoms.Atom`

        """
        self._new_types.append(atom)
        if not self._add_new_types:
            self._add_new_types = True

    def convert(self, return_reader=False):
        """Convert zmatrix to LAMMPS data chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            return an instance of :py:class:`~DATAReader`

        Returns
        -------
        `DATAReader` (only if `return_reader` is True)

        """
        pass
