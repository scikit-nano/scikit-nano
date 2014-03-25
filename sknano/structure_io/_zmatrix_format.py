# -*- coding: utf-8 -*-
"""
============================================================
ZMATRIX format (:mod:`sknano.structure_io._zmatrix_format`)
============================================================

.. currentmodule:: sknano.structure_io._zmatrix_format

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import os

from ..chemistry import Atom
from ..tools import get_fpath

from ._structure_data import StructureReader, StructureReaderError, \
    StructureWriter, StructureConverter, default_comment_line

__all__ = ['ZMATRIXDATA', 'ZMATRIXReader', 'ZMATRIXWriter',
           'ZMATRIX2DATAConverter']


class ZMATRIXReader(StructureReader):
    """Class for reading zmatrix chemical file format.

    Parameters
    ----------
    zmatrixfile : str
        zmatrix structure file

    """
    def __init__(self, fpath=None):
        super(ZMATRIXReader, self).__init__(fpath=fpath)

        if fpath is not None:
            self.read()

    def read(self):
        with open(self.fpath, 'r') as f:
            Natoms = int(f.readline().strip())
            self._comment_line = f.readline().strip()
            lines = f.readlines()
            for line in lines:
                s = line.strip().split()
                if len(s) != 0:
                    atom = \
                        Atom(s[0], x=float(s[1]), y=float(s[2]), z=float(s[3]))
                    self._structure_atoms.append(atom)
            if self._structure_atoms.Natoms != Natoms:
                error_msg = '`zmatrixfile` contained {} atoms '.format(
                    self._structure_atoms.Natoms) + 'but should contain ' + \
                    '{}'.format(Natoms)
                raise StructureReaderError(error_msg)


class ZMATRIXWriter(StructureWriter):
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
        atoms : :py:class:`~sknano.chemistry.Atoms`
            An :py:class:`~sknano.chemistry.Atoms` instance.
        comment_line : str, optional
            A string written to the first line of `zmatrix` file. If `None`,
            then it is set to the full path of the output `zmatrix` file.

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='zmatrix', outpath=outpath,
                              overwrite=True, add_fnum=False)
        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero_coords()

        with open(fpath, 'w') as f:
            f.write('{:d}\n'.format(atoms.Natoms))
            f.write('{}\n'.format(comment_line))
            for atom in atoms:
                f.write('{:>3s}{:15.8f}{:15.8f}{:15.8f}\n'.format(
                    atom.symbol, atom.x, atom.y, atom.z))


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
                                atoms=self._structure_atoms,
                                comment_line=self._comment_line)
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

        self._new_atomtypes = []
        self._add_new_atomtypes = False

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
        atom : `Atom`

        """
        self._new_atoms.append(atom)
        if not self._add_new_atoms:
            self._add_new_atoms = True

    def add_atomtype(self, atom=None):
        """Add new atom type to atom type dictionary.

        Parameters
        ----------
        atom : `Atom`

        """
        self._new_atomtypes.append(atom)
        if not self._add_new_atomtypes:
            self._add_new_atomtypes = True

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
        from ._lammps_data_format import DATAReader, DATAWriter

        zmatrixreader = ZMATRIXReader(fpath=self.infile)
        atoms = zmatrixreader.atoms
        comment_line = zmatrixreader.comment_line
        if self._add_new_atoms:
            atoms.extend(self._new_atoms)
        if self._add_new_atomtypes:
            atoms.add_atomtypes(self._new_atomtypes)

        if self._boxbounds is None:

            boxbounds = {'x': {'min': None, 'max': None},
                         'y': {'min': None, 'max': None},
                         'z': {'min': None, 'max': None}}

            for i, dim in enumerate(('x', 'y', 'z')):
                boxbounds[dim]['min'] = atoms.coords[:, i].min()
                boxbounds[dim]['max'] = atoms.coords[:, i].max()
        else:
            boxbounds = self._boxbounds

        DATAWriter.write(fname=self.outfile, atoms=atoms,
                         boxbounds=boxbounds, comment_line=comment_line,
                         pad_box=self._pad_box, xpad=self._xpad,
                         ypad=self._ypad, zpad=self._zpad)

        if return_reader:
            return DATAReader(fpath=self.outfile)
