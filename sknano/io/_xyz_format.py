# -*- coding: utf-8 -*-
"""
====================================================
XYZ format (:mod:`sknano.io._xyz_format`)
====================================================

.. currentmodule:: sknano.io._xyz_format

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import os

from .atoms import XYZAtom as Atom, XYZAtoms
from ..core import get_fpath

from ._base import StructureReader, StructureReaderError, \
    StructureWriter, StructureConverter, StructureFormat, default_comment_line

__all__ = ['XYZReader', 'XYZWriter', 'XYZData',
           'XYZ2DATAConverter', 'XYZFormat']


class XYZReader(StructureReader):
    """`StructureReader` class for reading `xyz` chemical file format.

    Parameters
    ----------
    fpath : str
        `xyz` structure file path.

    """
    def __init__(self, fpath=None, **kwargs):
        super(XYZReader, self).__init__(fpath=fpath, **kwargs)

        self._structure_atoms = XYZAtoms()

        if fpath is not None:
            self.read()

    def read(self):
        """Read `xyz` file."""
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
                error_msg = '`xyz` data contained {} atoms '.format(
                    self._structure_atoms.Natoms) + 'but should contain ' + \
                    '{}'.format(Natoms)
                raise StructureReaderError(error_msg)


class XYZWriter(StructureWriter):
    """`StructureWriter` class for writing `xyz` chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, atoms=None,
              comment_line=None, **kwargs):
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
            A string written to the first line of `xyz` file. If `None`,
            then it is set to the full path of the output `xyz` file.

        """
        if fpath is None:
            fpath = get_fpath(fname=fname, ext='xyz', outpath=outpath,
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


class XYZData(XYZReader):
    """Class for reading and writing `StructureData` in `xyz` format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None, **kwargs):
        super(XYZData, self).__init__(fpath=fpath, **kwargs)

    def write(self, xyzfile=None, **kwargs):
        """Write xyz file.

        Parameters
        ----------
        xyzfile : {None, str}, optional

        """
        try:
            kwargs.update(self._kwargs)
            if (xyzfile is None or xyzfile == '') and \
                    (self.fpath is None or self.fpath == ''):
                error_msg = '`xyzfile` must be a string at least 1 ' + \
                    'character long.'
                if xyzfile is None:
                    raise TypeError(error_msg)
                else:
                    raise ValueError(error_msg)
            else:
                xyzfile = self.fpath

            XYZWriter.write(fname=xyzfile, atoms=self._structure_atoms,
                            comment_line=self._comment_line, **kwargs)

        except (TypeError, ValueError) as e:
            print(e)


class XYZ2DATAConverter(StructureConverter):
    """
    `StructureConverter` class for converting `xyz` to `LAMMPS data` format.

    Parameters
    ----------
    xyzfile : str
    boxbounds : {None, dict}, optional
        dict specifying `min` and `max` box bounds in
        :math:`x,y,z` dimensions. If None, then determine bounds based
        atom coordinates.
    pad_box : bool, optional
        If True, after determining minimum simulation box bounds,
        expand :math:`\\pm x,\\pm y,\\pm z` dimensions of simulation box by
        `xpad`, `ypad`, `zpad` distance.
    xpad : float, optional
    ypad : float, optional
    zpad : float, optional

    """
    def __init__(self, xyzfile, **kwargs):

        self._xyzfile = xyzfile
        self._datafile = os.path.splitext(self._xyzfile)[0] + '.data'

        super(XYZ2DATAConverter, self).__init__(
            infile=self._xyzfile, outfile=self._datafile, **kwargs)

        self._new_atoms = []
        self._add_new_atoms = False

        self._new_atomtypes = []
        self._add_new_atomtypes = False

    @property
    def xyzfile(self):
        """`xyz` file name."""
        return self.infile

    @property
    def datafile(self):
        """`LAMMPS data` file name."""
        return self.outfile

    def add_atom(self, atom=None):
        """Add new `Atom` to `Atoms`.

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

    def convert(self, return_reader=False, **kwargs):
        """Convert `xyz` to `LAMMPS data` chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            return an instance of :class:`~DATAReader`

        Returns
        -------
        `DATAReader` (only if `return_reader` is True)

        """
        from .atoms import AtomsConverter
        from ._lammps_data_format import DATAReader, DATAWriter

        kwargs.update(self._kwargs)

        xyzreader = XYZReader(fpath=self.infile, **kwargs)
        atoms = AtomsConverter(atoms=xyzreader.atoms, to='lammps').atoms
        comment_line = xyzreader.comment_line

        if self._add_new_atoms:
            atoms.extend(self._new_atoms)
        if self._add_new_atomtypes:
            atoms.add_atomtypes(self._new_atomtypes)

        DATAWriter.write(fpath=self.outfile, atoms=atoms,
                         comment_line=comment_line, **kwargs)

        if return_reader:
            return DATAReader(fpath=self.outfile, **kwargs)


class XYZFormat(StructureFormat):
    """`StructureFormat` class defining properties for `xyz` format."""
    pass
