# -*- coding: utf-8 -*-
"""
====================================================
XYZ format (:mod:`sknano.io._xyz_format`)
====================================================

.. currentmodule:: sknano.io._xyz_format

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import os

from monty.io import zopen
from sknano.core import get_fpath

from ._base import Atom, StructureIO, StructureIOError, StructureConverter, \
    StructureFormatSpec, default_comment_line

__all__ = ['XYZReader', 'XYZWriter', 'XYZData', 'XYZFormatSpec', 'XYZIOError',
           'XYZ2DATAConverter']


class XYZReader(StructureIO):
    """`StructureIO` class for reading `xyz` chemical file format.

    Parameters
    ----------
    fpath : str
        `xyz` structure file path.

    """
    def __init__(self, fpath, **kwargs):
        super().__init__(fpath=fpath, **kwargs)

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read `xyz` file."""
        self.structure_data.clear()
        with zopen(self.fpath) as f:
            Natoms = int(f.readline().strip())
            self.comment_line = f.readline().strip()
            lines = f.readlines()
            for line in lines:
                s = line.strip().split()
                if len(s) != 0:
                    atom = \
                        Atom(element=s[0],
                             x=float(s[1]), y=float(s[2]), z=float(s[3]))
                    self.atoms.append(atom)
            if self.atoms.Natoms != Natoms:
                error_msg = '`xyz` data contained {} atoms '.format(
                    self.atoms.Natoms) + 'but should contain ' + \
                    '{}'.format(Natoms)
                raise XYZIOError(error_msg)

        if len(set(self.atoms.elements)) > 1:
            self.assign_unique_types()


class XYZWriter:
    """`StructureWriter` class for writing `xyz` chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, comment_line=None, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
            Output file name.
        outpath : str, optional
            Output file path.
        fpath : str, optional
            Full path (directory path + file name) to output data file.
        atoms : :class:`~sknano.core.atoms.Atoms`
            An :class:`~sknano.core.atoms.Atoms` instance.
        comment_line : str, optional
            A string written to the first line of `xyz` file. If `None`,
            then it is set to the full path of the output `xyz` file.

        """
        if structure is None and atoms is None:
            raise ValueError('Expected either `structure` or `atoms` object.')

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            fpath = get_fpath(fname=fname, ext='xyz', outpath=outpath,
                              overwrite=True, add_fnum=False)
        if comment_line is None:
            comment_line = default_comment_line

        atoms.rezero_coords()

        with zopen(fpath, 'wt') as f:
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
        super().__init__(fpath, **kwargs)

    def write(self, xyzfile=None, **kwargs):
        """Write xyz file.

        Parameters
        ----------
        xyzfile : {None, str}, optional

        """
        try:
            kwargs.update(self.kwargs)
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

            XYZWriter.write(fname=xyzfile, atoms=self.atoms,
                            comment_line=self.comment_line, **kwargs)

        except (TypeError, ValueError) as e:
            print(e)


class XYZIOError(StructureIOError):
    pass


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

        super().__init__(infile=self._xyzfile, outfile=self._datafile,
                         **kwargs)

        self._new_atoms = []
        self._add_new_atoms = False

        self._new_types = []
        self._add_new_types = False

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

    def add_type(self, atom=None):
        """Add new atom type to atom type dictionary.

        Parameters
        ----------
        atom : `Atom`

        """
        self._new_types.append(atom)
        if not self._add_new_types:
            self._add_new_types = True

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
        from ._lammps_data_format import DATAReader, DATAWriter

        kwargs.update(self.kwargs)

        xyzreader = XYZReader(self.infile, **kwargs)
        structure = xyzreader.structure

        if self._add_new_atoms:
            structure.atoms.extend(self._new_atoms)
        if self._add_new_types:
            structure.atoms.add_types(self._new_types)

        DATAWriter.write(fpath=self.outfile, atoms=structure.atoms,
                         comment_line=xyzreader.comment_line, **kwargs)

        if return_reader:
            return DATAReader(self.outfile, **kwargs)


class XYZFormatSpec(StructureFormatSpec):
    """`StructureFormatSpec` class defining properties for `xyz` format."""
    pass
