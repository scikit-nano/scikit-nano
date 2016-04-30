# -*- coding: utf-8 -*-
"""
====================================================
XYZ format (:mod:`sknano.io.xyz`)
====================================================

.. currentmodule:: sknano.io.xyz

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import os

from monty.io import zopen
from sknano.core import get_fpath
from sknano.core.atoms import StructureAtom as Atom

from .base import StructureData, StructureDataError, StructureDataConverter, \
    StructureDataFormatter, default_comment_line

__all__ = ['XYZ', 'XYZData', 'XYZReader', 'XYZWriter', 'XYZFormatter',
           'XYZConverter', 'XYZError', 'XYZIO', 'XYZIOReader', 'XYZIOWriter',
           'XYZIOFormatter', 'XYZIOConverter', 'XYZIOError', 'XYZFormatSpec',
           'XYZ2DATAConverter']


class XYZReader(StructureData):
    """`StructureData` class for reading `xyz` chemical file format.

    Parameters
    ----------
    fpath : str
        `xyz` structure file path.

    """
    def __init__(self, fpath, formatter=None, **kwargs):
        if formatter is None or not isinstance(formatter, XYZFormatter):
            formatter = XYZFormatter()
        super().__init__(fpath=fpath, formatter=formatter, **kwargs)

        if self.fpath is not None:
            self.read()

    def read(self):
        """Read `xyz` file."""
        self.structure.clear()
        with zopen(self.fpath) as stream:
            Natoms = int(stream.readline().strip())
            self.comment_line = stream.readline().strip()
            lines = stream.readlines()
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
                raise StructureDataError(error_msg)

        if len(set(self.atoms.elements)) > 1:
            self.assign_unique_types()

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(self.formatter.todict())
        return attr_dict

XYZIOReader = XYZReader


class XYZWriter:
    """`StructureWriter` class for writing `xyz` chemical file format."""

    @classmethod
    def write(cls, fname=None, outpath=None, fpath=None, structure=None,
              atoms=None, **kwargs):
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

        """
        if structure is None and atoms is None:
            raise ValueError('Expected either `structure` or `atoms` object.')

        if structure is not None and atoms is None:
            atoms = structure.atoms

        if fpath is None:
            fpath = get_fpath(fname=fname, ext='xyz', outpath=outpath,
                              overwrite=True, add_fnum=False)

        xyz = XYZData()
        xyz.write(xyzfile=fpath, atoms=atoms, **kwargs)

XYZIOWriter = XYZWriter


class XYZData(XYZReader):
    """Class for reading and writing `StructureData` in `xyz` format.

    Parameters
    ----------
    fpath : str, optional

    """
    def __init__(self, fpath=None, **kwargs):
        super().__init__(fpath, **kwargs)

    def write(self, xyzfile=None, atoms=None, comment_line=None, **kwargs):
        """Write xyz file.

        Parameters
        ----------
        xyzfile : {None, str}, optional
        comment_line : str, optional
            A string written to the first line of `xyz` file. If `None`,
            then it is set to the full path of the output `xyz` file.

        """
        try:
            kwargs.update(self.kwargs)

            if not xyzfile:
                if self.fpath is None:
                    error_msg = 'Invalid `xyzfile` {}'.format(xyzfile)
                    raise ValueError(error_msg)
                else:
                    xyzfile = self.fpath

            if comment_line is None:
                comment_line = default_comment_line

            if atoms is not None:
                self._atoms = atoms

            super()._update_atoms(**kwargs)

            try:
                with zopen(xyzfile, 'wt') as stream:
                    self._write_header(stream, comment_line)
                    self._write_atoms(stream)
            except OSError as e:
                print(e)

            self._atoms = self._atoms_copy

        except (TypeError, ValueError) as e:
            print(e)

    def _write_header(self, stream, comment_line):
        stream.write('{:d}\n'.format(self.atoms.Natoms))
        stream.write('{}\n'.format(comment_line))

    def _write_atoms(self, stream):
        sformat = self.formatter.format
        for atom in self.atoms:
            stream.write(sformat(atom))

XYZ = XYZIO = XYZData


class XYZFormatter(StructureDataFormatter):
    """`StructureDataFormatter` class defining properties for `xyz` format."""
    def __init__(self, format_string=None):
        if format_string is None:
            format_string = '{:>3s}{:15.8f}{:15.8f}{:15.8f}\n'
        self.format_string = format_string
        self.fmtstr = "format_string={format_string!r}"

    def __str__(self):
        strrep = super().__str__()
        items = ['format_string']
        values = [self.format_string]
        table = self._tabulate(list(zip(items, values)))
        strrep = '\n'.join((strrep, table))
        return strrep

    def format(self, atom):
        """Return :class:`~python:str` of atom attributes formatted for an \
            output stream."""
        return self.format_string.format(atom.symbol, atom.x, atom.y, atom.z)

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        attr_dict = super().todict()
        attr_dict.update(dict(format_string=self.format_string))
        return attr_dict

XYZFormatSpec = XYZIOFormatter = XYZFormatter


class XYZConverter(StructureDataConverter):
    """:class:`StructureDataConverter` class for converting `xyz` data."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._new_atoms = []
        self._add_new_atoms = False

        self._new_types = []
        self._add_new_types = False

    @property
    def xyzfile(self):
        """`xyz` file name."""
        return self.infile

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

XYZIOConverter = XYZConverter


class XYZ2DATAConverter(XYZConverter):
    """:class:`XYZConverter` for converting to `LAMMPS data` format.

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
        datafile = os.path.splitext(xyzfile)[0] + '.data'
        super().__init__(infile=xyzfile, outfile=datafile, **kwargs)

    @property
    def datafile(self):
        """`LAMMPS data` file name."""
        return self.outfile

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
        from .lammps_data import DATAReader, DATAWriter

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


class XYZError(StructureDataError):
    """Exception class for `XYZData` IO errors."""
    pass

XYZIOError = XYZDataError = XYZError
