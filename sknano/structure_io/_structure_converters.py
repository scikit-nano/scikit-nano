# -*- coding: utf-8 -*-
"""
===============================================================================
Classes for converting structure data (:mod:`sknano.structure_io._converters`)
===============================================================================

.. currentmodule:: sknano.structure_io._converters

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

import os
#import re

from abc import ABCMeta, abstractmethod

__all__ = ['StructureConverter',
           'FormatConverterError',
           'DATA2XYZConverter',
           'XYZ2DATAConverter',
           'StructureConverterError']


class StructureConverterError(Exception):
    """Base class for StructureConverter exceptions."""
    pass


class FormatConverterError(Exception):
    """Exception raised for errors in the file format."""
    pass


class StructureConverter(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for converting structure data."""
    def __init__(self, infile=None, outfile=None):
        self._infile = infile
        self._outfile = outfile

    @property
    def infile(self):
        return self._infile

    @property
    def outfile(self):
        return self._outfile

    @abstractmethod
    def convert(self):
        """Convert structure data from one format to another format."""
        return NotImplemented


class DATA2XYZConverter(StructureConverter):
    """
    Class for converting structure data from LAMMPS ``data`` to ``xyz`` format.

    .. versionadded:: 0.2.9

    Parameters
    ----------
    datafile : str

    """
    def __init__(self, datafile):
        self._datafile = datafile
        self._xyzfile = os.path.splitext(self._datafile)[0] + '.xyz'

        super(DATA2XYZConverter, self).__init__(
            infile=self._datafile, outfile=self._xyzfile)

    @property
    def datafile(self):
        """LAMMPS data file"""
        return self.infile

    @property
    def xyzfile(self):
        """xyz file name."""
        return self.outfile

    def convert(self, return_reader=False):
        """Convert LAMMPS ``data`` to ``xyz`` chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            return an instance of `XYZReader`

        Returns
        -------
        `XYZReader` (only if `return_reader` is True)

        """
        from ._lammps_data_structure_data import DATAReader
        from ._xyz_structure_data import XYZReader, XYZWriter

        datareader = DATAReader(fname=self._datafile)
        atoms = datareader.atoms
        comment_line = datareader.comment_line

        XYZWriter.write(fname=self._xyzfile, atoms=atoms,
                        comment_line=comment_line)

        if return_reader:
            return XYZReader(fname=self._xyzfile)


class XYZ2DATAConverter(StructureConverter):
    """
    Class for converting structure data from ``xyz`` to LAMMPS ``data`` format.

    Parameters
    ----------
    xyzfile : str
    boxbounds : {None, dict}, optional
        dict specifying ``min`` and ``max`` box bounds in
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
    def __init__(self, xyzfile, boxbounds=None, pad_box=True,
                 xpad=10., ypad=10., zpad=10.):

        self._xyzfile = xyzfile
        self._datafile = os.path.splitext(self._xyzfile)[0] + '.data'

        super(XYZ2DATAConverter, self).__init__(
            infile=self._xyzfile, outfile=self._datafile)

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
    def xyzfile(self):
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
        """Convert xyz to LAMMPS data chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            return an instance of :py:class:`~DATAReader`

        Returns
        -------
        `DATAReader` (only if `return_reader` is True)

        """
        from ._lammps_data_structure_data import DATAReader, DATAWriter
        from ._xyz_structure_data import XYZReader

        xyzreader = XYZReader(fname=self._xyzfile)
        atoms = xyzreader.atoms
        comment_line = xyzreader.comment_line
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

        DATAWriter.write(fname=self._datafile, atoms=atoms,
                         boxbounds=boxbounds, comment_line=comment_line,
                         pad_box=self._pad_box, xpad=self._xpad,
                         ypad=self._ypad, zpad=self._zpad)

        if return_reader:
            return DATAReader(fname=self._datafile)
