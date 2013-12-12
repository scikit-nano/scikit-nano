# -*- coding: utf-8 -*-
"""
=============================================================
Structure converters (:mod:`sknano.structure_io._converters`)
=============================================================

.. currentmodule:: sknano.structure_io._converters

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

import os
#import re

from abc import ABCMeta

#import numpy as np

#from ._readers import StructureReader
#from ._writers import StructureWriter

__all__ = ['StructureConverter',
           'FormatConverterError',
           'DATA2XYZConverter',
           'XYZ2DATAConverter']


class StructureConverter(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for converting structure data."""
    pass
#StructureReader.register(StructureConverter)
#StructureWriter.register(StructureConverter)


class FormatConverterError(Exception):
    """Exception raised for errors in the file format."""
    pass


class DATA2XYZConverter(StructureConverter):
    """
    Class for converting structure data from ``xyz`` to LAMMPS ``data`` format.

    Parameters
    ----------
    datafile : str

    """
    def __init__(self, datafile):
        pass


class XYZ2DATAConverter(StructureConverter):
    """
    Class for converting structure data from ``xyz`` to LAMMPS ``data`` format.

    Parameters
    ----------
    xyzfile : str
    boxbounds : {None, dict}, optional
        dict specifying box bounds
        ..versionadded:: 0.3.32
    pad_box : bool, optional
        pad simulation box bounds
        ..versionadded:: 0.3.32
    xpad : float, optional
    ypad : float, optional
    zpad : float, optional

    """
    def __init__(self, xyzfile, boxbounds=None, pad_box=True,
                 xpad=10., ypad=10., zpad=10.):

        self._xyzfile = xyzfile
        self._datafile = os.path.splitext(self._xyzfile)[0] + '.data'

        self._boxbounds = boxbounds
        self._pad_box = pad_box
        self._xpad = xpad
        self._ypad = ypad
        self._zpad = zpad
        self._boxpad = {'x': xpad, 'y': ypad, 'z': zpad}

        self._new_atoms = []
        self._add_new_atoms = False

        self._new_atomtypes = []
        self._add_new_atomtypes = False

    @property
    def datafile(self):
        """LAMMPS data file name."""
        return self._datafile

    def add_atom(self, atom=None):
        """Add new atom to atoms."""
        self._new_atoms.append(atom)
        self._add_new_atoms = True

    def add_atomtype(self, atom=None):
        """Add new atom type to atom type dictionary."""
        self._new_atomtypes.append(atom)
        self._add_new_atomtypes = True

    def convert(self, return_reader=False):
        """Convert xyz to LAMMPS data chemical file format.

        Parameters
        ----------
        return_reader : bool, optional
            return an instance of :py:class:`~DATAReader`
            ..versionadded:: 0.3.32

        Returns
        -------
        None or DATAReader

        """
        from ._readers import XYZReader, DATAReader
        from ._writers import DATAWriter
        xyzreader = XYZReader(self._xyzfile)
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

        if self._pad_box:
            for i, dim in enumerate(('x', 'y', 'z')):
                boxbounds[dim]['min'] = \
                    boxbounds[dim]['min'] - self._boxpad[dim]
                boxbounds[dim]['max'] = \
                    boxbounds[dim]['max'] + self._boxpad[dim]

        DATAWriter.write(fname=self._datafile, atoms=atoms,
                         boxbounds=boxbounds, comment_line=comment_line)

        if return_reader:
            return DATAReader(self._datafile)
