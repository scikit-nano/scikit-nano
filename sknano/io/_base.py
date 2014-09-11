# -*- coding: utf-8 -*-
"""
=============================================================================
Base classes for structure data (:mod:`sknano.io._base`)
=============================================================================

.. currentmodule:: sknano.io._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

from sknano.version import short_version as version

try:
    from sknano.core.atoms import KDTAtom as Atom, KDTAtoms as Atoms
except ImportError:
    from sknano.core.atoms import XAtom as Atom, XAtoms as Atoms

default_comment_line = \
    'Structure data generated using scikit-nano version {}'.format(version)
default_structure_format = 'xyz'
supported_structure_formats = ('xyz', 'data')

__all__ = ['Atom', 'Atoms',
           'StructureIO',
           'StructureReader',
           'StructureWriter',
           'StructureConverter',
           'StructureIOError',
           'StructureFormatSpec',
           'default_comment_line',
           'default_structure_format',
           'supported_structure_formats']


class StructureIO(object):
    """Base class defining common properties for structure data.

    Parameters
    ----------
    fpath : {None, str}, optional

    """
    def __init__(self, fpath=None, fname=None, **kwargs):
        self.atoms = Atoms()
        self.comment_line = default_comment_line
        if fpath is None and fname is not None:
            fpath = fname
        self.fpath = fpath
        self.kwargs = kwargs

    @property
    def atoms(self):
        """Return :class:`~sknano.core.atoms.Atoms` instance object."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        """Set :attr:`~StructureIO.atoms` attribute."""
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` instance.')
        self._atoms = value

    @atoms.deleter
    def atoms(self):
        del self._atoms

    @property
    def comment_line(self):
        """Comment line."""
        return self._comment_line

    @comment_line.setter
    def comment_line(self, value):
        """Set the comment line string.

        Parameters
        ----------
        value : str

        """
        if not isinstance(value, str):
            raise TypeError('Expected a string.')
        self._comment_line = value

    @comment_line.deleter
    def comment_line(self):
        del self._comment_line


class StructureReader(object):
    @classmethod
    def read(cls, fpath, structure_format=None, **kwargs):
        if fpath.endswith('.data') or structure_format == 'data':
            from ._lammps_data_format import DATAReader
            return DATAReader(fpath, **kwargs)
        else:
            from ._xyz_format import XYZReader
            return XYZReader.read(fpath)


class StructureWriter(object):
    @classmethod
    def write(cls, fname=None, atoms=None, structure_format=None, **kwargs):
        if structure_format == 'data':
            from ._lammps_data_format import DATAWriter
            DATAWriter.write(fname=fname, atoms=atoms, **kwargs)
        else:
            from ._xyz_format import XYZWriter
            XYZWriter.write(fname=fname, atoms=atoms, **kwargs)


class StructureConverter(object):
    """Abstract base class for converting structure data.

    Parameters
    ----------
    infile : str
    outfile : str

    """
    __metaclass__ = ABCMeta

    def __init__(self, infile=None, outfile=None, **kwargs):
        self.infile = infile
        self.outfile = outfile
        self.kwargs = kwargs

    @abstractmethod
    def convert(self):
        """Convert structure data from one format to another format."""
        raise NotImplementedError('Subclasses of `StructureConverter` need '
                                  'to implement the `convert` method.')


class StructureFormatSpec(object):
    """Base class for defining a format specification.

    Parameters
    ----------

    """
    def __init__(self, **kwargs):
        pass


class StructureIOError(Exception):
    """Base class for `StructureIO` errors."""
    pass
