# -*- coding: utf-8 -*-
"""
=============================================================================
Base classes for structure data (:mod:`sknano.io._base`)
=============================================================================

.. currentmodule:: sknano.io._base

"""
from __future__ import absolute_import, division, print_function
import six
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.utils.analysis import StructureAnalyzer
from sknano.version import short_version as version

default_comment_line = \
    'Structure data generated using scikit-nano version {}'.format(version)
default_structure_format = 'xyz'
supported_structure_formats = ('xyz', 'data', 'dump')

__all__ = ['Atom', 'Atoms',
           'StructureData',
           'StructureIO',
           'StructureReader',
           'StructureWriter',
           'StructureConverter',
           'StructureIOError',
           'StructureFormatSpec',
           'default_comment_line',
           'default_structure_format',
           'supported_structure_formats']


class StructureData(object):
    """Container class for structure data atoms."""
    def __init__(self):
        self.atoms = Atoms()

    @property
    def atoms(self):
        """Return :class:`~sknano.core.atoms.Atoms` instance object."""
        return self._atoms

    @atoms.setter
    def atoms(self, value):
        """Set :attr:`~StructureData.atoms` attribute."""
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` instance.')
        self._atoms = value

    @atoms.deleter
    def atoms(self):
        del self._atoms

    @property
    def bonds(self):
        """Return :attr:`~sknano.core.atoms.Atoms.bonds` attribute."""
        return self.atoms.bonds

    def analyze(self):
        """Analyze structure data."""
        self.atoms.update_attrs()
        #StructureAnalyzer(self)

    def clear(self):
        """Clear :class:`~sknano.core.atoms.Atoms` list."""
        self.atoms.clear()


class StructureIO(object):
    """Base class for structure data file input and output.

    Parameters
    ----------
    fpath : {None, str}, optional

    """
    def __init__(self, fpath=None, fname=None, **kwargs):
        self.structure_data = StructureData()
        self.comment_line = default_comment_line
        if fpath is None and fname is not None:
            fpath = fname
        self.fpath = fpath
        self.kwargs = kwargs

    @property
    def atoms(self):
        return self.structure_data.atoms

    @atoms.setter
    def atoms(self, value):
        """Set :attr:`~StructureData.atoms` attribute."""
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` instance.')
        self.structure_data.atoms = value

    @atoms.deleter
    def atoms(self):
        del self.structure_data.atoms

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
    """Structure data reader base class."""
    @classmethod
    def read(cls, fpath, structure_format=None, **kwargs):
        """Read structure data from file.

        Parameters
        ----------
        fpath : str
        structure_format : {None, str}
        """
        if (structure_format is None and fpath.endswith('.data')) or \
                structure_format == 'data':
            from ._lammps_data_format import DATAReader
            return DATAReader(fpath, **kwargs)
        elif (structure_format is None and fpath.endswith('.dump')) or \
                structure_format == 'dump':
            from ._lammps_dump_format import DUMPReader
            return DUMPReader(fpath, **kwargs)
        elif (structure_format is None and fpath.endswith('.xyz')) or \
                structure_format == 'xyz':
            from ._xyz_format import XYZReader
            return XYZReader.read(fpath)


class StructureWriter(object):
    """Structure data writer base class."""
    @classmethod
    def write(cls, fname=None, atoms=None, structure_format=None, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
        atoms : `~sknano.core.atoms.Atoms` instance
        structure_format : {None, str}, optional

        """
        if structure_format == 'data':
            from ._lammps_data_format import DATAWriter
            DATAWriter.write(fname=fname, atoms=atoms, **kwargs)
        elif structure_format == 'dump':
            from ._lammps_dump_format import DUMPWriter
            DUMPWriter.write(fname=fname, atoms=atoms, **kwargs)
        elif structure_format == 'xyz':
            from ._xyz_format import XYZWriter
            XYZWriter.write(fname=fname, atoms=atoms, **kwargs)


class StructureConverter(six.with_metaclass(ABCMeta, object)):
    """Abstract base class for converting structure data.

    Parameters
    ----------
    infile : str
    outfile : str

    """

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
