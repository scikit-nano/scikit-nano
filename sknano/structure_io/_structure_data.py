# -*- coding: utf-8 -*-
"""
=============================================================================
Base classes for structure data (:mod:`sknano.structure_io._structure_data`)
=============================================================================

.. currentmodule:: sknano.structure_io._structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from abc import ABCMeta, abstractmethod
from collections import OrderedDict

from ..chemistry import Atoms
from sknano.version import short_version as version

default_comment_line = \
    'Structure data generated using scikit-nano version {}'.format(version)
default_structure_format = 'xyz'
supported_structure_formats = ('xyz', 'data')

__all__ = ['StructureData',
           'StructureReader',
           'StructureWriter',
           'StructureConverter',
           'StructureFormat',
           'StructureDataError',
           'StructureReaderError',
           'StructureWriterError',
           'StructureConverterError',
           'StructureFormatError',
           'default_comment_line',
           'default_structure_format',
           'supported_structure_formats']


class StructureData(object):
    """Base class defining common properties for structure data.

    Parameters
    ----------
    fname : {None, str}, optional

    """
    def __init__(self, fname=None):
        self._structure_atoms = Atoms()
        self._comment_line = default_comment_line
        self._fname = fname
        self._properties = OrderedDict()

    @property
    def structure_atoms(self):
        """:py:class:`~sknano.chemistry.Atoms` instance."""
        return self._structure_atoms

    @property
    def atoms(self):
        """Alias for :py:attr:`~StructureData.structure_atoms`."""
        return self.structure_atoms

    @property
    def comment_line(self):
        """Comment line."""
        return self._comment_line

    @comment_line.setter
    def comment_line(self, value=str):
        """Set the comment line string.

        Parameters
        ----------
        value : str

        """
        self._comment_line = value

    @property
    def fname(self):
        """File name"""
        return self._fname

    @fname.setter
    def fname(self, value=str):
        """Set the file name string.

        Parameters
        ----------
        value : str

        """
        self._fname = value

    @fname.deleter
    def fname(self):
        """Delete file name string"""
        del self._fname

    @property
    def Natoms(self):
        """Number of atoms in :py:class:`~sknano.chemistry.Atoms`"""
        return self._structure_atoms.Natoms

    @property
    def properties(self):
        """OrderedDict of format properties."""
        return self._properties


class StructureReader(StructureData):
    __metaclass__ = ABCMeta
    """Abstract superclass for reading structure data.

    Parameters
    ----------
    fname : str
        structure file

    """
    def __init__(self, fname=None):
        super(StructureReader, self).__init__(fname=fname)

    @abstractmethod
    def read(self):
        """Read structure data from file"""
        return NotImplementedError('Subclasses of `StructureReader` need to '
                                   'implement the `read` method.')


class StructureWriter(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for writing structure data."""

    @abstractmethod
    def write(self):
        """Write structure data to file"""
        return NotImplementedError('Subclasses of `StructureWriter` need to '
                                   'implement the `write` method.')


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
        return NotImplementedError('Subclasses of `StructureConverter` need '
                                   'to implement the `convert` method.')


class StructureFormat(object):
    """Base class for containing structure format properties"""
    def __init__(self):
        self._properties = OrderedDict()

    @property
    def properties(self):
        """OrderedDict of format properties."""
        return self._properties


class StructureFormatter(object):
    """Base class for wrapping `Atom` attributes into a formatted string."""

    def format(self):
        pass


class StructureDataError(Exception):
    """Base class for `StructureData` exceptions."""
    pass


class StructureReaderError(StructureDataError):
    """Exception raised for `StructureReader` errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class StructureWriterError(StructureDataError):
    """Exception raised for `StructureWriter` errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class StructureConverterError(StructureDataError):
    """Exception raised for `StructureConverter` errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class StructureFormatError(StructureDataError):
    """Exception raised for structure format errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)
