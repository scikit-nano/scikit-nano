# -*- coding: utf-8 -*-
"""
=============================================================================
Base classes for structure data (:mod:`sknano.io._structure_io`)
=============================================================================

.. currentmodule:: sknano.io._structure_io

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
from collections import OrderedDict

from sknano.version import short_version as version

default_comment_line = \
    'Structure data generated using scikit-nano version {}'.format(version)
default_structure_format = 'xyz'
supported_structure_formats = ('xyz', 'data')

__all__ = ['StructureIO',
           'StructureReader',
           'StructureWriter',
           'StructureConverter',
           'StructureFormat',
           'StructureIOError',
           'StructureReaderError',
           'StructureWriterError',
           'StructureConverterError',
           'default_comment_line',
           'default_structure_format',
           'supported_structure_formats']


class StructureIO(object):
    """Base class defining common properties for structure data.

    Parameters
    ----------
    fpath : {None, str}, optional

    """
    def __init__(self, fpath=None, **kwargs):
        self._comment_line = default_comment_line
        self._fpath = fpath
        self._structure_atoms = None
        self._kwargs = kwargs

    @property
    def structure_atoms(self):
        """Return structure_atoms attribute"""
        return self._structure_atoms

    @property
    def atoms(self):
        """Alias for :attr:`~StructureIO.structure_atoms`."""
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
    def fpath(self):
        """File path"""
        return self._fpath

    @fpath.setter
    def fpath(self, value=str):
        """Set the file path string.

        Parameters
        ----------
        value : str

        """
        self._fpath = value

    @fpath.deleter
    def fpath(self):
        """Delete file path string"""
        del self._fpath

    @classmethod
    def write(cls, fname, outpath, atoms, structure_format, **kwargs):
        if structure_format == 'data':
            from ._lammps_data_format import DATAWriter
            DATAWriter.write(fname=fname, outpath=outpath, atoms=atoms)
        else:
            from ._xyz_format import XYZWriter
            XYZWriter.write(fname=fname, outpath=outpath, atoms=atoms)


class StructureReader(StructureIO):
    """Abstract base class for reading structure data.

    Parameters
    ----------
    fpath : str
        structure file

    """
    __metaclass__ = ABCMeta

    def __init__(self, fpath=None, **kwargs):
        super(StructureReader, self).__init__(fpath=fpath, **kwargs)

    @abstractmethod
    def read(self):
        """Read structure data from file"""
        return NotImplementedError('Subclasses of `StructureReader` need to '
                                   'implement the `read` method.')


class StructureWriter(object):
    """Abstract base class for writing structure data."""
    __metaclass__ = ABCMeta

    @abstractmethod
    def write(self):
        """Write structure data to file"""
        return NotImplementedError('Subclasses of `StructureWriter` need to '
                                   'implement the `write` method.')


class StructureConverter(object):
    """Abstract base class for converting structure data.

    Parameters
    ----------
    infile : str
    outfile : str

    """
    __metaclass__ = ABCMeta

    def __init__(self, infile=None, outfile=None, **kwargs):
        self._infile = infile
        self._outfile = outfile
        self._kwargs = kwargs

    @property
    def infile(self):
        """Return `infile`"""
        return self._infile

    @property
    def outfile(self):
        """Return `outfile`"""
        return self._outfile

    @abstractmethod
    def convert(self):
        """Convert structure data from one format to another format."""
        return NotImplementedError('Subclasses of `StructureConverter` need '
                                   'to implement the `convert` method.')


class StructureFormat(object):
    """Base class for containing structure format properties"""
    def __init__(self, **kwargs):
        self._properties = OrderedDict()
        self._kwargs = kwargs

    @property
    def properties(self):
        """:class:`~python:collections.OrderedDict` of format properties."""
        return self._properties


class StructureIOError(Exception):
    """Base class for `StructureIO` exceptions."""
    pass


class StructureReaderError(StructureIOError):
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


class StructureWriterError(StructureIOError):
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


class StructureConverterError(StructureIOError):
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
