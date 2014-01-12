# -*- coding: utf-8 -*-
"""
=========================================================================
Structure data (:mod:`sknano.structure_io._structure_data`)
=========================================================================

.. currentmodule:: sknano.structure_io._structure_data

"""
from __future__ import division, print_function, absolute_import

__docformat__ = 'restructuredtext'

from abc import ABCMeta, abstractmethod
from collections import OrderedDict

from ..chemistry import Atoms

supported_structure_formats = ('xyz', 'data')
default_structure_format = 'xyz'

__all__ = ['StructureData', 'StructureDataError',
           'StructureReader', 'StructureWriter',
           'StructureReaderError', 'StructureInputError',
           'StructureWriterError', 'StructureOutputError',
           'StructureFormatError', 'supported_structure_formats',
           'default_structure_format']


class StructureData(object):
    """Base class defining common properties for structure data.

    Parameters
    ----------
    fname : {None, str}, optional

    """
    def __init__(self, fname=None):
        self._atoms = Atoms()
        self._comment_line = None
        self._fname = fname
        self._properties = OrderedDict()

    @property
    def atoms(self):
        """:py:class:`Atoms` instance."""
        return self._atoms

    @property
    def comment_line(self):
        """Comment line."""
        return self._comment_line

    @comment_line.setter
    def comment_line(self, value=str):
        """comment line setter.

        Parameters
        ----------
        value : str

        """
        self._comment_line = value

    @property
    def fname(self):
        """fname."""
        return self._fname

    @fname.setter
    def fname(self, value=str):
        """fname setter.

        Parameters
        ----------
        value : str

        """
        self._fname = value

    @fname.deleter
    def fname(self):
        """fname deleter."""
        del self._fname

    @property
    def Natoms(self):
        """Number of atoms."""
        return self._atoms.Natoms

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
        """Read in structure data from file"""
        return NotImplemented


class StructureWriter(object):
    __metaclass__ = ABCMeta
    """Abstract superclass for writing structure data."""

    @abstractmethod
    def write(self):
        """Read in structure data from file"""
        return NotImplemented


class StructureFormatError(Exception):
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


class StructureDataError(Exception):
    """Base class for StructureData exceptions."""
    pass


class StructureReaderError(Exception):
    """Base class for StructureReader exceptions."""
    pass


class StructureInputError(StructureReaderError):
    """Exception raised for input file errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class StructureWriterError(Exception):
    """Base class for StructureWriter exceptions."""
    pass


class StructureOutputError(StructureWriterError):
    """Exception raised for file output errors.

    Parameters
    ----------
    msg : str
        Error message.

    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)
