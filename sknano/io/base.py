# -*- coding: utf-8 -*-
"""
=============================================================================
Base classes for structure data (:mod:`sknano.io.base`)
=============================================================================

.. currentmodule:: sknano.io.base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod
import copy
import os

from sknano.core import BaseClass, TabulateMixin, get_fpath
from sknano.core.atoms import MDAtoms as Atoms
from sknano.core.structures import StructureBase
# from sknano.utils.analysis import StructureAnalyzer
from sknano.version import version

default_comment_line = \
    'Structure data generated using scikit-nano version {}'.format(version)
default_structure_format = 'xyz'
supported_structure_formats = ('xyz', 'data', 'dump', 'pdb')

__all__ = ['StructureData',
           'StructureDataConverter',
           'StructureDataFormatter',
           'StructureDataError',
           'StructureDataMixin',
           'StructureDataReader',
           'StructureDataWriter',
           'StructureDataReaderMixin',
           'StructureDataWriterMixin',
           'StructureIO',
           'StructureIOConverter',
           'StructureIOFormatter',
           'StructureIOError',
           'StructureIOMixin',
           'StructureIOReader',
           'StructureIOWriter',
           'StructureIOReaderMixin',
           'StructureIOWriterMixin',
           'StructureReader',
           'StructureWriter',
           'StructureReaderMixin',
           'StructureWriterMixin',
           'StructureFormatSpec',
           'default_comment_line',
           'default_structure_format',
           'supported_structure_formats']


class StructureDataReaderMixin:
    """Mixin class for reading structure data."""
    def read(self, *args, fname=None, structure_format=None, **kwargs):
        """Read structure data."""
        if fname is None:
            if len(args) == 0:
                raise ValueError('`fname` is required')
            else:
                fname = list(args).pop()

        if fname.endswith(supported_structure_formats) and \
                structure_format is None:
            for ext in supported_structure_formats:
                if fname.endswith(ext):
                    structure_format = ext
                    break

        if not fname.endswith(structure_format):
            fname += '.' + structure_format

        if structure_format is None:
            raise ValueError(('Unknown structure format: ' +
                              '{}'.format(structure_format)))

        getattr(self, 'read_' + structure_format)(fname=fname, **kwargs)

    def read_data(self, *args, **kwargs):
        """Read LAMMPS Data file.

        Returns
        -------
        :class:`~sknano.io.DATAReader`

        """
        from sknano.io import DATAReader
        return DATAReader(*args, **kwargs)

    def read_dump(self, *args, **kwargs):
        """Read LAMMPS Dump file.

        Returns
        -------
        :class:`~sknano.io.DUMPReader`

        """
        from sknano.io import DUMPReader
        return DUMPReader(*args, **kwargs)

    def read_pdb(self, *args, **kwargs):
        """Read PDB file.

        Returns
        -------
        :class:`~sknano.io.PDBReader`

        """
        from sknano.io import PDBReader
        return PDBReader(*args, **kwargs)

    def read_xyz(self, *args, **kwargs):
        """Read XYZ file.

        Returns
        -------
        :class:`~sknano.io.XYZReader`

        """
        from sknano.io import XYZReader
        return XYZReader.read(*args, **kwargs)

StructureReaderMixin = StructureIOReaderMixin = StructureDataReaderMixin


class StructureDataWriterMixin:
    """Mixin class for saving structure data."""
    def _update_atoms(self, deepcopy=True, center_centroid=True,
                      center_com=False, region_bounds=None,
                      filter_condition=None, rotation_parameters=None,
                      **kwargs):
        atoms_copy = self._atoms[:]
        if deepcopy:
            atoms_copy = copy.deepcopy(atoms_copy)

        self._atoms_copy = atoms_copy
        atoms = self._atoms

        if any([kw in kwargs for kw
                in ('center_CM', 'center_center_of_mass')]):
            if 'center_CM' in kwargs:
                center_com = kwargs['center_CM']
                del kwargs['center_CM']
            else:
                center_com = kwargs['center_center_of_mass']
                del kwargs['center_center_of_mass']

        if region_bounds is not None:
            atoms.clip_bounds(region_bounds)

        if center_centroid:
            atoms.center_centroid()
        elif center_com:
            atoms.center_com()

        if filter_condition is not None:
            atoms.filter(filter_condition)
            # atoms = atoms.filtered(filter_condition)

        rotation_kwargs = ['rotation_angle', 'angle', 'rot_axis', 'axis',
                           'anchor_point', 'deg2rad', 'degrees', 'rot_point',
                           'from_vector', 'to_vector', 'transform_matrix']

        if rotation_parameters is None and \
                any([kw in kwargs for kw in rotation_kwargs]):
            rotation_parameters = {kw: kwargs[kw] for kw in rotation_kwargs
                                   if kw in kwargs}
            if 'rotation_angle' in rotation_parameters:
                rotation_parameters['angle'] = \
                    rotation_parameters['rotation_angle']
                del rotation_parameters['rotation_angle']
            if 'rot_axis' in rotation_parameters:
                rotation_parameters['axis'] = rotation_parameters['rot_axis']
                del rotation_parameters['rot_axis']
            if 'deg2rad' in rotation_parameters:
                rotation_parameters['degrees'] = rotation_parameters['deg2rad']
                del rotation_parameters['deg2rad']

            kwargs = {k: v for k, v in kwargs.items()
                      if k not in rotation_kwargs}

        if rotation_parameters is not None and \
                isinstance(rotation_parameters, dict):
            atoms.rotate(**rotation_parameters)

        atoms.rezero()
        self._atoms = atoms

    def _update_fpath(self, kwargs):
        fname = kwargs.get('fname', None)
        outpath = kwargs.get('outpath', None)
        structure_format = kwargs.get('structure_format', None)
        if fname.endswith(supported_structure_formats) and \
                structure_format is None:
            for ext in supported_structure_formats:
                if fname.endswith(ext):
                    structure_format = ext
                    break
        elif structure_format is None or \
            structure_format not in supported_structure_formats or \
            (not fname.endswith(supported_structure_formats) and
             structure_format not in supported_structure_formats):
            structure_format = default_structure_format

        self.structure_format = kwargs['structure_format'] = structure_format

        if not fname.endswith(structure_format):
            fname += '.' + structure_format
        self.fname = kwargs['fname'] = fname

        if outpath is not None:
            fpath = os.path.join(outpath, fname)
        else:
            fpath = os.path.join(os.getcwd(), fname)
        self.fpath = fpath

    def write(self, *args, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        outpath : str, optional
            Output path for structure data file.
        structure_format : {None, 'xyz', 'pdb', 'data', 'dump'}, optional
            chemical file format of saved structure data. Must be one of:

                - `xyz`
                - `pdb`
                - `data`
                - `dump`

            If `None`, then guess based on `fname` file extension or
            default to `xyz` format.
        center_centroid : bool, optional
            Center centroid on origin.
        center_com : bool, optional
            Center center-of-mass on origin.
        region_bounds : :class:`GeometricRegion`, optional
            An instance of a
            :class:`~sknano.core.geometric_regions.GeometricRegion` having
            a method
            :meth:`~sknano.core.geometric_regions.GeometricRegion.contains`
            to filter the atoms not contained within the
            `GeometricRegion`.
        filter_condition : array_like, optional

        """
        self._update_fpath(kwargs)
        # self._update_atoms(**kwargs)
        structure_format = kwargs.pop('structure_format')
        structure = kwargs.pop('structure', self)
        print('writing {}'.format(self.fname))
        getattr(self, 'write_' + structure_format)(structure=structure,
                                                   **kwargs)

    def write_data(self, **kwargs):
        """Write LAMMPS data file."""
        from sknano.io import DATAWriter
        DATAWriter.write(**kwargs)

    def write_dump(self, **kwargs):
        """Write LAMMPS dump file."""
        from sknano.io import DUMPWriter
        DUMPWriter.write(**kwargs)

    def write_pdb(self, **kwargs):
        """Write pdb file."""
        from sknano.io import PDBWriter
        PDBWriter.write(**kwargs)

    def write_xyz(self, **kwargs):
        """Write xyz file."""
        from sknano.io import XYZWriter
        XYZWriter.write(**kwargs)

StructureWriterMixin = StructureIOWriterMixin = StructureDataWriterMixin


class StructureDataMixin(StructureDataReaderMixin, StructureDataWriterMixin):
    """Mixin class providing I/O methods for structure data."""
    pass

StructureIOMixin = StructureDataMixin


class StructureData(StructureDataMixin, StructureBase, BaseClass):
    """Base class for structure data file input and output.

    Parameters
    ----------
    fpath : {None, str}, optional

    """
    def __init__(self, fpath=None, fname=None, formatter=None, **kwargs):
        super().__init__(**kwargs)

        self._atoms = Atoms()
        self.comment_line = default_comment_line
        if fpath is None and fname is not None:
            fpath = fname
        self.fpath = fpath
        self.formatter = formatter
        self.kwargs = kwargs
        fmtstr = "fpath={fpath!r}"
        if isinstance(formatter, StructureDataFormatter):
            fmtstr = ', '.join((fmtstr, formatter.fmtstr))
        self.fmtstr = fmtstr

    def __str__(self):
        # strrep = self._table_title_str()
        strrep = self._obj_mro_str()
        atoms = self.atoms
        # xtal_cell = self.crystal_cell

        if atoms.data:
            title = '.'.join((strrep, atoms.__class__.__name__))
            strrep = '\n'.join((title, str(atoms)))

        # formatter = self.formatter
        # if formatter is not None:
        #     strrep = '\n'.join((strrep, str(formatter)))
        return strrep

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self is other or self._atoms == other._atoms

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

    def __deepcopy__(self, memo):
        from copy import deepcopy
        cp = self.__class__(None)
        memo[id(self)] = cp
        for attr in dir(self):
            if not attr.startswith('_'):
                setattr(cp, attr, deepcopy(getattr(self, attr), memo))
        return cp

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(fpath=self.fpath, formatter=self.formatter)

StructureIO = StructureData


class StructureDataReader:
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
            from .lammps_data import DATAReader
            return DATAReader(fpath, **kwargs)
        elif (structure_format is None and fpath.endswith('.dump')) or \
                structure_format == 'dump':
            from .lammps_dump import DUMPReader
            return DUMPReader(fpath, **kwargs)
        elif (structure_format is None and fpath.endswith('.pdb')) or \
                structure_format == 'pdb':
            from .pdb import PDBReader
            return PDBReader.read(fpath)
        elif (structure_format is None and fpath.endswith('.xyz')) or \
                structure_format == 'xyz':
            from .xyz import XYZReader
            return XYZReader.read(fpath)
        else:
            raise StructureDataError("Unable to determine `structure_format`")

StructureReader = StructureIOReader = StructureDataReader


class StructureDataWriter:
    """Structure data writer base class."""
    @classmethod
    def write(cls, fname=None, structure_format=None, **kwargs):
        """Write structure data to file.

        Parameters
        ----------
        fname : str, optional
        structure_format : {None, str}, optional

        """
        if fname is None and structure_format is None:
            structure_format = default_structure_format

        if fname is None:
            fname = get_fpath(fname='structure_data', ext=structure_format,
                              add_fnum=True)

        if (structure_format is None and fname.endswith('.data')) or \
                structure_format == 'data':
            from .lammps_data import DATAWriter
            DATAWriter.write(fname=fname, **kwargs)
        elif (structure_format is None and fname.endswith('.dump')) or \
                structure_format == 'dump':
            from .lammps_dump import DUMPWriter
            DUMPWriter.write(fname=fname, **kwargs)
        elif (structure_format is None and fname.endswith('.pdb')) or \
                structure_format == 'pdb':
            from .pdb import PDBWriter
            PDBWriter.write(fname=fname, **kwargs)
        # elif (structure_format is None and fname.endswith('.xyz')) or \
        #        structure_format == 'xyz':
        else:
            from .xyz import XYZWriter
            XYZWriter.write(fname=fname, **kwargs)

StructureWriter = StructureIOWriter = StructureDataWriter


class StructureDataConverter(BaseClass, metaclass=ABCMeta):
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
        self.fmtstr = "infile={infile!r}, outfile={outfile!r}"

    @abstractmethod
    def convert(self):
        """Convert structure data from one format to another format."""
        raise NotImplementedError('Subclasses of `StructureConverter` need '
                                  'to implement the `convert` method.')

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(infile=self.infile, outfile=self.outfile)

StructureConverter = StructureIOConverter = StructureDataConverter


class StructureDataFormatter(TabulateMixin, BaseClass):
    """Base class for defining a format specification."""

    def __str__(self):
        strrep = self._table_title_str()
        objstr = self._obj_mro_str()
        return '\n'.join((strrep, objstr))

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict()

StructureFormatSpec = StructureIOFormatter = StructureDataFormatter


class StructureDataError(Exception):
    """Base class for `StructureData` errors."""
    pass

StructureIOError = StructureDataError
