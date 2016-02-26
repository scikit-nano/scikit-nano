# -*- coding: utf-8 -*-
"""
===============================================================================
Structure generator base module (:mod:`sknano.generators._base`)
===============================================================================

.. currentmodule:: sknano.generators._base

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

import copy
import os

import numpy as np

from sknano.core.atoms import MDAtom as Atom, MDAtoms as Atoms
from sknano.io import default_structure_format, supported_structure_formats

__all__ = ['Atom', 'Atoms', 'GeneratorBase', 'BaseGenerator', 'GeneratorMixin',
           'BulkGeneratorBase', 'BaseBulkGenerator',
           'STRUCTURE_GENERATORS']


#: Tuple of structure generator classes.
STRUCTURE_GENERATORS = ('FullereneGenerator',
                        'GrapheneGenerator',
                        'PrimitiveCellGrapheneGenerator',
                        'ConventionalCellGrapheneGenerator',
                        'BilayerGrapheneGenerator',
                        'UnrolledSWNTGenerator',
                        'SWNTGenerator',
                        'MWNTGenerator',
                        'AlphaQuartzGenerator',
                        'DiamondStructureGenerator',
                        'FCCStructureGenerator',
                        'CaesiumChlorideStructureGenerator',
                        'RocksaltStructureGenerator',
                        'ZincblendeStructureGenerator',
                        'MoS2Generator')


class GeneratorBase(metaclass=ABCMeta):
    """Base structure generator class."""
    def __init__(self, *args, autogen=True, **kwargs):

        super().__init__(*args, **kwargs)

        if autogen:
            self.generate()

    @abstractmethod
    def generate(self):
        """Generate structure data."""
        return NotImplementedError

    def finalize(self):
        """Finalize structure data by assigning unique ids and types to \
            structure atoms."""
        self.assign_unique_ids()
        self.assign_unique_types()

    def save(self, fname=None, outpath=None, structure_format=None,
             deepcopy=True, center_centroid=True, center_com=False,
             region_bounds=None, filter_condition=None,
             rotation_parameters=None, **kwargs):
        """Save structure data.

        Parameters
        ----------
        fname : {None, str}, optional
            file name string
        outpath : str, optional
            Output path for structure data file.
        structure_format : {None, str}, optional
            chemical file format of saved structure data. Must be one of:

                - `xyz`
                - `data`

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

        if not fname.endswith(structure_format):
            fname += '.' + structure_format
        self.fname = fname

        if outpath is not None:
            fpath = os.path.join(outpath, fname)
        else:
            fpath = os.path.join(os.getcwd(), fname)
        self.fpath = fpath

        # self._structure_format = structure_format

        if deepcopy:
            atoms = copy.deepcopy(self.atoms)
        else:
            atoms = self.atoms[:]

        if any([kw in kwargs for kw
                in ('center_CM', 'center_center_of_mass')]):
            if 'center_CM' in kwargs:
                center_com = kwargs['center_CM']
                del kwargs['center_CM']
            else:
                center_com = kwargs['center_center_of_mass']
                del kwargs['center_center_of_mass']

        if center_centroid:
            atoms.center_centroid()
        elif center_com:
            atoms.center_com()

        if region_bounds is not None:
            atoms.clip_bounds(region_bounds)

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
            self.rotate(**rotation_parameters)

        getattr(self, 'write_' + structure_format)(
            fname=fname, outpath=outpath, structure=self, atoms=atoms,
            **kwargs)

        # StructureWriter.write(fname=fname, outpath=outpath,
        #                       structure_format=structure_format,
        #                       structure=self.structure)

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

    def read_xyz(self, *args, **kwargs):
        """Read XYZ file.

        Returns
        -------
        :class:`~sknano.io.XYZReader`

        """
        from sknano.io import XYZReader
        return XYZReader.read(*args, **kwargs)

    def write_data(self, **kwargs):
        """Write LAMMPS DATA file."""
        from sknano.io import DATAWriter
        if 'structure' not in kwargs:
            kwargs['structure'] = self
        DATAWriter.write(**kwargs)

    def write_dump(self, **kwargs):
        """Write LAMMPS dump file."""
        from sknano.io import DUMPWriter
        if 'structure' not in kwargs:
            kwargs['structure'] = self
        DUMPWriter.write(**kwargs)

    def write_xyz(self, **kwargs):
        """Write XYZ file."""
        from sknano.io import XYZWriter
        if 'structure' not in kwargs:
            kwargs['structure'] = self
        XYZWriter.write(**kwargs)


BaseGenerator = GeneratorBase


class GeneratorMixin:
    """Mixin class with concrete implementation of \
        :meth:`~GeneratorBase.generate` method."""
    def generate(self):
        """Concrete implementation of :meth:`~GeneratorBase.generate` \
            method."""
        self.structure.clear()
        for atom in self.crystal_cell:
            self.atoms.append(Atom(**atom.todict()))
        self.finalize()


class BulkGeneratorBase(GeneratorMixin, GeneratorBase):
    """Base class for the *bulk structure generator* classes."""
    def save(self, fname=None, scaling_matrix=None, **kwargs):
        """Save structure data."""
        if fname is None:
            fname = self.__class__.__name__[:-len('Generator')]
            if fname.endswith('CC'):
                fname = '_'.join((fname, '-'.join(set(self.basis.symbols))))
        if scaling_matrix is None:
            scaling_matrix = self.scaling_matrix.A
        elif isinstance(scaling_matrix, np.matrix):
            scaling_matrix = scaling_matrix.A

        if fname is not None and scaling_matrix is not None:
            ext = None
            if fname.endswith(supported_structure_formats):
                fname, ext = os.path.splitext(fname)

            if isinstance(scaling_matrix, np.ndarray):
                if scaling_matrix.ndim == 2:
                    if scaling_matrix.shape == (1, 3):
                        scaling_matrix = scaling_matrix.flatten()
                    elif scaling_matrix.shape == (3, 3) and \
                            np.all(scaling_matrix -
                                   np.diag(np.diag(scaling_matrix)) == 0):
                        scaling_matrix = scaling_matrix.diagonal()
                    else:
                        scaling_matrix = scaling_matrix.flatten()

            fname = '_'.join((fname, 'x'.join(map(str, scaling_matrix))))
            if ext is not None:
                fname = fname + ext
        super().save(fname=fname, **kwargs)

BaseBulkGenerator = BulkGeneratorBase
