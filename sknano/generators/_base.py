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

# from abc import ABCMeta, abstractmethod

import copy
import os

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.io import default_structure_format, supported_structure_formats

__all__ = ['Atom', 'Atoms', 'GeneratorBase', 'BulkGeneratorBase',
           'STRUCTURE_GENERATORS']


#: Tuple of structure generator classes.
STRUCTURE_GENERATORS = ('FullereneGenerator',
                        'GrapheneGenerator',
                        'PrimitiveCellGrapheneGenerator',
                        'ConventionalCellGrapheneGenerator',
                        'BilayerGrapheneGenerator',
                        'UnrolledSWNTGenerator',
                        'SWNTGenerator',
                        'SWNTBundleGenerator',
                        'MWNTGenerator',
                        'MWNTBundleGenerator',
                        'AlphaQuartzGenerator',
                        'DiamondStructureGenerator',
                        'FCCStructureGenerator',
                        'GoldGenerator',
                        'CopperGenerator',
                        'CaesiumChlorideStructureGenerator',
                        'RocksaltStructureGenerator',
                        'ZincblendeStructureGenerator',
                        'MoS2Generator')


class GeneratorBase:
    """Base structure generator class."""
    def __init__(self, *args, autogen=True, **kwargs):
        super().__init__(*args, **kwargs)

        if autogen:
            self.generate()

    def generate(self):
        """Common :meth:`~GeneratorBase.generate` method structure \
            generators."""
        self.structure_data.clear()
        for atom in self.crystal_cell:
            self.atoms.append(Atom(**atom.todict()))

    def save(self, fname=None, outpath=None, structure_format=None,
             deepcopy=True, center_CM=True, region_bounds=None,
             filter_condition=None, rotation_parameters=None, **kwargs):
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
        center_CM : bool, optional
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

        if center_CM:
            atoms.center_CM()

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
            fname=fname, outpath=outpath, structure=self, **kwargs)

        # StructureWriter.write(fname=fname, outpath=outpath,
        #                       structure_format=structure_format,
        #                       structure=self.structure)


class BulkGeneratorBase(GeneratorBase):
    """Base class for the *bulk structure generator* classes."""

    def save(self, fname=None, scaling_matrix=None, **kwargs):
        if fname is not None:
            if scaling_matrix is not None:
                fname = '_'.join((fname,
                                  'x'.join(map(str, scaling_matrix))))

        super().save(fname=fname, **kwargs)
