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

from sknano.core.atoms import StructureAtom as Atom, StructureAtoms as Atoms
from sknano.io import StructureData, StructureWriter, \
    default_structure_format, supported_structure_formats

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


class GeneratorBase(metaclass=ABCMeta):
    """Base class for generator classes"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.structure_data = StructureData()

    def __getattr__(self, name):
        if name != '_structure_data':
            return getattr(self._structure_data, name)

    @property
    @abstractmethod
    def generate(self):
        """Generate the structure coordinates."""
        return NotImplementedError

    @property
    def structure_data(self):
        """Return :class:`~sknano.io.StructureData` instance."""
        return self._structure_data

    @structure_data.setter
    def structure_data(self, value):
        """Set :class:`~sknano.io.StructureData` instance."""
        if not isinstance(value, StructureData):
            raise TypeError('Expected a `StructureData` object.')
        self._structure_data = value

    @structure_data.deleter
    def structure_data(self):
        del self._structure_data

    @property
    def structure(self):
        """Alias to :attr:`~GeneratorBase.structure_data`."""
        return self.structure_data

    @structure.setter
    def structure(self, value):
        """Alias to :attr:`~GeneratorBase.structure_data`."""
        self.structure_data = value

    @structure.deleter
    def structure(self):
        del self.structure_data

    def save(self, fname=None, outpath=None, structure_format=None,
             deepcopy=True, center_CM=True, region_bounds=None,
             filter_condition=None, **kwargs):
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

        if kwargs:
            atoms.rotate(**kwargs)

        StructureWriter.write(fname=fname, outpath=outpath,
                              structure_format=structure_format, atoms=atoms)


class BulkGeneratorBase(GeneratorBase):
    """Base class for the *bulk structure generator* classes."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.generate()

    def generate(self):
        """Common :meth:`~GeneratorBase.generate` method for buk structure \
            generators."""
        self.structure_data.clear()
        for atom in self.unit_cell:
            self.atoms.append(Atom(atom.element, x=atom.x, y=atom.y, z=atom.z))

    def save(self, fname=None, scaling_matrix=None, **kwargs):
        if fname is not None:
            if scaling_matrix is not None:
                fname = '_'.join((fname,
                                  'x'.join(map(str, scaling_matrix))))

        super().save(fname=fname, **kwargs)
