# -*- coding: utf-8 -*-
"""
===============================================================================
Bulk structure generator (:mod:`sknano.generators._bulk_structure_generator`)
===============================================================================

.. currentmodule:: sknano.generators._bulk_structure_generator

.. todo::

   Add methods to perform fractional translation and cartesian translation
   before structure generation.

.. todo::

   Handle different units in output coordinates.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

__docformat__ = 'restructuredtext en'

# import numpy as np

# from sknano.core import pluralize
# from sknano.core.math import Point, Vector
from sknano.core.crystallography import AlphaQuartz, DiamondStructure, \
    Gold, Copper, FCCStructure, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure
from sknano.core.refdata import crystal_bonds as bonds
from ._base import BulkGeneratorBase

__all__ = ['AlphaQuartzGenerator', 'DiamondStructureGenerator',
           'FCCStructureGenerator', 'GoldGenerator', 'CopperGenerator',
           'CaesiumChlorideStructureGenerator',
           'RocksaltStructureGenerator', 'ZincblendeStructureGenerator']


class AlphaQuartzGenerator(AlphaQuartz, BulkGeneratorBase):
    """Class for generating bulk crystal structures.

    Parameters
    ----------
    """
    def __init__(self, a=bonds['alpha_quartz']['a'],
                 c=bonds['alpha_quartz']['c'],
                 scaling_matrix=None, autogen=True):
        super().__init__(a=a, c=c, scaling_matrix=scaling_matrix)
        self.generate()

    def save(self, fname='alpha_quartz', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class DiamondStructureGenerator(DiamondStructure, BulkGeneratorBase):
    def __init__(self, a=bonds['diamond'], scaling_matrix=None):
        super().__init__(a=a, scaling_matrix=scaling_matrix)
        self.generate()

    def save(self, fname='diamond', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class FCCStructureGenerator(FCCStructure, BulkGeneratorBase):
    def __init__(self, a, basis, coords, scaling_matrix=None):
        pass


class GoldGenerator(Gold, BulkGeneratorBase):

    def save(self, fname='gold', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class CopperGenerator(Copper, BulkGeneratorBase):
    # self.fname = 'copper'

    def save(self, fname='copper', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class CaesiumChlorideStructureGenerator(CaesiumChlorideStructure,
                                        BulkGeneratorBase):
    def __init__(self, a=bonds['caesium_chloride'], scaling_matrix=None):
        super().__init__(a=a, scaling_matrix=scaling_matrix)
        self.generate()

    def save(self, fname='caesium_chloride', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class RocksaltStructureGenerator(RocksaltStructure, BulkGeneratorBase):
    def __init__(self, a=bonds['rock_salt'], scaling_matrix=None):
        super().__init__(a=a, scaling_matrix=scaling_matrix)
        self.generate()

    def save(self, fname='rock_salt', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class ZincblendeStructureGenerator(ZincblendeStructure, BulkGeneratorBase):
    def __init__(self, a=bonds['zincblende'], scaling_matrix=None):
        super().__init__(a=a, scaling_matrix=scaling_matrix)
        self.generate()

    def save(self, fname='zincblende', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)
