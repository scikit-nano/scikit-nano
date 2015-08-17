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
    Iron, Gold, Copper, BCCStructure, FCCStructure, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2
from ._base import BulkGeneratorBase

__all__ = ['DiamondStructureGenerator',
           'BCCStructureGenerator', 'FCCStructureGenerator',
           'CaesiumChlorideStructureGenerator',
           'RocksaltStructureGenerator', 'ZincblendeStructureGenerator',
           'AlphaQuartzGenerator', 'IronGenerator',
           'GoldGenerator', 'CopperGenerator', 'MoS2Generator']


class AlphaQuartzGenerator(BulkGeneratorBase, AlphaQuartz):
    """:class:`AlphaQuartz` generator class.

    Parameters
    ----------
    a, c : :class:`~python:float`
    scaling_matrix : {None, :class:`~python:float`, :class:`~python:list`}
    """
    def save(self, fname='alpha_quartz', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class DiamondStructureGenerator(BulkGeneratorBase, DiamondStructure):
    """:class:`DiamondStructure` generator class.

    Parameters
    ----------
    a : class:`~python:float`

    """
    def save(self, fname='diamond', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class CaesiumChlorideStructureGenerator(BulkGeneratorBase,
                                        CaesiumChlorideStructure):
    """:class:`CaesiumChlorideStructure` generator class."""
    def save(self, fname='caesium_chloride', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class RocksaltStructureGenerator(BulkGeneratorBase, RocksaltStructure):
    """:class:`RocksaltStructure` generator class."""
    def save(self, fname='rock_salt', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class ZincblendeStructureGenerator(BulkGeneratorBase, ZincblendeStructure):
    """:class:`ZincblendeStructure` generator class."""
    def save(self, fname='zincblende', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class BCCStructureGenerator(BulkGeneratorBase, BCCStructure):
    """:class:`BCCStructure` generator class."""
    def save(self, fname='bcc_structure', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class FCCStructureGenerator(BulkGeneratorBase, FCCStructure):
    """:class:`FCCStructure` generator class."""
    def save(self, fname='fcc_structure', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class IronGenerator(BulkGeneratorBase, Iron):
    """:class:`Iron` generator class."""
    def save(self, fname='iron', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class GoldGenerator(BulkGeneratorBase, Gold):
    """:class:`Gold` generator class."""
    def save(self, fname='gold', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class CopperGenerator(BulkGeneratorBase, Copper):
    """:class:`Copper` generator class."""
    def save(self, fname='copper', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)


class MoS2Generator(BulkGeneratorBase, MoS2):
    """:class:`MoS2` generator class."""
    def save(self, fname='MoS2', **kwargs):
        super().save(fname=fname, scaling_matrix=self.scaling_matrix, **kwargs)
