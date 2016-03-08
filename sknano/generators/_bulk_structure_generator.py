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
from sknano.core.structures import AlphaQuartz, DiamondStructure, \
    Iron, Gold, Copper, BCCStructure, FCCStructure, CaesiumChlorideStructure, \
    RocksaltStructure, ZincblendeStructure, MoS2
from ._base import BulkGeneratorBase

__all__ = ['DiamondGenerator', 'DiamondStructureGenerator',
           'BCCGenerator', 'BCCStructureGenerator',
           'FCCGenerator', 'FCCStructureGenerator',
           'CsClGenerator', 'CaesiumChlorideGenerator',
           'CaesiumChlorideStructureGenerator', 'NaClGenerator',
           'RocksaltGenerator', 'RocksaltStructureGenerator',
           'ZincblendeGenerator', 'ZincblendeStructureGenerator',
           'AlphaQuartzGenerator', 'IronGenerator',
           'GoldGenerator', 'CopperGenerator', 'MoS2Generator']


class AlphaQuartzGenerator(BulkGeneratorBase, AlphaQuartz):
    """:class:`~sknano.core.structures.AlphaQuartz` generator class.

    Parameters
    ----------
    a, c : :class:`~python:float`
    scaling_matrix : {None, :class:`~python:float`, :class:`~python:list`}

    Examples
    --------
    >>> from sknano.generators import AlphaQuartzGenerator
    >>> alpha_quartz = AlphaQuartzGenerator(scaling_matrix=5)
    >>> alpha_quartz.save()

    .. image:: /images/alpha_quartz_5x5x5-1.png

    """
    pass


class DiamondGenerator(BulkGeneratorBase, DiamondStructure):
    """:class:`~sknano.core.structures.DiamondStructure` generator class.

    Parameters
    ----------
    a : class:`~python:float`

    Examples
    --------
    >>> from sknano.generators import DiamondGenerator
    >>> diamond = DiamondGenerator(scaling_matrix=5)
    >>> diamond.save()

    .. image:: /images/diamond_5x5x5-1.png

    """
    pass

DiamondStructureGenerator = DiamondGenerator


class CaesiumChlorideGenerator(BulkGeneratorBase, CaesiumChlorideStructure):
    """:class:`~sknano.core.structures.CaesiumChlorideStructure` generator \
        class.

    Examples
    --------
    >>> from sknano.generators import CaesiumChlorideGenerator
    >>> caesium_chloride = CaesiumChlorideGenerator(scaling_matrix=10)
    >>> caesium_chloride.save()

    .. image:: /images/caesium_chloride_10x10x10-1.png

    """
    pass

CaesiumChlorideStructureGenerator = CsClGenerator = CaesiumChlorideGenerator


class RocksaltGenerator(BulkGeneratorBase, RocksaltStructure):
    """:class:`~sknano.core.structures.RocksaltStructure` generator class.

    Examples
    --------
    >>> from sknano.generators import RocksaltGenerator
    >>> nacl = RocksaltGenerator(scaling_matrix=5)
    >>> nacl.save()

    .. image:: /images/rock_salt_5x5x5-1.png

    """
    pass

RocksaltStructureGenerator = NaClGenerator = RocksaltGenerator


class ZincblendeGenerator(BulkGeneratorBase, ZincblendeStructure):
    """:class:`~sknano.core.structures.ZincblendeStructure` generator class.

    Examples
    --------
    >>> from sknano.generators import ZincblendeGenerator
    >>> zincblende = ZincblendeGenerator(scaling_matrix=5)
    >>> zincblende.save()

    .. image:: /images/zincblende_5x5x5-1.png

    """
    pass

ZincblendeStructureGenerator = ZincblendeGenerator


class BCCGenerator(BulkGeneratorBase, BCCStructure):
    """:class:`~sknano.core.structures.BCCStructure` generator class.

    Examples
    --------
    >>> from sknano.generators import BCCGenerator
    >>> iron = BCCGenerator(basis='Fe', scaling_matrix=5)
    >>> iron.save()

    .. image:: /images/BCC_Fe_10x10x10-1.png

    """
    pass

BCCStructureGenerator = BCCGenerator


class FCCGenerator(BulkGeneratorBase, FCCStructure):
    """:class:`~sknano.core.structures.FCCStructure` generator class.

    Examples
    --------
    >>> from sknano.generators import FCCGenerator
    >>> gold = FCCGenerator(basis='Au', scaling_matrix=10)
    >>> gold.save()

    .. image:: /images/FCC_Au_5x5x5-1.png

    """
    pass

FCCStructureGenerator = FCCGenerator


class MoS2Generator(BulkGeneratorBase, MoS2):
    """:class:`~sknano.core.structures.MoS2` generator class.

    Examples
    --------
    >>> from sknano.generators import MoS2Generator
    >>> mos2 = MoS2Generator(scaling_matrix=5)
    >>> mos2.save()

    .. image:: /images/MoS2_5x5x5-1.png

    """
    pass


class IronGenerator(BulkGeneratorBase, Iron):
    """:class:`~sknano.core.structures.Iron` generator class."""
    pass


class GoldGenerator(BulkGeneratorBase, Gold):
    """:class:`~sknano.core.structures.Gold` generator class."""
    pass


class CopperGenerator(BulkGeneratorBase, Copper):
    """:class:`~sknano.core.structures.Copper` generator class."""
    pass
