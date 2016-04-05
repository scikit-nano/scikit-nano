# -*- coding: utf-8 -*-
"""
======================================================================
Structure generators (:mod:`sknano.generators`)
======================================================================

.. currentmodule:: sknano.generators

Contents
========

Base/mixin generator classes
-----------------------------

.. autosummary::
   :toctree: generated/

   GeneratorBase
   GeneratorMixin
   CrystalStructureGenerator
   NanoStructureGenerator
   GrapheneGeneratorBase
   NanotubeBundleGeneratorBase
   SWNTGeneratorBase
   MWNTGeneratorBase


Crystal structure generators
-----------------------------

.. autosummary::
   :toctree: generated/

   AlphaQuartzGenerator
   DiamondGenerator
   CaesiumChlorideGenerator
   RocksaltGenerator
   ZincblendeGenerator
   BCCGenerator
   FCCGenerator
   MoS2Generator


NanoStructure generator classes
--------------------------------

.. autosummary::
   :toctree: generated/

   FullereneGenerator
   GrapheneGenerator
   PrimitiveCellGrapheneGenerator
   ConventionalCellGrapheneGenerator
   BilayerGrapheneGenerator
   MWNTGenerator
   SWNTGenerator
   UnrolledSWNTGenerator

Compound structure generators
------------------------------

.. autosummary::
   :toctree: generated/

   LayeredStructureGenerator

Other
-----

.. autodata:: STRUCTURE_GENERATORS
   :annotation: = tuple of recognized generator classes.

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._base import *
from ._generator_configparser import *
from ._xtal_structure_generator import *
from ._mixins import *
from ._fullerene_generator import *
from ._graphene_generator import *
from ._bilayer_graphene_generator import *
from ._nanotube_bundle_generator import *
from ._swnt_generator import *
from ._mwnt_generator import *
from ._unrolled_swnt_generator import *
from ._layered_structure_generator import *
# from ._defect_generators import *

__all__ = [s for s in dir() if not s.startswith('_')]
