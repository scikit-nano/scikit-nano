# -*- coding: utf-8 -*-
"""
======================================================================
Structure generators (:mod:`sknano.generators`)
======================================================================

.. currentmodule:: sknano.generators

Contents
========

Nanostructure generators
------------------------

.. autosummary::
   :toctree: generated/

   GeneratorBase
   FullereneGenerator
   GrapheneGenerator
   PrimitiveCellGrapheneGenerator
   ConventionalCellGrapheneGenerator
   BilayerGrapheneGenerator
   MWNTGenerator
   MWNTBundleGenerator
   SWNTGenerator
   SWNTBundleGenerator
   UnrolledSWNTGenerator

Bulk structure generators
--------------------------

.. autosummary::
   :toctree: generated/

   BulkGeneratorBase
   AlphaQuartzGenerator
   GoldGenerator
   CopperGenerator
   MoS2Generator
   CaesiumChlorideStructureGenerator
   DiamondStructureGenerator
   BCCStructureGenerator
   FCCStructureGenerator
   RocksaltStructureGenerator
   ZincblendeStructureGenerator

Other
-----

.. autodata:: STRUCTURE_GENERATORS
   :annotation: = tuple of recognized generator classes.

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._base import *
from ._bulk_structure_generator import *
from ._mixins import *
from ._fullerene_generator import *
from ._graphene_generator import *
from ._bilayer_graphene_generator import *
from ._swnt_generator import *
from ._mwnt_generator import *
from ._nanotube_bundle_generator import *
from ._swnt_bundle_generator import *
from ._mwnt_bundle_generator import *
from ._unrolled_swnt_generator import *
# from ._defect_generators import *

__all__ = [s for s in dir() if not s.startswith('_')]
