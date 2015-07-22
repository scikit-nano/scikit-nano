# -*- coding: utf-8 -*-
"""
======================================================================
Crystallography code (:mod:`sknano.core.crystallography`)
======================================================================

.. currentmodule:: sknano.core.crystallography

Contents
========

2D lattices and structures
---------------------------
.. autosummary::
   :toctree: generated/

   Direct2DLatticeMixin
   Reciprocal2DLatticeMixin
   Crystal2DLattice
   Reciprocal2DLattice
   Crystal2DStructure

3D lattices and structures
---------------------------
.. autosummary::
   :toctree: generated/

   Direct3DLatticeMixin
   Reciprocal3DLatticeMixin
   UnitCellMixin
   Crystal3DLattice
   Reciprocal3DLattice
   Crystal3DStructure

   CaesiumChlorideStructure
   CsClStructure
   DiamondStructure
   RocksaltStructure
   NaClStructure
   SphaleriteStructure
   ZincblendeStructure
   BCCStructure
   FCCStructure
   HexagonalClosePackedStructure
   CubicClosePackedStructure

   AlphaQuartz
   Iron
   Copper
   Gold
   MoS2

Helper functions
----------------
.. autosummary::
   :toctree: generated/

   pymatgen_structure

"""
from __future__ import absolute_import, unicode_literals
__docformat__ = 'restructuredtext en'

from ._base import *
from ._extras import *
from ._mixins import *
from ._unit_cell import *

from ._2D_lattices import *
from ._2D_structures import *

from ._3D_lattices import *
from ._3D_structures import *

__all__ = [s for s in dir() if not s.startswith('_')]
