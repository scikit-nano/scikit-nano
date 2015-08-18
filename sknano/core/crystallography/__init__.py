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

   supercell_lattice_points
   pymatgen_structure

"""
from __future__ import absolute_import, unicode_literals
__docformat__ = 'restructuredtext en'

from ._extras import *

from ._xtal_cells import *
from ._xtal_lattices import *
from ._xtal_structures import *

__all__ = [s for s in dir() if not s.startswith('_')]
