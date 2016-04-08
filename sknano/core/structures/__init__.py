# -*- coding: utf-8 -*-
"""
===============================================================================
Core structure classes (:mod:`sknano.core.structures`)
===============================================================================

.. currentmodule:: sknano.core.structures

Contents
========

This package defines classes for creating
abstract object representations of nanostructures
including fullerenes, graphene, and nanotubes.

See the specific class doc pages below for detailed documentation on its use.

Base/Mixin structure classes
----------------------------

.. autosummary::
   :toctree: api/

   StructureBase
   StructureMixin
   CrystalStructureBase
   NanoStructureBase
   GrapheneBase
   GrapheneMixin
   SWNTBase
   SWNTMixin
   MWNTBase
   MWNTMixin
   NanotubeBundleBase
   NanotubeBundleMixin
   UnrolledSWNTBase
   UnrolledSWNTMixin

Crystal structure classes
--------------------------

.. autosummary::
   :toctree: api/

   Crystal2DStructure
   Crystal3DStructure
   CubicStructure
   BCCStructure
   FCCStructure
   CubicClosePackedStructure
   HexagonalClosePackedStructure
   HexagonalStructure
   CaesiumChlorideStructure
   CsClStructure
   DiamondStructure
   RocksaltStructure
   NaClStructure
   ZincblendeStructure
   SphaleriteStructure
   AlphaQuartz
   MoS2

Nanostructure Classes
----------------------

.. autosummary::
   :toctree: api/

   Fullerene
   Graphene
   PrimitiveCellGraphene
   HexagonalGraphene
   ConventionalCellGraphene
   RectangularGraphene
   BilayerGraphene
   NanotubeUnitCell
   SWNT
   MWNT
   UnrolledSWNT

Composite structure classes
----------------------------

.. autosummary::
   :toctree: api/

   Composition
   Compositions

Helper functions for nanotube properties
------------------------------------------

.. autosummary::
   :toctree: api/

   pymatgen_structure

Helper functions for nanotube properties
------------------------------------------

Nanotube compute functions:

.. autosummary::
   :toctree: api/

   compute_d
   compute_dR
   compute_N
   compute_t1
   compute_t2
   compute_Ch
   compute_chiral_angle
   compute_T
   compute_dt
   compute_rt
   compute_M
   compute_R
   compute_R_chiral_angle
   compute_symmetry_operation
   compute_psi
   compute_tau
   compute_Lz
   compute_electronic_type
   compute_Natoms_per_unit_cell
   compute_Natoms_per_tube
   compute_Natoms
   compute_unit_cell_mass
   compute_linear_mass_density
   compute_tube_mass
   compute_bundle_density

Helper functions for working with :math:`(n, m)` chirality data
-----------------------------------------------------------------

.. autosummary::
   :toctree: api/

   cmp_Ch
   filter_Ch
   filter_Ch_list
   generate_Ch_list
   generate_Ch_property_grid
   get_Ch_indices
   get_Ch_type
   map_Ch

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from .base import *
from .xtal_structures import *
from .compositions import *
from .nanotube_bundle import *
from .fullerenes import *
from .graphene import *
from .bilayer_graphene import *
from .swnt import *
from .mwnt import *
from .unrolled_swnt import *
from .extras import *
from .defects import *

__all__ = [s for s in dir() if not s.startswith('_')]
