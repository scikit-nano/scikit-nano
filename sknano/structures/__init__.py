# -*- coding: utf-8 -*-
"""
===============================================================================
Structure classes (:mod:`sknano.structures`)
===============================================================================

.. currentmodule:: sknano.structures

Contents
========

This package defines classes for creating
abstract object representations of nanostructures
including fullerenes, graphene, and nanotubes.

See the specific class doc pages below for detailed documentation on its use.


Nanostructure Classes
----------------------

.. autosummary::
   :toctree: generated/

   Fullerene
   Graphene
   PrimitiveCellGraphene
   HexagonalGraphene
   ConventionalCellGraphene
   RectangularGraphene
   BilayerGraphene
   SWNT
   SWNTBundle
   MWNT
   MWNTBundle

Base/Mixin Classes
------------------

.. autosummary::
   :toctree: generated/

   NanoStructureBase
   GrapheneBase
   NanotubeBundleBase
   GrapheneMixin
   SWNTMixin
   MWNTMixin
   NanotubeBundleMixin
   UnrolledSWNTMixin

Compute Functions
------------------

Nanotube compute functions:

.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

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

from ._base import *
from ._nanotube_bundle import *

from ._fullerenes import *

from ._graphene import *
from ._bilayer_graphene import *

from ._swnt import *
from ._swnt_bundle import *
from ._mwnt import *
from ._mwnt_bundle import *
from ._unrolled_swnt import *

from ._extras import *
from ._defects import *

__all__ = [s for s in dir() if not s.startswith('_')]
