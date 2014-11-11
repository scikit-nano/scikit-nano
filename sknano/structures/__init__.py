# -*- coding: utf-8 -*-
"""
===============================================================================
Abstract nanostructure data structures (:mod:`sknano.structures`)
===============================================================================

.. currentmodule:: sknano.structures

This package defines classes for creating
abstract object representations of nanostructures
including fullerenes, graphene, and nanotubes.

See the specific class doc pages below for detailed documentation on its use.

Contents
========

Classes for creating abstract object representations of nanostructures
--------------------------------------------------------------------------
.. autosummary::
   :toctree: generated/

   StructureBase
   NanotubeBundleBase

   MWNTMixin
   NanotubeBundleMixin
   UnrolledSWNTMixin

   Fullerene

   GraphenePrimitiveCell
   Graphene
   BilayerGraphene

   SWNT
   SWNTBundle
   MWNT
   MWNTBundle

Functions for computing nanostructure properties
--------------------------------------------------------------------------
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

Helper functions and data structures related to :math:`(n, m)` chirality data
-------------------------------------------------------------------------------
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
   chiral_type_name_mappings
   Ch_types
   filter_key_type_mappings

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from ._base import *
from ._mixins import *
from ._nanotube_bundle import *

from ._fullerenes import *

from ._bilayer_graphene import *
from ._graphene import *

from ._swnt import *
from ._swnt_bundle import *
from ._mwnt import *
from ._mwnt_bundle import *
from ._unrolled_swnt import *

from ._compute_funcs import *
from ._extras import *

__all__ = [s for s in dir() if not s.startswith('_')]
