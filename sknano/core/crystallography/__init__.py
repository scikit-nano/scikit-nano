# -*- coding: utf-8 -*-
"""
======================================================================
Crystallography code (:mod:`sknano.core.crystallography`)
======================================================================

.. currentmodule:: sknano.core.crystallography

Contents
========

Base crystallography classes
------------------------------

.. autosummary::
   :toctree: generated/

   LatticeBase
   ReciprocalLatticeBase
   Domain
   CrystalCell
   UnitCell
   SuperCell

2D crystal lattices
---------------------------

.. autosummary::
   :toctree: generated/

   Direct2DLatticeMixin
   Reciprocal2DLatticeMixin
   Crystal2DLattice
   Reciprocal2DLattice

3D crystal lattices
---------------------------

.. autosummary::
   :toctree: generated/

   Direct3DLatticeMixin
   Reciprocal3DLatticeMixin
   Crystal3DLattice
   Reciprocal3DLattice

Helper functions
----------------

.. autosummary::
   :toctree: generated/

   pbc_diff
   supercell_lattice_points

"""
from __future__ import absolute_import, unicode_literals
__docformat__ = 'restructuredtext en'

from ._extras import *
from ._xtal_cells import *
from ._xtal_lattices import *

__all__ = [s for s in dir() if not s.startswith('_')]
