# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of nature's building blocks (:mod:`sknano.core.atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms

Contents
========

Sub-package
--------------

* Mixin `Atom`/`Atoms` classes (:mod:`sknano.core.atoms.mixins`)

Base `Atom`/`Atoms` classes
---------------------------
The `Atom` class represents a single atom. The `Atoms` class is a container
class for `Atom` class instances. Sub-classes of `Atom` classes
add new atom attributes to the `Atom` class. Every `Atom` sub-class
has a corresponding container class that sub-classes the `Atoms` class.

.. autosummary::
   :toctree: api/

   Atom
   Atoms

`Atom`/`Atoms` sub-classes
---------------------------
.. autosummary::
   :toctree: api/

   ChargedAtom
   ChargedAtoms
   CNAtom
   CNAtoms
   DipoleAtom
   DipoleAtoms
   EnergyAtom
   EnergyAtoms
   ForceAtom
   ForceAtoms
   IDAtom
   IDAtoms
   ImageAtom
   ImageAtoms
   LatticeAtom
   LatticeAtoms
   NeighborAtom
   NeighborAtoms
   TypeAtom
   TypeAtoms
   VanDerWaalsAtom
   VanDerWaalsAtoms
   VelocityAtom
   VelocityAtoms
   XYZAtom
   XYZAtoms

Composite `Atom`/`Atoms` classes
--------------------------------
.. autosummary::
   :toctree: api/

   BasisAtom
   BasisAtoms
   MDAtom
   MDAtoms
   StructureAtom
   StructureAtoms

Classes for molecular dynamics simulations
------------------------------------------
.. autosummary::
   :toctree: api/

   Trajectory
   Snapshot

Helper functions for atom objects
---------------------------------

.. autosummary::
   :toctree: api/

   compute_angle
   compute_bond
   compute_dihedral
   compute_improper
   vdw_radius_from_basis

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from .atoms import *
from .xyz_atoms import *
from .velocity_atoms import *
from .force_atoms import *
from .cn_atoms import *
from .image_atoms import *
from .charged_atoms import *
from .energy_atoms import *
from .id_atoms import *
from .type_atoms import *
from .lattice_atoms import *
from .dipole_atoms import *
from .basis_atoms import *
from .md_atoms import *
from .trajectory import *
from .structure_atoms import *
from .neighbor_atoms import *
from .vdW_atoms import *
from .selections import *

from .mixins import *

__all__ = [s for s in dir() if not s.startswith('_')]
