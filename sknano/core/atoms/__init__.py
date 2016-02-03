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
   :toctree: generated/

   Atom
   Atoms

`Atom`/`Atoms` sub-classes
---------------------------
.. autosummary::
   :toctree: generated/

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
   :toctree: generated/

   BasisAtom
   BasisAtoms
   MDAtom
   MDAtoms
   StructureAtom
   StructureAtoms

`Bond`/`Bonds` classes
----------------------
.. autosummary::
   :toctree: generated/

   Bond
   Bonds

Classes for molecular dynamics simulations
------------------------------------------
.. autosummary::
   :toctree: generated/

   Trajectory
   Snapshot

Helper functions for atom objects
---------------------------------

.. autosummary::
   :toctree: generated/

   vdw_radius_from_basis

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._atoms import *
from ._xyz_atoms import *
from ._velocity_atoms import *
from ._force_atoms import *
from ._cn_atoms import *
from ._image_atoms import *
from ._charged_atoms import *
from ._energy_atoms import *
from ._id_atoms import *
from ._type_atoms import *
from ._lattice_atoms import *
from ._dipole_atoms import *
from ._basis_atoms import *
from ._bonds import *
from ._md_atoms import *
from ._trajectory import *
from ._structure_atoms import *
from ._neighbor_atoms import *
from ._vdW_atoms import *
from ._selections import *
from ._topology import *

__all__ = [s for s in dir() if not s.startswith('_')]
