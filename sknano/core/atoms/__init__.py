# -*- coding: utf-8 -*-
"""
===============================================================================
Class representations of nature's building blocks (:mod:`sknano.core.atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms

Contents
========

The `Atom` class represents a single atom. The `Atoms` class is a container
class for `Atom` class instances. Every new atom class has a corresponding
container class.

Base classes:

.. autosummary::
   :toctree: generated/

   Atom
   Atoms

`Atom`/`Atoms` sub-classes:

.. autosummary::
   :toctree: generated/

   XYZAtom
   XYZAtoms


Atoms with an eXtended set of attributes:

.. autosummary::
   :toctree: generated/

   XAtom
   XAtoms

Atoms for nearest-neighbor structure analysis:

.. autosummary::
   :toctree: generated/

   KDTAtom
   KDTAtoms

Atoms and data structures for POAV structure analysis:

.. autosummary::
   :toctree: generated/

   POAVAtom
   POAVAtoms
   POAVAtomMixin
   POAV
   POAV1
   POAV2
   POAVR

Class representation of atom bonds.

.. autosummary::
   :toctree: generated/

   Bond
   Bonds

Classes for molecular dynamics simulations:

.. autosummary::
   :toctree: generated/

   Trajectory
   Snapshot
   MDAtom
   MDAtoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

__all__ = ['StructureAtom', 'StructureAtoms']

from ._atom import *
from ._atoms import *
from ._xyz_atom import *
from ._xyz_atoms import *
from ._velocity_atom import *
from ._velocity_atoms import *
from ._force_atom import *
from ._force_atoms import *
from ._cn_atom import *
from ._cn_atoms import *
from ._image_atom import *
from ._image_atoms import *
from ._charged_atom import *
from ._charged_atoms import *
from ._energy_atom import *
from ._energy_atoms import *
from ._id_atom import *
from ._id_atoms import *
from ._type_atom import *
from ._type_atoms import *

from ._basis_atom import *
from ._basis_atoms import *

from ._extended_atom import *
from ._extended_atoms import *
from ._kdtree_atom import *
from ._kdtree_atoms import *
from ._poav_atom import *
from ._poav_atoms import *

from ._bond import *
from ._bonds import *

from ._md_atom import *
from ._md_atoms import *
from ._trajectory import *

from ._structure_atom import *
from ._structure_atoms import *

from ._neighbor_atom import *
from ._neighbor_atoms import *

__all__ += [s for s in dir() if not s.startswith('_')]
