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

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

#from ._base import *
from ._atom import *
from ._atoms import *
from ._extended_atom import *
from ._extended_atoms import *
from ._kdtree_atom import *
from ._kdtree_atoms import *
from ._poav_atom import *
from ._poav_atoms import *
from ._structure_atoms import *

from ._bond import *
from ._bonds import *

from ._trajectory import *

__all__ = [s for s in dir() if not s.startswith('_')]
