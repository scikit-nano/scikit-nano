# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for LAMMPS structure data (:mod:`sknano.core.atoms._lammps_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._lammps_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

# from operator import attrgetter

from ._cn_atoms import CNAtom, CNAtoms
from ._id_atoms import IDAtom, IDAtoms
from ._xyz_atoms import XYZAtom, XYZAtoms
from ._charged_atoms import ChargedAtom, ChargedAtoms
from ._velocity_atoms import VelocityAtom, VelocityAtoms
from ._image_atoms import ImageAtom, ImageAtoms
from ._type_atoms import TypeAtom, TypeAtoms
from ._kdtree_atoms import KDTreeAtomMixin, KDTreeAtomsMixin
from ._poav_atoms import POAVAtomMixin, POAVAtomsMixin
from ._neighbor_atoms import NeighborAtomMixin, NeighborAtomsMixin

# from ._bonds import Bonds

__all__ = ['LAMMPSAtom', 'LAMMPSAtoms']


class LAMMPSAtom(NeighborAtomMixin, POAVAtomMixin, KDTreeAtomMixin, CNAtom,
                 VelocityAtom, ImageAtom, XYZAtom, ChargedAtom, TypeAtom,
                 IDAtom):
    """An Atom class for LAMMPS structure data."""
    pass


class LAMMPSAtoms(NeighborAtomsMixin, POAVAtomsMixin, KDTreeAtomsMixin,
                  CNAtoms, VelocityAtoms, ImageAtoms, XYZAtoms,
                  ChargedAtoms, TypeAtoms, IDAtoms):
    """An Atoms class for LAMMPS structure data."""
    @property
    def __atom_class__(self):
        return LAMMPSAtom
