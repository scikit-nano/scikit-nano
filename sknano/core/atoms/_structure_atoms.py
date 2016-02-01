# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for structure analysis (:mod:`sknano.core.atoms._structure_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._structure_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

# from operator import attrgetter

from ._cn_atoms import CNAtom, CNAtoms
from ._id_atoms import IDAtom, IDAtoms
from ._charged_atoms import ChargedAtom, ChargedAtoms
from ._image_atoms import ImageAtom, ImageAtoms
from ._type_atoms import TypeAtom, TypeAtoms
from ._neighbor_atoms import NeighborAtom, NeighborAtoms
from ._lattice_atoms import LatticeAtom, LatticeAtoms
from ._xyz_atoms import XYZAtom, XYZAtoms
from ._velocity_atoms import VelocityAtom, VelocityAtoms
from ._topology import AtomTopologyMixin, AtomsTopologyMixin
from .mixins import AtomAdapterMixin, AtomsAdapterMixin, \
    POAVAtomMixin, POAVAtomsMixin, RingAtomMixin, RingAtomsMixin

__all__ = ['StructureAtom', 'StructureAtoms']


class StructureAtom(AtomAdapterMixin, POAVAtomMixin, RingAtomMixin,
                    AtomTopologyMixin, NeighborAtom,
                    CNAtom, VelocityAtom, ImageAtom, LatticeAtom,
                    XYZAtom, ChargedAtom, TypeAtom, IDAtom):
    """An `Atom` class for structure analysis."""
    pass


class StructureAtoms(AtomsAdapterMixin, POAVAtomsMixin, RingAtomsMixin,
                     AtomsTopologyMixin, NeighborAtoms,
                     CNAtoms, VelocityAtoms, ImageAtoms,
                     LatticeAtoms, XYZAtoms, ChargedAtoms, TypeAtoms, IDAtoms):
    """An `Atoms` sub-class for structure analysis."""

    @property
    def __atom_class__(self):
        return StructureAtom

    # def sort(self, key=attrgetter('CN', 'v', 'i', 'r', 'q', 'type', 'mol',
    #                               'id', 'mass', 'Z', 'element'),
    #          reverse=False):
    #     super().sort(key=key, reverse=reverse)
