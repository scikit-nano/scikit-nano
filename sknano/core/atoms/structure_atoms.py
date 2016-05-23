# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for structure analysis (:mod:`sknano.core.atoms.structure_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.structure_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals

__docformat__ = 'restructuredtext en'

# from operator import attrgetter

# from .cn_atoms import CNAtom, CNAtoms
from .id_atoms import IDAtom, IDAtoms
from .charged_atoms import ChargedAtom, ChargedAtoms
from .image_atoms import ImageAtom, ImageAtoms
from .type_atoms import TypeAtom, TypeAtoms
from .neighbor_atoms import NeighborAtom, NeighborAtoms
from .lattice_atoms import LatticeAtom, LatticeAtoms
from .xyz_atoms import XYZAtom, XYZAtoms
# from .velocity_atoms import VelocityAtom, VelocityAtoms
from .mixins import AtomAdapterMixin, AtomsAdapterMixin, \
    AtomTopologyMixin, AtomsTopologyMixin, AtomTransformationsMixin, \
    AtomsTransformationsMixin, BoundingRegionsMixin, \
    POAVAtomMixin, POAVAtomsMixin, RingAtomMixin, RingAtomsMixin


__all__ = ['StructureAtom', 'StructureAtoms']


class StructureAtom(AtomAdapterMixin, AtomTopologyMixin,
                    AtomTransformationsMixin, POAVAtomMixin, RingAtomMixin,
                    NeighborAtom, ImageAtom, LatticeAtom, XYZAtom,
                    ChargedAtom, TypeAtom, IDAtom):
    """An :class:`Atom` class for structure analysis."""
    @property
    def __atoms_class__(self):
        return StructureAtoms


class StructureAtoms(AtomsAdapterMixin, AtomsTopologyMixin,
                     AtomsTransformationsMixin, BoundingRegionsMixin,
                     POAVAtomsMixin, RingAtomsMixin, NeighborAtoms, ImageAtoms,
                     LatticeAtoms, XYZAtoms, ChargedAtoms, TypeAtoms, IDAtoms):
    """An :class:`Atoms` sub-class for structure analysis."""

    @property
    def __atom_class__(self):
        return StructureAtom
