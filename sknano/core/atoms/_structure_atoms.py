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

from operator import attrgetter

from ._cn_atoms import CNAtom, CNAtoms
from ._id_atoms import IDAtom, IDAtoms
from ._charged_atoms import ChargedAtom, ChargedAtoms
from ._image_atoms import ImageAtom, ImageAtoms
from ._type_atoms import TypeAtom, TypeAtoms
from ._kdtree_atoms import KDTreeAtomMixin, KDTreeAtomsMixin
from ._poav_atoms import POAVAtomMixin, POAVAtomsMixin
from ._neighbor_atoms import NeighborAtomMixin, NeighborAtomsMixin
from ._periodic_atoms import PBCAtom, PBCAtoms
from ._lattice_atoms import LatticeAtom, LatticeAtoms
from ._xyz_atoms import XYZAtom, XYZAtoms
from ._velocity_atoms import VelocityAtom, VelocityAtoms

from ._bonds import Bonds

__all__ = ['StructureAtom', 'StructureAtoms']


class StructureAtom(NeighborAtomMixin, POAVAtomMixin, KDTreeAtomMixin,
                    CNAtom, VelocityAtom, ImageAtom, PBCAtom, LatticeAtom,
                    XYZAtom, ChargedAtom, TypeAtom, IDAtom):
    """An `Atom` class for structure analysis.

    Parameters
    ----------
    element : {str, int}, optional
        A string representation of the element symbol or an integer specifying
        an element atomic number.
    id : int, optional
        atom ID
    mol : int, optional
        molecule ID
    type : int, optional
        atom type
    x, y, z : float, optional
        :math:`x, y, z` components of `StructureAtom` position vector relative
        to origin.
    vx, vy, vz : float, optional
        :math:`v_x, v_y, v_z` components of `StructureAtom` velocity.

    """
    def __init__(self, *args, NN=None, **kwargs):
        super().__init__(*args, **kwargs)

        self._neighbors = None
        if NN is not None:
            self.NN = NN

        self._POAV1 = None
        self._POAV2 = None
        self._POAVR = None
        # self.fmtstr = super().fmtstr + ", NN={NN!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('NN')
        # attrs.extend(['NN', 'bonds'])
        return attrs

    @CNAtom.CN.getter
    def CN(self):
        """`StructureAtom` coordination number."""
        try:
            return self.NN.Natoms
        except AttributeError:
            return super().CN


class StructureAtoms(NeighborAtomsMixin, POAVAtomsMixin, KDTreeAtomsMixin,
                     CNAtoms, VelocityAtoms, ImageAtoms, PBCAtoms,
                     LatticeAtoms, XYZAtoms, ChargedAtoms, TypeAtoms, IDAtoms):
    """An `Atoms` sub-class for structure analysis.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.StructureAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `StructureAtoms`}, optional
        if not `None`, then a list of `StructureAtom` instance objects or an
        existing `StructureAtoms` instance object.
    kNN : :class:`~python:int`
        Number of nearest neighbors to return when querying the kd-tree.
    NNrc : :class:`~python:float`
        Nearest neighbor radius cutoff.

    """
    def __init__(self, atoms=None, kNN=16, NNrc=2.0, **kwargs):

        super().__init__(atoms, **kwargs)
        self.kNN = kNN
        self.NNrc = NNrc
        self.bonds = atoms.bonds if hasattr(atoms, 'bonds') else Bonds()

    @property
    def __atom_class__(self):
        return StructureAtom

    def sort(self, key=attrgetter('CN', 'v', 'i', 'r', 'q', 'type', 'mol',
                                  'id', 'mass', 'Z', 'element'),
             reverse=False):
        super().sort(key=key, reverse=reverse)

    def compute_rdf(self):
        pass

    @property
    def volume(self):
        """Volume of region containing atoms."""
        try:
            return self._volume
        except AttributeError:
            return self.bounds.volume

    @volume.setter
    def volume(self, value):
        self._volume = float(value)
