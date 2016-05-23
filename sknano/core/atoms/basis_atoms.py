# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for crystal structures (:mod:`sknano.core.atoms.basis_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.basis_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from .id_atoms import IDAtom, IDAtoms
from .lattice_atoms import LatticeAtom, LatticeAtoms
from .xyz_atoms import XYZAtom, XYZAtoms
from .mixins import AtomTransformationsMixin, AtomsTransformationsMixin, \
    BoundingRegionsMixin

__all__ = ['BasisAtom', 'BasisAtoms']


class BasisAtom(AtomTransformationsMixin, LatticeAtom, XYZAtom, IDAtom):
    """An `Atom` sub-class for a crystal structure basis atom.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : :class:`~python:float`
    """
    @property
    def __atoms_class__(self):
        return BasisAtoms


class BasisAtoms(AtomsTransformationsMixin, BoundingRegionsMixin,
                 LatticeAtoms, XYZAtoms, IDAtoms):
    """An `Atoms` sub-class for crystal structure basis atoms.

    Sub-class of `Atoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.BasisAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `BasisAtoms`}, optional
        if not `None`, then a list of `BasisAtom` instance objects or an
        existing `BasisAtoms` instance object.

    """
    @property
    def __atom_class__(self):
        return BasisAtom
