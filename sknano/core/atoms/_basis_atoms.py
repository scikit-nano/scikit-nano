# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for crystal structures (:mod:`sknano.core.atoms._basis_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._basis_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._lattice_atoms import LatticeAtom, LatticeAtoms
from ._periodic_atoms import PBCAtom, PBCAtoms
from ._xyz_atoms import XYZAtom, XYZAtoms

__all__ = ['BasisAtom', 'BasisAtoms']


class BasisAtom(PBCAtom, LatticeAtom, XYZAtom):
    """An abstract object representation of a crystal structure basis atom.

    Parameters
    ----------
    lattice : :class:`~sknano.core.crystallography.Crystal3DLattice`
    xs, ys, zs : float
    """
    pass


class BasisAtoms(PBCAtoms, LatticeAtoms, XYZAtoms):
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
