# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atoms`)
===============================================================================

An `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._kdtree_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from ._extended_atoms import XAtoms
from ._mixins import NNAtomsMixin
from ._bonds import Bonds
from ._kdtree_atom import KDTAtom

__all__ = ['KDTAtoms']


class KDTAtoms(XAtoms, NNAtomsMixin):
    """An `Atoms` sub-class for KDTree analysis.

    Sub-class of `XAtoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.KDTAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `KDTAtoms`}, optional
        if not `None`, then a list of `KDTAtom` instance objects or an
        existing `KDTAtoms` instance object.

    """
    def __init__(self, kNN=3, NNrc=2.0, atoms=None):

        super().__init__(atoms=atoms)

        self.kNN = kNN
        self.NNrc = NNrc

        try:
            self.bonds = atoms.bonds
        except AttributeError:
            self.bonds = Bonds()

    @property
    def __atom_class__(self):
        return KDTAtom
