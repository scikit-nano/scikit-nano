# -*- coding: utf-8 -*-
"""
===============================================================================
Atom classes for NN analysis (:mod:`sknano.core.atoms._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numbers

from ._atoms import Atom, Atoms
from .mixins import NeighborAtomMixin, NeighborAtomsMixin

__all__ = ['NeighborAtom', 'NeighborAtoms']


class NeighborAtom(NeighborAtomMixin, Atom):
    """An `Atom` class for neighbor analysis."""
    def __init__(self, *args, neighbors=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._neighbors = neighbors
        self.nn_adjacency_map = {}
        # self.nn_adjacency_list = []
        # self.fmtstr = super().fmtstr + ", neighbors={neighbors!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('neighbors')
        return attrs

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(neighbors=self.neighbors))
        return super_dict


class NeighborAtoms(NeighborAtomsMixin, Atoms):

    def __init__(self, *args, kNN=16, NNrc=2.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.kNN = kNN
        self.NNrc = NNrc

        self.idx = []
        self.nn_idx = []

        self._nn_adjacency_matrix = None
        self._nn_adjacency_map = None
        self._nn_adjacency_list = None
        self._nn_vectors = None
        self._nn_seed_list = None

    @property
    def __atom_class__(self):
        return NeighborAtom
