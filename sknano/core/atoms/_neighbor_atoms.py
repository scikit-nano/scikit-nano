# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin `Atoms` class for NN analysis (:mod:`sknano.core.atoms._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import numpy as np

from sknano.core import ordinal_form

__all__ = ['NeighborAtomsMixin']


class NeighborAtomsMixin:
    """Mixin :class:`~sknano.core.atoms.Atoms` class for NN analysis."""

    def update_neighbors(self, cutoffs=None):
        super().update_attrs()
        if cutoffs is not None:
            for n, cutoff in enumerate(cutoffs, start=1):
                try:
                    NNd, NNi = self.query_atom_tree(k=self.kNN, rc=cutoff)
                    for j, atom in enumerate(self):
                        neighbors = self.__class__(
                            [self[NNi[j][k]] for k, d in
                             enumerate(NNd[j]) if d <= cutoff],
                            casttype=False, **self.kwargs)
                        neighbors.distances = \
                            [NNd[j][k] for k, d in enumerate(NNd[j])
                             if d <= cutoff]
                        atom.set_neighbors(n, neighbors)
                except ValueError:
                    pass

    @property
    def first_neighbors(self):
        return [atom.first_neighbors for atom in self]

    @property
    def second_neighbors(self):
        return [atom.second_neighbors for atom in self]

    @property
    def third_neighbors(self):
        return [atom.third_neighbors for atom in self]

    def get_neighbors(self, n):
        return [getattr(atom, '_{}_neighbors'.format(ordinal_form(n)))
                for atom in self]
