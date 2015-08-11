# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin Atom classes for NN analysis (:mod:`sknano.core.atoms._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from sknano.core import ordinal_form

__all__ = ['NeighborAtomMixin', 'NeighborAtomsMixin']


class NeighborAtomMixin:
    """Mixin class for neighbor analysis."""
    @property
    def first_neighbors(self):
        return self._first_neighbors

    @first_neighbors.setter
    def first_neighbors(self, value):
        self._first_neighbors = value

    @property
    def second_neighbors(self):
        return self._second_neighbors

    @second_neighbors.setter
    def second_neighbors(self, value):
        self._second_neighbors = value

    @property
    def third_neighbors(self):
        return self._third_neighbors

    @third_neighbors.setter
    def third_neighbors(self, value):
        self._third_neighbors = value

    def get_neighbors(self, n):
        if n == 1:
            return self.first_neighbors

        elif n == 2:
            return self.second_neighbors

        elif n == 3:
            return self.third_neighbors

        else:
            return getattr(self, '_{}_neighbors'.format(ordinal_form(n)))

    def set_neighbors(self, n, neighbors):
        if n == 1:
            self._first_neighbors = neighbors
        elif n == 2:
            try:
                for neighbor in self.first_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
            except TypeError:
                pass
            self._second_neighbors = neighbors
        elif n == 3:
            try:
                for neighbor in self.first_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
                for neighbor in self.second_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
            except TypeError:
                pass
            self._third_neighbors = neighbors
        else:
            try:
                for neighbor in self.first_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
                for neighbor in self.second_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
                for neighbor in self.third_neighbors:
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
            except TypeError:
                pass
            for i in range(4, n):
                for neighbor in self.get_neighbors(i):
                    if neighbor in neighbors:
                        neighbors.remove(neighbor)
            setattr(self, '_{}_neighbors'.format(ordinal_form(n)), neighbors)


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
