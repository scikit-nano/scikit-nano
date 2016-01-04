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

# from functools import total_ordering
# from operator import attrgetter
# import numbers

import numpy as np

from sknano.core import ordinal_form
from ._atoms import Atoms
from ._cn_atoms import CNAtom
from ._kdtree_atoms import KDTreeAtomsMixin
from ._bonds import Bond, Bonds

__all__ = ['NeighborAtomMixin', 'NeighborAtomsMixin']


class NeighborAtomMixin:
    """An `Atom` class for neighbor analysis."""
    def __init__(self, *args, neighbors=None, **kwargs):
        super().__init__(*args, **kwargs)

        self._neighbors = neighbors
        # self.fmtstr = super().fmtstr + ", neighbors={neighbors!r}"

    def __dir__(self):
        attrs = super().__dir__()
        attrs.append('neighbors')
        return attrs

    @CNAtom.CN.getter
    def CN(self):
        """`NeighborAtom` coordination number."""
        try:
            return self.neighbors.Natoms
        except AttributeError:
            return super().CN

    @property
    def NN(self):
        """Nearest-neighbor `Atoms`."""
        try:
            return self.neighbors
        except AttributeError:
            return None

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        self.neighbors = value

    @NN.deleter
    def NN(self):
        del self.neighbors

    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        try:
            return Bonds([Bond(self, nn) for nn in self.NN])
        except (AttributeError, TypeError):
            return Bonds()

    @property
    def neighbors(self):
        """Neighbor atoms."""
        return self._neighbors

    @neighbors.setter
    def neighbors(self, value):
        if not isinstance(value, Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._neighbors = value

    @neighbors.deleter
    def neighbors(self):
        del self._neighbors

    @property
    def neighbor_distances(self):
        """Neighbor atom distances."""
        return self.neighbors.distances

    # @neighbor_distances.setter
    # def neighbor_distances(self, value):
    #     self._neighbor_distances = np.asarray(value)

    @property
    def first_neighbors(self):
        return self.get_nth_nearest_neighbors(1)

    @first_neighbors.setter
    def first_neighbors(self, value):
        self.set_nth_nearest_neighbors(1, value)

    @property
    def second_neighbors(self):
        return self.get_nth_nearest_neighbors(2)

    @second_neighbors.setter
    def second_neighbors(self, value):
        self.set_nth_nearest_neighbors(2, value)

    @property
    def third_neighbors(self):
        return self.get_nth_nearest_neighbors(3)

    @third_neighbors.setter
    def third_neighbors(self, value):
        self.set_nth_nearest_neighbors(3, value)

    def get_nth_nearest_neighbors(self, n, exclusive=True):
        """Get the `n`th set of neighbors."""
        if exclusive:
            return getattr(self, '_{}_neighbors'.format(ordinal_form(n)))
        else:
            neighbors = self.__class__(**self.kwargs)
            for i in range(1, n + 1):
                neighbors.extend(
                    getattr(self, '_{}_neighbors'.format(ordinal_form(n))))
            return neighbors

    def set_nth_nearest_neighbors(self, n, neighbors):
        # Since we're setting the nth set of neighbors, we will remove
        # all 1..n-1 neighbors from the nth set.
        # distances = distances
        # neighbors = neighbors
        for i in range(1, n):
            for neighbor in self.get_nth_nearest_neighbors(i):
                if neighbor in neighbors:
                    # print(neighbors.index(neighbor))
                    # distances.remove(neighbors.index(neighbor))
                    neighbors.remove(neighbor)
        # neighbors.distances = distances
        # print(type(neighbors))

        setattr(self, '_{}_neighbors'.format(ordinal_form(n)), neighbors)

    def todict(self):
        super_dict = super().todict()
        super_dict.update(dict(neighbors=self.neighbors))
        return super_dict


class NeighborAtomsMixin(KDTreeAtomsMixin):
    """Atoms class for KDTree neighbor analysis.

    Parameters
    ----------
    atoms : {None, sequence, `NeighborAtoms`}, optional
        if not `None`, then a list of `StructureAtom` instance objects or an
        existing `NeighborAtoms` instance object.
    kNN : :class:`~python:int`
        Number of nearest neighbors to return when querying the kd-tree.
    NNrc : :class:`~python:float`
        Nearest neighbor radius cutoff.

    """
    def __init__(self, *args, kNN=16, NNrc=2.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.kNN = kNN
        self.NNrc = NNrc
        self.bonds = Bonds()
        # self.bonds = atoms.bonds if hasattr(atoms, 'bonds') else Bonds()

    @property
    def neighbors(self):
        """Return array of neighbor atoms."""
        return np.asarray([atom.neighbors for atom in self])

    @property
    def distances(self):
        """Neighbor atoms distances."""
        return self._distances

    @distances.setter
    def distances(self, value):
        self._distances = np.asarray(value)

    @property
    def neighbor_distances(self):
        distances = []
        # [distances.extend(atom.neighbor_distances.tolist()) for atom in self]
        [distances.extend(atom.neighbors.distances.tolist()) for atom in self]
        return np.asarray(distances)

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `KDTAtom`."""
        # self._update_neighbors()
        return np.asarray([atom.NN for atom in self])

    @property
    def first_neighbors(self):
        return [atom.first_neighbors for atom in self]

    @property
    def second_neighbors(self):
        return [atom.second_neighbors for atom in self]

    @property
    def third_neighbors(self):
        return [atom.third_neighbors for atom in self]

    def get_nth_nearest_neighbors(self, n, exclusive=True):
        return [getattr(atom, '_{}_neighbors'.format(ordinal_form(n)))
                for atom in self]

    def update_attrs(self, **kwargs):
        """Update :class:`KDTAtom`\ s attributes."""
        self.update_neighbors(**kwargs)
        self.update_bonds()

    def update_neighbors(self, cutoffs=None):
        """Update :attr:`NeighborAtom.neighbors`."""
        if cutoffs is None:
            cutoffs = []
        cutoffs.append(self.NNrc)

        for n, cutoff in enumerate(sorted(cutoffs), start=1):
            try:
                NNd, NNi = self.query_atom_tree(k=self.kNN, rc=cutoff)
                for j, atom in enumerate(self):
                    neighbors = \
                        self.__class__([self[NNi[j][k]] for k, d in
                                        enumerate(NNd[j]) if d <= cutoff],
                                       casttype=False, **self.kwargs)
                    distances = \
                        [NNd[j][k] for k, d in enumerate(NNd[j])
                         if d <= cutoff]
                    neighbors.distances = distances
                    if np.allclose(cutoff, self.NNrc):
                        atom.neighbors = neighbors
                    atom.set_nth_nearest_neighbors(n, neighbors[:])
            except ValueError:
                pass

    def update_bonds(self):
        """Update :attr:`KDTAtom.bonds`."""
        self.bonds.clear()
        [self.bonds.extend(atom.bonds) for atom in self]
