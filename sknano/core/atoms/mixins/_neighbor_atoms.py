# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for NN analysis (:mod:`sknano.core.atoms.mixins._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.mixins._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
# import numbers

import numpy as np

from sknano.core import ordinal_form
from sknano.core.math import Vector, Vectors
from .._atoms import Atoms
from .._cn_atoms import CNAtom
from ._kdtree_atoms import KDTreeAtomsMixin

__all__ = ['NeighborAtomMixin', 'NeighborAtomsMixin']


class NeighborAtomMixin:
    """Mixin `Atom` class for neighbor analysis."""

    @CNAtom.CN.getter
    def CN(self):
        """`NeighborAtom` coordination number."""
        try:
            return self.neighbors.Natoms
        except AttributeError:
            return super().CN

    @property
    def Nneighbors(self):
        """Number of neighbors."""
        return self.neighbors.Natoms

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
        """First nearest neighbors."""
        return self.get_nth_nearest_neighbors(1)

    @first_neighbors.setter
    def first_neighbors(self, value):
        self.set_nth_nearest_neighbors(1, value)

    @property
    def second_neighbors(self):
        """Second nearest neighbors"""
        return self.get_nth_nearest_neighbors(2)

    @second_neighbors.setter
    def second_neighbors(self, value):
        self.set_nth_nearest_neighbors(2, value)

    @property
    def third_neighbors(self):
        """Third nearest neighbors"""
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
        """Set nth nearest neighbors"""
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


class NeighborAtomsMixin(KDTreeAtomsMixin):
    """Mixin `Atoms` class for neighbor analysis."""

    @property
    def neighbors_analyzed(self):
        """Return `True` if neighbors have been analyzed."""
        return self._neighbors_analyzed

    @neighbors_analyzed.setter
    def neighbors_analyzed(self, value):
        if not isinstance(value, bool):
            raise ValueError('Expected a boolean value.')
        self._neighbors_analyzed = self.kwargs['neighbors_analyzed'] = value

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
        """Neighbor distances"""
        distances = []
        # [distances.extend(atom.neighbor_distances.tolist()) for atom in self]
        [distances.extend(atom.neighbors.distances.tolist()) for atom in self]
        return np.asarray(distances)

    @property
    def NNN(self):
        """Number of first nearest-neighbors."""
        return len(self.nn_idx)

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `KDTAtom`."""
        # self._update_neighbors()
        return np.asarray([atom.NN for atom in self])

    @property
    def first_neighbors(self):
        """First neighbors."""
        return [atom.first_neighbors for atom in self]

    @property
    def second_neighbors(self):
        """Second neighbors."""
        return [atom.second_neighbors for atom in self]

    @property
    def third_neighbors(self):
        """Third neighbors."""
        return [atom.third_neighbors for atom in self]

    def get_nth_nearest_neighbors(self, n, exclusive=True):
        """Return `n`th nearest neighbors."""
        return [getattr(atom, '_{}_neighbors'.format(ordinal_form(n)))
                for atom in self]

    def update_attrs(self, **kwargs):
        """Update :class:`NeighborAtom`\ s attributes."""
        self.update_neighbors(**kwargs)
        # self.update_bonds()

    def update_neighbors(self, cutoffs=None, **kwargs):
        """Update :attr:`NeighborAtom.neighbors`."""
        if 'cutoff' in kwargs and cutoffs is None:
            cutoffs = kwargs['cutoff']
            del kwargs['cutoff']
        if cutoffs is None:
            cutoffs = []
        elif not isinstance(cutoffs, list):
            cutoffs = [cutoffs]
        if not any([np.allclose(self.NNrc, cutoff) for cutoff in cutoffs]):
            cutoffs.append(self.NNrc)

        if not self._neighbors_analyzed:
            self.neighbors_analyzed = True

        for n, cutoff in enumerate(sorted(cutoffs), start=1):
            try:
                NNd, NNi = self.query_atom_tree(k=self.kNN, rc=cutoff)
                for i, atom in enumerate(self):
                    neighbors = \
                        self.__class__([self[NNi[i][j]] for j, d in
                                        enumerate(NNd[i]) if d <= cutoff],
                                       update_item_class=False, **self.kwargs)
                    distances = \
                        [NNd[i][j] for j, d in enumerate(NNd[i])
                         if d <= cutoff]
                    neighbors.distances = distances
                    if np.allclose(cutoff, self.NNrc):
                        atom.neighbors = neighbors
                    atom.set_nth_nearest_neighbors(n, neighbors[:])

            except ValueError:
                pass

    def update_neighbor_lists(self):
        """Update neighbor lists"""
        self._update_nn_lists()
        self._update_nn_seed()
        self._update_nn_vectors()
        self._update_nn_adjacency_matrix()

    @property
    def nn_adjacency_matrix(self):
        """Return nearest-neighbor adjacency matrix."""
        if self._nn_adjacency_matrix is None:
            self._update_nn_adjacency_matrix()
        return self._nn_adjacency_matrix

    @property
    def nn_adjacency_map(self):
        """Return nearest-neighbor adjacency map"""
        if self._nn_adjacency_map is None:
            self._update_nn_adjacency_map()
        return self._nn_adjacency_map

    @property
    def nn_adjacency_list(self):
        """Return nearest-neighbor adjacency list"""
        if self._nn_adjacency_list is None:
            self._update_nn_adjacency_list()
        return self._nn_adjacency_list

    @property
    def nn_seed(self):
        """Return nearest-neighbor seed list"""
        if self._nn_seed is None:
            self._update_nn_seed()
        return self._nn_seed

    @property
    def nn_vectors(self):
        """Return nearest-neighbor vectors."""
        if self._nn_vectors is None:
            self._update_nn_vectors()
        return self._nn_vectors

    def _update_nn_lists(self, NNi=None, cutoff=None):
        if NNi is None:
            if cutoff is None:
                cutoff = self.NNrc
            _, NNi = self.query_atom_tree(k=self.kNN, rc=cutoff)
        Natoms = self.Natoms
        idx = []
        nn_idx = []
        for i, nn_indices in enumerate(NNi):
            for ni in nn_indices[:]:
                if ni < Natoms:
                    idx.append(i)
                    nn_idx.append(ni)
        self.idx = np.asarray(idx, dtype=int)
        self.nn_idx = np.asarray(nn_idx, dtype=int)

    def _update_nn_seed(self):
        n = self.Natoms
        nnn = self.NNN
        idx = self.idx
        nn_seed = (n + 1) * [0]
        for k in range(n):
            nn_seed[k] = -1
        nn_seed[n] = nnn
        nn_seed[idx[0]] = 0
        for k in range(1, nnn):
            if idx[k] != idx[k - 1]:
                nn_seed[idx[k]] = k
        self._nn_seed = np.asarray(nn_seed, dtype=int)

    def _update_nn_adjacency_matrix(self):
        Natoms = self.Natoms
        nn_seed = self.nn_seed
        nn_idx = self.nn_idx
        nn_adjacency_matrix = np.zeros((Natoms, Natoms), dtype=int)
        # self._update_nn_lists()

        def _update_nn_adjacency_map(nn_map, indices, nindices):
            for i in indices:
                for ni in nindices:
                    if ni not in nn_map:
                        nn_map[ni] = nn_map[i] + 1

        for i, atom in enumerate(self):
            nn_map = atom.nn_adjacency_map = {}
            nn_map[i] = 0
            indices = [i]
            nindices = [nn_idx[si] for si in range(nn_seed[i], nn_seed[i+1])]
            _update_nn_adjacency_map(nn_map, indices, nindices)
            while len(nindices) > 0:
                indices = nindices[:]
                nindices = []
                for ni in indices[:]:
                    nindices.extend(
                        [nn_idx[nsi] for nsi in
                         range(nn_seed[ni], nn_seed[ni+1])
                         if nn_idx[nsi] not in nn_map])
                    _update_nn_adjacency_map(nn_map, indices, nindices)
            [nn_adjacency_matrix[i].__setitem__(j, nn_map[j])
             for j in range(Natoms) if j in nn_map]
        self._nn_adjacency_matrix = nn_adjacency_matrix

    def _update_nn_adjacency_map(self):
        Natoms = self.Natoms
        nn_idx = self.nn_idx
        nn_seed = self.nn_seed
        nn_adjacency_map = OrderedDict()
        for i in range(Natoms):
            nn_adjacency_map[i] = \
                [nn_idx[si] for si in range(nn_seed[i], nn_seed[i+1])]
        self._nn_adjacency_map = nn_adjacency_map

    def _update_nn_adjacency_list(self):
        Natoms = self.Natoms
        nn_idx = self.nn_idx
        nn_seed = self.nn_seed
        nn_adjacency_list = []
        for i in range(Natoms):
            nn_adjacency_list.append(
                [nn_idx[si] for si in range(nn_seed[i], nn_seed[i+1])])
        self._nn_adjacency_list = nn_adjacency_list

    def _update_nn_vectors(self):
        from sknano.core.crystallography import pbc_diff
        idx = self.idx
        nn_idx = self.nn_idx
        lattice = self.lattice
        pbc = np.any(self.pbc)
        nn_vectors = Vectors()
        for i, ni in zip(idx, nn_idx):
            if pbc:
                fdiff = pbc_diff(self[ni].rs, self[i].rs)
                r = Vector(lattice.fractional_to_cartesian(fdiff))
            else:
                r = self[ni].r - self[i].r
            nn_vectors.append(r)
        self._nn_vectors = nn_vectors

    # def update_bonds(self):
    #     """Update :attr:`KDTAtom.bonds`."""
    #     self.bonds.clear()
    #     [self.bonds.extend(atom.bonds) for atom in self]
