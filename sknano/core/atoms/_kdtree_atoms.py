# -*- coding: utf-8 -*-
"""
===============================================================================
Atoms class for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atoms`)
===============================================================================

An `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._kdtree_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers

import numpy as np

try:
    from scipy.spatial import KDTree
except ImportError:
    raise ImportError('Install scipy version >= 0.13.0 to allow '
                      'nearest-neighbor queries between atoms.')

from ._extended_atoms import XAtoms
from ._neighbor_atoms import NeighborAtoms

__all__ = ['KDTAtoms']


class KDTAtoms(XAtoms):
    """An `Atoms` class for KDTree analysis."""

    def __init__(self, **kwargs):

        super(KDTAtoms, self).__init__(**kwargs)

        self._kNN = 3
        self._NN_cutoff = 0

    @property
    def atom_tree(self):
        """Return the :py:class:`~scipy:scipy.spatial.KDTree` of coords."""
        try:
            return KDTree(self.coords)
        except ValueError:
            return None

    @property
    def coordination_numbers(self):
        """Return array of `KDTAtom` coordination numbers."""
        self._update_coordination_numbers()
        return np.asarray([atom.CN for atom in self])

    def _update_coordination_numbers(self):
        """Update `KDTAtom` coordination numbers."""
        try:
            NN_d, NN_i = \
                self.query_atom_tree(n=self.kNN,
                                     cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self):
                atom.CN = 0
                for d in NN_d[i]:
                    if d < self.NN_cutoff:
                        atom.CN += 1
        except ValueError:
            pass

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `KDTAtom`."""
        self._update_nearest_neighbors()
        return np.asarray([atom.NN for atom in self])

    def _update_nearest_neighbors(self):
        """Update `KDTAtom` nearest-neighbors."""
        try:
            NN_d, NN_i = self.query_atom_tree(n=self.kNN,
                                              cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self):
                atom.NN = NeighborAtoms()
                for j, d in enumerate(NN_d[i]):
                    if d < self.NN_cutoff:
                        atom.NN.append(self[NN_i[i][j]])
        except ValueError:
            pass

    def query_atom_tree(self, n=6, eps=0, p=2, cutoff_radius=np.inf):
        """Query atom tree for nearest neighbors distances and indices.

        Parameters
        ----------
        n : integer
            The number of nearest neighbors to return.
        eps : nonnegative float
            Return approximate nearest neighbors; the kth returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real kth nearest neighbor.
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        cutoff_radius : nonnegative float
            Return only neighbors within this distance. This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns
        -------
        NN_d : array of floats
            The distances to the nearest neighbors.
        NN_i : array of integers
            The locations of the neighbors in self.atom_tree.data. NN_i is the
            same shape as NN_d.

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            d, i = atom_tree.query(self.coords,
                                   k=n+1, eps=eps, p=p,
                                   distance_upper_bound=cutoff_radius)
            return d[:, 1:], i[:, 1:]

    @property
    def kNN(self):
        """Number of nearest neighbors to return when querying the kd-tree."""
        return self._kNN

    @kNN.setter
    def kNN(self, value):
        """Set maximum number of nearest neighbors to return when querying
        the kd-tree."""
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected an integer >= 0')
        self._kNN = int(value)

    @property
    def NN_cutoff(self):
        """Only return neighbors within this distance when querying the
        kd-tree."""
        return self._NN_cutoff

    @NN_cutoff.setter
    def NN_cutoff(self, value):
        """Set the cutoff distance to check for neighest neighbors."""
        if not (isinstance(value, numbers.Number) and value >= 0):
            raise TypeError('Expected a real number greater >= 0')
        self._NN_cutoff = value

    def select(self, **kwargs):
        pass

    def select_within(self, volume):
        pass
