# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for KDTree analysis (:mod:`sknano.core.atoms._kdtree_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._kdtree_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

import numbers
import warnings
warnings.filterwarnings('ignore', "Mean of empty slice.")
warnings.filterwarnings('ignore',
                        'invalid value encountered in double_scalars')

import numpy as np

try:
    from scipy.spatial import KDTree
    # from scipy.spatial import cKDTree as KDTree
except ImportError:
    raise ImportError('Install scipy version >= 0.13.0 to allow '
                      'nearest-neighbor queries between atoms.')

from sknano.core import dedupe, flatten
from sknano.core.analysis import PeriodicKDTree

__all__ = ['KDTreeAtomsMixin']


class KDTreeAtomsMixin:
    """Mixin Atoms class for KDTree analysis."""

    @property
    def kNN(self):
        """Max number of nearest-neighbors to return from kd-tree search."""
        return self._kNN

    @kNN.setter
    def kNN(self, value):
        if not isinstance(value, numbers.Number):
            raise TypeError('Expected an integer >= 0')
        self._kNN = self.kwargs['kNN'] = int(value)

    @property
    def NNrc(self):
        """Nearest neighbor radius cutoff."""
        return self._NNrc

    @NNrc.setter
    def NNrc(self, value):
        if not (isinstance(value, numbers.Number) and value >= 0):
            raise TypeError('Expected a real number greater >= 0')
        self._NNrc = self.kwargs['NNrc'] = value

    @property
    def atom_tree(self):
        """:class:`~scipy:scipy.spatial.KDTree` of :attr:`~XAtoms.coords.`"""
        try:
            if np.any(self.pbc):
                boxsize = np.asarray(self.lattice.lengths)
                return PeriodicKDTree(np.asarray(self.coords), boxsize=boxsize,
                                      lattice=self.lattice)
            return KDTree(np.asarray(self.coords))
        except ValueError as e:
            print(e)
            return None

    def query_atom_tree(self, k=16, eps=0, p=2, rc=np.inf):
        """Query atom tree for nearest neighbors distances and indices.

        Parameters
        ----------
        k : integer
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
        rc : nonnegative float
            Radius cutoff. Return only neighbors within this distance.
            This is used to prune tree searches, so if you are doing a series
            of nearest-neighbor queries, it may help to supply the distance to
            the nearest neighbor of the most recent point.

        Returns
        -------
        d : array of floats
            The distances to the nearest neighbors, sorted by distance.
        i : array of integers
            The locations of the neighbors in self.atom_tree.data. `i` is the
            same shape as `d`.

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            d, i = atom_tree.query(np.asarray(self.coords), k=k+1, eps=eps,
                                   p=p, distance_upper_bound=rc)
            return d[:, 1:], i[:, 1:]

    def query_ball_point(self, pts, r, p=2.0, eps=0):
        """Find all `Atoms` within distance `r` of point(s) `pts`.

        Parameters
        ----------
        pts : :class:`~sknano.core.math.Point`
            The :class:`~sknano.core.math.Point` or
            :class:`~sknano.core.math.Points` to search for neighbors of.
        r : positive :class:`~python:float`
            The radius of :class:`~sknano.core.atoms.KDTAtoms` to return
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        eps : nonnegative :class:`~python:float`, optional
            Approximate search.

        Returns
        -------
        :class:`~sknano.core.atoms.Atoms`

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            NNi = \
                list(dedupe(flatten(
                    atom_tree.query_ball_point(pts, r, p=p, eps=eps))))

        return self.__class__(atoms=np.asarray(self)[NNi].tolist(),
                              **self.kwargs)

    def query_ball_tree(self, other, r, p=2.0, eps=0):
        """Find all pairs of `Atoms` whose distance is at more `r`.

        Parameters
        ----------
        other : :class:`~scipy:scipy.spatial.KDTree`
            The tree containing points to search against
        r : positive :class:`~python:float`
            The radius of :class:`~sknano.core.atoms.KDTAtoms` to return
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        eps : nonnegative :class:`~python:float`, optional
            Approximate search.

        Returns
        -------
        :class:`~sknano.core.atoms.Atoms`

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            # NNi = \
            #     atom_tree.query_ball_tree(other.atom_tree, r, p=p, eps=eps)
            NNi = \
                other.atom_tree.query_ball_tree(atom_tree, r, p=p, eps=eps)
            NNi = list(dedupe(flatten(NNi)))

        return self.__class__(atoms=np.asarray(self)[NNi].tolist(),
                              **self.kwargs)

    def query_pairs(self, r, p=2.0, eps=0):
        """Find all pairs of points within a distance `r`.

        Parameters
        ----------
        r : positive float
        p : float, optional
        eps : float, optional

        Returns
        -------
        results : set

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            return atom_tree.query_pairs(r, p=p, eps=eps)

    def count_neighbors(self, other, r, p=2.0):
        """Count how many nearby neighbor pairs can be formed.

        Count the number of pairs (x1, x2) that can be formed, with
        `x1` drawn from `self` and `x2` drawn from `other`, and
        where ``distance(x1, x2, p) <= r``.

        Parameters
        ----------
        other : `KDTree`
        r : float or one-dimensional array of floats
        p : float, 1<=p<=infinity, optional

        Returns
        -------
        result : int or 1-D array of ints

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            return atom_tree.count_neighbors(other, r, p=p)

    def count_neighbors_in_self(self, r, p=2.0):
        """Count number of neighbor pairs for each atom in self.

        Parameters
        ----------
        r : float or one-dimensional array of floats
        p : float, 1<=p<=infinity, optional

        Returns
        -------
        result : 1-D array of ints

        """
        r = np.asarray(r)
        ids = self.ids
        return np.asarray([KDTree([atom.r]).count_neighbors(
                           self.filtered(ids != atom.id).atom_tree,
                           r, p=p) for atom in self])
