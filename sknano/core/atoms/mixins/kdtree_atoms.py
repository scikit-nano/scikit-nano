# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin classes for KDTree analysis (:mod:`sknano.core.atoms.kdtree_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms.kdtree_atoms

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import ABCMeta, abstractmethod

# import warnings
# warnings.filterwarnings('ignore', "Mean of empty slice.")
# warnings.filterwarnings('ignore',
#                         'invalid value encountered in double_scalars')

import numpy as np

from sknano.core import dedupe, flatten

__all__ = ['KDTreeAtomsMixin']


class KDTreeAtomsMixin(metaclass=ABCMeta):
    """Mixin Atoms class for KDTree analysis."""

    @property
    @abstractmethod
    def atom_tree(self):
        """:class:`~scipy:scipy.spatial.KDTree` of atom coordinates."""
        raise NotImplementedError

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
