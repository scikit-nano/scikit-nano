# -*- coding: utf-8 -*-
"""
==========================================================================
Periodic Table of Elements (:mod:`sknano.core.refdata._periodic_kdtree`)
==========================================================================

.. currentmodule:: sknano.core.refdata._periodic_kdtree

.. versionadded:: 0.3.22

A wrapper around scipy.spatial.kdtree to implement periodic boundary
conditions

Original source code written by Patrick Varilly, 6 Jul 2012.

Modified and included with scikit-nano v0.3.22.

"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals
__docformat__ = 'restructuredtext en'

import itertools
import heapq

import numpy as np

try:
    from scipy.spatial import KDTree, cKDTree
except ImportError:
    raise ImportError('Install scipy version >= 0.13.0 to allow '
                      'nearest-neighbor queries between atoms.')

__all__ = ['PeriodicKDTree', 'PeriodicCKDTree']


def _generate_images(x, boxsize, distance_upper_bound):
    # Map x onto the canonical unit cell, then produce the relevant
    # mirror images
    xs = x - np.where(boxsize > 0, np.floor(x / boxsize) * boxsize, 0)
    m = len(x)
    images = [xs]
    for i in range(m):
        # print('boxsize[i]: {}'.format(boxsize[i]))
        if boxsize[i] > 0:
            disp = np.zeros(m)
            disp[i] = boxsize[i]

            if distance_upper_bound == np.inf:
                images = list(
                    itertools.chain.from_iterable(
                        (xi + disp, xi, xi - disp) for xi in images))
            else:
                extra_xi = []

                # Point near lower boundary, include image on upper side
                if abs(xs[i]) < distance_upper_bound:
                    extra_xi.extend(xi + disp for xi in images)

                # Point near upper boundary, include image on lower side
                if abs(boxsize[i] - xs[i]) < distance_upper_bound:
                    extra_xi.extend(xi - disp for xi in images)

                images.extend(extra_xi)

    return images


class PeriodicKDTree(KDTree):
    """
    kd-tree for quick nearest-neighbor lookup with periodic boundaries

    See scipy.spatial.kdtree for details on kd-trees.

    Searches with periodic boundaries are implemented by mapping all
    initial data points to one canonical periodic image, building an
    ordinary kd-tree with these points, then querying this kd-tree multiple
    times, if necessary, with all the relevant periodic images of the
    query point.

    Note that to ensure that no two distinct images of the same point
    appear in the results, it is essential to restrict the maximum
    distance between a query point and a data point to half the smallest
    box dimension.
    """

    def __init__(self, data, boxsize=None, lattice=None, center=None,
                 leafsize=10):
        """Construct a kd-tree.

        Parameters
        ----------
        data : array_like, shape (n,k)
            The data points to be indexed. This array is not copied, and
            so modifying this data will result in bogus results.
        boxsize : array_like, shape (k,)
            Size of the periodic box along each spatial dimension.  A
            negative or zero size for dimension k means that space is not
            periodic along k.
        leafsize : positive int
            The number of points at which the algorithm switches over to
            brute-force.

        """

        # Map all points to canonical periodic image
        data = np.asarray(data)
        boxsize = self.boxsize = np.asarray(boxsize)
        wrapped_data = data - \
            np.where(boxsize > 0.0, (np.floor(data / boxsize) * boxsize), 0.0)

        # Calculate maximum distance_upper_bound
        self.max_distance_upper_bound = \
            np.amin(np.where(boxsize > 0, 0.5 * boxsize, np.inf))

        # Set up underlying kd-tree
        super().__init__(wrapped_data, leafsize)

    # The following name is a kludge to override KDTree's private method
    def _KDTree__query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        # This is the internal query method, which guarantees that x
        # is a single point, not an array of points
        #
        # A slight complication: k could be "None", which means "return
        # all neighbors within the given distance_upper_bound".

        # Cap distance_upper_bound
        distance_upper_bound = np.amin([distance_upper_bound,
                                        self.max_distance_upper_bound])

        # Run queries over all relevant images of x
        hits_list = []
        for xi in _generate_images(x, self.boxsize, distance_upper_bound):
            hits_list.append(
                super()._KDTree__query(xi, k, eps, p, distance_upper_bound)
            )

        # Now merge results
        if k is None:
            return list(heapq.merge(*hits_list))
        elif k > 1:
            return heapq.nsmallest(k, itertools.chain(*hits_list))
        elif k == 1:
            return [min(itertools.chain(*hits_list))]
        else:
            raise ValueError("Invalid k in periodic_kdtree._KDTree__query")

    # The following name is a kludge to override KDTree's private method
    def _KDTree__query_ball_point(self, x, r, p=2., eps=0):
        # This is the internal query method, which guarantees that x
        # is a single point, not an array of points

        # Cap r
        r = np.amin([r, self.max_distance_upper_bound])

        # Run queries over all relevant images of x
        results = []
        for xi in _generate_images(x, self.boxsize, r):
            results.extend(
                super()._KDTree__query_ball_point(xi, r, p, eps))
        return results

    # def query_ball_tree(self, other, r, p=2., eps=0):
    #     raise NotImplementedError()

    # def query_pairs(self, r, p=2., eps=0):
    #     raise NotImplementedError()

    # def count_neighbors(self, other, r, p=2.):
    #     raise NotImplementedError()

    # def sparse_distance_matrix(self, other, max_distance, p=2.):
    #     raise NotImplementedError()


class PeriodicCKDTree(cKDTree):
    """
    Cython kd-tree for quick nearest-neighbor lookup with periodic boundaries

    See scipy.spatial.ckdtree for details on kd-trees.

    Searches with periodic boundaries are implemented by mapping all
    initial data points to one canonical periodic image, building an
    ordinary kd-tree with these points, then querying this kd-tree multiple
    times, if necessary, with all the relevant periodic images of the
    query point.

    Note that to ensure that no two distinct images of the same point
    appear in the results, it is essential to restrict the maximum
    distance between a query point and a data point to half the smallest
    box dimension.
    """

    def __init__(self, data, boxsize=None, leafsize=10):
        """Construct a kd-tree.

        Parameters
        ----------
        Construct a kd-tree.

        Parameters
        ----------
        data : array_like, shape (n,k)
            The n data points of dimension k to be indexed. This array is
            not copied unless this is necessary to produce a contiguous
            array of doubles, and so modifying this data will result in
            bogus results.
        boxsize : array_like, shape (k,)
            Size of the periodic box along each spatial dimension.  A
            negative or zero size for dimension k means that space is not
            periodic along k.
        leafsize : positive int
            The number of points at which the algorithm switches over to
            brute-force.

        """

        # Map all points to canonical periodic image
        boxsize = self.boxsize = np.asarray(boxsize)
        data = np.asarray(data)
        wrapped_data = data - \
            np.where(boxsize > 0.0, np.floor(data / boxsize) * boxsize, 0.0)

        # Calculate maximum distance_upper_bound
        self.max_distance_upper_bound = \
            np.amin(np.where(boxsize > 0, 0.5 * boxsize, np.inf))

        # Set up underlying kd-tree
        super().__init__(wrapped_data, leafsize)

    # Ideally, KDTree and cKDTree would expose identical query and __query
    # interfaces.  But they don't, and cKDTree.__query is also inaccessible
    # from Python.  We do our best here to cope.
    def __query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        # This is the internal query method, which guarantees that x
        # is a single point, not an array of points
        #
        # A slight complication: k could be "None", which means "return
        # all neighbors within the given distance_upper_bound".

        # Cap distance_upper_bound
        distance_upper_bound = np.amin([distance_upper_bound,
                                        self.max_distance_upper_bound])

        # Run queries over all relevant images of x
        hits_list = []
        for xi in _generate_images(x, self.boxsize, distance_upper_bound):
            d, i = \
                super().query(xi, k, eps, p, distance_upper_bound)
            if k > 1:
                hits_list.append(list(zip(d, i)))
            else:
                hits_list.append([(d, i)])

        # Now merge results
        if k > 1:
            return heapq.nsmallest(k, itertools.chain(*hits_list))
        elif k == 1:
            return [min(itertools.chain(*hits_list))]
        else:
            raise ValueError("Invalid k in periodic_kdtree._KDTree__query")

    def query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        """
        Query the kd-tree for nearest neighbors

        Parameters
        ----------
        x : array_like, last dimension self.m
            An array of points to query.
        k : integer
            The number of nearest neighbors to return.
        eps : non-negative float
            Return approximate nearest neighbors; the kth returned value
            is guaranteed to be no further than (1+eps) times the
            distance to the real k-th nearest neighbor.
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        distance_upper_bound : nonnegative float
            Return only neighbors within this distance.  This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns
        -------
        d : array of floats
            The distances to the nearest neighbors.
            If x has shape tuple+(self.m,), then d has shape tuple+(k,).
            Missing neighbors are indicated with infinite distances.
        i : ndarray of ints
            The locations of the neighbors in self.data.
            If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
            Missing neighbors are indicated with self.n.

        """
        x = np.asarray(x)
        if np.shape(x)[-1] != self.m:
            raise ValueError("x must consist of vectors of length "
                             "{:d} but has shape {!s}".format(self.m,
                                                              np.shape(x)))
        if p < 1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        retshape = np.shape(x)[:-1]
        if retshape != ():
            if k > 1:
                dd = np.empty(retshape+(k,), dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(retshape+(k,), dtype=np.int)
                ii.fill(self.n)
            elif k == 1:
                dd = np.empty(retshape, dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(retshape, dtype=np.int)
                ii.fill(self.n)
            else:
                raise ValueError("k must be >= 1, or None")
            for c in np.ndindex(retshape):
                hits = self.__query(x[c], k=k, eps=eps, p=p,
                                    distance_upper_bound=distance_upper_bound)
                if k > 1:
                    for j in range(len(hits)):
                        dd[c+(j,)], ii[c+(j,)] = hits[j]
                elif k == 1:
                    if len(hits) > 0:
                        dd[c], ii[c] = hits[0]
                    else:
                        dd[c] = np.inf
                        ii[c] = self.n
            return dd, ii
        else:
            hits = self.__query(x, k=k, eps=eps, p=p,
                                distance_upper_bound=distance_upper_bound)
            if k == 1:
                if len(hits) > 0:
                    return hits[0]
                else:
                    return np.inf, self.n
            elif k > 1:
                dd = np.empty(k, dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(k, dtype=np.int)
                ii.fill(self.n)
                for j in range(len(hits)):
                    dd[j], ii[j] = hits[j]
                return dd, ii
            else:
                raise ValueError("k must be >= 1, or None")

    # Ideally, KDTree and cKDTree would expose identical __query_ball_point
    # interfaces.  But they don't, and cKDTree.__query_ball_point is also
    # inaccessible from Python.  We do our best here to cope.
    def __query_ball_point(self, x, r, p=2., eps=0):
        # This is the internal query method, which guarantees that x
        # is a single point, not an array of points

        # Cap r
        r = min(r, self.max_distance_upper_bound)

        # Run queries over all relevant images of x
        results = []
        for xi in _generate_images(x, self.boxsize, r):
            results.extend(super().query_ball_point(xi, r, p, eps))
        return results

    def query_ball_point(self, x, r, p=2., eps=0):
        """
        Find all points within distance r of point(s) x.

        Parameters
        ----------
        x : array_like, shape tuple + (self.m,)
            The point or points to search for neighbors of.
        r : positive float
            The radius of points to return.
        p : float, optional
            Which Minkowski p-norm to use.  Should be in the range [1, inf].
        eps : nonnegative float, optional
            Approximate search. Branches of the tree are not explored if their
            nearest points are further than ``r / (1 + eps)``, and branches are
            added in bulk if their furthest points are nearer than
            ``r * (1 + eps)``.

        Returns
        -------
        results : list or array of lists
            If `x` is a single point, returns a list of the indices of the
            neighbors of `x`. If `x` is an array of points, returns an object
            array of shape tuple containing lists of neighbors.

        Notes
        -----
        If you have many points whose neighbors you want to find, you may
        save substantial amounts of time by putting them in a
        PeriodicCKDTree and using query_ball_tree.
        """
        x = np.asarray(x).astype(np.float)
        if x.shape[-1] != self.m:
            raise ValueError("Searching for a %d-dimensional point in a "
                             "%d-dimensional KDTree" % (x.shape[-1], self.m))
        if len(x.shape) == 1:
            return self.__query_ball_point(x, r, p, eps)
        else:
            retshape = x.shape[:-1]
            result = np.empty(retshape, dtype=np.object)
            for c in np.ndindex(retshape):
                result[c] = self.__query_ball_point(x[c], r, p, eps)
            return result

    # def query_ball_tree(self, other, r, p=2., eps=0):
    #     raise NotImplementedError()

    # def query_pairs(self, r, p=2., eps=0):
    #     raise NotImplementedError()

    # def count_neighbors(self, other, r, p=2.):
    #     raise NotImplementedError()

    # def sparse_distance_matrix(self, other, max_distance, p=2.):
    #     raise NotImplementedError()
