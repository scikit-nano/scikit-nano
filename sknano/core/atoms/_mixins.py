# -*- coding: utf-8 -*-
"""
===============================================================================
Mixin Atom/Atoms classes (:mod:`sknano.core.atoms._mixins`)
===============================================================================

Mixin `Atom` and `Atoms` classes.

.. currentmodule:: sknano.core.atoms._mixins

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
except ImportError:
    raise ImportError('Install scipy version >= 0.13.0 to allow '
                      'nearest-neighbor queries between atoms.')

import sknano.core.atoms
from ._bond import Bond
from ._bonds import Bonds

__all__ = ['NNAtomMixin', 'NNAtomsMixin']


class NNAtomMixin:
    """Mixin :class:`~sknano.core.atoms.Atom` class for nearest-neighbor \
        :class:`~scipy:scipy.spatial.KDTree` analysis."""

    @property
    def NN(self):
        """Nearest-neighbor `Atoms`."""
        try:
            return self._NN
        except AttributeError:
            return None

    @NN.setter
    def NN(self, value):
        """Set nearest-neighbor `Atoms`."""
        if not isinstance(value, sknano.core.atoms.Atoms):
            raise TypeError('Expected an `Atoms` object.')
        self._NN = value

    @property
    def bonds(self):
        """Return atom `Bonds` instance."""
        try:
            return Bonds(bonds=[Bond(self, nn) for nn in self.NN])
        except (AttributeError, TypeError):
            return Bonds()


class NNAtomsMixin:
    """Mixin :class:`~sknano.core.atoms.Atoms` class for nearest-neighbor \
        :class:`~scipy:scipy.spatial.KDTree` analysis."""

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
        self._kNN = self.kwargs['kNN'] = int(value)

    @property
    def NNrc(self):
        """Only return neighbors within this distance when querying the
        kd-tree."""
        return self._NNrc

    @NNrc.setter
    def NNrc(self, value):
        """Set the cutoff distance to check for neighest neighbors."""
        if not (isinstance(value, numbers.Number) and value >= 0):
            raise TypeError('Expected a real number greater >= 0')
        self._NNrc = self.kwargs['NNrc'] = value

    @property
    def atom_tree(self):
        """:class:`~scipy:scipy.spatial.KDTree` of :attr:`~XAtoms.coords.`"""
        try:
            return KDTree(self.coords)
        except ValueError:
            return None

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `KDTAtom`."""
        # self._update_nearest_neighbors()
        return np.asarray([atom.NN for atom in self])

    def query_atom_tree(self, k=3, eps=0, p=2, rc=np.inf):
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
        rc : nonnegative float
            Radius cutoff. Return only neighbors within this distance.
            This is used to prune tree searches, so if you are doing a series
            of nearest-neighbor queries, it may help to supply the distance to
            the nearest neighbor of the most recent point.

        Returns
        -------
        d : array of floats
            The distances to the nearest neighbors.
        i : array of integers
            The locations of the neighbors in self.atom_tree.data. `i` is the
            same shape as `d`.

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            d, i = atom_tree.query(self.coords, k=k+1, eps=eps, p=p,
                                   distance_upper_bound=rc)
            return d[:, 1:], i[:, 1:]

    def query_ball_point(self, pts, r, p=2.0, eps=0):
        """Find all `Atoms` within distance `r` of point(s) `pts`.

        Parameters
        ----------
        pts : `Point`
            The `Point` or `Points` to search for neighbors of.
        r : positive float
            The radius of `KDTAtoms` to return
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use.
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        eps : nonnegative float, optional
            Approximate search.

        Returns
        -------
        list or array of lists
            `KDTAtoms`

        """
        atom_tree = self.atom_tree
        if atom_tree is not None:
            NNi = atom_tree.query_ball_point(pts, r, p=p, eps=eps)

        return self.__class__(atoms=np.asarray(self)[NNi].tolist(),
                              **self.kwargs)

    def update_attrs(self):
        """Update :class:`KDTAtom`\ s attributes."""
        self.update_nearest_neighbors()
        self.update_bonds()

    def update_nearest_neighbors(self):
        """Update :attr:`KDTAtom.NN`."""
        try:
            NNd, NNi = self.query_atom_tree(k=self.kNN, rc=self.NNrc)
            for j, atom in enumerate(self):
                NN = self.__class__(**self.kwargs)
                for k, d in enumerate(NNd[j]):
                    if d <= self.NNrc:
                        NN.append(self[NNi[j][k]])
                atom.NN = NN
        except ValueError:
            pass

    def update_bonds(self):
        """Update :attr:`KDTAtom.bonds`."""
        self.bonds.clear()
        [self.bonds.extend(atom.bonds) for atom in self]
