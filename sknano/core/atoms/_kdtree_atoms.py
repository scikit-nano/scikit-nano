# -*- coding: utf-8 -*-
"""
===============================================================================
Container class for `KDTAtom`\ s (:mod:`sknano.core.atoms._kdtree_atoms`)
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
except ImportError:
    raise ImportError('Install scipy version >= 0.13.0 to allow '
                      'nearest-neighbor queries between atoms.')

# from ._extended_atoms import XAtoms
from ._cn_atoms import CNAtoms
from ._id_atoms import IDAtoms
from ._xyz_atoms import XYZAtoms
from ._bonds import Bonds
from ._kdtree_atom import KDTAtom
from ._neighbor_atoms import NeighborAtoms

__all__ = ['KDTAtoms']


class KDTAtoms(CNAtoms, XYZAtoms, IDAtoms):
    """An `Atoms` sub-class for KDTree analysis.

    Sub-class of `XAtoms` class, and a container class for lists of
    :class:`~sknano.core.atoms.KDTAtom` instances.

    Parameters
    ----------
    atoms : {None, sequence, `KDTAtoms`}, optional
        if not `None`, then a list of `KDTAtom` instance objects or an
        existing `KDTAtoms` instance object.
    kNN : :class:`~python:int`
        Number of nearest neighbors to return when querying the kd-tree.
    NNrc : :class:`~python:float`
        Nearest neighbor radius cutoff.

    """
    def __init__(self, atoms=None, kNN=16, NNrc=2.0, **kwargs):

        super().__init__(atoms, **kwargs)
        self.kNN = kNN
        self.NNrc = NNrc
        self.bonds = atoms.bonds if hasattr(atoms, 'bonds') else Bonds()

    @property
    def __atom_class__(self):
        return KDTAtom

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
            return KDTree(self.coords)
        except ValueError:
            return None

    @property
    def neighbors(self):
        return np.asarray([atom.neighbors for atom in self])

    @property
    def neighbor_distances(self):
        distances = []
        [distances.extend(atom.neighbor_distances.tolist()) for atom in self]
        return distances

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `KDTAtom`."""
        # self._update_neighbors()
        return np.asarray([atom.NN for atom in self])

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
        self.update_neighbors()
        self.update_bonds()

    def update_neighbors(self):
        """Update :attr:`KDTAtom.NN`."""
        try:
            NNd, NNi = self.query_atom_tree(k=self.kNN, rc=self.NNrc)
            for j, atom in enumerate(self):
                # atom.neighbors = self.__class__(**self.kwargs)

                # atom.neighbors = NeighborAtoms()
                # [atom.neighbors.append(self[NNi[j][k]])
                #  for k, d in enumerate(NNd[j]) if d <= self.NNrc]

                atom.neighbors = \
                    NeighborAtoms([self[NNi[j][k]] for k, d in
                                   enumerate(NNd[j]) if d <= self.NNrc],
                                  casttype=False)

                # atom.neighbors = \
                #     self.__class__([self[NNi[j][k]] for k, d in
                #                     enumerate(NNd[j]) if d <= self.NNrc],
                #                    casttype=False, **self.kwargs)

                atom.neighbor_distances = \
                    [NNd[j][k] for k, d in enumerate(NNd[j]) if d <= self.NNrc]
                atom.NN = atom.neighbors
        except ValueError:
            pass

    def update_bonds(self):
        """Update :attr:`KDTAtom.bonds`."""
        self.bonds.clear()
        [self.bonds.extend(atom.bonds) for atom in self]

    def neighbor_counts(self, r):
        return np.asarray([KDTree([atom.r]).count_neighbors(
                           self.filtered(self.ids != atom.id).atom_tree,
                           np.asarray(r)) for atom in self])
