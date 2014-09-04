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

from sknano.core.math import Vector, vector as vec
#from ._atom_bonds import AtomBonds
from ._bond import Bond
from ._bonds import Bonds
from ._extended_atoms import XAtoms
from ._neighbor_atoms import NeighborAtoms

__all__ = ['KDTAtoms']


class KDTAtoms(XAtoms):
    """An `Atoms` class for KDTree analysis."""

    def __init__(self, **kwargs):

        super(KDTAtoms, self).__init__(**kwargs)

        self._kNN = 3
        self._NNrc = 2.0

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
            NNd, NNi = self.query_atom_tree(k=self.kNN, rc=self.NNrc)
            for i, atom in enumerate(self):
                atom.CN = 0
                for d in NNd[i]:
                    if d < self.NNrc:
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
            NNd, NNi = self.query_atom_tree(k=self.kNN, rc=self.NNrc)
            for j, atom in enumerate(self):
                atom.NN = NeighborAtoms()
                atom.bonds = Bonds()
                for k, d in enumerate(NNd[j]):
                    if d < self.NNrc:
                        atom.NN.append(self[NNi[j][k]])
                        atom.bonds.append(Bond(atom, atom.NN[-1]))
        except ValueError:
            pass

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
    def NNrc(self):
        """Only return neighbors within this distance when querying the
        kd-tree."""
        return self._NNrc

    @NNrc.setter
    def NNrc(self, value):
        """Set the cutoff distance to check for neighest neighbors."""
        if not (isinstance(value, numbers.Number) and value >= 0):
            raise TypeError('Expected a real number greater >= 0')
        self._NNrc = value

    @property
    def bonds(self):
        self._update_nearest_neighbors()
        #return np.asarray([atom.bonds for atom in self])
        bonds = Bonds()
        [bonds.extend(atom.bonds) for atom in self]
        return bonds

    @property
    def pyramidalization_angles(self):
        self._update_pyramidalization_angles()
        angles = []
        [angles.append(atom.pyramidalization_angle) for atom in self if
         atom.poav is not None]
        return np.asarray(angles)

    def _update_pyramidalization_angles(self):
        self._update_nearest_neighbors()
        for atom in self:
            if atom.bonds.Nbonds == 3:
                b1, b2, b3 = atom.bonds
                v21 = Vector(b2.vector - b1.vector, p0=b1.vector.p)
                v31 = Vector(b3.vector - b1.vector, p0=b1.vector.p)
                #poav = vec.cross(b2.vector - b1.vector, b3.vector - b1.vector)
                poav = vec.cross(v21, v31)
                atom.poav = poav.unit_vector
                atom.sigma_bond_angle = vec.angle(atom.poav, b1.vector)
                #print('angle(poav, b1): {}'.format(np.degrees(vec.angle(
                #    atom.poav, b1.vector))))
                #print('angle(poav, b2): {}'.format(np.degrees(vec.angle(
                #    atom.poav, b2.vector))))
                #print('angle(poav, b3): {}'.format(np.degrees(vec.angle(
                #    atom.poav, b3.vector))))
                #print(np.degrees(atom.sigma_bond_angle))
                if atom.sigma_bond_angle < np.pi / 2:
                    atom.sigma_bond_angle = np.pi - atom.sigma_bond_angle
                    atom.poav = -atom.poav
                    #print('angle(poav, b1): {}'.format(np.degrees(vec.angle(
                    #    atom.poav, b1.vector))))
                    #print('angle(poav, b2): {}'.format(np.degrees(vec.angle(
                    #    atom.poav, b2.vector))))
                    #print('angle(poav, b3): {}'.format(np.degrees(vec.angle(
                    #    atom.poav, b3.vector))))

                atom.pyramidalization_angle = atom.sigma_bond_angle - np.pi / 2
                #print('atom.pyramidalization_angle: {}'.format(np.degrees(
                #    atom.pyramidalization_angle)))

    @property
    def poav_misalignment_angles(self):
        self._update_poav_misalignment_angles()
        angles = []
        [angles.extend(angle) for angle in
         [atom.poav_misalignment_angles for atom in self
          if len(atom.poav_misalignment_angles) > 0]]
        return np.asarray(angles)

    def _update_poav_misalignment_angles(self):
        self._update_pyramidalization_angles()
        for atom in self:
            if atom.bonds.Nbonds == 3:
                atom.poav_misalignment_angles = []
                for i, NN in enumerate(atom.NN):
                    bond = atom.bonds[i]
                    if NN.poav is not None:
                        nvec = vec.cross(bond.vector, atom.poav)
                        misalignment_angle = \
                            np.pi / 2 - vec.angle(NN.poav, nvec)
                        atom.poav_misalignment_angles.append(
                            misalignment_angle)
                            #vec.angle(atom.poav, NN.poav))
