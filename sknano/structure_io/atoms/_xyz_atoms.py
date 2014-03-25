# -*- coding: utf-8 -*-
"""
==================================================================
Class for XYZ atoms (:mod:`sknano.structure_io.atoms._xyz_atoms`)
==================================================================

.. currentmodule:: sknano.structure_io.atoms._xyz_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import OrderedDict

import numpy as np

try:
    from scipy.spatial import KDTree
    has_kdtree = True
except ImportError:
    print('Install scipy version >= 0.13.0 to allow '
          'nearest-neighbor queries between atoms.')
    has_kdtree = False

from ...tools import xyz_axes
from ._atoms import Atoms

__all__ = ['XYZAtoms']


class XYZAtoms(Atoms):
    """`Atoms` class for creating collection of `XYZAtom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `XYZAtoms`}, optional
        if not `None`, then a list of `XYZAtom` instance objects or an
        existing `XYZAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list
    use_kdtree : bool, optional
        use :py:class:`~scipy:scipy.spatial.KDTree` to perform
        nearest-neighbor analysis.

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False,
                 use_kdtree=True):

        super(XYZAtoms, self).__init__(
            atoms=atoms, copylist=copylist, deepcopy=deepcopy,
            update_property_lists=False)

        #self._property_lists['atomID'] = self._atom_ids = []
        #self._property_lists['CN'] = self._coordination_numbers = []
        #self._property_lists['NN'] = self._nearest_neighbors = []

        self._atom_tree = None
        if use_kdtree and has_kdtree is False:
            use_kdtree = False
        self._use_kdtree = use_kdtree

        self._NN_number = 6
        self._NN_cutoff = np.inf

        if atoms is not None:
            for atom in self._atoms:
                for p, plist in self._property_lists.iteritems():
                    plist.append(getattr(atom, p))

            if use_kdtree:
                self._atom_tree = self.atom_tree

    @property
    def atom_tree(self):
        """Return the :py:class:`~scipy:scipy.spatial.KDTree` of coords."""
        if self._use_kdtree and len(self._atoms) != 0:
            self._atom_tree = KDTree(self.coords)
        return self._atom_tree

    @property
    def coordination_numbers(self):
        """Return array of `XYZAtom` coordination numbers."""
        self.update_coordination_numbers()
        coordination_numbers = []
        for atom in self._atoms:
            coordination_numbers.append(atom.CN)
        self._coordination_numbers = coordination_numbers[:]
        return np.asarray(self._coordination_numbers)

    def update_coordination_numbers(self):
        """Update `XYZAtom` coordination numbers."""
        if self._use_kdtree:
            NN_d, NN_i = \
                self.query_atom_tree(n=self.NN_number,
                                     cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self._atoms):
                CN = 0
                for d in NN_d[i]:
                    if d < self.NN_cutoff:
                        CN += 1
                atom.CN = CN

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `XYZAtom`."""
        self.update_nearest_neighbors()
        nearest_neighbors = []
        for atom in self._atoms:
            nearest_neighbors.append(atom.NN)
        self._nearest_neighbors = nearest_neighbors[:]
        return np.asarray(self._nearest_neighbors)

    def update_nearest_neighbors(self):
        """Update `XYZAtom` nearest-neighbors."""
        if self._use_kdtree:
            NN_d, NN_i = self.query_atom_tree(n=self.NN_number,
                                              cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self._atoms):
                NN_atoms = []
                for j, d in enumerate(NN_d[i]):
                    if d < self.NN_cutoff:
                        NN_atoms.append(self._atoms[NN_i[i][j]])
                atom.NN = NN_atoms

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
        NN_d = NN_i = None
        if atom_tree is not None:
            d, i = atom_tree.query(self.coords, k=n+1, eps=eps, p=p,
                                   distance_upper_bound=cutoff_radius)
            NN_d, NN_i = d[:, 1:], i[:, 1:]
        return NN_d, NN_i

    def query_coordination_numbers(self, n=6, rc=np.inf):
        """Query and update `XYZAtom` coordination numbers.

        Parameters
        ----------
        n : int, optional
        rc : nonnegative float, optional

        """
        if self._use_kdtree:
            NN_d, NN_i = self.query_atom_tree(n=n, cutoff_radius=rc)
            for i, atom in enumerate(self._atoms):
                CN = 0
                for d in NN_d[i]:
                    if d < rc:
                        CN += 1
                atom.CN = CN
        coordination_numbers = []
        for atom in self._atoms:
            coordination_numbers.append(atom.CN)
        self._coordination_numbers = coordination_numbers[:]
        return np.asarray(self._coordination_numbers)

    def query_nearest_neighbors(self, n=6, rc=np.inf):
        """Query and update `XYZAtom` nearest neighbors.

        Parameters
        ----------
        n : int, optional
        rc : nonnegative float, optional

        """
        if self._use_kdtree:
            NN_d, NN_i = self.query_atom_tree(n=n, cutoff_radius=rc)
            for i, atom in enumerate(self._atoms):
                NN_atoms = []
                for j, d in enumerate(NN_d[i]):
                    if d < rc:
                        NN_atoms.append(self._atoms[NN_i[i][j]])
                atom.NN = NN_atoms
        nearest_neighbors = []
        for atom in self._atoms:
            nearest_neighbors.append(atom.NN)
        self._nearest_neighbors = nearest_neighbors[:]
        return np.asarray(self._nearest_neighbors)

    @property
    def NN_number(self):
        return self._NN_number

    @NN_number.setter
    def NN_number(self, value):
        self._NN_number = int(value)

    @property
    def NN_cutoff(self):
        return self._NN_cutoff

    @NN_cutoff.setter
    def NN_cutoff(self, value):
        self._NN_cutoff = value

    def assign_unique_ids(self, starting_id=1):
        """Assign unique ID to each `XYZAtom` in `Atoms`."""
        for i, atom in enumerate(self._atoms, start=starting_id):
            atom.atomID = i

    def filter_atoms(self, filter_atom_ids, invert=False):
        """Filter `XYZAtoms`.

        Parameters
        ----------
        filter_atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `XYZAtoms`

        """
        #if invert:
        #    filtered_atoms = \
        #        XYZAtoms([atom for atom in self.atoms if
        #               atom.atomID not in filter_atom_ids])

        #else:
        #    filtered_atoms = \
        #        XYZAtoms([atom for atom in self.atoms if
        #               atom.atomID in filter_atom_ids])

        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        filtered_atoms = \
            XYZAtoms(np.asarray(self._atoms)[filter_indices].tolist())
        return filtered_atoms

    def get_atom(self, atomID=None, index=None):
        try:
            return self._atoms[atomID - 1]
        except (TypeError, IndexError):
            try:
                return self._atoms[index]
            except (TypeError, IndexError):
                return None

    def get_filtered_atom_ids(self, filter_atom_ids, invert=False):
        """Return atom ids filtered by list of `filter_atom_ids`.

        Parameters
        ----------
        invert : bool, optional

        Returns
        -------
        ndarray

        """
        atom_ids = self.atom_ids
        filter_indices = \
            np.in1d(atom_ids, filter_atom_ids, invert=invert).nonzero()
        return atom_ids[filter_indices]

    def get_filtered_coords(self, filter_atom_ids, components=None,
                            as_dict=False, invert=False):
        """Return filtered coordinates filtered by filter_atom_ids.

        Parameters
        ----------
        filter_atom_ids : array_like
        components : {None, sequence}, optional
        as_dict : bool, optional
        invert : bool, optional

        Returns
        -------
        filtered_coords : :py:class:`python:~collections.OrderedDict` or
        ndarray

        """
        atom_ids = self.atom_ids
        coords = self.coords
        filter_indices = \
            np.in1d(atom_ids, filter_atom_ids, invert=invert).nonzero()
        filtered_coords = coords[filter_indices]

        if components is None or components == 'r':
            components = ('x', 'y', 'z')
        elif isinstance(components, str):
            components = (components,)

        if as_dict:
            return OrderedDict(zip(
                components, [filtered_coords[:, xyz_axes.index(component)]
                             for component in components]))
        else:
            filtered_coords
