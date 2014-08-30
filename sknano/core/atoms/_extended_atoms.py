# -*- coding: utf-8 -*-
"""
===============================================================================
Extended Atoms class feature set (:mod:`sknano.core.atoms._extended_atoms`)
===============================================================================

An "eXtended" `Atoms` class for structure analysis.

.. currentmodule:: sknano.core.atoms._extended_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import OrderedDict
from functools import total_ordering
from operator import attrgetter
import numbers

import numpy as np

try:
    from scipy.spatial import KDTree
    has_kdtree = True
except ImportError:
    print('Install scipy version >= 0.13.0 to allow '
          'nearest-neighbor queries between atoms.')
    has_kdtree = False

from sknano.core import xyz
from ._atoms import Atoms

__all__ = ['XAtoms']


@total_ordering
class XAtoms(Atoms):
    """An eXtended `Atoms` class for structure analysis.

    Parameters
    ----------
    atoms : {None, sequence, `XAtoms`}, optional
        if not `None`, then a list of `XAtom` instance objects or an
        existing `XAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list
    use_kdtree : bool, optional
        use :py:class:`~scipy:scipy.spatial.KDTree` to perform
        nearest-neighbor analysis.

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False,
                 use_kdtree=True, update_tree=True):

        super(XAtoms, self).__init__(atoms=atoms,
                                     copylist=copylist,
                                     deepcopy=deepcopy)

        self._atomtypes = {}

        self._atom_tree = None
        if use_kdtree and has_kdtree is False:
            use_kdtree = False
        self._use_kdtree = use_kdtree

        self._kNN = 3
        self._NN_cutoff = 0

    def __eq__(self, other):
        return self[:] == other[:]

    def __lt__(self, other):
        return self[:] < other[:]

    def sort(self, key=None, reverse=False):
        if key is None:
            self._data.sort(key=attrgetter('element', 'Z', 'atomtype',
                                           'moleculeID', 'atomID'),
                            reverse=reverse)
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def atomtypes(self):
        """Return the atom type dictionary."""
        self._update_atomtypes()
        return self._atomtypes

    def _update_atomtypes(self):
        for atom in self:
            if atom.atomtype not in self._atomtypes:
                self._atomtypes[atom.atomtype] = {}
                self._atomtypes[atom.atomtype]['mass'] = atom.m
                self._atomtypes[atom.atomtype]['q'] = atom.q

    @property
    def atom_ids(self):
        """Return array of `XAtom` IDs."""
        return np.asarray([atom.atomID for atom in self])

    @property
    def atom_tree(self):
        """Return the :py:class:`~scipy:scipy.spatial.KDTree` of coords."""
        if self._use_kdtree and len(self) != 0:
            return KDTree(self.coords)

    @property
    def charges(self):
        """Return array of `XAtom` charges."""
        return np.asarray([atom.q for atom in self])

    @property
    def coordination_numbers(self):
        """Return array of `XAtom` coordination numbers."""
        self._update_coordination_numbers()
        return np.asarray([atom.CN for atom in self])

    def _update_coordination_numbers(self):
        """Update `XAtom` coordination numbers."""
        if self._use_kdtree:
            NN_d, NN_i = \
                self.query_atom_tree(n=self.kNN,
                                     cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self):
                for d in NN_d[i]:
                    if d < self.NN_cutoff:
                        atom.CN += 1

    @property
    def nearest_neighbors(self):
        """Return array of nearest-neighbor atoms for each `XAtom`."""
        self._update_nearest_neighbors()
        return np.asarray([atom.NN for atom in self])

    def _update_nearest_neighbors(self):
        """Update `XAtom` nearest-neighbors."""
        if self._use_kdtree:
            NN_d, NN_i = self.query_atom_tree(n=self.kNN,
                                              cutoff_radius=self.NN_cutoff)
            for i, atom in enumerate(self):
                for j, d in enumerate(NN_d[i]):
                    if d < self.NN_cutoff:
                        atom.NN.append(self[NN_i[i][j]])

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
    def Ntypes(self):
        """Number of atom types in `XAtoms`."""
        return len(self.atomtypes.keys())

    @property
    def q(self):
        """Return the total net charge of `XAtoms`."""
        return self.charges.sum()

    @property
    def velocities(self):
        """Return array of `XAtom` velocities."""
        return np.asarray([atom.v for atom in self])

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

    def add_atomtype(self, atom):
        """Add atom type to :attr:`~XAtoms.atomtypes`.

        Parameters
        ----------
        atom : :class:`~sknano.core.atoms.XAtom`
            A :class:`~sknano.core.atoms.XAtom` instance.

        """
        if atom.atomtype not in self._atomtypes:
            self._atomtypes[atom.atomtype] = {}
            self._atomtypes[atom.atomtype]['mass'] = atom.m
            self._atomtypes[atom.atomtype]['q'] = atom.q

    def add_atomtypes(self, atoms=None):
        """Add atomtype for each atom in atoms to atomtypes dictionary.
        Parameters
        ----------
        atoms : sequence
            a list of `XAtom` object instances

        """
        try:
            [self.add_atomtype(atom) for atom in atoms]
        except TypeError:
            print('Expected an iterable sequence of `XAtom` objects.')

    def assign_unique_ids(self, starting_id=1):
        """Assign unique ID to each `XAtom` in `XAtoms`."""
        for i, atom in enumerate(self, start=starting_id):
            atom.atomID = i

    def filter_atoms(self, filter_atom_ids, invert=False):
        """Filter `XAtoms`.

        Parameters
        ----------
        filter_atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `XAtoms`

        """
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        return XAtoms(np.asarray(self)[filter_indices].tolist())

    def get_atom(self, atomID=None, index=None):
        try:
            return self[atomID - 1]
        except (TypeError, IndexError):
            try:
                return self[index]
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
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        return self.atom_ids[filter_indices]

    def get_filtered_coords(self, filter_atom_ids, components=None,
                            asdict=False, invert=False):
        """Return filtered coordinates filtered by filter_atom_ids.

        Parameters
        ----------
        filter_atom_ids : array_like
        components : {None, sequence}, optional
        asdict : bool, optional
        invert : bool, optional

        Returns
        -------
        filtered_coords : :py:class:`python:~collections.OrderedDict` or
        ndarray

        """
        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        filtered_coords = self.coords[filter_indices]

        if components is None or components == 'r':
            components = ('x', 'y', 'z')
        elif isinstance(components, (str, unicode)):
            components = (components,)

        if asdict:
            return OrderedDict(zip(components,
                                   [filtered_coords[:, xyz.index(component)]
                                    for component in components]))
        else:
            filtered_coords

    def select(self, **kwargs):
        pass

    def select_within(self, volume):
        pass
