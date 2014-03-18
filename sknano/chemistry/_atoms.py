# -*- coding: utf-8 -*-
"""
===================================================================
Abstract representation of Atoms (:mod:`sknano.chemistry._atoms`)
===================================================================

.. currentmodule:: sknano.chemistry._atoms

.. todo::

   Apply more consistent use of upper and lower case vars for `Atoms`
   attributes.

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
import math

from collections import OrderedDict, MutableSequence

import numpy as np

try:
    from scipy.spatial import KDTree
    has_kdtree = True
except ImportError:
    print('Install scipy version >= 0.13.0 to allow '
          'nearest-neighbor queries between `Atoms`.')
    has_kdtree = False

from ..tools import xyz_axes

__all__ = ['Atoms']


class Atoms(MutableSequence):

    """Class for creating abstract represention of a collection of `Atoms`.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.
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

        self._atoms = []
        self._atom_ids = []
        self._molecule_ids = []
        self._charges = []
        self._coords = []
        self._masses = []
        self._velocities = []
        self._forces = []
        self._symbols = []
        self._coordination_numbers = []
        self._nearest_neighbors = []

        self._property_lists = \
            {'m': self._masses, 'q': self._charges, 'r': self._coords,
             'v': self._velocities, 'atomID': self._atom_ids,
             'moleculeID': self._molecule_ids, 'symbol': self._symbols,
             'CN': self._coordination_numbers, 'NN': self._nearest_neighbors}

        self._atomtypes = {}

        self._atom_tree = None

        if use_kdtree and has_kdtree is False:
            use_kdtree = False
        self._use_kdtree = use_kdtree
        self._NN_number = 6
        self._NN_cutoff = np.inf

        if atoms is not None:
            try:
                if copylist and not deepcopy:
                    self._atoms.extend(atoms.atoms[:])
                elif deepcopy:
                    self._atoms.extend(copy.deepcopy(atoms.atoms))
                else:
                    self._atoms.extend(atoms.atoms)
            except AttributeError:
                if isinstance(atoms, list):
                    if copylist and not deepcopy:
                        self._atoms = atoms[:]
                    elif deepcopy:
                        self._atoms = copy.deepcopy(atoms)
                    else:
                        self._atoms = atoms

                    for atom in self._atoms:
                        for p, plist in self._property_lists.iteritems():
                            plist.append(getattr(atom, p))
                else:
                    raise TypeError('`Atoms(atoms={!r})` '.format(atoms) +
                                    'is not a valid `Atoms` constructor '
                                    'argument.\n atoms must be `None`, a list '
                                    'of `Atom` objects, or an `Atoms` object '
                                    'instance.')

            if use_kdtree:
                self._atom_tree = self.atom_tree

    def __str__(self):
        """Return string representation of `Atoms`."""
        atoms_str = ''
        for atom in self._atoms:
            atoms_str += str(atom)
        return atoms_str

    @property
    def atoms(self):
        """Return the list of `Atom` objects"""
        return self._atoms

    @property
    def atomlist(self):
        """Return the list of `Atom` objects

        .. deprecated:: 0.1.0
           Use :py:attr:`atoms` attribute.

        """
        return self.atoms

    @property
    def atomtypes(self):
        """Return the atom type dictionary."""
        for atom in self._atoms:
            if atom.atomtype not in self._atomtypes:
                self._atomtypes[atom.atomtype] = {}
                self._atomtypes[atom.atomtype]['mass'] = atom.m
                self._atomtypes[atom.atomtype]['q'] = atom.q
        return self._atomtypes

    @property
    def atom_ids(self):
        """Return array of `Atom` IDs."""
        atom_ids = []
        for atom in self._atoms:
            atom_ids.append(atom.atomID)
        self._atom_ids = atom_ids[:]
        return np.asarray(self._atom_ids)

    @property
    def atom_tree(self):
        """Return the :py:class:`~scipy:scipy.spatial.KDTree` of coords."""
        if self._use_kdtree and len(self._atoms) != 0:
            self._atom_tree = KDTree(self.coords)
        return self._atom_tree

    @property
    def CM(self):
        """Center-of-Mass coordinates of `Atoms`.

        Returns
        -------
        ndarray
            3-element ndarray specifying center-of-mass coordinates of `Atoms`.

        """
        masses = np.asarray([self.masses])
        coords = np.asarray(self.coords)
        MxR = masses.T * coords
        return np.sum(MxR, axis=0) / np.sum(masses)

    @property
    def M(self):
        """Total mass of `Atoms`."""
        return math.fsum(self.masses)

    @property
    def Natoms(self):
        """Number of atoms in `Atoms`."""
        return len(self._atoms)

    @property
    def NN_cutoff(self):
        return self._NN_cutoff

    @NN_cutoff.setter
    def NN_cutoff(self, value):
        self._NN_cutoff = value

    @property
    def NN_number(self):
        return self._NN_number

    @NN_number.setter
    def NN_number(self, value):
        self._NN_number = int(value)

    @property
    def Ntypes(self):
        """Number of atom types in `Atoms`."""
        return len(self.atomtypes.keys())

    @property
    def charges(self):
        """Return list of `Atom` charges."""
        #self._charges_array = np.asarray(self._charges)
        charges = []
        for atom in self._atoms:
            charges.append(atom.q)
        self._charges = charges[:]
        return np.asarray(self._charges)

    @property
    def coords(self):
        """Return array of `Atom` coordinates."""
        coords = []
        for atom in self._atoms:
            coords.append(atom.r)
        self._coords = coords[:]
        return np.asarray(self._coords)

    @property
    def coordination_numbers(self):
        """Return array of `Atom` coordination numbers."""
        self.update_coordination_numbers()
        coordination_numbers = []
        for atom in self._atoms:
            coordination_numbers.append(atom.CN)
        self._coordination_numbers = coordination_numbers[:]
        return np.asarray(self._coordination_numbers)

    def update_coordination_numbers(self):
        """Update `Atom` coordination numbers."""
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
        """Return array of nearest-neighbor atoms for each `Atom`."""
        self.update_nearest_neighbors()
        nearest_neighbors = []
        for atom in self._atoms:
            nearest_neighbors.append(atom.NN)
        self._nearest_neighbors = nearest_neighbors[:]
        return np.asarray(self._nearest_neighbors)

    def update_nearest_neighbors(self):
        """Update `Atom` nearest-neighbors."""
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
        """Query and update `Atom` coordination numbers.

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
        """Query and update `Atom` nearest neighbors.

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
    def masses(self):
        """Return the list of `Atom` masses."""
        #self._masses_array = np.asarray(self._masses)
        masses = []
        for atom in self._atoms:
            masses.append(atom.m)
        self._masses = masses[:]
        return self._masses

    @property
    def m(self):
        """Total mass of `Atoms`.

        .. deprecated:: 0.1.0
           Use :py:attr:`M` attribute instead.

        """
        return self.M

    @property
    def q(self):
        """Return the total net charge of `Atoms`."""
        return np.asarray(self.charges).sum()

    @property
    def r(self):
        """Alias for :py:attr:`~sknano.chemistry.Atoms.CM`"""
        return self.CM

    @property
    def symbols(self):
        """Return array of `Atom` symbols."""
        symbols = []
        for atom in self._atoms:
            symbols.append(atom.symbol)
        self._symbols = symbols[:]
        return np.asarray(self._symbols)

    @property
    def velocities(self):
        """Return array of `Atom` velocities."""
        velocities = []
        for atom in self._atoms:
            velocities.append(atom.v)
        self._velocities = velocities[:]
        return np.asarray(self._velocities)

    def add_atomtype(self, atom):
        """Add atom type to :py:attr:`~sknano.chemistry.Atoms.atomtypes` dict.

        Parameters
        ----------
        atom : `Atom`
            an intance of an `Atom` object

        """
        if atom.atomtype not in self._atomtypes:
            self._atomtypes[atom.atomtype] = {}
            self._atomtypes[atom.atomtype]['mass'] = atom.m
            self._atomtypes[atom.atomtype]['q'] = atom.q

    def add_atomtypes(self, atomtypes=[]):
        """Add atom type in atomtypes list to atomtype dictionary.
        Parameters
        ----------
        atomtypes : sequence
            a list of `Atom` object instances

        """
        for atom in atomtypes:
            self.add_atomtype(atom)

    def assign_unique_ids(self, starting_id=1):
        """Assign unique ID to each `Atom` in `Atoms`."""
        for i, atom in enumerate(self._atoms, start=starting_id):
            atom.atomID = i

    def center_CM(self, r_indices=[0, 1, 2]):
        """Center atoms on CM coordinates

        Parameters
        ----------
        r_indices : sequence, optional
            list of component indices of CM position vector to translate
            and center on the origin.

        """
        dr = -self.CM
        self.translate(dr, r_indices=r_indices)

    def clip_bounds(self, min_limit=None, max_limit=None,
                    abs_limit=None, r_indices=[0, 1, 2]):
        """Remove atoms outside the given limits along given dimension.

        Parameters
        ----------
        min_limit : {None, float}, optional
        max_limit : {None, float}, optional
        abs_limit : {None, float}, optional
        r_indices : sequence, optional

        """
        for atom in self._atoms[:]:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if i in r_indices:
                    if abs_limit is not None:
                        if abs(ri) <= abs_limit:
                            continue
                        else:
                            del self._atoms[self._atoms.index(atom)]
                            break
                    else:
                        if min_limit is not None:
                            if ri >= min_limit:
                                continue
                            else:
                                del self._atoms[self._atoms.index(atom)]
                                break
                        if max_limit is not None:
                            if ri <= max_limit:
                                continue
                            else:
                                del self._atoms[self._atoms.index(atom)]
                                break

    def filter_atoms(self, filter_atom_ids, invert=False):
        """Filter `Atoms`.

        Parameters
        ----------
        filter_atom_ids : array_like
        invert : bool, optional

        Returns
        -------
        filtered_atoms : `Atoms`

        """
        #if invert:
        #    filtered_atoms = \
        #        Atoms([atom for atom in self.atoms if
        #               atom.atomID not in filter_atom_ids])

        #else:
        #    filtered_atoms = \
        #        Atoms([atom for atom in self.atoms if
        #               atom.atomID in filter_atom_ids])

        filter_indices = \
            np.in1d(self.atom_ids, filter_atom_ids, invert=invert).nonzero()
        filtered_atoms = \
            Atoms(np.asarray(self._atoms)[filter_indices].tolist())
        return filtered_atoms

    def fix_minus_zero_coords(self, epsilon=1.0e-10):
        """Set really really small negative coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon zero so we don't end up with -0.00000
        coordinates in structure data output.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any component of an
            atom's position

        """
        for atom in self._atoms:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if ri < 0 and abs(ri) < epsilon:
                    r[i] = 0.0
            atom.x, atom.y, atom.z = r

    def getatomsattr(self, asarray=False, as_array=False):
        pass

    def get_atom(self, atomID=None, index=None):
        try:
            return self._atoms[atomID - 1]
        except (TypeError, IndexError):
            try:
                return self._atoms[index]
            except (TypeError, IndexError):
                return None

    def get_atoms(self, asarray=False, as_array=False):
        """Return list of `Atoms`.

        Parameters
        ----------
        asarray, as_array : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray or as_array:
            return np.asarray(self._atoms)
        else:
            return self._atoms

    def get_coords(self, components=None, as_dict=False):
        """Return atom coords.

        Parameters
        ----------
        components : {None, sequence}, optional
        as_dict : bool, optional

        Returns
        -------
        coords : :py:class:`python:~collections.OrderedDict` or ndarray

        """
        coords = self.coords
        if as_dict:
            if components is None or components == 'r':
                components = ('x', 'y', 'z')
            elif isinstance(components, str):
                components = (components,)

            return OrderedDict(zip(
                components, [coords[:, xyz_axes.index(component)] for
                             component in components]))
        else:
            return coords

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

    def rezero_coords(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        for atom in self._atoms:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if abs(ri) < epsilon:
                    r[i] = 0.0
            atom.x, atom.y, atom.z = r

    def rezero_xyz(self, epsilon=1.0e-10):
        return self.rezero_coords(epsilon=epsilon)

    def rotate(self, R_matrix):
        """Rotate atom coordinates using rotation matrix `R_matrix`.

        Parameters
        ----------
        R_matrix : array_like
            3x3 array representation of *3D rotation matrix*.

        """
        for atom in self._atoms:
            atom.r = np.dot(R_matrix, atom.r.T).T

    def translate(self, dr, r_indices=[0, 1, 2]):
        """Translate atom coordinates.

        Parameters
        ----------
        dr : array_like
            array representation of displacement vector to translate
            components of atoms position vector through
        r_indices : sequence, optional
            list of component indices of position vector that
            `dr` vector acts upon

        """
        for atom in self._atoms:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if i in r_indices:
                    r[i] += dr[i]
            atom.x, atom.y, atom.z = r

    def __delitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Delete list element `self.atoms[index]` and delete all elements
        from atom properties lists `self.masses[index]`,
        `self.charges[index]`, and `self.coords[index]`

        Parameters
        ----------
        index : int
            index of target list element

        """
        del self._atoms[index]
        for plist in self._property_lists.itervalues():
            del plist[index]

    def __getitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Get `Atom` object instance at list element `self.atoms[index]`

        Parameters
        ----------
        index : int
            index of target list element

        Returns
        -------
        `Atom`
            `Atom` instance at target `self.atoms[index]`

        """
        return self._atoms[index]

    def __len__(self):
        """Concrete implementation of @abstractmethod.

        Returns
        -------
        int
            length of `self.atoms` list.

        """
        return len(self._atoms)

    def __setitem__(self, index, atom):
        """Concrete implementation of @abstractmethod.

        set target list element `self.atoms[index] = atom`

        Also set element of all atom properties lists (`self.masses[index]`,
        `self.charges[index]`, and `self.coords[index]`) to atom instance
        properties (`atom.m`, `atom.q`, `atom.r`), respectively.

        Parameters
        ----------
        index : int
            index of target list element
        atom : `Atom`
            `Atom` instance object to set target list element to.

        """
        self._atoms[index] = atom
        for p, plist in self._property_lists.iteritems():
            plist[index] = getattr(atom, p)

    def insert(self, index, atom):
        """Concrete implementation of @abstractmethod.

        insert `Atom` instance at target list `index`

        Also insert `Atom` instance properties at the given target list index
        for all `Atom` properties in `self._property_lists.keys()`
        into their respective target lists of `Atoms` properties
        `self._property_lists.values()`.

        Parameters
        ----------
        index : int
            index of target list element
        atom : `Atom`
            `Atom` object instance to set target list element to

        """
        self._atoms.insert(index, atom)
        for p, plist in self._property_lists.iteritems():
            plist.insert(index, getattr(atom, p))
