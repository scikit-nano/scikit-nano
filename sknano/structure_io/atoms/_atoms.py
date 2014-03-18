# -*- coding: utf-8 -*-
"""
==============================================================================
Base class for structure data atoms (:mod:`sknano.structure_io.atoms._atoms`)
==============================================================================

.. currentmodule:: sknano.structure_io.atoms._atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
import math

from collections import OrderedDict, MutableSequence

import numpy as np

from ...tools import xyz_axes

__all__ = ['Atoms']


class Atoms(MutableSequence):
    """Base class for collection of `Atom` objects.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects or an
        existing `Atoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False,
                 update_property_lists=True):
        self._atoms = []

        self._coords = []
        self._masses = []
        self._symbols = []

        self._property_lists = \
            {'m': self._masses, 'r': self._coords, 'symbol': self._symbols}

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

                    if update_property_lists:
                        for atom in self._atoms:
                            for p, plist in self._property_lists.iteritems():
                                plist.append(getattr(atom, p))
                else:
                    raise TypeError('`Atoms(atoms={!r})` '.format(atoms) +
                                    'is not a valid `Atoms` constructor '
                                    'argument.\n atoms must be `None`, a list '
                                    'of `Atom` objects, or an `Atoms` object '
                                    'instance.')

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
    def coords(self):
        """Return array of `Atom` coordinates."""
        coords = []
        for atom in self._atoms:
            coords.append(atom.r)
        self._coords = coords[:]
        return np.asarray(self._coords)

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
    def symbols(self):
        """Return array of `Atom` symbols."""
        symbols = []
        for atom in self._atoms:
            symbols.append(atom.symbol)
        self._symbols = symbols[:]
        return np.asarray(self._symbols)

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
