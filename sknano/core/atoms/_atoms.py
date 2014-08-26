# -*- coding: utf-8 -*-
"""
==============================================================================
Base class for structure data atoms (:mod:`sknano.core.atoms._atoms`)
==============================================================================

.. currentmodule:: sknano.core.atoms._atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import copy
import math

from collections import OrderedDict, MutableSequence

import numpy as np

from sknano.core import xyz
from sknano.core.math import Vector, transformation_matrix

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
        self._data = []

        self._coords = []
        self._masses = []
        self._symbols = []

        self._property_lists = \
            {'m': self._masses, 'r': self._coords, 'symbol': self._symbols}

        if atoms is not None:
            try:
                if copylist and not deepcopy:
                    self._data.extend(atoms.data[:])
                elif deepcopy:
                    self._data.extend(copy.deepcopy(atoms.data))
                else:
                    self._data.extend(atoms.data)
            except AttributeError:
                if isinstance(atoms, list):
                    if copylist and not deepcopy:
                        self._data = atoms[:]
                    elif deepcopy:
                        self._data = copy.deepcopy(atoms)
                    else:
                        self._data = atoms

                    if update_property_lists:
                        for atom in self._data:
                            for p, plist in self._property_lists.iteritems():
                                plist.append(getattr(atom, p))
                else:
                    raise TypeError('`atoms={!r}` '.format(atoms) +
                                    'is not a valid `Atoms` constructor '
                                    'argument.\n atoms must be `None`, a list '
                                    'of `Atom` objects, or an `Atoms` object '
                                    'instance.')

    def __repr__(self):
        """Return `repr` string of `Atoms`."""
        return "Atoms(atoms={!r})".format(self._data)

    def __delitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Delete list element `self[index]` and delete all elements
        from atom properties lists `self.masses[index]`,
        `self.charges[index]`, and `self.coords[index]`

        Parameters
        ----------
        index : int
            index of target list element

        """
        del self._data[index]
        for plist in self._property_lists.itervalues():
            del plist[index]

    def __getitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Get `Atom` object instance at list element `self[index]`

        Parameters
        ----------
        index : int
            index of target list element

        Returns
        -------
        `Atom`
            `Atom` instance at target `self[index]`

        """
        return self._data[index]

    def __len__(self):
        """Concrete implementation of @abstractmethod.

        Returns
        -------
        int
            length of `self` list.

        """
        return len(self._data)

    def __setitem__(self, index, atom):
        """Concrete implementation of @abstractmethod.

        set target list element `self[index] = atom`

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
        self._data[index] = atom
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
        self._data.insert(index, atom)
        for p, plist in self._property_lists.iteritems():
            plist.insert(index, getattr(atom, p))

    @property
    def data(self):
        """Return the list of `Atom` objects"""
        return self._data

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
        return Vector(np.sum(MxR, axis=0) / np.sum(masses))

    @property
    def M(self):
        """Total mass of `Atoms`."""
        return math.fsum(self.masses)

    @property
    def Natoms(self):
        """Number of atoms in `Atoms`."""
        return len(self._data)

    @property
    def coords(self):
        """Return array of `Atom` coordinates."""
        coords = [atom.r for atom in self._data]
        self._coords = coords[:]
        return np.asarray(self._coords)

    @property
    def masses(self):
        """Return the list of `Atom` masses."""
        #self._masses_array = np.asarray(self._masses)
        masses = [atom.m for atom in self._data]
        self._masses = masses[:]
        return self._masses

    @property
    def symbols(self):
        """Return array of `Atom` symbols."""
        symbols = [atom.symbol for atom in self._data]
        self._symbols = symbols[:]
        return np.asarray(self._symbols)

    def center_CM(self):
        """Center atoms on CM coordinates."""
        dr = -self.CM
        self.translate(dr)

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
        for atom in self._data[:]:
            r = atom.r.tolist()
            for i, ri in enumerate(r[:]):
                if i in r_indices:
                    if abs_limit is not None:
                        if abs(ri) <= abs_limit:
                            continue
                        else:
                            del self._data[self._data.index(atom)]
                            break
                    else:
                        if min_limit is not None:
                            if ri >= min_limit:
                                continue
                            else:
                                del self._data[self._data.index(atom)]
                                break
                        if max_limit is not None:
                            if ri <= max_limit:
                                continue
                            else:
                                del self._data[self._data.index(atom)]
                                break

    def get_atoms(self, asarray=False):
        """Return list of `Atoms`.

        Parameters
        ----------
        asarray : bool, optional

        Returns
        -------
        sequence or ndarray

        """
        if asarray:
            return np.asarray(self._data)
        else:
            return self._data

    def get_coords(self, asdict=False):
        """Return atom coords.

        Parameters
        ----------
        asdict : bool, optional

        Returns
        -------
        coords : :py:class:`python:~collections.OrderedDict` or ndarray

        """
        coords = self.coords
        if asdict:
            return OrderedDict(zip(xyz, self.coords))
        else:
            return coords

    def rezero_coords(self, epsilon=1.0e-10):
        """Alias for :meth:`Atoms.rezero`."""
        self.rezero(epsilon=epsilon)

    def rezero_xyz(self, epsilon=1.0e-10):
        self.rezero(epsilon=epsilon)

    def rezero(self, epsilon=1.0e-10):
        """Set really really small coordinates to zero.

        Set all coordinates with absolute value less than
        epsilon to zero.

        Parameters
        ----------
        epsilon : float
            smallest allowed absolute value of any :math:`x,y,z` component.

        """
        [atom.rezero(epsilon=epsilon) for atom in self._data]

    def rotate(self, angle=None, rot_axis=None, anchor_point=None,
               deg2rad=False, transform_matrix=None):
        """Rotate atom coordinates about arbitrary axis.

        Parameters
        ----------
        angle : float

        """
        if transform_matrix is None:
            transform_matrix = \
                transformation_matrix(angle, rot_axis=rot_axis,
                                      anchor_point=anchor_point,
                                      deg2rad=deg2rad)
        [atom.rotate(transform_matrix=transform_matrix) for atom in self._data]

    def translate(self, t):
        """Translate atom coordinates.

        Parameters
        ----------
        t : array_like
            3-elment array of :math:`x,y,z` components of translation vector
        """
        [atom.translate(t) for atom in self._data]
