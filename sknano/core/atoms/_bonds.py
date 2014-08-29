# -*- coding: utf-8 -*-
"""
===============================================================================
Atom class with extended feature set (:mod:`sknano.core.atoms._bonds`)
===============================================================================

An "eXtended" `Atom` class for structure analysis.

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence
from operator import attrgetter

import numpy as np

from sknano.core.math import Vector, transformation_matrix

__all__ = ['Bonds', 'BondMixin']


class Bonds(MutableSequence):
    """Base class for collection of atom `Bonds`.

    Parameters
    ----------
    atoms : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        self._data = []
        self._coords = []
        self._masses = []
        self._symbols = []

        if atoms is not None:
            try:
                if copylist and not deepcopy:
                    self._data.extend(atoms[:])
                elif deepcopy:
                    self._data.extend(copy.deepcopy(atoms))
                else:
                    self._data.extend(atoms)
            except AttributeError:
                raise TypeError('`atoms={!r}` '.format(atoms) +
                                'is not a valid `Bonds` constructor '
                                'argument.\n atoms must be `None`, a list '
                                'of `Atom` objects, or an `Atoms` object '
                                'instance.')

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return `repr` string of `Bonds`."""
        return "Bonds(atoms={!r})".format(self._data)

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

    def __len__(self):
        """Concrete implementation of @abstractmethod.

        Returns
        -------
        int
            length of `self` list.

        """
        return len(self._data)

    def insert(self, index, atom):
        """Concrete implementation of @abstractmethod.

        insert `Atom` instance at target list `index`

        Parameters
        ----------
        index : int
            index of target list element
        atom : `Atom`
            `Atom` object instance to set target list element to

        """
        self._data.insert(index, atom)

    @property
    def data(self):
        """Return the list of `Atom` objects"""
        return self._data

    def sort(self, key=None, reverse=False):

        if key is None:
            self._data.sort(key=attrgetter('element', 'Z', 'z'),
                            reverse=reverse)
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def Nbonds(self):
        """Number of bonds in `Bonds`."""
        return len(self._data)

    @property
    def bonds(self):
        """Return array of `Atom` coordinates."""
        coords = [atom.r for atom in self._data]
        self._coords = coords[:]
        return np.asarray(self._coords)

    @property
    def lengths(self):
        """Return the list of `Atom` masses."""
        masses = [atom.m for atom in self._data]
        self._masses = masses[:]
        return self._masses

    @property
    def symbols(self):
        """Return array of `Atom` symbols."""
        symbols = [atom.symbol for atom in self._data]
        self._symbols = symbols[:]
        return np.asarray(self._symbols)

    def clip_bounds(self, region, center_before_clipping=False):
        """Remove atoms outside the given limits along given dimension.

        Parameters
        ----------
        region : :class:`~sknano.utils.geometric_shapes.`GeometricRegion`

        """
        CM0 = None
        if center_before_clipping:
            CM0 = self.CM
            self.translate(-CM0)

        bonds = self.get_bonds(asarray=True)
        limits = region.limits
        self._data = \
            bonds[np.logical_and(
                np.logical_and(
                    self.x <= limits['x']['max'],
                    np.logical_and(
                        self.y <= limits['y']['max'],
                        self.z <= limits['z']['max'])),
                np.logical_and(
                    self.x >= limits['x']['min'],
                    np.logical_and(
                        self.y >= limits['y']['min'],
                        self.z >= limits['z']['min'])))].tolist()

        #for dim, limits in region.limits.iteritems():
        #    bonds = bonds[np.where(getattr(self, dim) <= limits['max'])]
        #    bonds = bonds[np.where(getattr(self, dim) >= limits['min'])]
        #    self = bonds.tolist()

        if CM0 is not None:
            self.translate(CM0)

    def get_bonds(self, asarray=False):
        """Return list of `Bonds`.

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


class BondMixin(object):

    def bond_length(self):
        return self._bond_length
