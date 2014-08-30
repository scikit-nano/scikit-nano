# -*- coding: utf-8 -*-
"""
===============================================================================
Classes for atom bonds (:mod:`sknano.core.atoms._bonds`)
===============================================================================

Classes for `Atom` bonds

.. currentmodule:: sknano.core.atoms._bonds

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence
from operator import attrgetter
import copy

#import numpy as np

from sknano.core.math import Vector, vector as vec

__all__ = ['Bonds', 'BondMixin']


class Bonds(MutableSequence):
    """Base class for collection of atom `Bonds`.

    Parameters
    ----------
    bonds : {None, sequence, `Atoms`}, optional
        if not `None`, then a list of `Atom` instance objects
    copylist : bool, optional
        perform shallow copy of bonds list
    deepcopy : bool, optional
        perform deepcopy of bonds list

    """

    def __init__(self, bonds=None, copylist=True, deepcopy=False):
        self._data = []
        self._vectors = []

        if bonds is not None:
            try:
                if copylist and not deepcopy:
                    self._data.extend(bonds[:])
                elif deepcopy:
                    self._data.extend(copy.deepcopy(bonds))
                else:
                    self._data.extend(bonds)
            except AttributeError:
                raise TypeError('`bonds={!r}` '.format(bonds) +
                                'is not a valid `Bonds` constructor '
                                'argument.\n bonds must be `None`, a list '
                                'of `Atom` objects, or an `Bonds` object '
                                'instance.')

    def __str__(self):
        return repr(self)

    def __repr__(self):
        """Return `repr` string of `Bonds`."""
        return "Bonds(bonds={!r})".format(self._data)

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
        """Return the list of `Bonds` objects"""
        return self._data

    def sort(self, key=None, reverse=False):

        if key is None:
            self._data.sort(key=attrgetter('length'), reverse=reverse)
        else:
            self._data.sort(key=key, reverse=reverse)

    @property
    def bond_vectors(self):
        return self._vectors

    @bond_vectors.setter
    def bonds_vectors(self, value):
        self._vectors = value


class BondMixin(object):

    def bond_length(self):
        return self._bond_length

    def compute_sigma_bond_angles(self):

        for atom in enumerate(self):
            for nnatom in enumerate(atom.NN):
                nnv = Vector(nnatom.r, p0=atom.r.p0)
                atom.bonds.vectors.append(nnv)

    def compute_pyramidalization_angle(self):
        pass
