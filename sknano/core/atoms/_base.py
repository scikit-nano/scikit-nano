# -*- coding: utf-8 -*-
"""
==============================================================================
Base classes for atoms package (:mod:`sknano.core.atoms._base`)
==============================================================================

.. currentmodule:: sknano.core.atoms._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence, OrderedDict
from functools import total_ordering

import copy

__all__ = ['AtomList', 'BondList', 'BondDict']


@total_ordering
class AtomList(MutableSequence):
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

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        self._data = []

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
                                'is not a valid `AtomList` constructor '
                                'argument.\n `atoms` must be `None`, a list '
                                'of `Atom` objects, or an `AtomList` '
                                'instance.')

    def __eq__(self, other):
        return self[:] == other[:]

    def __lt__(self, other):
        return self[:] < other[:]

    def __delitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Delete list element `self[index]`.

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

    def clear(self):
        self._data = []


@total_ordering
class BondList(MutableSequence):
    """Base class for collection of atom `Bonds`.

    Parameters
    ----------
    bonds : {None, sequence, `Bonds`}, optional
        if not `None`, then a list of `Bond` instance objects
    copylist : bool, optional
        perform shallow copy of bonds list
    deepcopy : bool, optional
        perform deepcopy of bonds list

    """
    def __init__(self, bonds=None, copylist=True, deepcopy=False):
        self._data = []

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
                                'of `Bond` objects, or a `Bonds` object '
                                'instance.')

    def __eq__(self, other):
        return self[:] == other[:]

    def __lt__(self, other):
        return self[:] < other[:]

    def __delitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Delete list element `self[index]`.

        Parameters
        ----------
        index : int
            index of target list element

        """
        del self._data[index]

    def __getitem__(self, index):
        """Concrete implementation of @abstractmethod.

        Get `Bond` object instance at list element `self[index]`

        Parameters
        ----------
        index : int
            index of target list element

        Returns
        -------
        `Bond`
            `Bond` instance at target `self[index]`

        """
        return self._data[index]

    def __setitem__(self, index, bond):
        """Concrete implementation of @abstractmethod.

        set target list element `self[index] = bond`

        Parameters
        ----------
        index : int
            index of target list element
        bond : `Bond`
            `Bond` instance object to set target list element to.

        """
        self._data[index] = bond

    def __len__(self):
        """Concrete implementation of @abstractmethod.

        Returns
        -------
        int
            length of `self` list.

        """
        return len(self._data)

    def insert(self, index, bond):
        """Concrete implementation of @abstractmethod.

        insert `Bond` instance at target list `index`

        Parameters
        ----------
        index : int
            index of target list element
        bond : `Bond`
            `Bond` object instance to set target list element to

        """
        self._data.insert(index, bond)

    @property
    def data(self):
        """Return the list of `Bonds` objects"""
        return self._data

    def clear(self):
        self._data = []


class BondDict(OrderedDict):
    pass
