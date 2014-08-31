# -*- coding: utf-8 -*-
"""
===============================================================================
Class container for neighbor atoms (:mod:`sknano.core.atoms._neighbor_atoms`)
===============================================================================

.. currentmodule:: sknano.core.atoms._neighbor_atoms

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence
from functools import total_ordering

#import numpy as np

#from sknano.core.math import Vector

__all__ = ['NeighborAtoms']


@total_ordering
class NeighborAtoms(MutableSequence):
    """An eXtended `Atoms` class for structure analysis.

    Parameters
    ----------
    atoms : {None, sequence, `NeighborAtoms`}, optional
        if not `None`, then a list of `XAtom` instance objects or an
        existing `NeighborAtoms` instance object.
    copylist : bool, optional
        perform shallow copy of atoms list
    deepcopy : bool, optional
        perform deepcopy of atoms list
    use_kdtree : bool, optional
        use :py:class:`~scipy:scipy.spatial.KDTree` to perform
        nearest-neighbor analysis.

    """

    def __init__(self, atoms=None, copylist=True, deepcopy=False):
        self._data = []

    def __str__(self):
        """Return a nice string representation of `NeighborAtoms`."""
        return "NeighborAtoms(atoms={!s})".format(self._data)

    def __repr__(self):
        """Return the canonical string representation of `NeighborAtoms`."""
        return "NeighborAtoms(atoms={!r})".format(self._data)

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

    def __eq__(self, other):
        return self[:] == other[:]

    def __lt__(self, other):
        return self[:] < other[:]

    @property
    def data(self):
        """Return the list of `NeighborAtom` objects"""
        return self._data
