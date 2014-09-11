# -*- coding: utf-8 -*-
"""
==============================================================================
Base classes for atoms package (:mod:`sknano.core.atoms._base`)
==============================================================================

.. currentmodule:: sknano.core.atoms._base

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence
from functools import total_ordering

import copy

__all__ = ['UserList', 'AtomList', 'BondList']


@total_ordering
class UserList(MutableSequence):
    """Base class around list objects.

    Modified implementation of :class:`~python:collections.UserList`.

    Parameters
    ----------
    initlist : {None, sequence, UserList}, optional
        if not `None`, then a list or an instance of `UserList`
    copylist : bool, optional
        perform shallow copy of list items
    deepcopy : bool, optional
        perform deepcopy of list items
    """
    def __init__(self, initlist=None, copylist=True, deepcopy=False):
        self.data = []

        if initlist is not None:
            if type(initlist) == type(self.data):
                if copylist and not deepcopy:
                    self.data.extend(initlist[:])
                elif deepcopy:
                    self.data.extend(copy.deepcopy(initlist))
                else:
                    self.data.extend(initlist)
            elif isinstance(initlist, UserList):
                if copylist and not deepcopy:
                    self.data.extend(initlist.data[:])
                elif deepcopy:
                    self.data.extend(copy.deepcopy(initlist.data))
                else:
                    self.data.extend(initlist.data)
            else:
                self.data = list(initlist)

    def __repr__(self):
        return repr(self.data)

    def __lt__(self, other):
        return self.data < self.__cast(other)

    def __le__(self, other):
        return self.data <= self.__cast(other)

    def __eq__(self, other):
        return self.data == self.__cast(other)

    def __ne__(self, other):
        return self.data != self.__cast(other)

    def __gt__(self, other):
        return self.data > self.__cast(other)

    def __ge__(self, other):
        return self.data >= self.__cast(other)

    def __cast(self, other):
        return other.data if isinstance(other, UserList) else other

    def __contains__(self, item):
        return item in self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, i):
        return self.data[i]

    def __setitem__(self, i, item):
        self.data[i] = item

    def __delitem__(self, i):
        del self.data[i]

    def __add__(self, other):
        if isinstance(other, UserList):
            return self.__class__(self.data + other.data)
        elif isinstance(other, type(self.data)):
            return self.__class__(self.data + other)
        return self.__class__(self.data + list(other))

    def __radd__(self, other):
        if isinstance(other, UserList):
            return self.__class__(other.data + self.data)
        elif isinstance(other, type(self.data)):
            return self.__class__(other + self.data)
        return self.__class__(list(other) + self.data)

    def __iadd__(self, other):
        if isinstance(other, UserList):
            self.data += other.data
        elif isinstance(other, type(self.data)):
            self.data += other
        else:
            self.data += list(other)
        return self

    def __mul__(self, n):
        return self.__class__(self.data*n)

    __rmul__ = __mul__

    def __imul__(self, n):
        self.data *= n
        return self

    def append(self, item):
        self.data.append(item)

    def insert(self, i, item):
        self.data.insert(i, item)

    def pop(self, i=-1):
        return self.data.pop(i)

    def remove(self, item):
        self.data.remove(item)

    def clear(self):
        #self.data.clear()
        del self.data[:]

    def copy(self):
        return self.__class__(self)

    def count(self, item):
        return self.data.count(item)

    def index(self, item, *args):
        return self.data.index(item, *args)

    def reverse(self):
        self.data.reverse()

    def sort(self, *args, **kwds):
        self.data.sort(*args, **kwds)

    def extend(self, other):
        if isinstance(other, UserList):
            self.data.extend(other.data)
        else:
            self.data.extend(other)


class AtomList(UserList):
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
        super(AtomList, self).__init__(initlist=atoms,
                                       copylist=copylist,
                                       deepcopy=deepcopy)


class BondList(UserList):
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
        super(BondList, self).__init__(initlist=bonds,
                                       copylist=copylist,
                                       deepcopy=deepcopy)
