# -*- coding: utf-8 -*-
"""
==============================================================================
Custom list class (:mod:`sknano.core._user_list`)
==============================================================================

.. currentmodule:: sknano.core._user_list

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

from collections import MutableSequence
from functools import total_ordering

import copy

__all__ = ['UserList']


@total_ordering
class UserList(MutableSequence):
    """Modified implementation of :class:`~python:collections.UserList`.

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
