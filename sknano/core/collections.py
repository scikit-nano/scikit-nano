# -*- coding: utf-8 -*-
"""
==============================================================================
Custom container datatypes (:mod:`sknano.core.collections`)
==============================================================================

.. currentmodule:: sknano.core.collections

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

from abc import abstractmethod

try:
    from collections.abc import Mapping, MutableSequence, Set
except ImportError:
    from collections import Mapping, MutableSequence, Set

from .meta import BaseClass

__all__ = ['ListBasedSet', 'UserList', 'frozendict']


# class AttrDict(dict):
#     def __init__(self, *args, **kwargs):
#         super(AttrDict, self).__init__(*args, **kwargs)
#         self.__dict__ = self


class ListBasedSet(Set):
    """Alternate set implementation favoring space over speed and not
    requiring the set elements to be hashable.

    Parameters
    ----------
    iterable : :class:`~python:collections.Iterable`

    """
    def __init__(self, iterable):
        self.elements = elements = []
        for value in iterable:
            if value not in elements:
                elements.append(value)

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, value):
        return value in self.elements

    def __len__(self):
        return len(self.elements)


class UserList(MutableSequence, BaseClass):
    """Modified implementation of :class:`~python:collections.UserList`.

    Sub-class of both :class:`~python:MutableSequence` and
    :class:`~sknano.core.BaseClass` and is itself an
    abstract class requiring concrete implementations of
    :attr:`~sknano.core.UserList.__item_class__` and
    :meth:`~sknano.core.BaseClass.todict`.

    Constructor takes additional variable arguments and
    variable keyword arguments.

    Parameters
    ----------
    initlist : {None, sequence, :class:`UserList`}, optional
        if not `None`, then a list or an instance of :class:`UserList`

    """
    def __init__(self, *args, initlist=None, **kwargs):
        self.data = []
        if initlist is not None:
            if type(initlist) == type(self.data):
                self.data[:] = initlist
            elif isinstance(initlist, UserList):
                self.data[:] = initlist.data[:]
            else:
                self.data = list(initlist)
        self.kwargs = kwargs
        super().__init__(*args, **kwargs)
        self.fmtstr = "{initlist!r}"

    @property
    @abstractmethod
    def __item_class__(self):
        """The list item's class."""
        return NotImplementedError

    def __cast(self, other):
        return other.data if isinstance(other, UserList) else other

    def __cast_item(self, item):
        if not isinstance(item, self.__item_class__):
            item = self.__item_class__(item)
        return item

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

    def __contains__(self, item):
        return item in self.data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, i):
        item = self.data[i]
        if isinstance(i, slice) and isinstance(item, list):
            item = self.__class__(item, **self.kwargs)
        elif not isinstance(item, self.__item_class__):
            item = self.__item_class__(item)
        return item

    def __setitem__(self, i, item):
        if not isinstance(item, (self.__class__, self.__item_class__)):
            if isinstance(i, slice):
                item = self.__class__(item, **self.kwargs)
            else:
                try:
                    item = self.__item_class__(**item)
                except TypeError:
                    item = self.__item_class__(item)
        self.data[i] = item

    def __delitem__(self, i):
        del self.data[i]

    def __add__(self, other):
        if isinstance(other, UserList):
            return self.__class__(self.data + other.data, **self.kwargs)
        elif isinstance(other, type(self.data)):
            return self.__class__(self.data + other, **self.kwargs)
        return self.__class__(self.data + list(other), **self.kwargs)

    def __radd__(self, other):
        if isinstance(other, UserList):
            return self.__class__(other.data + self.data, **self.kwargs)
        elif isinstance(other, type(self.data)):
            return self.__class__(other + self.data, **self.kwargs)
        return self.__class__(list(other) + self.data, **self.kwargs)

    def __iadd__(self, other):
        if isinstance(other, UserList):
            self.data += other.data
        elif isinstance(other, type(self.data)):
            self.data += other
        else:
            self.data += list(other)
        return self

    def __mul__(self, n):
        return self.__class__(self.data * n, **self.kwargs)

    __rmul__ = __mul__

    def __imul__(self, n):
        self.data *= n
        return self

    def append(self, item):
        self.data.append(self.__cast_item(item))

    def insert(self, i, item):
        self.data.insert(i, self.__cast_item(item))

    def pop(self, i=-1):
        return self.data.pop(i)

    def remove(self, item):
        self.data.remove(item)

    def clear(self):
        try:
            self.data.clear()
        except AttributeError:
            del self.data[:]

    def copy(self):
        return self.__class__(self.data[:], **self.kwargs)

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

    def todict(self):
        """Return :class:`~python:dict` of constructor parameters."""
        return dict(initlist=self.data)


class frozendict(Mapping):
    """A simple implementation of a read-only *frozen* `dict`."""
    def __init__(self, data):
        self._data = data

    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)
