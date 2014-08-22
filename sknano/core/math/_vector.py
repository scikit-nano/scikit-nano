# -*- coding: utf-8 -*-
"""
==================================================================
Custom NumPy Vector class (:mod:`sknano.core.math._vector`)
==================================================================

.. currentmodule:: sknano.core.math._vector

"""
from __future__ import absolute_import, division, print_function
__docformat__ = 'restructuredtext en'

import numbers
import numpy as np

from ._point import Point

__all__ = ['Vector', 'cross', 'dot']


class Vector(np.ndarray):
    """Abstract object representation of a vector in :math:`R^n`

    Parameters
    ----------
    v : array_like, optional
        components of vector
    nd : {None, int}, optional
    p : array_like, optional
        Terminating `Point` of vector.
        :math:`x, y` coordinates of point in :math:`R^2` space.
        :math:`x, y, z` coordinates of point in :math:`R^3` space.
    p0 : array_like, optional
        Origin `Point` of vector in :math:`R^n` space.
    dtype : data-type, optional
    copy : bool, optional

    Notes
    -----
    .. todo::

       add new methods for coordinate transformations

    """
    __array_priority__ = 10.0
    _verbose = False

    def __new__(cls, v=None, nd=None, p=None, p0=None, dtype=None, copy=True,
                verbose=False):
        if Vector._verbose:
            print('In __new__\n'
                  'cls: {}\n'.format(cls) +
                  'v: {}\n'.format(v) +
                  'type(v): {}\n'.format(type(v)))

        if isinstance(v, Vector):
            if dtype is None:
                intype = v.dtype
            else:
                intype = np.dtype(dtype)

            vec = v.view(cls)
            if intype != v.dtype:
                return vec.astype(intype)

            if copy:
                return vec.copy()
            else:
                return vec

        dtype = np.dtype(dtype)

        if isinstance(v, (tuple, list, np.ndarray)):
            try:
                for i, coord in enumerate(v[:]):
                    if coord is None:
                        v[i] = 0.0
            except TypeError:
                v = np.zeros(len(v), dtype=dtype)
            else:
                v = np.asarray(v, dtype=dtype)
            nd = len(v)

            if p0 is None:
                p0 = Point(nd=nd, dtype=dtype)
            else:
                p0 = Point(p0, nd=nd, dtype=dtype)
            p = p0 + v
        else:
            if nd is None or not isinstance(nd, numbers.Number):
                nd = 3
            if p is None:
                p = Point(nd=nd, dtype=dtype)
            else:
                p = Point(p, nd=nd, dtype=dtype)

            if p0 is None:
                p0 = Point(nd=nd, dtype=dtype)
            else:
                p0 = Point(p0, nd=nd, dtype=dtype)

            v = p - p0

        arr = np.array(v, dtype=dtype, copy=copy).view(cls)
        vec = np.ndarray.__new__(cls, arr.shape, arr.dtype, buffer=arr)
        #vec = \
        #    super(Vector, cls).__new__(cls, arr.shape, arr.dtype, buffer=arr)

        vec.nd = nd
        vec._p = p
        vec._p0 = p0
        if nd == 2:
            vec.x, vec.y = vec
        elif nd == 3:
            vec.x, vec.y, vec.z = vec

        #vec.verbose = verbose

        return vec

    def __array_finalize__(self, obj):
        if obj is None:
            return None

        self.nd = len(obj)
        #self.verbose = getattr(obj, 'verbose', False)
        if Vector._verbose:
            print('In __array_finalize__\n'
                  'self: {}\n'.format(self) +
                  'type(self): {}\n'.format(type(self)) +
                  'obj: {}\n'.format(obj) +
                  'type(obj): {}\n'.format(type(obj)) +
                  'obj.shape: {}\n'.format(obj.shape))

        if self.nd == 2:
            self.x, self.y = obj
        elif self.nd == 3:
            self.x, self.y, self.z = obj

        self._p = getattr(obj, 'p', None)
        self._p0 = getattr(obj, 'p0', None)

    def __array_prepare__(self, obj, context=None):
        if Vector._verbose:
            print('In __array_prepare__\n'
                  'self: {}\n'.format(self) +
                  'type(self): {}\n'.format(type(self)) +
                  'obj: {}\n'.format(obj) +
                  'type(obj): {}\n'.format(type(obj)) +
                  'context: {}\n'.format(context))

        if self.__array_priority__ >= Vector.__array_priority__:
            ret = obj if isinstance(obj, type(self)) else obj.view(type(self))
        else:
            ret = obj.view(Vector)

        if context is None:
            return ret
        return super(Vector, self).__array_prepare__(obj, context)

    def __array_wrap__(self, obj, context=None):
        if Vector._verbose:
            print('In __array_wrap__\n'
                  'self: {}\n'.format(self) +
                  'type(self): {}\n'.format(type(self)) +
                  'type(self): {}\n'.format(type(self)) +
                  'self.p: {}\n'.format(self.p) +
                  'type(self.p): {}\n'.format(type(self.p)) +
                  'self.p0: {}\n'.format(self.p0) +
                  'type(self.p0): {}\n'.format(type(self.p0)) +
                  'obj: {}\n'.format(obj) +
                  'type(obj): {}\n'.format(type(obj)) +
                  'obj.p: {}\n'.format(obj.p) +
                  'type(obj.p): {}\n'.format(type(obj.p)) +
                  'obj.p0: {}\n'.format(obj.p0) +
                  'type(obj.p0): {}\n'.format(type(obj.p0)) +
                  'context: {}\n'.format(context))

        if obj.shape == ():
            return obj[()]
        else:
            #return np.ndarray.__array_wrap__(self, obj, context)
            res = super(Vector, self).__array_wrap__(obj, context)
            #res._p0 = Point(self.p0)
            #res._p = Point(self[:] + res._p0)
            res = Vector(res.__array__(), p0=self.p0)
            if Vector._verbose:
                print('In __array_wrap__\n'
                      'self: {}\n'.format(self) +
                      'type(self): {}\n'.format(type(self)) +
                      'self.p: {}\n'.format(self.p) +
                      'type(self.p): {}\n'.format(type(self.p)) +
                      'self.p0: {}\n'.format(self.p0) +
                      'type(self.p0): {}\n'.format(type(self.p0)) +
                      'obj: {}\n'.format(obj) +
                      'type(obj): {}\n'.format(type(obj)) +
                      'obj.p: {}\n'.format(obj.p) +
                      'type(obj.p): {}\n'.format(type(obj.p)) +
                      'obj.p0: {}\n'.format(obj.p0) +
                      'type(obj.p0): {}\n'.format(type(obj.p0)) +
                      'res: {}\n'.format(res) +
                      'type(res): {}\n'.format(type(res)) +
                      'res.p: {}\n'.format(res.p) +
                      'type(res.p): {}\n'.format(type(res.p)) +
                      'res.p0: {}\n'.format(res.p0) +
                      'type(res.p0): {}\n'.format(type(res.p0)) +
                      'context: {}\n'.format(context))
            #return super(Vector, self).__array_wrap__(obj, context)
            return res

    def __str__(self):
        return repr(self)

    def __repr__(self):
        reprstr = \
            "Vector({v!r}, nd={nd!r}, p={p!r}, p0={p0!r}, dtype={dtype!r})"
        parameters = dict(v=self.__array__(),
                          nd=getattr(self, 'nd', None),
                          p=getattr(self, 'p', None),
                          p0=getattr(self, 'p0', None),
                          dtype=self.dtype)
        return reprstr.format(**parameters)

    def __getattr__(self, name):
        #if Vector._verbose and name.find('ufunc') > 0:
        #    print('In __getattr__\n'
        #          'self: {}\n'.format(self.__array__()) +
        #          'name: {}\n'.format(name))
        try:
            nd = len(self)
            if nd == 2 and name in ('x', 'y'):
                if name == 'x':
                    return self[0]
                else:
                    return self[1]
            elif nd == 3 and name in ('x', 'y', 'z'):
                if name == 'x':
                    return self[0]
                elif name == 'y':
                    return self[1]
                else:
                    return self[2]
        except TypeError:
            pass
        return super(Vector, self).__getattribute__(name)

    def __setattr__(self, name, value):
        #nd = len(self)
        nd = getattr(self, 'nd', None)

        #if Vector._verbose:
        #    print('In __setattr__\n'
        #          'self: {}\n'.format(self) +
        #          'type(self): {}\n'.format(type(self)) +
        #          'name: {}\n'.format(name) +
        #          'value: {}\n'.format(value))
        if nd is not None and nd == 2 and name in ('x', 'y'):
            if name == 'x':
                self[0] = value
                try:
                    self._p.x = self._p0.x + value
                except AttributeError:
                    pass
            else:
                self[1] = value
                try:
                    self._p.y = self._p0.y + value
                except AttributeError:
                    pass
        elif nd is not None and nd == 3 and name in ('x', 'y', 'z'):
            if name == 'x':
                self[0] = value
                try:
                    self._p.x = self._p0.x + value
                except AttributeError:
                    pass
            elif name == 'y':
                self[1] = value
                try:
                    self._p.y = self._p0.y + value
                except AttributeError:
                    pass
            else:
                self[2] = value
                try:
                    self._p.z = self._p0.z + value
                except AttributeError:
                    pass
        else:
            super(Vector, self).__setattr__(name, value)
            #np.ndarray.__setattr__(self, name, value)

    #def __getitem__(self, key):
    #    super(Vector, self).__getitem__(key)

    #def __setitem__(self, key, value):
    #    super(Vector, self).__setitem__(key, value)

    def __add__(self, other):
        if Vector._verbose:
            print('in __add__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        return super(Vector, self).__add__(other)

    def __mul__(self, other):
        if Vector._verbose:
            print('in __mul__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        if isinstance(other, (np.ndarray, list, tuple)):
            return np.multiply(np.asarray(self), np.asarray(other))
        return NotImplemented
        #return super(Vector, self).__mul__(other)

    def __sub__(self, other):
        if Vector._verbose:
            print('in __sub__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        return super(Vector, self).__sub__(np.asarray(other))

    def __radd__(self, other):
        if Vector._verbose:
            print('in __radd__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        return self.__add__(other)

    def __rmul__(self, other):
        if Vector._verbose:
            print('in __rmul__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        return np.multiply(np.asarray(other), np.asarray(self))

    def __imul__(self, other):
        if Vector._verbose:
            print('in __imul__\n'
                  'self: {}\n'.format(self) +
                  'other: {}\n'.format(other))
        self[:] = self * other
        return self

    @property
    def p(self):
        """Terminating :class:`Point` of vector."""
        return self._p

    @p.setter
    def p(self, value=np.ndarray):
        """Set new terminating :class:`Point` of vector."""
        self._p[:] = value
        self[:] = self._p - self._p0

    @property
    def p0(self):
        """Origin :class:`Point` of vector."""
        return self._p0

    @p0.setter
    def p0(self, value=np.ndarray):
        """Set new origin :class:`Point` of vector."""
        self._p0[:] = value
        self[:] = self._p - self._p0

    def dot(self, other, out=None):
        """Computes dot product of two vectors."""
        return self.__array__().dot(other.__array__())

    def rezero_components(self, epsilon=1.0e-10):
        """Re-zero `Vector` coordinates near zero.

        Set `Vector` coordinates with absolute value less than `epsilon` to
        zero.

        Parameters
        ----------
        epsilon : float, optional
            Smallest allowed absolute value of any :math:`x,y,z` coordinate.

        """
        self[np.where(np.abs(self) <= epsilon)] = 0.0


def cross(v1, v2):
    """Vector cross product of two vectors.

    Parameters
    ----------
    a, b : `Vector`

    Returns
    -------
    np.number or `Vector`

    """
    val = np.cross(np.asarray(v1), np.asarray(v2))
    if val.shape == ():
        return val[()]
    else:
        return Vector(val, p0=v1.p0)


def dot(v1, v2):
    """Dot product of two vectors.

    Parameters
    ----------
    a, b : `Vector`

    Returns
    -------
    val

    """
    return np.dot(np.asarray(v1), np.asarray(v2))