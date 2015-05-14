# -*- coding: utf-8 -*-
"""
====================================================================
Custom NumPy Quaternion class (:mod:`sknano.core.math._quaternion`)
====================================================================

.. currentmodule:: sknano.core.math._quaternion

"""
from __future__ import absolute_import, division, print_function
from __future__ import unicode_literals
__docformat__ = 'restructuredtext en'

# import copy
import numbers
import warnings

import numpy as np
np.seterr(all='warn')

# from ._point import Point
# from ._transforms import rotate, transformation_matrix
from ._vector import Vector

__all__ = ['Quaternion']


class Quaternion(np.ndarray):
    """Abstract object representation of a quaternion.

    Parameters
    ----------
    dtype : data-type, optional
    copy : bool, optional

    Examples
    --------

    """
    __array_priority__ = 15.0
    _verbosity = 0

    def __new__(cls, H, dtype=None, copy=True):

        if isinstance(H, Quaternion):
            intype = H.dtype
            if dtype is None:
                dtype = intype
            if intype == dtype and not copy:
                return H
            return H.astype(dtype)

        if isinstance(H, np.ndarray):
            if dtype is None:
                intype = H.dtype
            else:
                intype = np.dtype(dtype)

            Hview = H.view(cls)

            if intype != H.dtype:
                return Hview.astype(intype)

            if copy:
                return Hview.copy()
            else:
                return Hview

        arr = np.array(H, dtype=dtype, copy=copy)
        if len(arr) != 4:
            raise ValueError("Expected array_like input: [w, x, y, z]")
        ret = np.ndarray.__new__(cls, arr.shape, arr.dtype, buffer=arr)

        return ret

    def __array_finalize__(self, obj):

        if obj is None:
            return None

        if len(self) == 4:
            return
        else:
            raise ValueError("Expected array_like input: [w, x, y, z]")

    def __str__(self):
        fmtstr = "Quaternion:\nw = {w!r}\nx = {x!r}\ny = {y!r}\nz = {z!r}"
        return fmtstr.format(**dict(w=self.w, x=self.x, y=self.y, z=self.z))

    def __repr__(self):
        return "Quaternion({!r})".format(self.tolist())

    def tolist(self):
        return np.around(self.__array__(), decimals=6).tolist()

    def __eq__(self, other):
        if isinstance(other, Quaternion) and \
            (self is other or np.allclose(self.__array__(),
                                          other.__array__())):
            return True
        return False

    def __lt__(self, other):
        if isinstance(other, Quaternion):
            return self.norm < other.norm
        return False

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not (self < other or self == other)

    def __ge__(self, other):
        return not (self < other)

    def __ne__(self, other):
        return not (self == other)

    def __add__(self, other):
        if isinstance(other, Quaternion):
            return super().__add__(other)

        if np.isscalar(other):
            return Quaternion.from_components(w=self.w + other,
                                              x=self.x, y=self.y, z=self.z)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            w = self.w * other.w - self.v.dot(other.v)
            v = self.v.cross(other.v) + self.w * other.v + other.w * self.v
            return Quaternion.from_components(w=w, x=v.x, y=v.y, z=v.z)

        if np.isscalar(other):
            return Quaternion(self.__array__() * other)

        return NotImplemented

    def __truediv__(self, other):
        if np.isscalar(other):
            return Quaternion(self.__array__() / other)
        return NotImplemented

    def __floordiv__(self, other):
        if np.isscalar(other):
            return Quaternion(self.__array__() // other)
        return NotImplemented

    def __mod__(self, other):
        return NotImplemented

    def __divmod__(self, other):
        return NotImplemented

    def __pow__(self, other):
        return NotImplemented

    def __lshift__(self, other):
        return NotImplemented

    def __rshift__(self, other):
        return NotImplemented

    def __and__(self, other):
        return NotImplemented

    def __xor__(self, other):
        return NotImplemented

    def __or__(self, other):
        return NotImplemented

    __radd__ = __add__
    __rmul__ = __mul__

    def __rtruediv__(self, other):
        return NotImplemented

    def __rfloordiv__(self, other):
        return NotImplemented

    @classmethod
    def from_components(cls, w, x, y, z):
        return cls([w, x, y, z])

    @property
    def w(self):
        """Real component of `Quaternion`."""
        return self.real

    @property
    def v(self):
        """Vector of imaginary components of `Quaternion`."""
        return Vector(self.imag)

    @property
    def real(self):
        """Real part :math:`w` of `Quaternion`."""
        return self[0]

    @property
    def imag(self):
        """:class:`~python:list` of `Quaternion` imaginary components \
            :math:`x,y,z`, with imaginary part \
            :math:`x\\mathbf{i}+y\\mathbf{j}+z\\mathbf{k}`."""
        return self.__array__()[1:].tolist()

    @property
    def x(self):
        return self[1]

    @property
    def y(self):
        return self[2]

    @property
    def z(self):
        return self[3]

    @property
    def axis(self):
        """Rotation axis."""
        return Vector([self.x, self.y, self.z]) / np.sin(np.arccos(self.w))

    @property
    def angle(self):
        """Rotation angle."""
        return 2 * np.arccos(self.w)

    @property
    def norm(self):
        return np.sqrt((self.__array__() ** 2).sum())

    @property
    def conjugate(self):
        return Quaternion([self.w, -self.x, -self.y, -self.z])

    @property
    def inverse(self):
        return self.conjugate / self.norm

    @property
    def unit_quaternion(self):
        return self / self.norm
